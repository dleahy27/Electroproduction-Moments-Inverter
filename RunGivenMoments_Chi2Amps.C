// Proton-only: Fit partial-wave magnitudes + phases to measured H^alpha_{L M}
// using explicit chi2 minimisation (ROOT Minuit2).
//
// This macro implements the same CG-sum structure as E2S0AmpLoader::CGMatrixReflectivity
// and the same amplitude/phase term construction as E2S0AmpLoader::BruTermCircle,
// but evaluates everything numerically (no BruFit Setup required).
//
// Outputs a ROOT file containing:
//   - log_val = log10(chi2)
//   - predicted H_alpha_L_M (for alpha=0..8, L=0..2, M=0..L)
//   - fitted amplitudes a/b_{T/L}_{l}_{m} and phases a/bphi_{T/L}_{l}_{m}
//
// You can run the macro either with a Q2 value (mapped to the nearest SDME table bin)
// or with a ROOT file containing the observed moments.

// General ROOT Stuff
#include "TBenchmark.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"

// Minimiser stuff
#include "Math/Factory.h"
#include "Math/IFunction.h"
#include "Math/Minimizer.h"
#include <Math/SpecFuncMathMore.h>

// Multiprocessing stuff
#include "ROOT/TProcessExecutor.hxx"
#include "TFileMerger.h"
#include "TSystem.h"
#include <thread>

// Gneral C++ stuff
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chi2_amp_fit {

// Fast sin/cos helper (uses ::sincos when available)
inline void FastSinCos(double x, double& s, double& c) {
#if defined(__GLIBC__) || defined(__APPLE__)
  ::sincos(x, &s, &c);
#else
  s = std::sin(x);
  c = std::cos(x);
#endif
}


// -------------------------
// Physics / indexing helpers
// -------------------------
inline double ClebschGordan(int l1, int l2, int l3, int m1, int m2, int m3) {
  // Matches E2S0AmpLoader::ClebschGordan
  using ROOT::Math::wigner_3j;
  return (TMath::Power(-1.0, l1 - l2 + m3) * TMath::Sqrt(2.0 * l3 + 1.0) *
          wigner_3j(2 * l1, 2 * l2, 2 * l3, 2 * m1, 2 * m2, -2 * m3));
}

// inline double WrapToPi(double phi) {
//   constexpr double kTwoPi = 2.0 * TMath::Pi();
//   double x = std::fmod(phi + TMath::Pi(), kTwoPi);
//   if (x < 0) x += kTwoPi;
//   return x - TMath::Pi();
// }

// -------------------------
// Observed moments container
// -------------------------
struct ObservedMoment {
  int alpha = 0;
  int L = 0;
  int M = 0;
  double value = 0.0;
  double sigma = 1.0;
  bool isMixed04 = false;
  std::string name;
};

enum class ObservedMomentsSource : uint8_t { kQ2Table = 0, kRootFile };

struct ObservedInput {
  ObservedMomentsSource source = ObservedMomentsSource::kQ2Table;
  int q2bin = 1;
  std::string momentsFile;
  std::string momentsTree = "syntheticMoments";
};

// -------------------------
// Parameter bookkeeping
// -------------------------
struct ParDef {
  std::string name;
  double init = 0.0;
  double step = 1e-2;
  double low = -1.0;
  double high = 1.0;
  bool fixed = false;
  bool isPhase = false;
};

struct FitConfig {
  // Partial Wave options
  int lmax = 1;
  int mmax = 1;
  bool useNegRef = true;   // nRefl==2 in loader
  bool onlyEven = false;
  bool negm = true;

  double epsR4 = 1.0;
  double rh04MixCoeff = -1.0;

  ObservedInput observedInput;

  // Minimiser options
  unsigned nStarts = 10000;
  unsigned maxCalls = 50000;
  unsigned maxIters = 50000;
  double tolerance = 1e-8;
  int strategy = 2;
  int printLevel = 0;
  bool runHesse = true;
  uint32_t randomSeed = 0;

  // Normalisation options
  double H0H4NormTarget = 2.0;
  double H0H4NormSigma  = 1e-3;

  // Optional per-FCN-call trace. Keep disabled by default for speed.
  unsigned recordEvery = 0;

  // Multiprocessing configs
  unsigned int cores = 1;
  bool verbose = true; // Print out each each tenth of progress
};

// -------------------------
// Model term list (precomputed)
// -------------------------
enum class TrigKind : uint8_t { kNone = 0, kCos, kSin };

struct Term {
  double coeff = 0.0;
  int idxMag1 = -1;
  int idxPhi1 = -1;
  int idxMag2 = -1;
  int idxPhi2 = -1;
  bool useTrig = false;
  TrigKind trig = TrigKind::kNone;
};

struct MomentModel {
  int alpha = 0;
  int L = 0;
  int M = 0;
  std::string name;
  std::vector<Term> terms;
};

// -------------------------
// Build parameter list: a/b, T/L, l,m magnitudes + phases
// -------------------------
static std::vector<std::pair<int,int>> EnumerateWaves(int lmax, int mmax, bool negm, bool onlyEven) {
  std::vector<std::pair<int,int>> waves;
  for (int l = 0; l <= lmax; ++l) {
    if (onlyEven && (l % 2 == 1)) continue;
    for (int m = -l; m <= l; ++m) {
      if (std::abs(m) > mmax) continue;
      if (!negm && m < 0) continue;
      waves.emplace_back(l, m);
    }
  }
  return waves;
}

static std::string MagName(char refl, char orient, int l, int m) {
  std::ostringstream os;
  std::string _m;
  std::string _l = std::to_string(l);
  if (m==-1)
  {
    _m="m1";
  } else
  {
    _m = std::to_string(m);
  }
  os << refl << '_' << orient << '_' << _l << '_' << _m;
  return os.str();
}

static std::string PhiName(char refl, char orient, int l, int m) {
  std::ostringstream os;
  std::string _m;
  std::string _l = std::to_string(l);
  if (m==-1)
  {
    _m="m1";
  } else
  {
    _m = std::to_string(m);
  }
  os << refl << "phi_" << orient << '_' << _l << '_' << _m;
  return os.str();
}

static std::vector<ParDef> BuildAmplitudePhaseParameters(const FitConfig& cfg) {
  // Use the same naming conventions as ElectroTwoSpin0AmpLoader:
  // magnitudes: a_T_l_m, a_L_l_m, b_T_l_m, b_L_l_m
  // phases:     aphi_T_l_m, aphi_L_l_m, bphi_T_l_m, bphi_L_l_m

  constexpr double kPi = TMath::Pi();

  std::vector<ParDef> pars;
  auto waves = EnumerateWaves(cfg.lmax, cfg.mmax, cfg.negm, cfg.onlyEven);

  auto addMag = [&](char refl, char orient, int l, int m) {
    ParDef p;
    p.name = MagName(refl, orient, l, m);
    p.init = 0.0;
    p.step = 1e-3;
    p.low = 0.0;
    p.high = 1.0;
    p.fixed = false;
    p.isPhase = false;
    pars.push_back(p);
  };

  auto addPhi = [&](char refl, char orient, int l, int m) {
    ParDef p;
    p.name = PhiName(refl, orient, l, m);
    p.init = 0.0;
    p.step = 1e-3;
    p.low = -kPi;
    p.high = kPi;
    p.fixed = false;
    p.isPhase = true;
    pars.push_back(p);
  };

  for (auto [l, m] : waves) {
    // magnitudes
    addMag('a', 'T', l, m);
    addMag('a', 'L', l, m);
    addMag('b', 'T', l, m);
    addMag('b', 'L', l, m);

    // phases
    addPhi('a', 'T', l, m);
    addPhi('a', 'L', l, m);
    addPhi('b', 'T', l, m);
    addPhi('b', 'L', l, m);
  }

  // Apply the same proton-only fixing/zeroing as RunGivenMoments.C
  auto fixTo = [&](const std::string& name, double val) {
    for (auto& p : pars) {
      if (p.name == name) {
        p.init = val;
        p.fixed = true;
        p.low = val;
        p.high = val;
        p.step = 0.0;
        return;
      }
    }
    throw std::runtime_error("Parameter not found to fix: " + name);
  };

  // Fix reference phases: D+1 real (same as original macro)
  fixTo("aphi_T_1_1", 0.0);
  fixTo("bphi_T_1_1", 0.0);
  // fixTo("aphi_L_1_1", 0.0);
  // fixTo("bphi_L_1_1", 0.0);

  // 0 by construction
  fixTo("a_L_0_0", 0.0);
  fixTo("a_L_1_0", 0.0);
  // fixTo("a_L_1_1", 0.0);
  // fixTo("a_L_1_m1", 0.0);
  // fixTo("b_L_1_1", 0.0);
  // fixTo("b_L_1_m1", 0.0);
  // fixTo("b_T_1_1", 0.0);
  // fixTo("b_T_1_m1", 0.0);
  // fixTo("b_T_1_0", 0.0);
  fixTo("aphi_L_0_0", 0.0);
  fixTo("aphi_L_1_0", 0.0);
  // fixTo("aphi_L_1_1", 0.0);
  // fixTo("aphi_L_1_m1", 0.0);
  // fixTo("bphi_L_1_1", 0.0);
  // fixTo("bphi_L_1_m1", 0.0);
  // fixTo("bphi_T_1_1", 0.0);
  // fixTo("bphi_T_1_m1", 0.0);
  // fixTo("bphi_T_1_0", 0.0);

  // 0 by hindsight (remove D-1 etc) - same list as RunGivenMoments.C
  fixTo("b_T_0_0", 0.0);
  fixTo("a_T_0_0", 0.0);
  fixTo("b_L_0_0", 0.0);
  fixTo("bphi_T_0_0", 0.0);
  fixTo("aphi_T_0_0", 0.0);
  fixTo("bphi_L_0_0", 0.0);
  return pars;
}

// -------------------------
// Observed moments: proton-only, from SDME table (4 Q^2 bins)
// q2bin: 0-><Q^2>=0.82, 1->1.19, 2->1.66, 3->3.06 GeV^2
// Values are (value ± stat ± syst) from your table; we combine stat+syst in quadrature.
// IMPORTANT: These are RH moments (as used by the fit). Comment out any add(...) lines you don't want.
// -------------------------
struct ValErr2 { double v{0}, stat{0}, syst{0}; };
static inline double Comb(const ValErr2& x) { return TMath::Sqrt(x.stat*x.stat + x.syst*x.syst); }

struct ProtonSDMEsTable {
  ValErr2 r00_04;      // r^04_00
  ValErr2 re_r10_04;   // Re r^04_10
  ValErr2 r1m1_04;     // r^04_1-1

  ValErr2 r1m1_1;      // r^1_1-1
  ValErr2 re_r10_1;    // Re r^1_10
  ValErr2 im_r10_2;    // Im r^2_10
  ValErr2 r00_1;       // r^1_00
  ValErr2 im_r10_3;    // Im r^3_10
  ValErr2 r00_8;       // r^8_00
  ValErr2 r11_5;       // r^5_11
  ValErr2 r1m1_5;      // r^5_1-1
  ValErr2 im_r1m1_6;   // Im r^6_1-1
  ValErr2 im_r1m1_7;   // Im r^7_1-1
  ValErr2 r11_8;       // r^8_11
  ValErr2 r1m1_8;      // r^8_1-1
  ValErr2 r11_1;       // r^1_11
  ValErr2 im_r1m1_3;   // Im r^3_1-1
  ValErr2 im_r1m1_2;   // Im r^2_1-1
  ValErr2 re_r10_5;    // Re r^5_10
  ValErr2 im_r10_6;    // Im r^6_10
  ValErr2 im_r10_7;    // Im r^7_10
  ValErr2 re_r10_8;    // Re r^8_10
  ValErr2 r00_5;       // r^5_00
};

static ProtonSDMEsTable GetProtonSDMEs_TableQ2Bin(int q2bin) {
  ProtonSDMEsTable p;
  switch (q2bin) {

    case 0: // <Q^2>=0.82
      p.r00_04     = { 0.349, 0.026, 0.061 };
      p.r1m1_1     = { 0.283, 0.023, 0.049 };
      p.im_r1m1_2  = { -0.294, 0.019, 0.038 };
      p.re_r10_5   = { 0.151, 0.028, 0.026 };
      p.im_r10_6   = { -0.149, 0.015, 0.010 };
      p.im_r10_7   = { 0.079, 0.068, 0.011 };
      p.re_r10_8   = { 0.040, 0.043, 0.011 };
      p.re_r10_04  = { 0.028, 0.028, 0.020 };
      p.re_r10_1   = { -0.037, 0.044, 0.032 };
      p.im_r10_2   = { 0.023, 0.019, 0.007 };
      p.r00_5      = { 0.121, 0.038, 0.039 };
      p.r00_1      = { -0.054, 0.039, 0.013 };
      p.im_r10_3   = { 0.002, 0.041, 0.008 };
      p.r00_8      = { 0.022, 0.079, 0.026 };
      p.r11_5      = { -0.015, 0.010, 0.007 };
      p.r1m1_5     = { 0.009, 0.011, 0.019 };
      p.im_r1m1_6  = { -0.011, 0.010, 0.013 };
      p.im_r1m1_7  = { -0.003, 0.078, 0.021 };
      p.r11_8      = { 0.019, 0.053, 0.007 };
      p.r1m1_8     = { 0.013, 0.062, 0.008 };
      p.r1m1_04    = { -0.024, 0.013, 0.021 };
      p.r11_1      = { -0.039, 0.017, 0.018 };
      p.im_r1m1_3  = { 0.021, 0.051, 0.010 };
      break;

    case 1: // <Q^2>=1.19
      p.r00_04     = { 0.368, 0.018, 0.011 };
      p.r1m1_1     = { 0.262, 0.018, 0.024 };
      p.im_r1m1_2  = { -0.255, 0.016, 0.022 };
      p.re_r10_5   = { 0.171, 0.007, 0.000 };
      p.im_r10_6   = { -0.167, 0.007, 0.003 };
      p.im_r10_7   = { 0.092, 0.038, 0.010 };
      p.re_r10_8   = { 0.020, 0.031, 0.008 };
      p.re_r10_04  = { 0.029, 0.007, 0.003 };
      p.re_r10_1   = { -0.043, 0.012, 0.006 };
      p.im_r10_2   = { 0.022, 0.012, 0.018 };
      p.r00_5      = { 0.094, 0.017, 0.017 };
      p.r00_1      = { 0.011, 0.032, 0.018 };
      p.im_r10_3   = { -0.041, 0.026, 0.005 };
      p.r00_8      = { 0.040, 0.084, 0.014 };
      p.r11_5      = { -0.011, 0.006, 0.006 };
      p.r1m1_5     = { 0.008, 0.007, 0.006 };
      p.im_r1m1_6  = { 0.002, 0.007, 0.007 };
      p.im_r1m1_7  = { 0.023, 0.056, 0.013 };
      p.r11_8      = { 0.056, 0.045, 0.004 };
      p.r1m1_8     = { 0.072, 0.053, 0.011 };
      p.r1m1_04    = { -0.014, 0.010, 0.010 };
      p.r11_1      = { -0.034, 0.013, 0.013 };
      p.im_r1m1_3  = { 0.000, 0.033, 0.004 };
      break;

    case 2: // <Q^2>=1.66
      p.r00_04     = { 0.397, 0.017, 0.018 };
      p.r1m1_1     = { 0.274, 0.019, 0.024 };
      p.im_r1m1_2  = { -0.239, 0.017, 0.011 };
      p.re_r10_5   = { 0.161, 0.006, 0.004 };
      p.im_r10_6   = { -0.167, 0.006, 0.005 };
      p.im_r10_7   = { 0.039, 0.036, 0.004 };
      p.re_r10_8   = { 0.074, 0.034, 0.002 };
      p.re_r10_04  = { 0.035, 0.007, 0.011 };
      p.re_r10_1   = { -0.036, 0.012, 0.012 };
      p.im_r10_2   = { 0.005, 0.012, 0.024 };
      p.r00_5      = { 0.057, 0.015, 0.019 };
      p.r00_1      = { 0.007, 0.031, 0.009 };
      p.im_r10_3   = { -0.074, 0.025, 0.005 };
      p.r00_8      = { 0.054, 0.086, 0.011 };
      p.r11_5      = { -0.008, 0.006, 0.011 };
      p.r1m1_5     = { -0.013, 0.007, 0.003 };
      p.im_r1m1_6  = { 0.002, 0.007, 0.004 };
      p.im_r1m1_7  = { -0.005, 0.055, 0.010 };
      p.r11_8      = { 0.051, 0.044, 0.006 };
      p.r1m1_8     = { -0.018, 0.054, 0.004 };
      p.r1m1_04    = { -0.019, 0.010, 0.003 };
      p.r11_1      = { -0.023, 0.013, 0.008 };
      p.im_r1m1_3  = { -0.031, 0.032, 0.007 };
      break;

    case 3: // <Q^2>=3.06
      p.r00_04     = { 0.454, 0.014, 0.011 };
      p.r1m1_1     = { 0.204, 0.017, 0.012 };
      p.im_r1m1_2  = { -0.197, 0.017, 0.012 };
      p.re_r10_5   = { 0.141, 0.006, 0.008 };
      p.im_r10_6   = { -0.156, 0.006, 0.010 };
      p.im_r10_7   = { 0.187, 0.034, 0.018 };
      p.re_r10_8   = { 0.098, 0.032, 0.005 };
      p.re_r10_04  = { 0.026, 0.007, 0.003 };
      p.re_r10_1   = { -0.009, 0.013, 0.010 };
      p.im_r10_2   = { 0.022, 0.013, 0.008 };
      p.r00_5      = { 0.151, 0.015, 0.007 };
      p.r00_1      = { 0.037, 0.034, 0.002 };
      p.im_r10_3   = { 0.048, 0.024, 0.006 };
      p.r00_8      = { 0.010, 0.085, 0.016 };
      p.r11_5      = { -0.021, 0.006, 0.016 };
      p.r1m1_5     = { 0.020, 0.007, 0.008 };
      p.im_r1m1_6  = { -0.010, 0.007, 0.007 };
      p.im_r1m1_7  = { -0.109, 0.047, 0.004 };
      p.r11_8      = { -0.002, 0.035, 0.005 };
      p.r1m1_8     = { 0.004, 0.045, 0.014 };
      p.r1m1_04    = { 0.001, 0.009, 0.007 };
      p.r11_1      = { -0.018, 0.012, 0.010 };
      p.im_r1m1_3  = { -0.026, 0.028, 0.005 };
      break;

    default:
      ::Error("GetProtonSDMEs_TableQ2Bin", "Invalid q2bin=%d (expected 0..3). Using bin 0.", q2bin);
      return GetProtonSDMEs_TableQ2Bin(0);
  }
  return p;
}

static int GetClosestQ2Bin(double q2Value) {
  static constexpr double kQ2Centers[] = {0.82, 1.19, 1.66, 3.06};
  int best = 0;
  double bestDiff = std::abs(q2Value - kQ2Centers[0]);
  for (int i = 1; i < 4; ++i) {
    const double diff = std::abs(q2Value - kQ2Centers[i]);
    if (diff < bestDiff) {
      best = i;
      bestDiff = diff;
    }
  }
  return best;
}

static std::vector<ObservedMoment> BuildObservedProtonMomentsFromQ2Bin(int q2bin) {
  const auto p = GetProtonSDMEs_TableQ2Bin(q2bin);

  std::vector<ObservedMoment> obs;
  obs.reserve(64);

  auto add = [&](int alpha, int L, int M, double val, double sig) {
    ObservedMoment m;
    m.alpha = alpha;
    m.L = L;
    m.M = M;
    m.value = val;
    m.sigma = sig;
    std::ostringstream os;
    os << "RH_" << alpha << "_" << L << "_" << M;
    m.name = os.str();
    obs.push_back(m);
  };

  auto add04 = [&](int L, int M, double val, double sig) {
    ObservedMoment m;
    m.alpha = 0;
    m.L = L;
    m.M = M;
    m.value = val;
    m.sigma = sig;
    m.isMixed04 = true;
    std::ostringstream os;
    os << "RH04_" << L << "_" << M;
    m.name = os.str();
    obs.push_back(m);
  };

  const double s_r00_04    = Comb(p.r00_04);
  const double s_re_r10_04 = Comb(p.re_r10_04);
  const double s_r1m1_04   = Comb(p.r1m1_04);
  const double s_r1m1_1    = Comb(p.r1m1_1);
  const double s_re_r10_1  = Comb(p.re_r10_1);
  const double s_r00_1     = Comb(p.r00_1);
  const double s_r11_1     = Comb(p.r11_1);
  const double s_im_r10_2  = Comb(p.im_r10_2);
  const double s_im_r1m1_2 = Comb(p.im_r1m1_2);
  const double s_im_r10_3  = Comb(p.im_r10_3);
  const double s_im_r1m1_3 = Comb(p.im_r1m1_3);
  const double s_r00_5     = Comb(p.r00_5);
  const double s_r11_5     = Comb(p.r11_5);
  const double s_re_r10_5  = Comb(p.re_r10_5);
  const double s_r1m1_5    = Comb(p.r1m1_5);
  const double s_im_r10_6  = Comb(p.im_r10_6);
  const double s_im_r1m1_6 = Comb(p.im_r1m1_6);
  const double s_im_r10_7  = Comb(p.im_r10_7);
  const double s_im_r1m1_7 = Comb(p.im_r1m1_7);
  const double s_r00_8     = Comb(p.r00_8);
  const double s_r11_8     = Comb(p.r11_8);
  const double s_re_r10_8  = Comb(p.re_r10_8);
  const double s_r1m1_8    = Comb(p.r1m1_8);

  const double f21 = TMath::Sqrt(12.0) / 5.0;
  const double f22 = TMath::Sqrt(6.0)  / 5.0;

  add04(2, 0, 0.2 * (3.0 * p.r00_04.v - 1.0),  TMath::Abs(0.6) * s_r00_04);
  add04(2, 1, f21 * p.re_r10_04.v,              TMath::Abs(f21) * s_re_r10_04);
  add04(2, 2, -f22 * p.r1m1_04.v,               TMath::Abs(f22) * s_r1m1_04);

  add(1, 0, 0, -(2.0 * p.r11_1.v + p.r00_1.v),
      TMath::Sqrt((2 * s_r11_1) * (2 * s_r11_1) + s_r00_1 * s_r00_1));
  add(1, 2, 0, 0.4 * (p.r11_1.v - p.r00_1.v),
      0.4 * TMath::Sqrt(s_r11_1 * s_r11_1 + s_r00_1 * s_r00_1));
  add(1, 2, 1, -f21 * p.re_r10_1.v,             TMath::Abs(f21) * s_re_r10_1);
  add(1, 2, 2,  f22 * p.r1m1_1.v,               TMath::Abs(f22) * s_r1m1_1);

  add(2, 2, 1, -f21 * p.im_r10_2.v,             TMath::Abs(f21) * s_im_r10_2);
  add(2, 2, 2,  f22 * p.im_r1m1_2.v,            TMath::Abs(f22) * s_im_r1m1_2);

  add(3, 2, 1, -f21 * p.im_r10_3.v,             TMath::Abs(f21) * s_im_r10_3);
  add(3, 2, 2,  f22 * p.im_r1m1_3.v,            TMath::Abs(f22) * s_im_r1m1_3);

  add(5, 0, 0, -(2.0 * p.r11_5.v + p.r00_5.v),
      TMath::Sqrt((2 * s_r11_5) * (2 * s_r11_5) + s_r00_5 * s_r00_5));
  add(5, 2, 0, 0.4 * (p.r11_5.v - p.r00_5.v),
      0.4 * TMath::Sqrt(s_r11_5 * s_r11_5 + s_r00_5 * s_r00_5));
  add(5, 2, 1, -f21 * p.re_r10_5.v,             TMath::Abs(f21) * s_re_r10_5);
  add(5, 2, 2,  f22 * p.r1m1_5.v,               TMath::Abs(f22) * s_r1m1_5);

  add(6, 2, 1, -f21 * p.im_r10_6.v,             TMath::Abs(f21) * s_im_r10_6);
  add(6, 2, 2,  f22 * p.im_r1m1_6.v,            TMath::Abs(f22) * s_im_r1m1_6);

  add(7, 2, 1, -f21 * p.im_r10_7.v,             TMath::Abs(f21) * s_im_r10_7);
  add(7, 2, 2,  f22 * p.im_r1m1_7.v,            TMath::Abs(f22) * s_im_r1m1_7);

  add(8, 0, 0, -(2.0 * p.r11_8.v + p.r00_8.v),
      TMath::Sqrt((2 * s_r11_8) * (2 * s_r11_8) + s_r00_8 * s_r00_8));
  add(8, 2, 0, 0.4 * (p.r11_8.v - p.r00_8.v),
      0.4 * TMath::Sqrt(s_r11_8 * s_r11_8 + s_r00_8 * s_r00_8));
  add(8, 2, 1, -f21 * p.re_r10_8.v,             TMath::Abs(f21) * s_re_r10_8);
  add(8, 2, 2,  f22 * p.r1m1_8.v,               TMath::Abs(f22) * s_r1m1_8);

  return obs;
}

static std::vector<ObservedMoment> BuildObservedProtonMomentsFromFile(const std::string& inFile,
                                                                      const std::string& treeName = "syntheticMoments") {
  constexpr double kSigma = 1e-3;

  std::unique_ptr<TFile> fin(TFile::Open(inFile.c_str(), "READ"));
  if (!fin || fin->IsZombie()) {
    throw std::runtime_error("BuildObservedProtonMomentsFromFile: failed to open file " + inFile);
  }

  TTree* t = dynamic_cast<TTree*>(fin->Get(treeName.c_str()));
  if (!t) {
    throw std::runtime_error("BuildObservedProtonMomentsFromFile: could not find tree '" + treeName + "'");
  }
  if (t->GetEntries() < 1) {
    throw std::runtime_error("BuildObservedProtonMomentsFromFile: tree '" + treeName + "' is empty");
  }

  auto readBranch = [&](const char* bname) -> double {
    if (!t->GetBranch(bname)) {
      throw std::runtime_error(std::string("BuildObservedProtonMomentsFromFile: missing branch '") + bname + "'");
    }
    double x = 0.0;
    t->SetBranchAddress(bname, &x);
    t->GetEntry(0);
    t->ResetBranchAddresses();
    return x;
  };

  std::vector<ObservedMoment> obs;
  obs.reserve(24);

  auto add = [&](int alpha, int L, int M, const char* branchName) {
    ObservedMoment m;
    m.alpha = alpha;
    m.L = L;
    m.M = M;
    m.value = readBranch(branchName);
    m.sigma = kSigma;
    m.name = branchName;
    obs.push_back(m);
  };

  auto add04 = [&](int L, int M, const char* branchName) {
    ObservedMoment m;
    m.alpha = 0;
    m.L = L;
    m.M = M;
    m.value = readBranch(branchName);
    m.sigma = kSigma;
    m.isMixed04 = true;
    m.name = branchName;
    obs.push_back(m);
  };

  add04(2, 0, "RH04_2_0");
  add04(2, 1, "RH04_2_1");
  add04(2, 2, "RH04_2_2");

  add(1, 0, 0, "RH_1_0_0");
  add(1, 2, 0, "RH_1_2_0");
  add(1, 2, 1, "RH_1_2_1");
  add(1, 2, 2, "RH_1_2_2");

  add(2, 2, 1, "RH_2_2_1");
  add(2, 2, 2, "RH_2_2_2");

  add(3, 2, 1, "RH_3_2_1");
  add(3, 2, 2, "RH_3_2_2");

  add(5, 0, 0, "RH_5_0_0");
  add(5, 2, 0, "RH_5_2_0");
  add(5, 2, 1, "RH_5_2_1");
  add(5, 2, 2, "RH_5_2_2");

  add(6, 2, 1, "RH_6_2_1");
  add(6, 2, 2, "RH_6_2_2");

  add(7, 2, 1, "RH_7_2_1");
  add(7, 2, 2, "RH_7_2_2");

  add(8, 0, 0, "RH_8_0_0");
  add(8, 2, 0, "RH_8_2_0");
  add(8, 2, 1, "RH_8_2_1");
  add(8, 2, 2, "RH_8_2_2");

  return obs;
}

static std::vector<ObservedMoment> BuildObservedMoments(const FitConfig& cfg) {
  if (cfg.observedInput.source == ObservedMomentsSource::kRootFile) {
    return BuildObservedProtonMomentsFromFile(cfg.observedInput.momentsFile, cfg.observedInput.momentsTree);
  }
  return BuildObservedProtonMomentsFromQ2Bin(cfg.observedInput.q2bin);
}

// -------------------------
// Term evaluation: matches BruTermCircle numerically
// -------------------------
static void BruTermCircleResolved(int reflsign, double factor, int l, int m, int lpr, int mpr,
                                  int alpha, bool negm, bool orientswap,
                                  double& outCoeff,
                                  char& refl1, char& orient1, int& l1, int& m1,
                                  char& refl2, char& orient2, int& l2, int& m2,
                                  bool& useTrig, TrigKind& trig) {
  // Implements the selection logic from E2S0AmpLoader::BruTermCircle, but instead of
  // returning a string, returns (overall coefficient, wave selectors, trig kind).

  if (!negm && (m < 0 || mpr < 0)) {
    outCoeff = 0.0;
    useTrig = false;
    trig = TrigKind::kNone;
    return;
  }

  refl1 = (reflsign == -1) ? 'b' : 'a';

  // primed reflectivity
  char reflpr = 'a';
  if (reflsign == -1 && alpha <= 4) reflpr = 'b';
  if (reflsign == 1 && alpha > 4) reflpr = 'b';

  // orientation selection
  orient1 = 'T';
  char orientpr = 'T';
  if (alpha == 4) {
    orient1 = orientpr = 'L';
  }
  if (alpha >= 5) {
    orient1 = 'L';
    orientpr = 'T';
  }
  if (alpha >= 5 && orientswap) {
    orient1 = 'T';
    orientpr = 'L';
  }

  // reflfactor
  int reflfactor = 1;
  if (alpha == 1 || alpha == 2) reflfactor = reflsign;

  // trig selection
  useTrig = !(l == lpr && m == mpr);
  trig = TrigKind::kCos;
  if (alpha == 3 || alpha == 7 || alpha == 8) trig = TrigKind::kSin;

  // extra minus for alpha==8
  if (alpha == 8) factor *= -1;

  outCoeff = reflfactor * factor;

  // wave selectors
  l1 = l; m1 = m;
  l2 = lpr; m2 = mpr;
  refl2 = reflpr;
  orient2 = orientpr;
}

// -------------------------
// Build numeric model for all desired H(alpha,L,M) as a list of terms
// -------------------------
static std::vector<MomentModel> BuildMomentModels(const FitConfig& cfg,
                                                  const std::map<std::string,int>& nameToIndex) {
  std::vector<MomentModel> models;

  // Build for alpha=0..8, L=0..2*lmax, M=0..L
  for (int alpha = 0; alpha <= 8; ++alpha) {
    for (int L = 0; L <= 2 * cfg.lmax; ++L) {
      for (int M = 0; M <= L; ++M) {
        MomentModel mm;
        mm.alpha = alpha;
        mm.L = L;
        mm.M = M;
        {
          // Internal model name corresponds to the *pure* (unscaled) moment H^alpha_{L M}.
          // Experimental inputs are RH^alpha_{L M} = R^alpha * H^alpha_{L M}.
          // We keep H_ names here so that:
          //   - R can be formed from H_4_0_0 and H_0_0_0
          //   - observed RH_* values can be mapped to the underlying H_* model index
          std::ostringstream os;
          os << "H_" << alpha << "_" << L << "_" << M;
          mm.name = os.str();
        }

        // Respect the loader behaviour: alpha 2,3,6,7 are zero for M=0
        if ((alpha == 2 || alpha == 3 || alpha == 6 || alpha == 7) && M == 0) {
          models.push_back(std::move(mm));
          continue;
        }

        for (int il = 0; il <= cfg.lmax; ++il) {
          if (cfg.onlyEven && (il % 2 == 1)) continue;
          for (int im = -il; im <= il; ++im) {
            if (std::abs(im) > cfg.mmax) continue;
            if (!cfg.negm && im < 0) continue;

            for (int ilpr = 0; ilpr <= cfg.lmax; ++ilpr) {
              if (cfg.onlyEven && (ilpr % 2 == 1)) continue;
              for (int impr = -ilpr; impr <= ilpr; ++impr) {
                if (std::abs(impr) > cfg.mmax) continue;
                if (!cfg.negm && impr < 0) continue;

                const double CM = ClebschGordan(ilpr, L, il, impr, M, im);
                const double C0 = ClebschGordan(ilpr, L, il, 0, 0, 0);
                const double factor = TMath::Sqrt((2.0 * ilpr + 1.0) / (2.0 * il + 1.0));
                double ccfactor = CM * C0 * factor;
                if (ccfactor == 0.0) continue;

                int mmprimesign = 1;
                if (std::abs((im - impr) % 2) == 1) mmprimesign = -1;

                int mprimesign = 1;
                if (std::abs(impr % 2) == 1) mprimesign = -1;

                int msign = 1;
                if (std::abs(im % 2) == 1) msign = -1;

                auto emitBru = [&](int reflsign, double f, int ll, int mm1, int llpr, int mm2,
                                   int alphaTerm, bool orSwap, bool asTrigSignFlip) {
                  // Resolve BruTermCircle selection
                  double baseCoeff = 0.0;
                  char r1='a', o1='T', r2='a', o2='T';
                  int l1=0,m1=0,l2=0,m2=0;
                  bool useTrig=true;
                  TrigKind trig=TrigKind::kNone;
                  BruTermCircleResolved(reflsign, f, ll, mm1, llpr, mm2, alphaTerm, cfg.negm, orSwap,
                                        baseCoeff, r1, o1, l1, m1, r2, o2, l2, m2, useTrig, trig);
                  if (baseCoeff == 0.0) return;

                  // Map wave selectors to parameter indices
                  const std::string mag1 = MagName(r1, o1, l1, m1);
                  const std::string phi1 = PhiName(r1, o1, l1, m1);
                  const std::string mag2 = MagName(r2, o2, l2, m2);
                  const std::string phi2 = PhiName(r2, o2, l2, m2);

                  auto itMag1 = nameToIndex.find(mag1);
                  auto itPhi1 = nameToIndex.find(phi1);
                  auto itMag2 = nameToIndex.find(mag2);
                  auto itPhi2 = nameToIndex.find(phi2);
                  if (itMag1 == nameToIndex.end() || itMag2 == nameToIndex.end()) {
                    throw std::runtime_error("Missing magnitude parameter in map: " + mag1 + " or " + mag2);
                  }
                  if (itPhi1 == nameToIndex.end() || itPhi2 == nameToIndex.end()) {
                    throw std::runtime_error("Missing phase parameter in map: " + phi1 + " or " + phi2);
                  }

                  Term t;
                  t.coeff = baseCoeff;
                  //if (asTrigSignFlip) t.coeff *= -1.0; // used to emulate the "+-" in alpha==3 strings

                  t.idxMag1 = itMag1->second;
                  t.idxPhi1 = itPhi1->second;
                  t.idxMag2 = itMag2->second;
                  t.idxPhi2 = itPhi2->second;
                  t.useTrig = useTrig;
                  t.trig = trig;
                  mm.terms.push_back(t);
                };

                // Replicate alpha-specific blocks from CGMatrixReflectivity
                if (alpha == 0) {
                  int refl = +1;
                  emitBru(refl, ccfactor, il, im, ilpr, impr, 0, false, false);
                  emitBru(refl, mmprimesign * ccfactor, il, -im, ilpr, -impr, 0, false, false);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, ccfactor, il, im, ilpr, impr, 0, false, false);
                    emitBru(refl, mmprimesign * ccfactor, il, -im, ilpr, -impr, 0, false, false);
                  }
                } else if (alpha == 1) {
                  int refl = +1;
                  emitBru(refl, msign*ccfactor, il, -im, ilpr, impr, 1, false, false);
                  emitBru(refl, mprimesign * ccfactor, il, im, ilpr, -impr, 1, false, false);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, msign*ccfactor, il, -im, ilpr, impr, 1, false, false);
                    emitBru(refl, mprimesign * ccfactor, il, im, ilpr, -impr, 1, false, false);
                  }
                } else if (alpha == 2) {
                  int refl = +1;
                  emitBru(refl, msign*ccfactor, il, -im, ilpr, impr, 2, false, false);
                  emitBru(refl, -1 * mprimesign * ccfactor, il, im, ilpr, -impr, 2, false, false);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, msign*ccfactor, il, -im, ilpr, impr, 2, false, false);
                    emitBru(refl, -1 * mprimesign * ccfactor, il, im, ilpr, -impr, 2, false, false);
                  }
                } else if (alpha == 3) {
                  int refl = +1;
                  double f = -1.0 * ccfactor; // ccfactor *= -1 in code
                  emitBru(refl, f, il, im, ilpr, impr, 3, false, false);
                  // emulate "+-" by flipping sign of the second term
                  emitBru(refl, -1*mmprimesign * f, il, -im, ilpr, -impr, 3, false, true);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, f, il, im, ilpr, impr, 3, false, false);
                    emitBru(refl, -1*mmprimesign * f, il, -im, ilpr, -impr, 3, false, true);
                  }
                } else if (alpha == 4) {
                  int refl = +1;
                  double f = -1.0 * ccfactor; // ccfactor*=-1
                  emitBru(refl, 2 * f, il, im, ilpr, impr, 4, false, false);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, 2 * f, il, im, ilpr, impr, 4, false, false);
                  }
                } else if (alpha == 5) {
                  int refl = +1;
                  double f = (-1.0 / TMath::Sqrt(2.0)) * ccfactor;
                  emitBru(refl, f, il, im, ilpr, impr, 5, false, false);
                  emitBru(refl, f, il, im, ilpr, impr, 5, true, false);
                  emitBru(refl, mmprimesign * f, il, -im, ilpr, -impr, 5, false, false);
                  emitBru(refl, mmprimesign * f, il, -im, ilpr, -impr, 5, true, false);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, f, il, im, ilpr, impr, 5, false, false);
                    emitBru(refl, f, il, im, ilpr, impr, 5, true, false);
                    emitBru(refl, mmprimesign * f, il, -im, ilpr, -impr, 5, false, false);
                    emitBru(refl, mmprimesign * f, il, -im, ilpr, -impr, 5, true, false);
                  }
                } else if (alpha == 6) {
                  int refl = +1;
                  double f = (-1.0 / TMath::Sqrt(2.0)) * ccfactor;
                  emitBru(refl, f, il, im, ilpr, impr, 6, false, false);
                  emitBru(refl, -1 * f, il, im, ilpr, impr, 6, true, false);
                  emitBru(refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 6, false, false);
                  emitBru(refl, mmprimesign * f, il, -im, ilpr, -impr, 6, true, false);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, f, il, im, ilpr, impr, 6, false, false);
                    emitBru(refl, -1 * f, il, im, ilpr, impr, 6, true, false);
                    emitBru(refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 6, false, false);
                    emitBru(refl, mmprimesign * f, il, -im, ilpr, -impr, 6, true, false);
                  }
                } else if (alpha == 7) {
                  int refl = +1;
                  double f = (-1.0 / TMath::Sqrt(2.0)) * ccfactor;
                  emitBru(refl, f, il, im, ilpr, impr, 7, false, false);
                  emitBru(refl, f, il, im, ilpr, impr, 7, true, false);
                  emitBru(refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 7, false, false);
                  emitBru(refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 7, true, false);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, f, il, im, ilpr, impr, 7, false, false);
                    emitBru(refl, f, il, im, ilpr, impr, 7, true, false);
                    emitBru(refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 7, false, false);
                    emitBru(refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 7, true, false);
                  }
                } else if (alpha == 8) {
                  int refl = +1;
                  double f = (-1.0 / TMath::Sqrt(2.0)) * ccfactor;
                  emitBru(refl, f, il, im, ilpr, impr, 8, false, false);
                  emitBru(refl, -1 * f, il, im, ilpr, impr, 8, true, false);
                  emitBru(refl, mmprimesign * f, il, -im, ilpr, -impr, 8, false, false);
                  emitBru(refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 8, true, false);
                  if (cfg.useNegRef) {
                    refl = -1;
                    emitBru(refl, f, il, im, ilpr, impr, 8, false, false);
                    emitBru(refl, -1 * f, il, im, ilpr, impr, 8, true, false);
                    emitBru(refl, mmprimesign * f, il, -im, ilpr, -impr, 8, false, false);
                    emitBru(refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 8, true, false);
                  }
                } else {
                }
              }
            }
          }
        }

        models.push_back(std::move(mm));
      }
    }
  }

  return models;
}

// -------------------------
// Evaluator / gradient for chi2
// -------------------------
struct EvalContext {
  FitConfig cfg;
  std::vector<ParDef> fullPars;
  std::vector<int> freeToFull;
  std::vector<int> fullToFree;

  std::vector<ObservedMoment> observed;
  std::vector<MomentModel> modelsRec;

  std::map<std::string, size_t> modelIndexByName;
  std::vector<int> observedModelIdx;
  std::vector<int> observedModelIdx0;
  std::vector<int> observedModelIdx4;

  std::vector<int> magFullIdx;

  int idxH0_00 = -1;
  int idxH4_00 = -1;

  mutable unsigned callCount = 0;
  mutable TTree* iterTree = nullptr;
  mutable double iter_log_val = 0.0;
  mutable std::vector<double> iter_parVals;
  mutable std::vector<double> iter_Hrec;
};

static void FillFullFromFree(const EvalContext& ctx, const double* x, std::vector<double>& v) {
  if (v.size() != ctx.fullPars.size()) v.assign(ctx.fullPars.size(), 0.0);
  for (size_t i = 0; i < ctx.fullPars.size(); ++i) v[i] = ctx.fullPars[i].init;
  for (size_t i = 0; i < ctx.freeToFull.size(); ++i) {
    const int fullIdx = ctx.freeToFull[i];
    v[(size_t)fullIdx] = x[i];
  }
}


static inline void EnsureSize(std::vector<double>& v, size_t n) {
  if (v.size() != n) v.assign(n, 0.0);
}


// Full-parameter derivative version (derivatives w.r.t ALL full parameters, including fixed).
static double EvalMomentAndDerivFull(const MomentModel& mm,
                                     const std::vector<double>& fullVals,
                                     std::vector<double>& dHdFull) {
  std::fill(dHdFull.begin(), dHdFull.end(), 0.0);
  double H = 0.0;

  for (const auto& t : mm.terms) {
    const double m1 = fullVals[t.idxMag1];
    const double p1 = fullVals[t.idxPhi1];
    const double m2 = fullVals[t.idxMag2];
    const double p2 = fullVals[t.idxPhi2];

    // if (!t.useTrig) {
    //   const double val = t.coeff * (m1 * m2);
    //   H += val;
    //   dHdFull[(size_t)t.idxMag1] += t.coeff * m2;
    //   dHdFull[(size_t)t.idxMag2] += t.coeff * m1;
    //   continue;
    // }

    const double dphi = p1 - p2;
    double s = 0.0, c = 1.0;
    FastSinCos(dphi, s, c);

    double trigv = 0.0;
    double dtrig_dphi = 0.0;
    if (t.trig == TrigKind::kCos) {
      trigv = c;
      dtrig_dphi = -s;
    } else if (t.trig == TrigKind::kSin) {
      trigv = s;
      dtrig_dphi = c;
    }

    const double val = t.coeff * m1 * m2 * trigv;
    H += val;

    dHdFull[(size_t)t.idxMag1] += t.coeff * (m2 * trigv);
    dHdFull[(size_t)t.idxMag2] += t.coeff * (m1 * trigv);

    const double common = t.coeff * m1 * m2 * dtrig_dphi;
    dHdFull[(size_t)t.idxPhi1] += common;
    dHdFull[(size_t)t.idxPhi2] -= common;
  }

  return H;
}


static double EvalMomentOnly(const MomentModel& mm,
                             const std::vector<double>& fullVals) {
  const double* fv = fullVals.data();
  double H = 0.0;

  for (const auto& t : mm.terms) {
    const double m1 = fv[t.idxMag1];
    const double p1 = fv[t.idxPhi1];
    const double m2 = fv[t.idxMag2];
    const double p2 = fv[t.idxPhi2];

    // if (!t.useTrig) {
    //   H += t.coeff * (m1 * m2);
    //   continue;
    // }

    const double dphi = p1 - p2;
    double s = 0.0, c = 1.0;
    FastSinCos(dphi, s, c);

    if (t.trig == TrigKind::kCos) {
      H += t.coeff * m1 * m2 * c;
    } else if (t.trig == TrigKind::kSin) {
      H += t.coeff * m1 * m2 * s;
    }
  }

  return H;
}

static void EvalAllMoments(const std::vector<MomentModel>& models,
                           const std::vector<double>& fullVals,
                           std::vector<double>& values) {
  EnsureSize(values, models.size());
  for (size_t i = 0; i < models.size(); ++i) {
    values[i] = EvalMomentOnly(models[i], fullVals);
  }
}

class Chi2Function final : public ROOT::Math::IMultiGradFunction {
public:
  explicit Chi2Function(std::shared_ptr<EvalContext> ctx) : ctx_(std::move(ctx)) {
    if (!ctx_) throw std::runtime_error("Chi2Function: null context");
  }

  unsigned int NDim() const override { return static_cast<unsigned>(ctx_->freeToFull.size()); }

  ROOT::Math::IBaseFunctionMultiDim* Clone() const override {
    return new Chi2Function(ctx_);
  }

  
  double DoEval(const double* x) const override {
  FillFullFromFree(*ctx_, x, fullVals_);
  const auto& fullVals = fullVals_;

  ctx_->callCount++;

  const bool writeTrace = (ctx_->iterTree && ctx_->cfg.recordEvery > 0 &&
                           (ctx_->callCount % ctx_->cfg.recordEvery == 0));

  EnsureSize(buf_momRec_, ctx_->modelsRec.size());
  if (buf_momSeen_.size() != ctx_->modelsRec.size()) buf_momSeen_.assign(ctx_->modelsRec.size(), 0);
  std::fill(buf_momSeen_.begin(), buf_momSeen_.end(), static_cast<unsigned char>(0));

  auto GetMoment = [&](int idx) -> double {
    if (idx < 0) return 0.0;
    const size_t uidx = static_cast<size_t>(idx);
    if (!buf_momSeen_[uidx]) {
      buf_momRec_[uidx] = EvalMomentOnly(ctx_->modelsRec[uidx], fullVals);
      buf_momSeen_[uidx] = 1;
    }
    return buf_momRec_[uidx];
  };

  const double H0 = (ctx_->idxH0_00 >= 0) ? GetMoment(ctx_->idxH0_00) : 0.0;
  const double H4 = (ctx_->idxH4_00 >= 0) ? GetMoment(ctx_->idxH4_00) : 0.0;

  double R = -H4/H0;

  const double eps = ctx_->cfg.epsR4;
  const double invDen = 1.0 / (1.0 + eps * R);
  const double sqrtR = std::sqrt(R);

  auto ScaleForAlpha = [&](int alpha) -> double {
    if (alpha >= 0 && alpha <= 3) return invDen;
    if (alpha == 4) return eps * R * invDen;
    if (alpha >= 5 && alpha <= 8) return sqrtR * invDen;
    return 1.0;
  };

  double chi2 = 0.0;

  const double rNorm = (H0 - H4 - ctx_->cfg.H0H4NormTarget) / ctx_->cfg.H0H4NormSigma;
  chi2 += rNorm * rNorm;

  for (size_t i = 0; i < ctx_->observed.size(); ++i) {
    const auto& ob = ctx_->observed[i];
    double RH = 0.0;

    if (ob.isMixed04) {
      const int idx0 = ctx_->observedModelIdx0[i];
      const int idx4 = ctx_->observedModelIdx4[i];
      if (idx0 < 0 || idx4 < 0) continue;
      RH = GetMoment(idx0) + ctx_->cfg.rh04MixCoeff * (eps * R) * GetMoment(idx4);
    } else {
      const int midx = ctx_->observedModelIdx[i];
      if (midx < 0) continue;
      RH = ScaleForAlpha(ob.alpha) * GetMoment(midx);
    }

    const double r = (RH - ob.value) / ob.sigma;
    chi2 += r * r;
  }

  if (writeTrace) {
    ctx_->iter_log_val = (chi2 > 0.0) ? std::log10(chi2) : -999.0;

    for (size_t j = 0; j < ctx_->fullPars.size(); ++j) {
      double v = fullVals[j];
      ctx_->iter_parVals[j] = v;
    }

    for (size_t i = 0; i < ctx_->modelsRec.size(); ++i) {
      if (!buf_momSeen_[i]) {
        buf_momRec_[i] = EvalMomentOnly(ctx_->modelsRec[i], fullVals);
        buf_momSeen_[i] = 1;
      }
      ctx_->iter_Hrec[i] = buf_momRec_[i];
    }
    ctx_->iterTree->Fill();
  }

  return chi2 / NDim();
}

void Gradient(const double* x, double* grad) const override {
  std::fill(grad, grad + NDim(), 0.0);

  FillFullFromFree(*ctx_, x, fullVals_);
  const auto& fullVals = fullVals_;
  const size_t nFull = ctx_->fullPars.size();

  auto mapFullDerivToFree = [&](const std::vector<double>& dFull, unsigned iFree) -> double {
    const int fullIdx = ctx_->freeToFull[iFree];
    return dFull[(size_t)fullIdx];
  };

  EnsureSize(buf_fullA_, nFull);
  EnsureSize(buf_fullB_, nFull);
  EnsureSize(buf_fullC_, nFull);
  EnsureSize(buf_fullD_, nFull);
  EnsureSize(buf_fullE_, nFull);
  EnsureSize(buf_fullF_, nFull);

  auto& dH0_full = buf_fullA_;
  auto& dH4_full = buf_fullB_;
  auto& dRraw_full = buf_fullC_;
  auto& dR_full = buf_fullD_;
  auto& dH_full = buf_fullE_;
  auto& dH2_full = buf_fullF_;

  double H0 = 0.0;
  double H4 = 0.0;
  if (ctx_->idxH0_00 >= 0) H0 = EvalMomentAndDerivFull(ctx_->modelsRec[(size_t)ctx_->idxH0_00], fullVals, dH0_full);
  if (ctx_->idxH4_00 >= 0) H4 = EvalMomentAndDerivFull(ctx_->modelsRec[(size_t)ctx_->idxH4_00], fullVals, dH4_full);

  constexpr double kEpsH0 = 1e-12;
  constexpr double kEpsR = 1e-12;

  const double H0abs = std::abs(H0);
  const bool denomClamped = !(H0abs > kEpsH0);
  const double H0safe = denomClamped ? (H0 >= 0.0 ? kEpsH0 : -kEpsH0) : H0;

  const double Rraw = -H4 / H0safe;
  const double signRraw = (Rraw >= 0.0 ? 1.0 : -1.0);
  double R = std::abs(Rraw);
  bool Rclamped = false;
  if (R < kEpsR) {
    R = kEpsR;
    Rclamped = true;
  }
  const double sqrtR = std::sqrt(R);

  const double eps = ctx_->cfg.epsR4;
  const double den = 1.0 + eps * R;
  const double invDen = 1.0 / den;
  const double invDen2 = invDen * invDen;

  std::fill(dRraw_full.begin(), dRraw_full.end(), 0.0);
  if (!denomClamped) {
    const double H0sq = H0safe * H0safe;
    for (size_t i = 0; i < nFull; ++i) {
      dRraw_full[i] = -((dH4_full[i] * H0safe) - (H4 * dH0_full[i])) / H0sq;
    }
  } else {
    for (size_t i = 0; i < nFull; ++i) {
      dRraw_full[i] = -(dH4_full[i]) / H0safe;
    }
  }

  std::fill(dR_full.begin(), dR_full.end(), 0.0);
  if (!Rclamped) {
    for (size_t i = 0; i < nFull; ++i) dR_full[i] = signRraw * dRraw_full[i];
  }

  EnsureSize(buf_dR_, NDim());
  auto& dR_free = buf_dR_;
  for (unsigned i = 0; i < NDim(); ++i) dR_free[i] = mapFullDerivToFree(dR_full, i);

  auto ScaleForAlpha = [&](int alpha) -> double {
    if (alpha >= 0 && alpha <= 3) return invDen;
    if (alpha == 4) return eps * R * invDen;
    if (alpha >= 5 && alpha <= 8) return sqrtR * invDen;
    return 1.0;
  };

  auto dScaleForAlpha_dR = [&](int alpha) -> double {
    if (alpha >= 0 && alpha <= 3) return -eps * invDen2;
    if (alpha == 4) return eps * invDen2;
    if (alpha >= 5 && alpha <= 8) return 0.5 * invDen / sqrtR - eps * sqrtR * invDen2;
    return 0.0;
  };

  const double normScale = 2.0 * (H0 - H4 - ctx_->cfg.H0H4NormTarget) /
                           (ctx_->cfg.H0H4NormSigma * ctx_->cfg.H0H4NormSigma);
  for (unsigned i = 0; i < NDim(); ++i) {
    grad[i] += normScale * (mapFullDerivToFree(dH0_full, i) - mapFullDerivToFree(dH4_full, i));
  }

  for (size_t iobs = 0; iobs < ctx_->observed.size(); ++iobs) {
    const auto& ob = ctx_->observed[iobs];

    if (ob.isMixed04) {
      const int idx0 = ctx_->observedModelIdx0[iobs];
      const int idx4 = ctx_->observedModelIdx4[iobs];
      if (idx0 < 0 || idx4 < 0) continue;

      const double H0_LM = EvalMomentAndDerivFull(ctx_->modelsRec[(size_t)idx0], fullVals, dH_full);
      const double H4_LM = EvalMomentAndDerivFull(ctx_->modelsRec[(size_t)idx4], fullVals, dH2_full);
      const double k = ctx_->cfg.rh04MixCoeff * eps;
      const double RH = H0_LM + k * R * H4_LM;
      const double r = (RH - ob.value) / ob.sigma;
      const double scale = 2.0 * r / ob.sigma;

      for (unsigned i = 0; i < NDim(); ++i) {
        const double dH0_i = mapFullDerivToFree(dH_full, i);
        const double dH4_i = mapFullDerivToFree(dH2_full, i);
        grad[i] += scale * (dH0_i + k * (H4_LM * dR_free[i] + R * dH4_i));
      }
      continue;
    }

    const int midx = ctx_->observedModelIdx[iobs];
    if (midx < 0) continue;

    const double H = EvalMomentAndDerivFull(ctx_->modelsRec[(size_t)midx], fullVals, dH_full);
    const double scaleAlpha = ScaleForAlpha(ob.alpha);
    const double dScale_dR = dScaleForAlpha_dR(ob.alpha);
    const double RH = scaleAlpha * H;
    const double r = (RH - ob.value) / ob.sigma;
    const double scale = 2.0 * r / ob.sigma;

    for (unsigned i = 0; i < NDim(); ++i) {
      const double dH_i = mapFullDerivToFree(dH_full, i);
      grad[i] += scale * (scaleAlpha * dH_i + H * dScale_dR * dR_free[i]);
    }
  }
}

double DoDerivative(const double* x, unsigned int icoord) const override {
  if (icoord >= NDim()) return 0.0;
  EnsureSize(buf_gradTmp_, NDim());
  std::fill(buf_gradTmp_.begin(), buf_gradTmp_.end(), 0.0);
  Gradient(x, buf_gradTmp_.data());
  return buf_gradTmp_[icoord];
}

private:
  std::shared_ptr<EvalContext> ctx_;
  mutable std::vector<double> fullVals_;
  mutable std::vector<double> buf_momRec_;
  mutable std::vector<unsigned char> buf_momSeen_;
  mutable std::vector<double> buf_dR_;
  mutable std::vector<double> buf_fullA_;
  mutable std::vector<double> buf_fullB_;
  mutable std::vector<double> buf_fullC_;
  mutable std::vector<double> buf_fullD_;
  mutable std::vector<double> buf_fullE_;
  mutable std::vector<double> buf_fullF_;
  mutable std::vector<double> buf_gradTmp_;
};

// -------------------------
// Utility: build context (param maps, models)
// -------------------------
static std::shared_ptr<EvalContext> BuildContext(const FitConfig& cfg) {
  auto ctx = std::make_shared<EvalContext>();
  ctx->cfg = cfg;
  ctx->fullPars = BuildAmplitudePhaseParameters(cfg);

  std::map<std::string, int> nameToIndex;
  for (int i = 0; i < (int)ctx->fullPars.size(); ++i) {
    nameToIndex[ctx->fullPars[i].name] = i;
    if (!ctx->fullPars[i].isPhase) ctx->magFullIdx.push_back(i);
  }

  ctx->fullToFree.assign(ctx->fullPars.size(), -1);
  for (int i = 0; i < (int)ctx->fullPars.size(); ++i) {
    if (ctx->fullPars[i].fixed) continue;
    ctx->fullToFree[i] = (int)ctx->freeToFull.size();
    ctx->freeToFull.push_back(i);
  }

  ctx->observed = BuildObservedMoments(cfg);
  ctx->modelsRec = BuildMomentModels(cfg, nameToIndex);

  for (size_t i = 0; i < ctx->modelsRec.size(); ++i) {
    ctx->modelIndexByName[ctx->modelsRec[i].name] = i;
  }

  auto it0 = ctx->modelIndexByName.find("H_0_0_0");
  if (it0 != ctx->modelIndexByName.end()) ctx->idxH0_00 = static_cast<int>(it0->second);
  auto it4 = ctx->modelIndexByName.find("H_4_0_0");
  if (it4 != ctx->modelIndexByName.end()) ctx->idxH4_00 = static_cast<int>(it4->second);

  ctx->observedModelIdx.assign(ctx->observed.size(), -1);
  ctx->observedModelIdx0.assign(ctx->observed.size(), -1);
  ctx->observedModelIdx4.assign(ctx->observed.size(), -1);

  for (size_t i = 0; i < ctx->observed.size(); ++i) {
    const auto& ob = ctx->observed[i];

    if (ob.isMixed04) {
      std::ostringstream os0, os4;
      os0 << "H_0_" << ob.L << "_" << ob.M;
      os4 << "H_4_" << ob.L << "_" << ob.M;

      auto itObs0 = ctx->modelIndexByName.find(os0.str());
      auto itObs4 = ctx->modelIndexByName.find(os4.str());
      if (itObs0 != ctx->modelIndexByName.end()) ctx->observedModelIdx0[i] = static_cast<int>(itObs0->second);
      if (itObs4 != ctx->modelIndexByName.end()) ctx->observedModelIdx4[i] = static_cast<int>(itObs4->second);
      continue;
    }

    std::string key = ob.name;
    if (key.rfind("RH_", 0) == 0) key = "H_" + key.substr(3);
    auto it = ctx->modelIndexByName.find(key);
    if (it != ctx->modelIndexByName.end()) ctx->observedModelIdx[i] = static_cast<int>(it->second);
  }

  return ctx;
}

// -------------------------
// Main runner
// -------------------------
static void MakeBranchesForPars(TTree* t, const std::vector<ParDef>& pars, std::vector<double>& storage) {
  storage.assign(pars.size(), 0.0);
  for (size_t i = 0; i < pars.size(); ++i) {
    t->Branch(pars[i].name.c_str(), &storage[i]);
  }
}

static void MakeBranchesForMoments(TTree* t,
                                  const std::vector<MomentModel>& models,
                                  std::vector<double>& storage) {
  storage.assign(models.size(), 0.0);
  for (size_t i = 0; i < models.size(); ++i) {
    t->Branch(models[i].name.c_str(), &storage[i]);
  }
}
} // namespace chi2_amp_fit

// =====================================================================================
// User-facing macro entry points
// =====================================================================================
namespace {

void RunGivenMoments_Chi2Amps_Impl(const chi2_amp_fit::FitConfig& cfg, const char* outFile) {
  using namespace chi2_amp_fit;

  auto ctx = BuildContext(cfg);

  TRandom3 rng(cfg.randomSeed);
  if (cfg.randomSeed == 0) rng.SetSeed(0);

  std::unique_ptr<TFile> fout(TFile::Open(outFile, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error(std::string("Failed to open output file: ") + outFile);
  }

  TTree* t = new TTree("fitResults", "Chi2 fit results (per start)");
  double log_val = 0.0;
  t->Branch("log_val", &log_val);

  std::vector<double> parVals;
  MakeBranchesForPars(t, ctx->fullPars, parVals);
  std::vector<double> momRecVals;
  MakeBranchesForMoments(t, ctx->modelsRec, momRecVals);

  Chi2Function fcn(ctx);

  TBenchmark bench;
  bench.Start("fit");

  std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

  min->SetMaxFunctionCalls(cfg.maxCalls);
  min->SetMaxIterations(cfg.maxIters);
  min->SetTolerance(cfg.tolerance);
  min->SetStrategy(cfg.strategy);
  min->SetPrintLevel(cfg.printLevel);
  min->SetFunction(fcn);

  for (unsigned iStart = 0; iStart < cfg.nStarts; ++iStart) {
    ctx->callCount = 0;
    if (iStart % 100 == 0) std::cout << iStart << std::endl;

    for (unsigned i = 0; i < fcn.NDim(); ++i) {
      const int fullIdx = ctx->freeToFull[i];
      const auto& p = ctx->fullPars[fullIdx];

      double start = p.isPhase ? rng.Gaus(0.0, 0.1) : rng.Gaus(0.5, 0.1);
      if (start < p.low) start = p.low;
      if (start > p.high) start = p.high;

      min->SetLimitedVariable(i, p.name.c_str(), start, (p.step > 0 ? p.step : 1e-3), p.low, p.high);
    }

    min->Minimize();
    if (cfg.runHesse) min->Hesse();

    const double chi2 = min->MinValue();
    log_val = (chi2 > 0.0) ? std::log10(chi2) : -999.0;

    FillFullFromFree(*ctx, min->X(), parVals);

    EvalAllMoments(ctx->modelsRec, parVals, momRecVals);
    t->Fill();
  }

  bench.Stop("fit");
  bench.Print("fit");

  fout->Write();
  fout->Close();

  std::cout << "Saved: " << outFile << std::endl;
}

} // namespace

// Standard constructor with q2 value that then picks from HERMES data
void RunGivenMoments_Chi2Amps_Setup(double q2Value,
                              unsigned nStarts = 10000,
                              const char* outFile = "resultsGivenMoments_chi2_amps.root",
                              uint32_t seed = 0,
                              double epsR4 = 1.0) {
  chi2_amp_fit::FitConfig cfg;
  cfg.nStarts = nStarts;
  cfg.randomSeed = seed;
  cfg.epsR4 = epsR4;
  cfg.observedInput.source = chi2_amp_fit::ObservedMomentsSource::kQ2Table;
  cfg.observedInput.q2bin = chi2_amp_fit::GetClosestQ2Bin(q2Value);
  RunGivenMoments_Chi2Amps_Impl(cfg, outFile);
}

// Constructor when you know the Q2 bin (maybe a better implementation) so for example
// Eventually move to file + q2 value and take the SDMEs out of this
// i.e. move the HERMEs etc data into root files then just load a root file + Q2
// update generator to have some Q2 value then...
void RunGivenMoments_Chi2Amps_Setup(int q2bin,
                                    unsigned nStarts = 10000,
                                    const char* outFile = "resultsGivenMoments_chi2_amps.root",
                                    uint32_t seed = 0,
                                    double epsR4 = 1.0)
{
  chi2_amp_fit::FitConfig cfg;
  cfg.nStarts = nStarts;
  cfg.randomSeed = seed;
  cfg.epsR4 = epsR4;
  cfg.observedInput.source = chi2_amp_fit::ObservedMomentsSource::kQ2Table;
  cfg.observedInput.q2bin = q2bin;
  RunGivenMoments_Chi2Amps_Impl(cfg, outFile);
}

// Constructor with moments file (i.e. when using generated fixed amplitudes)
void RunGivenMoments_Chi2Amps_Setup(const char* momentsFile,
                              unsigned nStarts = 10000,
                              const char* outFile = "resultsGivenMoments_chi2_amps.root",
                              uint32_t seed = 0,
                              double epsR4 = 1.0,
                              const char* treeName = "syntheticMoments") {
  chi2_amp_fit::FitConfig cfg;
  cfg.nStarts = nStarts;
  cfg.randomSeed = seed;
  cfg.epsR4 = epsR4;
  cfg.observedInput.source = chi2_amp_fit::ObservedMomentsSource::kRootFile;
  cfg.observedInput.momentsFile = momentsFile ? momentsFile : "";
  cfg.observedInput.momentsTree = treeName ? treeName : "syntheticMoments";
  RunGivenMoments_Chi2Amps_Impl(cfg, outFile);
}

// Below handles all of the multiprocessing stuff
namespace {
static std::string MakePartFileName(const char* outFile, unsigned workerId) {
  std::string base = outFile ? std::string(outFile) : std::string("resultsGivenMoments_chi2_amps.root");
  const std::string ext = ".root";
  if (base.size() >= ext.size() && base.compare(base.size() - ext.size(), ext.size(), ext) == 0) {
    base.erase(base.size() - ext.size());
  }
  std::ostringstream os;
  os << base << ".part_" << workerId << ".root";
  return os.str();
}
}

void RunGivenMoments_Chi2Amps(double q2Value,
                              unsigned nStarts = 10000,
                              const char* outFile = "resultsGivenMoments_MP_chi2_amps.root",
                              uint32_t seed = 0,
                              double epsR4 = 1.0,
                              unsigned int nCores = 1) {
  if (nStarts == 0) {
    throw std::runtime_error("RunGivenMoments_Chi2Amps_MP: nStarts must be > 0");
  }

  unsigned nWorkers = nCores;
  if (nWorkers == 0) {
    nWorkers = std::max(1u, std::thread::hardware_concurrency());
  }
  if (nWorkers > nStarts) nWorkers = nStarts;

  std::vector<unsigned> workerIds(nWorkers);
  for (unsigned i = 0; i < nWorkers; ++i) workerIds[i] = i;

  ROOT::TProcessExecutor pool(nWorkers);

  auto partFiles = pool.Map([=](unsigned workerId) {
    const unsigned baseStarts = nStarts / nWorkers;
    const unsigned extra = nStarts % nWorkers;
    const unsigned myStarts = baseStarts + (workerId < extra ? 1u : 0u);

    const auto partFile = MakePartFileName(outFile, workerId);
    const uint32_t workerSeed =
        (seed == 0) ? (0x9e3779b9u + 100003u * workerId)
                    : (seed + 100003u * workerId);

    RunGivenMoments_Chi2Amps_Setup(q2Value,
                                        myStarts,
                                  partFile.c_str(),
                                        workerSeed,
                                        epsR4);
    return partFile;
  }, workerIds);

  TFileMerger merger(/*isLocal=*/kTRUE, /*histoOneGo=*/kFALSE);
  merger.OutputFile(outFile, "RECREATE");
  for (const auto& partFile : partFiles) {
    merger.AddFile(partFile.c_str());
  }

  const Bool_t ok = merger.Merge();
  if (!ok) {
    throw std::runtime_error(std::string("RunGivenMoments_Chi2Amps_MP: merge failed for output file: ") + outFile);
  }


  for (const auto& partFile : partFiles) {
    gSystem->Unlink(partFile.c_str());
  }

  std::cout << "Merged " << partFiles.size() << " worker files into " << outFile << std::endl;
}