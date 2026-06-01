#include "TBenchmark.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"

#include "Math/Factory.h"
#include "Math/IFunction.h"
#include "Math/Minimizer.h"
#include <Math/SpecFuncMathMore.h>

#include "ROOT/TProcessExecutor.hxx"
#include "TFileMerger.h"
#include "TSystem.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace chi2_amp_fit_opt {

constexpr double kPi = TMath::Pi();
constexpr double kTwoPi = 2.0 * TMath::Pi();
constexpr double kSqrt2 = 1.4142135623730950488016887242097;
constexpr double kInvSqrt2 = 1.0 / kSqrt2;
constexpr double kSqrt6Over5 = 0.48989794855663561963945681494118;
constexpr double kSqrt12Over5 = 0.69282032302755092063339055356909;
constexpr int kNumBins = 18;

inline void FastSinCos(double x, double& s, double& c) {
#if defined(__GLIBC__) || defined(__APPLE__)
  ::sincos(x, &s, &c);
#else
  s = std::sin(x);
  c = std::cos(x);
#endif
}

inline double ClebschGordan(int l1, int l2, int l3, int m1, int m2, int m3) {
  using ROOT::Math::wigner_3j;
  const double sign = ((l1 - l2 + m3) & 1) ? -1.0 : 1.0;
  return sign * std::sqrt(2.0 * l3 + 1.0) *
         wigner_3j(2 * l1, 2 * l2, 2 * l3, 2 * m1, 2 * m2, -2 * m3);
}

struct ObservedMoment {
  int alpha = 0;
  int L = 0;
  int M = 0;
  double value = 0.0;
  double sigma = 1.0;
  bool isMixed04 = false;
  std::string name;
};

struct ParDef {
  std::string name;
  double init = 0.0;
  double step = 1e-3;
  double low = -1.0;
  double high = 1.0;
  bool fixed = false;
  bool isPhase = false;
};

struct FitConfig {
  int lmax = 2;
  int mmax = 2;
  bool useNegRef = true;
  bool onlyEven = false;
  bool negm = true;
  bool enforceLongitudinalParity = false;

  double epsilon = 1.0;

  unsigned nStarts = 10000;
  unsigned maxCalls = 50000;
  unsigned maxIters = 50000;
  double tolerance = 1e-6;
  int strategy = 2;
  int printLevel = 0;
  bool runHesse = true;
  uint32_t randomSeed = 0;


  bool verbose = true;
  bool useNumericalGradient = true;

  bool useMCMCPreScan = false;
  unsigned mcmcSteps = 2000;
  double mcmcTemperature = 1.0;
  double mcmcProposalMagSigma = 0.03;
  double mcmcProposalPhaseSigma = 0.10;

  double magnitudeStartMean = 0.5;
  double magnitudeStartSigma = 0.5;
  double magnitudeMin = 0.0;
  double magnitudeMax = 5.0;
  bool photoProduction = false;

  std::string momentsFile;
  std::string momentsTree;
  int bin;
};

enum class TrigKind : uint8_t { kCos = 0, kSin = 1 };

struct PhasePair {
  int idxPhi1 = -1;
  int idxPhi2 = -1;
};

struct Term {
  double coeff = 0.0;
  int idxMag1 = -1;
  int idxMag2 = -1;
  int phasePairIdx = -1;
  TrigKind trig = TrigKind::kCos;
  bool ignorePhase = false;
};

struct MomentModel {
  int alpha = 0;
  int L = 0;
  int M = 0;
  std::string name;
  std::vector<Term> terms;
};

static double ReadArrayBranchElement(TTree* t, const char* branchName, int bin, bool isErr) {
  TBranch* branch = t->GetBranch(branchName);
  if (!branch) return isErr ? 0.001 : 0.0;

  TLeaf* leaf = branch->GetLeaf(branchName);
  if (!leaf && branch->GetListOfLeaves() && branch->GetListOfLeaves()->GetEntries() > 0) {
    leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
  }
  if (!leaf) return isErr ? 0.001 : 0.0;

  t->GetEntry(0);
  const int nData = leaf->GetNdata();
  if (nData <= 1) return leaf->GetValue(0);
  if (bin < 0 || bin >= nData) throw std::runtime_error(std::string("Bin out of range for branch ") + branchName);
  return leaf->GetValue(bin);
}

// Input file handling below here

static std::vector<ObservedMoment> BuildObservedMoments(const std::string& inFile,
                                                                           const std::string& treeName,
                                                                           int bin,
                                                                           const FitConfig& cfg) {
  std::unique_ptr<TFile> fin(TFile::Open(inFile.c_str(), "READ"));
  if (!fin || fin->IsZombie()) throw std::runtime_error("Failed to open file " + inFile);
  TTree* t = dynamic_cast<TTree*>(fin->Get(treeName.c_str()));
  if (!t) throw std::runtime_error("Could not find tree '" + treeName + "'");
  if (t->GetEntries() < 1) throw std::runtime_error("Tree '" + treeName + "' is empty");
  if (bin < 0 || bin >= kNumBins) throw std::runtime_error("Bin out of range");

  std::vector<ObservedMoment> obs;
  obs.reserve(40);
  if (cfg.photoProduction) obs.reserve(18); // Photoproduction has max 12 elements
  auto add = [&](int alpha, int L, int M, const char* valName, const char* errName, bool mixed04) {
    ObservedMoment m;
    m.alpha = alpha; m.L = L; m.M = M; m.value = ReadArrayBranchElement(t, valName, bin, false); m.sigma = ReadArrayBranchElement(t, errName, bin, true); m.isMixed04 = mixed04; m.name = valName;
    obs.push_back(std::move(m));
  };

  if (cfg.photoProduction) {
    add(0,0,0,"RH_0_0_0","RH_0_0_0_err",false);
    add(0,1,0,"RH_0_1_0","RH_0_1_0_err",false);
    add(0,1,1,"RH_0_1_1","RH_0_1_1_err",false);
    add(0,2,0,"RH_0_2_0","RH_0_2_0_err",false);
    add(0,2,1,"RH_0_2_1","RH_0_2_1_err",false);
    add(0,2,2,"RH_0_2_2","RH_0_2_2_err",false);

    add(1,0,0,"RH_1_0_0","RH_1_0_0_err",false);
    add(1,1,0,"RH_1_1_0","RH_1_1_0_err",false);
    add(1,1,1,"RH_1_1_1","RH_1_1_1_err",false);
    add(1,2,0,"RH_1_2_0","RH_1_2_0_err",false);
    add(1,2,1,"RH_1_2_1","RH_1_2_1_err",false);
    add(1,2,2,"RH_1_2_2","RH_1_2_2_err",false);

    add(2,1,1,"RH_2_1_1","RH_2_1_1_err",false);
    add(2,2,1,"RH_2_2_1","RH_2_2_1_err",false);
    add(2,2,2,"RH_2_2_2","RH_2_2_2_err",false);

    // Need to add options for circular vs linear pol gluex in linear
    // Cant have these values zeroed out as they effect the plot,  not 0 in practice
    //add(3,1,1,"RH_3_1_1","RH_3_1_1_err",false);
    //add(3,2,1,"RH_3_2_1","RH_3_2_1_err",false);
    //add(3,2,2,"RH_3_2_2","RH_3_2_2_err",false);
    return obs;
  }

  // Electroproduction / full moment set. RH04_0_0 is read directly from the input file.
  add(0,0,0,"RH04_0_0","RH04_0_0_err",true);
  add(0,1,0,"RH04_1_0","RH04_1_0_err",true);
  add(0,1,1,"RH04_1_1","RH04_1_1_err",true);
  add(0,2,0,"RH04_2_0","RH04_2_0_err",true);
  add(0,2,1,"RH04_2_1","RH04_2_1_err",true);
  add(0,2,2,"RH04_2_2","RH04_2_2_err",true);
  add(0,3,0,"RH04_3_0","RH04_3_0_err",true);
  add(0,3,1,"RH04_3_1","RH04_3_1_err",true);
  add(0,3,2,"RH04_3_2","RH04_3_2_err",true);
  add(0,3,3,"RH04_3_3","RH04_3_3_err",true);
  add(0,4,0,"RH04_4_0","RH04_4_0_err",true);
  add(0,4,1,"RH04_4_1","RH04_4_1_err",true);
  add(0,4,2,"RH04_4_2","RH04_4_2_err",true);
  add(0,4,3,"RH04_4_3","RH04_4_3_err",true);
  add(0,4,4,"RH04_4_4","RH04_4_4_err",true);

  add(1,0,0,"RH_1_0_0","RH_1_0_0_err",false);
  add(1,1,0,"RH_1_1_0","RH_1_1_0_err",false);
  add(1,1,1,"RH_1_1_1","RH_1_1_1_err",false);
  add(1,2,0,"RH_1_2_0","RH_1_2_0_err",false);
  add(1,2,1,"RH_1_2_1","RH_1_2_1_err",false);
  add(1,2,2,"RH_1_2_2","RH_1_2_2_err",false);
  add(1,3,0,"RH_1_3_0","RH_1_3_0_err",false);
  add(1,3,1,"RH_1_3_1","RH_1_3_1_err",false);
  add(1,3,2,"RH_1_3_2","RH_1_3_2_err",false);
  add(1,3,3,"RH_1_3_3","RH_1_3_3_err",false);
  add(1,4,0,"RH_1_4_0","RH_1_4_0_err",false);
  add(1,4,1,"RH_1_4_1","RH_1_4_1_err",false);
  add(1,4,2,"RH_1_4_2","RH_1_4_2_err",false);
  add(1,4,3,"RH_1_4_3","RH_1_4_3_err",false);
  add(1,4,4,"RH_1_4_4","RH_1_4_4_err",false);

  add(2,1,1,"RH_2_1_1","RH_2_1_1_err",false);
  add(2,2,1,"RH_2_2_1","RH_2_2_1_err",false);
  add(2,2,2,"RH_2_2_2","RH_2_2_2_err",false);
  add(2,3,1,"RH_2_3_1","RH_2_3_1_err",false);
  add(2,3,2,"RH_2_3_2","RH_2_3_2_err",false);
  add(2,3,3,"RH_2_3_3","RH_2_3_3_err",false);
  add(2,4,1,"RH_2_4_1","RH_2_4_1_err",false);
  add(2,4,2,"RH_2_4_2","RH_2_4_2_err",false);
  add(2,4,3,"RH_2_4_3","RH_2_4_3_err",false);
  add(2,4,4,"RH_2_4_4","RH_2_4_4_err",false);

  add(3,1,1,"RH_3_1_1","RH_3_1_1_err",false);
  add(3,2,1,"RH_3_2_1","RH_3_2_1_err",false);
  add(3,2,2,"RH_3_2_2","RH_3_2_2_err",false);
  add(3,3,1,"RH_3_3_1","RH_3_3_1_err",false);
  add(3,3,2,"RH_3_3_2","RH_3_3_2_err",false);
  add(3,3,3,"RH_3_3_3","RH_3_3_3_err",false);
  add(3,4,1,"RH_3_4_1","RH_3_4_1_err",false);
  add(3,4,2,"RH_3_4_2","RH_3_4_2_err",false);
  add(3,4,3,"RH_3_4_3","RH_3_4_3_err",false);
  add(3,4,4,"RH_3_4_4","RH_3_4_4_err",false);

  add(5,0,0,"RH_5_0_0","RH_5_0_0_err",false);
  add(5,1,0,"RH_5_1_0","RH_5_1_0_err",false);
  add(5,1,1,"RH_5_1_1","RH_5_1_1_err",false);
  add(5,2,0,"RH_5_2_0","RH_5_2_0_err",false);
  add(5,2,1,"RH_5_2_1","RH_5_2_1_err",false);
  add(5,2,2,"RH_5_2_2","RH_5_2_2_err",false);
  add(5,3,0,"RH_5_3_0","RH_5_3_0_err",false);
  add(5,3,1,"RH_5_3_1","RH_5_3_1_err",false);
  add(5,3,2,"RH_5_3_2","RH_5_3_2_err",false);
  add(5,3,3,"RH_5_3_3","RH_5_3_3_err",false);
  add(5,4,0,"RH_5_4_0","RH_5_4_0_err",false);
  add(5,4,1,"RH_5_4_1","RH_5_4_1_err",false);
  add(5,4,2,"RH_5_4_2","RH_5_4_2_err",false);
  add(5,4,3,"RH_5_4_3","RH_5_4_3_err",false);
  add(5,4,4,"RH_5_4_4","RH_5_4_4_err",false);

  add(6,1,1,"RH_6_1_1","RH_6_1_1_err",false);
  add(6,2,1,"RH_6_2_1","RH_6_2_1_err",false);
  add(6,2,2,"RH_6_2_2","RH_6_2_2_err",false);
  add(6,3,1,"RH_6_3_1","RH_6_3_1_err",false);
  add(6,3,2,"RH_6_3_2","RH_6_3_2_err",false);
  add(6,3,3,"RH_6_3_3","RH_6_3_3_err",false);
  add(6,4,1,"RH_6_4_1","RH_6_4_1_err",false);
  add(6,4,2,"RH_6_4_2","RH_6_4_2_err",false);
  add(6,4,3,"RH_6_4_3","RH_6_4_3_err",false);
  add(6,4,4,"RH_6_4_4","RH_6_4_4_err",false);

  add(7,1,1,"RH_7_1_1","RH_7_1_1_err",false);
  add(7,2,1,"RH_7_2_1","RH_7_2_1_err",false);
  add(7,2,2,"RH_7_2_2","RH_7_2_2_err",false);
  add(7,3,1,"RH_7_3_1","RH_7_3_1_err",false);
  add(7,3,2,"RH_7_3_2","RH_7_3_2_err",false);
  add(7,3,3,"RH_7_3_3","RH_7_3_3_err",false);
  add(7,4,1,"RH_7_4_1","RH_7_4_1_err",false);
  add(7,4,2,"RH_7_4_2","RH_7_4_2_err",false);
  add(7,4,3,"RH_7_4_3","RH_7_4_3_err",false);
  add(7,4,4,"RH_7_4_4","RH_7_4_4_err",false);

  add(8,0,0,"RH_8_0_0","RH_8_0_0_err",false);
  add(8,1,0,"RH_8_1_0","RH_8_1_0_err",false);
  add(8,1,1,"RH_8_1_1","RH_8_1_1_err",false);
  add(8,2,0,"RH_8_2_0","RH_8_2_0_err",false);
  add(8,2,1,"RH_8_2_1","RH_8_2_1_err",false);
  add(8,2,2,"RH_8_2_2","RH_8_2_2_err",false);
  add(8,3,0,"RH_8_3_0","RH_8_3_0_err",false);
  add(8,3,1,"RH_8_3_1","RH_8_3_1_err",false);
  add(8,3,2,"RH_8_3_2","RH_8_3_2_err",false);
  add(8,3,3,"RH_8_3_3","RH_8_3_3_err",false);
  add(8,4,0,"RH_8_4_0","RH_8_4_0_err",false);
  add(8,4,1,"RH_8_4_1","RH_8_4_1_err",false);
  add(8,4,2,"RH_8_4_2","RH_8_4_2_err",false);
  add(8,4,3,"RH_8_4_3","RH_8_4_3_err",false);
  add(8,4,4,"RH_8_4_4","RH_8_4_4_err",false);

  return obs;
}


// Numerical Calculation stuff below here

static std::vector<std::pair<int, int>> EnumerateWaves(int lmax, int mmax, bool negm, bool onlyEven) {
  std::vector<std::pair<int, int>> waves;
  waves.reserve((lmax + 1) * (2 * mmax + 1));
  for (int l = 0; l <= lmax; ++l) {
    if (onlyEven && (l & 1)) continue;
    for (int m = -l; m <= l; ++m) {
      if (std::abs(m) > mmax) continue;
      if (!negm && m < 0) continue;
      waves.emplace_back(l, m);
    }
  }
  return waves;
}

static std::string MString(int m) { return (m < 0) ? "m"+std::to_string(-m) : std::to_string(m); }
static std::string MagName(char refl, char orient, int l, int m) { return std::string(1, refl) + '_' + orient + '_' + std::to_string(l) + '_' + MString(m); }
static std::string PhiName(char refl, char orient, int l, int m) { return std::string(1, refl) + "phi_" + orient + '_' + std::to_string(l) + '_' + MString(m); }

struct ParamLabel {
  char refl = '\0';
  char orient = '\0';
  int l = -1;
  int m = 0;
  bool isPhase = false;
  bool valid = false;
};

static ParamLabel ParseParamLabel(const std::string& name) {
  ParamLabel out;
  if (name.size() < 5) return out;
  out.refl = name[0];
  out.isPhase = (name.find("phi_") != std::string::npos);
  out.orient = (name.find("_L_") != std::string::npos || name.find("phi_L_") != std::string::npos) ? 'L' : 'T';
  const size_t last = name.find_last_of('_');
  if (last == std::string::npos || last + 1 >= name.size()) return out;
  const size_t prev = name.find_last_of('_', last - 1);
  if (prev == std::string::npos || prev + 1 >= last) return out;
  out.l = std::stoi(name.substr(prev + 1, last - prev - 1));
  const std::string mstr = name.substr(last + 1);
  out.m = (!mstr.empty() && mstr[0] == 'm')
          ? -std::stoi(mstr.substr(1))
          : std::stoi(mstr);
  out.valid = true;
  return out;
}

static long long MakeParamKey(char refl, char orient, int l, int m, bool isPhase) {
  const long long reflBit = (refl == 'b');
  const long long orientBit = (orient == 'L');
  const long long phaseBit = isPhase;
  const long long mEnc = static_cast<long long>(m + 32);
  return ( ((((phaseBit << 1) | reflBit) << 1) | orientBit) << 12) | (static_cast<long long>(l) << 6) | mEnc;
}

static int ReflectivitySign(char refl) {
  return (refl == 'b') ? -1 : +1;
}

static int LongitudinalParitySign(char refl, int absM) {
  const int mParity = (absM & 1) ? -1 : +1;
  return ReflectivitySign(refl) * mParity;
}

static bool SkipLongitudinalNegativeM(const FitConfig& cfg, char orient, int m) {
  return cfg.enforceLongitudinalParity && orient == 'L' && m < 0;
}

static std::vector<ParDef> BuildAmplitudePhaseParameters(const FitConfig& cfg) {
  std::vector<ParDef> pars;
  const auto waves = EnumerateWaves(cfg.lmax, cfg.mmax, cfg.negm, cfg.onlyEven);
  pars.reserve(waves.size() * 8);

  auto add = [&](char refl, char orient, int l, int m, bool isPhase) {
    if (SkipLongitudinalNegativeM(cfg, orient, m)) return;

    ParDef p;
    p.name = isPhase ? PhiName(refl, orient, l, m) : MagName(refl, orient, l, m);
    p.init = 0.0;
    p.step = isPhase ? 0.3 : 0.05;
    if (cfg.photoProduction)
    {
      p.low = 0.0; // set to 0 for photoproduction due to ambiguity i.e. one appears in both sides
    } else
    {
      p.low = isPhase ? -kPi : 0.0;
    }
    p.high = isPhase ? kPi : cfg.magnitudeMax;
    p.isPhase = isPhase;

    // For longitudinal m=0, negative reflectivity is odd under the parity relation
    // A_m = eps (-1)^m A_-m and therefore vanishes.
    if (cfg.enforceLongitudinalParity && orient == 'L' && refl == 'b' && m == 0) {
      p.init = 0.0;
      p.fixed = true;
      p.low = 0.0;
      p.high = 0.0;
      p.step = 0.0;
    }

    pars.push_back(std::move(p));
  };

  for (const auto& w : waves) {
    const int l = w.first;
    const int m = w.second;
    add('a', 'T', l, m, false); add('a', 'L', l, m, false); add('b', 'T', l, m, false); add('b', 'L', l, m, false);
    add('a', 'T', l, m, true);  add('a', 'L', l, m, true);  add('b', 'T', l, m, true);  add('b', 'L', l, m, true);
  }

  auto fixTo = [&](const std::string& name, double value) {
    auto it = std::find_if(pars.begin(), pars.end(), [&](const ParDef& p) { return p.name == name; });
    if (it == pars.end()) throw std::runtime_error("Parameter not found to fix: " + name);
    it->init = value; it->fixed = true; it->low = value; it->high = value; it->step = 0.0;
  };

  // Fix one transverse reference phase per reflectivity to remove the global phase ambiguity.
  const int refL = cfg.lmax;
  const int refM = std::min(cfg.mmax, refL);
  fixTo(PhiName('a', 'T', 2, 2), 0.0);
  fixTo(PhiName('b', 'T', 2, 2), 0.0);

  // P-Wave Transveres
  //fixTo("a_L_1_0", 0.0); fixTo("aphi_L_1_0", 0.0); //fixTo("b_L_1_0", 0.0); fixTo("bphi_L_1_0", 0.0);
  // fixTo("a_L_1_1", 0.0); fixTo("aphi_L_1_1", 0.0); fixTo("b_L_1_1", 0.0); fixTo("bphi_L_1_1", 0.0);
  // fixTo("a_L_1_m1", 0.0); fixTo("aphi_L_1_m1", 0.0); fixTo("b_L_1_m1", 0.0); fixTo("bphi_L_1_m1", 0.0);

  // P-Wave Longitudinal
  //fixTo("a_L_1_0", 0.0); fixTo("aphi_L_1_0", 0.0);
  //fixTo("b_L_1_0", 0.0); fixTo("bphi_L_1_0", 0.0);
  //fixTo("a_L_1_1", 0.0); fixTo("aphi_L_1_1", 0.0); fixTo("b_L_1_1", 0.0); fixTo("bphi_L_1_1", 0.0);
  //fixTo("a_L_1_m1", 0.0); fixTo("aphi_L_1_m1", 0.0); fixTo("b_L_1_m1", 0.0); fixTo("bphi_L_1_m1", 0.0);

  // S-Wave
 // fixTo("b_T_0_0", 0.0); fixTo("a_T_0_0", 0.0); fixTo("a_L_0_0", 0.0); fixTo("b_L_0_0", 0.0);
 // fixTo("bphi_T_0_0", 0.0); fixTo("aphi_T_0_0", 0.0); fixTo("aphi_L_0_0", 0.0); fixTo("bphi_L_0_0", 0.0);

  if (cfg.photoProduction) {
    for (auto& p : pars) {
      const auto label = ParseParamLabel(p.name);
      if (!label.valid || label.orient != 'L') continue;
      p.init = 0.0;
      p.fixed = true;
      p.low = 0.0;
      p.high = 0.0;
      p.step = 0.0;
    }
  }

  return pars;
}

struct BruSelection {
  double coeff = 0.0;
  char refl1 = 'a';
  char orient1 = 'T';
  int l1 = 0;
  int m1 = 0;
  char refl2 = 'a';
  char orient2 = 'T';
  int l2 = 0;
  int m2 = 0;
  TrigKind trig = TrigKind::kCos;
};

static bool ResolveBruSelection(int reflsign, double factor, int l, int m, int lpr, int mpr,
                                int alpha, bool negm, bool orientSwap, BruSelection& out) {
  if (!negm && (m < 0 || mpr < 0)) return false;
  out.refl1 = (reflsign == -1) ? 'b' : 'a';
  out.refl2 = (reflsign == -1) ? 'b' : 'a';

  out.orient1 = 'T';
  out.orient2 = 'T';
  if (alpha == 4) out.orient1 = out.orient2 = 'L';
  if (alpha >= 5) { out.orient1 = 'L'; out.orient2 = 'T'; }
  if (alpha >= 5 && orientSwap) { out.orient1 = 'T'; out.orient2 = 'L'; }

  int reflfactor = 1;
  if (alpha == 1 || alpha == 2) reflfactor = reflsign;
  if (alpha == 8) factor *= -1.0;

  out.coeff = reflfactor * factor;
  out.l1 = l; out.m1 = m; out.l2 = lpr; out.m2 = mpr;
  out.trig = (alpha == 3 || alpha == 7 || alpha == 8) ? TrigKind::kSin : TrigKind::kCos;
  return out.coeff != 0.0;
}

struct EvalContext {
  FitConfig cfg;
  std::vector<ParDef> fullPars;
  std::vector<int> freeToFull;
  std::vector<int> fullToFree;

  std::vector<ObservedMoment> observed;
  std::vector<MomentModel> modelsRec;
  std::unordered_map<std::string, size_t> modelIndexByName;
  std::vector<int> observedModelIdx;
  std::vector<int> observedModelIdx0;
  std::vector<int> observedModelIdx4;

  std::vector<PhasePair> phasePairs;
  int idxH0_00 = -1;
  int idxH4_00 = -1;

  mutable unsigned callCount = 0;
  mutable TTree* iterTree = nullptr;
  mutable double iter_log_val = 0.0;
  mutable std::vector<double> iter_parVals;
  mutable std::vector<double> iter_Hrec;
};

struct MCMCResult {
  std::vector<double> xStart;
  std::vector<double> xBest;
  double chi2Start = 1e300;
  double chi2Best = 1e300;
  double chi2Last = 1e300;
  unsigned accepted = 0;
  unsigned proposed = 0;
};

static std::vector<MomentModel> BuildMomentModels(const FitConfig& cfg,
                                                  const std::unordered_map<long long, int>& paramIndex,
                                                  std::vector<PhasePair>& phasePairs) {
  std::vector<MomentModel> models;
  const int alphaMax = cfg.photoProduction ? 3 : 8;
  models.reserve((alphaMax + 1) * (2 * cfg.lmax + 1) * (2 * cfg.lmax + 2) / 2);

  std::unordered_map<long long, int> phasePairLookup;
  const auto waves = EnumerateWaves(cfg.lmax, cfg.mmax, cfg.negm, cfg.onlyEven);

  struct CGKey {
    int lpr, L, l, mpr, M, m;
    bool operator==(const CGKey& o) const { return lpr == o.lpr && L == o.L && l == o.l && mpr == o.mpr && M == o.M && m == o.m; }
  };
  struct CGKeyHash {
    size_t operator()(const CGKey& k) const {
      size_t h = 1469598103934665603ull;
      auto mix = [&](int v) { h ^= static_cast<size_t>(v + 32); h *= 1099511628211ull; };
      mix(k.lpr); mix(k.L); mix(k.l); mix(k.mpr); mix(k.M); mix(k.m);
      return h;
    }
  };
  std::unordered_map<CGKey, double, CGKeyHash> cgCache;
  cgCache.reserve(2048);

  auto getCG = [&](int lpr, int L, int l, int mpr, int M, int m) {
    const CGKey key{lpr, L, l, mpr, M, m};
    auto it = cgCache.find(key);
    if (it != cgCache.end()) return it->second;
    const double v = ClebschGordan(lpr, L, l, mpr, M, m);
    cgCache.emplace(key, v);
    return v;
  };

  auto getPhasePairIdx = [&](int idxPhi1, int idxPhi2) {
    const long long key = (static_cast<long long>(idxPhi1) << 32) | static_cast<unsigned int>(idxPhi2);
    auto it = phasePairLookup.find(key);
    if (it != phasePairLookup.end()) return it->second;
    const int idx = static_cast<int>(phasePairs.size());
    phasePairs.push_back({idxPhi1, idxPhi2});
    phasePairLookup.emplace(key, idx);
    return idx;
  };

  auto paramIdx = [&](char refl, char orient, int l, int m, bool isPhase) -> int {
    const auto it = paramIndex.find(MakeParamKey(refl, orient, l, m, isPhase));
    if (it == paramIndex.end()) {
      std::ostringstream os;
      os << "Missing parameter index for " << (isPhase ? "phi" : "mag") << ' ' << refl << ' ' << orient << ' ' << l << ' ' << m;
      throw std::runtime_error(os.str());
    }
    return it->second;
  };

  auto emit = [&](MomentModel& mm, int reflsign, double factor, int l, int m, int lpr, int mpr, int alpha, bool orientSwap) {
    BruSelection sel;
    if (!ResolveBruSelection(reflsign, factor, l, m, lpr, mpr, alpha, cfg.negm, orientSwap, sel)) return;

    // Longitudinal parity relates the negative-m amplitude to the positive-m one,
    // A_m = eps (-1)^m A_-m.  We therefore fit only m >= 0 longitudinal
    // amplitudes and fold negative-m terms onto the positive-m parameter with the
    // corresponding sign in the bilinear coefficient.
    double coeff = sel.coeff;
    auto foldLongitudinalM = [&](char refl, char orient, int& waveM) {
      if (!cfg.enforceLongitudinalParity || orient != 'L' || waveM >= 0) return;
      const int absM = -waveM;
      coeff *= LongitudinalParitySign(refl, absM);
      waveM = absM;
    };
    foldLongitudinalM(sel.refl1, sel.orient1, sel.m1);
    foldLongitudinalM(sel.refl2, sel.orient2, sel.m2);

    Term t;
    t.coeff = coeff;
    t.idxMag1 = paramIdx(sel.refl1, sel.orient1, sel.l1, sel.m1, false);
    t.idxMag2 = paramIdx(sel.refl2, sel.orient2, sel.l2, sel.m2, false);
    const int idxPhi1 = paramIdx(sel.refl1, sel.orient1, sel.l1, sel.m1, true);
    const int idxPhi2 = paramIdx(sel.refl2, sel.orient2, sel.l2, sel.m2, true);
    t.phasePairIdx = getPhasePairIdx(idxPhi1, idxPhi2);
    t.trig = sel.trig;
    t.ignorePhase = (sel.l1 == sel.l2 && sel.m1 == sel.m2);
    mm.terms.push_back(t);
  };

  for (int alpha = 0; alpha <= alphaMax; ++alpha) {
    for (int L = 0; L <= 2 * cfg.lmax; ++L) {
      for (int M = 0; M <= L; ++M)
      {
        MomentModel mm;
        mm.alpha = alpha; mm.L = L; mm.M = M;
        mm.name = std::string("H_") + std::to_string(alpha) + "_" + std::to_string(L) + "_" + std::to_string(M);
        if ((alpha == 2 || alpha == 3 || alpha == 6 || alpha == 7) && M == 0) {
          models.push_back(std::move(mm));
          continue;
        }

        for (const auto& w1 : waves)
        {
          const int il = w1.first;
          const int im = w1.second;
          for (const auto& w2 : waves)
          {
            const int ilpr = w2.first;
            const int impr = w2.second;

            const double CM = getCG(ilpr, L, il, impr, M, im);
            if (CM == 0.0) continue;
            const double C0 = getCG(ilpr, L, il, 0, 0, 0);
            if (C0 == 0.0) continue;
            const double ccfactor = CM * C0 * std::sqrt((2.0 * ilpr + 1.0) / (2.0 * il + 1.0));
            if (ccfactor == 0.0) continue;

            const int mmprimesign = ((std::abs(im - impr) & 1) ? -1 : 1);
            const int mprimesign = ((std::abs(impr) & 1) ? -1 : 1);
            const int msign = ((std::abs(im) & 1) ? -1 : 1);

            if (alpha == 0) {
              emit(mm, +1, ccfactor, il, im, ilpr, impr, 0, false);
              emit(mm, +1, mmprimesign * ccfactor, il, -im, ilpr, -impr, 0, false);
              if (cfg.useNegRef) { emit(mm, -1, ccfactor, il, im, ilpr, impr, 0, false); emit(mm, -1, mmprimesign * ccfactor, il, -im, ilpr, -impr, 0, false); }
            } else if (alpha == 1) {
              emit(mm, +1, msign * ccfactor, il, -im, ilpr, impr, 1, false);
              emit(mm, +1, mprimesign * ccfactor, il, im, ilpr, -impr, 1, false);
              if (cfg.useNegRef) { emit(mm, -1, msign * ccfactor, il, -im, ilpr, impr, 1, false); emit(mm, -1, mprimesign * ccfactor, il, im, ilpr, -impr, 1, false); }
            } else if (alpha == 2) {
              emit(mm, +1, msign * ccfactor, il, -im, ilpr, impr, 2, false);
              emit(mm, +1, -mprimesign * ccfactor, il, im, ilpr, -impr, 2, false);
              if (cfg.useNegRef) { emit(mm, -1, msign * ccfactor, il, -im, ilpr, impr, 2, false); emit(mm, -1, -mprimesign * ccfactor, il, im, ilpr, -impr, 2, false); }
            } else if (alpha == 3) {
              const double f = -ccfactor;
              emit(mm, +1, f, il, im, ilpr, impr, 3, false);
              emit(mm, +1, -mmprimesign * f, il, -im, ilpr, -impr, 3, false);
              if (cfg.useNegRef) { emit(mm, -1, f, il, im, ilpr, impr, 3, false); emit(mm, -1, -mmprimesign * f, il, -im, ilpr, -impr, 3, false); }
            } else if (alpha == 4) {
              const double f = 2.0 * ccfactor;
              emit(mm, +1, f, il, im, ilpr, impr, 4, false);
              if (cfg.useNegRef) emit(mm, -1, f, il, im, ilpr, impr, 4, false);
            }  else if (alpha == 5) {
              int refl = +1;
              double f = (-1.0 / TMath::Sqrt(2.0)) * ccfactor;
              emit(mm, refl, f, il, im, ilpr, impr, 5, false);
              emit(mm, refl, f, il, im, ilpr, impr, 5, true);
              emit(mm, refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 5, false);
              emit(mm, refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 5, true);
              if (cfg.useNegRef) {
                refl = -1;
                emit(mm, refl, f, il, im, ilpr, impr, 5, false);
                emit(mm, refl, f, il, im, ilpr, impr, 5, true);
                emit(mm, refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 5, false);
                emit(mm, refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 5, true);
              }
            } else if (alpha == 6) {
              int refl = +1;
              double f = (1.0 / TMath::Sqrt(2.0)) * ccfactor;
              emit(mm, refl, f, il, im, ilpr, impr, 6, false);
              emit(mm, refl, f, il, im, ilpr, impr, 6, true);
              emit(mm, refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 6, false);
              emit(mm, refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 6, true);
              if (cfg.useNegRef) {
                refl = -1;
                emit(mm, refl, f, il, im, ilpr, impr, 6, false);
                emit(mm, refl, f, il, im, ilpr, impr, 6, true);
                emit(mm, refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 6, false);
                emit(mm, refl, -1 * mmprimesign * f, il, -im, ilpr, -impr, 6, true);
              }
            } else if (alpha == 7) {
              int refl = +1;
              double f = (-1.0 / TMath::Sqrt(2.0)) * ccfactor;
              emit(mm, refl, f, il, im, ilpr, impr, 7, false);
              emit(mm, refl, f, il, im, ilpr, impr, 7, true);
              emit(mm, refl, mmprimesign * f, il, -im, ilpr, -impr, 7, false);
              emit(mm, refl, mmprimesign * f, il, -im, ilpr, -impr, 7, true);
              if (cfg.useNegRef) {
                refl = -1;
                emit(mm, refl, f, il, im, ilpr, impr, 7, false);
                emit(mm, refl, f, il, im, ilpr, impr, 7, true);
                emit(mm, refl, mmprimesign * f, il, -im, ilpr, -impr, 7, false);
                emit(mm, refl, mmprimesign * f, il, -im, ilpr, -impr, 7, true);
              }
            } else if (alpha == 8) {
              int refl = +1;
              double f = (-1.0 / TMath::Sqrt(2.0)) * ccfactor;
              emit(mm, refl, f, il, im, ilpr, impr, 8, false);
              emit(mm, refl, -1 * f, il, im, ilpr, impr, 8, true);
              emit(mm, refl, -1*mmprimesign * f, il, -im, ilpr, -impr, 8, false);
              emit(mm, refl,  mmprimesign * f, il, -im, ilpr, -impr, 8, true);
              if (cfg.useNegRef) {
                refl = -1;
                emit(mm, refl, f, il, im, ilpr, impr, 8, false);
                emit(mm, refl, -1 * f, il, im, ilpr, impr, 8, true);
                emit(mm, refl, -1*mmprimesign * f, il, -im, ilpr, -impr, 8, false);
                emit(mm, refl, mmprimesign * f, il, -im, ilpr, -impr, 8, true);
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

static double SampleGammaMT(TRandom3& rng, double shape, double scale = 1.0) {
  if (!(shape > 0.0) || !(scale > 0.0) || !std::isfinite(shape) || !std::isfinite(scale)) return 0.0;

  if (shape < 1.0) {
    const double u = std::max(rng.Rndm(), 1e-300);
    return SampleGammaMT(rng, shape + 1.0, scale) * std::pow(u, 1.0 / shape);
  }

  const double d = shape - 1.0 / 3.0;
  const double c = 1.0 / std::sqrt(9.0 * d);
  for (;;) {
    const double x = rng.Gaus(0.0, 1.0);
    double v = 1.0 + c * x;
    if (v <= 0.0) continue;
    v = v * v * v;
    const double u = std::max(rng.Rndm(), 1e-300);
    if (u < 1.0 - 0.0331 * x * x * x * x) return scale * d * v;
    if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return scale * d * v;
  }
}

static bool FillFullFromFree(const EvalContext& ctx,
                             const double* x,
                             std::vector<double>& fullVals) {
  if (fullVals.size() != ctx.fullPars.size()) fullVals.resize(ctx.fullPars.size());
  for (size_t i = 0; i < ctx.fullPars.size(); ++i) fullVals[i] = ctx.fullPars[i].init;

  for (size_t i = 0; i < ctx.freeToFull.size(); ++i) {
    const int fullIdx = ctx.freeToFull[i];
    const auto& p = ctx.fullPars[static_cast<size_t>(fullIdx)];
    const double value = x[i];
    if (!std::isfinite(value)) return false;
    if (!p.isPhase && (value < p.low || value > p.high)) return false;
    fullVals[static_cast<size_t>(fullIdx)] = value;
  }
  return true;
}

static inline void EnsureSize(std::vector<double>& v, size_t n, double fill = 0.0) {
  if (v.size() != n) v.assign(n, fill);
}

static void BuildPhasePairTrigCache(const EvalContext& ctx,
                                    const std::vector<double>& fullVals,
                                    std::vector<double>& pairSin,
                                    std::vector<double>& pairCos) {
  EnsureSize(pairSin, ctx.phasePairs.size());
  EnsureSize(pairCos, ctx.phasePairs.size());
  for (size_t i = 0; i < ctx.phasePairs.size(); ++i) {
    const auto& pp = ctx.phasePairs[i];
    FastSinCos(fullVals[pp.idxPhi1] - fullVals[pp.idxPhi2], pairSin[i], pairCos[i]);
  }
}

static double EvalMomentOnly(const MomentModel& mm,
                             const std::vector<double>& fullVals,
                             const std::vector<double>& pairSin,
                             const std::vector<double>& pairCos) {
  double H = 0.0;
  for (const auto& t : mm.terms) {
    // Had these as we do not need to calculate these but stopped working...
    //if (t.ignorePhase && (mm.alpha==3 || mm.alpha==7 || mm.alpha==8)) continue; // No imaginary parts for these
    //const double trig = t.ignorePhase ? 1.0 : ((t.trig == TrigKind::kCos) ? pairCos[t.phasePairIdx] : pairSin[t.phasePairIdx]);
    const double trig = (t.trig == TrigKind::kCos) ? pairCos[t.phasePairIdx] : pairSin[t.phasePairIdx];
    H += t.coeff * fullVals[t.idxMag1] * fullVals[t.idxMag2] * trig;
  }
  return H;
}

static double EvalMomentAndDerivFull(const EvalContext& ctx,
                                     const MomentModel& mm,
                                     const std::vector<double>& fullVals,
                                     const std::vector<double>& pairSin,
                                     const std::vector<double>& pairCos,
                                     std::vector<double>& dHdFull) {
  EnsureSize(dHdFull, ctx.fullPars.size());
  std::fill(dHdFull.begin(), dHdFull.end(), 0.0);
  double H = 0.0;
  for (const auto& t : mm.terms) {
    const auto& pp = ctx.phasePairs[t.phasePairIdx];
    const double m1 = fullVals[t.idxMag1]; // Magnitudes
    const double m2 = fullVals[t.idxMag2];
    //if (t.ignorePhase && (mm.alpha==3 || mm.alpha==7 || mm.alpha==8)) continue; // No imaginary parts for these
    const double trig = (t.trig == TrigKind::kCos) ? pairCos[t.phasePairIdx] : pairSin[t.phasePairIdx];
    const double dtrig =(t.trig == TrigKind::kCos) ? -pairSin[t.phasePairIdx] : pairCos[t.phasePairIdx];
    const double val = t.coeff * m1 * m2 * trig;
    H += val;
    dHdFull[t.idxMag1] += t.coeff * m2 * trig;
    dHdFull[t.idxMag2] += t.coeff * m1 * trig;
    //if (!t.ignorePhase) {
    const double common = t.coeff * m1 * m2 * dtrig;
    dHdFull[pp.idxPhi1] += common;
    dHdFull[pp.idxPhi2] -= common;
    //}
  }
  return H;
}

static void EvalAllMoments(const EvalContext& ctx,
                           const std::vector<double>& fullVals,
                           std::vector<double>& values) {
  EnsureSize(values, ctx.modelsRec.size());
  std::vector<double> pairSin, pairCos;
  BuildPhasePairTrigCache(ctx, fullVals, pairSin, pairCos);
  for (size_t i = 0; i < ctx.modelsRec.size(); ++i) values[i] = EvalMomentOnly(ctx.modelsRec[i], fullVals, pairSin, pairCos);
}

static double WrapToRange(double x, double low, double high) {
  if (!(high > low)) return low;
  const double width = high - low;
  while (x < low) x += width;
  while (x > high) x -= width;
  if (x < low) x = low;
  if (x > high) x = high;
  return x;
}


static double EvaluateChi2AtPoint(const ROOT::Math::IBaseFunctionMultiDim& fcn,
                                  const std::vector<double>& x) {
  if (x.empty()) return 1e300;
  const double chi2 = fcn(x.data());
  return (std::isfinite(chi2) ? chi2 : 1e300);
}

static void ProposeMCMCStep(const EvalContext& ctx,
                            const FitConfig& cfg,
                            TRandom3& rng,
                            const std::vector<double>& current,
                            std::vector<double>& proposal) {
  proposal = current;
  if (proposal.empty()) return;
  const unsigned iFree = rng.Integer(proposal.size());
  const int fullIdx = ctx.freeToFull[iFree];
  const auto& p = ctx.fullPars[static_cast<size_t>(fullIdx)];
  const double step = p.isPhase
                        ? std::max(cfg.mcmcProposalPhaseSigma, (p.step > 0.0 ? p.step : 0.0))
                        : std::max(cfg.mcmcProposalMagSigma, (p.step > 0.0 ? p.step : 0.0));
  const double trial = proposal[iFree] + rng.Gaus(0.0, step);
  proposal[iFree] = p.isPhase ? WrapToRange(trial, p.low, p.high)
                              : std::min(std::max(trial, p.low), p.high);
}

static MCMCResult RunMCMCPreScan(const EvalContext& ctx,
                                 const FitConfig& cfg,
                                 const ROOT::Math::IBaseFunctionMultiDim& fcn,
                                 TRandom3& rng,
                                 const std::vector<double>& xSeed) {
  MCMCResult out;
  out.xStart = xSeed;
  out.xBest = xSeed;
  if (xSeed.empty() || cfg.mcmcSteps == 0) {
    out.chi2Start = EvaluateChi2AtPoint(fcn, xSeed);
    out.chi2Best = out.chi2Start;
    out.chi2Last = out.chi2Start;
    return out;
  }

  std::vector<double> current = xSeed;
  std::vector<double> proposal = xSeed;
  double currentChi2 = EvaluateChi2AtPoint(fcn, current);
  if (!(currentChi2 < 1e299)) {
    out.chi2Start = currentChi2;
    out.chi2Best = currentChi2;
    out.chi2Last = currentChi2;
    return out;
  }

  out.chi2Start = currentChi2;
  out.chi2Best = currentChi2;
  out.chi2Last = currentChi2;

  const double temperature = (cfg.mcmcTemperature > 0.0) ? cfg.mcmcTemperature : 1.0;

  for (unsigned istep = 0; istep < cfg.mcmcSteps; ++istep) {
    ProposeMCMCStep(ctx, cfg, rng, current, proposal);
    const double propChi2 = EvaluateChi2AtPoint(fcn, proposal);
    ++out.proposed;
    bool accept = false;
    if (propChi2 < currentChi2) {
      accept = true;
    } else if (propChi2 < 1e299) {
      const double delta = (propChi2 - currentChi2) / temperature;
      accept = (delta <= 0.0) || (rng.Uniform() < std::exp(-0.5 * delta));
    }

    if (accept) {
      current.swap(proposal);
      currentChi2 = propChi2;
      ++out.accepted;
      if (currentChi2 < out.chi2Best) {
        out.chi2Best = currentChi2;
        out.xBest = current;
      }
    }
  }

  out.chi2Last = currentChi2;
  return out;
}

class Chi2Function final : public ROOT::Math::IMultiGradFunction {
public:
  explicit Chi2Function(std::shared_ptr<EvalContext> ctx) : ctx_(std::move(ctx)) {
    if (!ctx_) throw std::runtime_error("Chi2Function: null context");
  }

  unsigned int NDim() const override { return static_cast<unsigned>(ctx_->freeToFull.size()); }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override { return new Chi2Function(ctx_); }

  double DoEval(const double* x) const override {
    if (!FillFullFromFree(*ctx_, x, fullVals_)) return 1e300;
    BuildPhasePairTrigCache(*ctx_, fullVals_, pairSin_, pairCos_);
    ++ctx_->callCount;

    EnsureSize(buf_momRec_, ctx_->modelsRec.size());
    if (buf_momSeen_.size() != ctx_->modelsRec.size()) buf_momSeen_.assign(ctx_->modelsRec.size(), 0);
    std::fill(buf_momSeen_.begin(), buf_momSeen_.end(), static_cast<unsigned char>(0));

    auto getMomentRaw = [&](int idx) -> double {
      if (idx < 0) throw std::runtime_error("Observed moment is not mapped to a model moment");
      const size_t uidx = static_cast<size_t>(idx);
      if (!buf_momSeen_[uidx]) {
        buf_momRec_[uidx] = EvalMomentOnly(ctx_->modelsRec[uidx], fullVals_, pairSin_, pairCos_);
        buf_momSeen_[uidx] = 1;
      }
      return buf_momRec_[uidx];
    };

    double chi2 = 0.0;
    for (size_t i = 0; i < ctx_->observed.size(); ++i) {
      const auto& ob = ctx_->observed[i];
      double model = 0.0;
      if (ob.isMixed04) {
        const double H0 = getMomentRaw(ctx_->observedModelIdx0[i]);
        const double H4 = getMomentRaw(ctx_->observedModelIdx4[i]);
        model = H0 + ctx_->cfg.epsilon * H4;
      } else {
        model = getMomentRaw(ctx_->observedModelIdx[i]);
      }

      const double sigma = ob.sigma;
      const double r = (ob.value - model) / sigma;
      chi2 += r * r;
    }

    return chi2;
  }

  void Gradient(const double* x, double* grad) const override {
    std::fill(grad, grad + NDim(), 0.0);
    if (!FillFullFromFree(*ctx_, x, fullVals_)) return;
    BuildPhasePairTrigCache(*ctx_, fullVals_, pairSin_, pairCos_);

    const size_t nFull = ctx_->fullPars.size();
    EnsureSize(buf_fullA_, nFull);
    EnsureSize(buf_fullB_, nFull);

    for (size_t iobs = 0; iobs < ctx_->observed.size(); ++iobs) {
      const auto& ob = ctx_->observed[iobs];
      const double sigma = (std::isfinite(ob.sigma) && std::abs(ob.sigma) > 0.0) ? std::abs(ob.sigma) : 1.0;
      double model = 0.0;

      if (ob.isMixed04) {
        const int idx0 = ctx_->observedModelIdx0[iobs];
        const int idx4 = ctx_->observedModelIdx4[iobs];
        if (idx0 < 0 || idx4 < 0) throw std::runtime_error("RH04 moment is not mapped to H0/H4 model moments");

        const double H0 = EvalMomentAndDerivFull(*ctx_, ctx_->modelsRec[idx0], fullVals_, pairSin_, pairCos_, buf_fullA_);
        const double H4 = EvalMomentAndDerivFull(*ctx_, ctx_->modelsRec[idx4], fullVals_, pairSin_, pairCos_, buf_fullB_);
        model = H0 + ctx_->cfg.epsilon * H4;
        const double pull = (ob.value - model) / sigma;

        for (unsigned i = 0; i < NDim(); ++i) {
          const int fullIdx = ctx_->freeToFull[i];
          const double dModel = buf_fullA_[static_cast<size_t>(fullIdx)]
                              + ctx_->cfg.epsilon * buf_fullB_[static_cast<size_t>(fullIdx)];
          grad[i] += -2.0 * pull * dModel / sigma;
        }
      } else {
        const int midx = ctx_->observedModelIdx[iobs];
        if (midx < 0) throw std::runtime_error("Observed moment is not mapped to a model moment");

        model = EvalMomentAndDerivFull(*ctx_, ctx_->modelsRec[midx], fullVals_, pairSin_, pairCos_, buf_fullA_);
        const double pull = (ob.value - model) / sigma;

        for (unsigned i = 0; i < NDim(); ++i) {
          const int fullIdx = ctx_->freeToFull[i];
          const double dModel = buf_fullA_[static_cast<size_t>(fullIdx)];
          grad[i] += -2.0 * pull * dModel / sigma;
        }
      }
    }
  }

  double DoDerivative(const double* x, unsigned int icoord) const override {
    EnsureSize(buf_gradTmp_, NDim());
    Gradient(x, buf_gradTmp_.data());
    return (icoord < NDim()) ? buf_gradTmp_[icoord] : 0.0;
  }

private:
  std::shared_ptr<EvalContext> ctx_;
  mutable std::vector<double> fullVals_;
  mutable std::vector<double> pairSin_;
  mutable std::vector<double> pairCos_;
  mutable std::vector<double> buf_momRec_;
  mutable std::vector<unsigned char> buf_momSeen_;
  mutable std::vector<double> buf_fullA_, buf_fullB_;
  mutable std::vector<double> buf_gradTmp_;
};

class Chi2FunctionNoGrad final : public ROOT::Math::IBaseFunctionMultiDim {
public:
  explicit Chi2FunctionNoGrad(std::shared_ptr<EvalContext> ctx) : fcn_(std::move(ctx)) {}

  unsigned int NDim() const override { return fcn_.NDim(); }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override { return new Chi2FunctionNoGrad(*this); }
  double DoEval(const double* x) const override { return fcn_.DoEval(x); }

private:
  Chi2Function fcn_;
};

static std::shared_ptr<EvalContext> BuildContext(const FitConfig& cfg) {
  auto ctx = std::make_shared<EvalContext>();
  ctx->cfg = cfg;
  ctx->fullPars = BuildAmplitudePhaseParameters(cfg);

  std::unordered_map<long long, int> paramIndex;
  paramIndex.reserve(ctx->fullPars.size() * 2);
  for (int i = 0; i < static_cast<int>(ctx->fullPars.size()); ++i) {
    const auto& p = ctx->fullPars[static_cast<size_t>(i)];
    const auto label = ParseParamLabel(p.name);
    if (!label.valid) throw std::runtime_error("Could not parse parameter label: " + p.name);
    paramIndex.emplace(MakeParamKey(label.refl, label.orient, label.l, label.m, p.isPhase), i);
  }

  if (cfg.verbose) {
    std::cout << "Production mode: "
              << (cfg.photoProduction ? "photoproduction (alpha <= 3, L fixed to 0)" : "electroproduction/full")
              << std::endl;
    std::cout << "Magnitude parameters are fitted directly with no unit-normalisation constraint. Longitudinal negative-m waves are folded onto m >= 0 using parity." << std::endl;
  }

  ctx->fullToFree.assign(ctx->fullPars.size(), -1);
  for (int i = 0; i < static_cast<int>(ctx->fullPars.size()); ++i) {
    if (ctx->fullPars[static_cast<size_t>(i)].fixed) continue;
    ctx->fullToFree[static_cast<size_t>(i)] = static_cast<int>(ctx->freeToFull.size());
    ctx->freeToFull.push_back(i);
  }

  ctx->observed = BuildObservedMoments(cfg.momentsFile, cfg.momentsTree, cfg.bin, cfg);
  ctx->modelsRec = BuildMomentModels(cfg, paramIndex, ctx->phasePairs);
  ctx->modelIndexByName.reserve(ctx->modelsRec.size() * 2);
  for (size_t i = 0; i < ctx->modelsRec.size(); ++i) ctx->modelIndexByName.emplace(ctx->modelsRec[i].name, i);

  auto it0 = ctx->modelIndexByName.find("H_0_0_0");
  auto it4 = ctx->modelIndexByName.find("H_4_0_0");
  if (it0 != ctx->modelIndexByName.end()) ctx->idxH0_00 = static_cast<int>(it0->second);
  if (it4 != ctx->modelIndexByName.end()) ctx->idxH4_00 = static_cast<int>(it4->second);

  ctx->observedModelIdx.assign(ctx->observed.size(), -1);
  ctx->observedModelIdx0.assign(ctx->observed.size(), -1);
  ctx->observedModelIdx4.assign(ctx->observed.size(), -1);

  for (size_t i = 0; i < ctx->observed.size(); ++i) {
    const auto& ob = ctx->observed[i];
    if (ob.isMixed04) {
      const std::string key0 = std::string("H_0_") + std::to_string(ob.L) + "_" + std::to_string(ob.M);
      const std::string key4 = std::string("H_4_") + std::to_string(ob.L) + "_" + std::to_string(ob.M);
      auto itObs0 = ctx->modelIndexByName.find(key0);
      auto itObs4 = ctx->modelIndexByName.find(key4);
      if (itObs0 != ctx->modelIndexByName.end()) ctx->observedModelIdx0[i] = static_cast<int>(itObs0->second);
      if (itObs4 != ctx->modelIndexByName.end()) ctx->observedModelIdx4[i] = static_cast<int>(itObs4->second);
      if (ctx->observedModelIdx0[i] < 0 || ctx->observedModelIdx4[i] < 0) {
        throw std::runtime_error("Failed to map " + ob.name + " to H0/H4 model moments");
      }
    } else {
      std::string key = ob.name;
      if (key.rfind("RH_", 0) == 0) key = "H_" + key.substr(3);
      auto it = ctx->modelIndexByName.find(key);
      if (it != ctx->modelIndexByName.end()) ctx->observedModelIdx[i] = static_cast<int>(it->second);
      if (ctx->observedModelIdx[i] < 0) throw std::runtime_error("Failed to map observed moment " + ob.name);
    }
  }

  return ctx;
}

static void MakeBranchesForPars(TTree* t, const std::vector<ParDef>& pars, std::vector<double>& storage) {
  storage.assign(pars.size(), 0.0);
  for (size_t i = 0; i < pars.size(); ++i) t->Branch(pars[i].name.c_str(), &storage[i]);
}

static void MakeBranchesForMoments(TTree* t, const std::vector<MomentModel>& models, std::vector<double>& storage) {
  storage.assign(models.size(), 0.0);
  for (size_t i = 0; i < models.size(); ++i) t->Branch(models[i].name.c_str(), &storage[i]);
}

static void MakeBranchesForObservedModel(TTree* t, const std::vector<ObservedMoment>& observed, std::vector<double>& storage) {
  storage.assign(observed.size(), 0.0);
  for (size_t i = 0; i < observed.size(); ++i) t->Branch(observed[i].name.c_str(), &storage[i]);
}

static void FillObservedModelValues(const EvalContext& ctx,
                                    const std::vector<double>& rawMoments,
                                    std::vector<double>& modelVals) {
  modelVals.assign(ctx.observed.size(), 0.0);

  for (size_t i = 0; i < ctx.observed.size(); ++i) {
    const auto& ob = ctx.observed[i];
    double model = 0.0;
    if (ob.isMixed04) {
      const int idx0 = ctx.observedModelIdx0[i];
      const int idx4 = ctx.observedModelIdx4[i];
      if (idx0 < 0 || idx4 < 0) throw std::runtime_error("RH04 output moment is not mapped to H0/H4");
      model = rawMoments[static_cast<size_t>(idx0)]
            + ctx.cfg.epsilon * rawMoments[static_cast<size_t>(idx4)];
    } else {
      const int midx = ctx.observedModelIdx[i];
      if (midx < 0) throw std::runtime_error("Output moment is not mapped to a model moment");
      model = rawMoments[static_cast<size_t>(midx)];
    }
    modelVals[i] = model;
  }
}

static void BuildRandomStartPoint(const EvalContext& ctx,
                                  TRandom3& rng,
                                  std::vector<double>& xStart) {
  xStart.assign(ctx.freeToFull.size(), 0.0);

  for (unsigned i = 0; i < ctx.freeToFull.size(); ++i) {
    const int fullIdx = ctx.freeToFull[i];
    const auto& p = ctx.fullPars[static_cast<size_t>(fullIdx)];
    if (p.isPhase) {
      xStart[i] = rng.Uniform(p.low, p.high);
    } else {
      double value = rng.Gaus(ctx.cfg.magnitudeStartMean, ctx.cfg.magnitudeStartSigma);
      while (value < p.low || value > p.high) {
        value = rng.Gaus(ctx.cfg.magnitudeStartMean, ctx.cfg.magnitudeStartSigma);
      }
      xStart[i] = value;
    }
  }
}

} // namespace chi2_amp_fit_opt

namespace {

void RunGivenMoments_Chi2Amps_Impl(const chi2_amp_fit_opt::FitConfig& cfg, const char* outFile) {
  using namespace chi2_amp_fit_opt;
  auto ctx = BuildContext(cfg);

  TRandom3 rng(cfg.randomSeed);
  if (cfg.randomSeed == 0) rng.SetSeed(0);

  std::unique_ptr<TFile> fout(TFile::Open(outFile, "RECREATE"));
  if (!fout || fout->IsZombie()) throw std::runtime_error(std::string("Failed to open output file: ") + outFile);

  TTree* t = new TTree("fitResults", "Optimized chi2 fit results (per start)");
  double log_val = 0.0;
  double start_log_val = -999.0;
  double prescan_log_val = -999.0;
  double mcmc_acceptance = 0.0;

  t->Branch("log_val", &log_val);
  if (cfg.useMCMCPreScan){
    t->Branch("start_log_val", &start_log_val);
    t->Branch("prescan_log_val", &prescan_log_val);
    t->Branch("mcmc_acceptance", &mcmc_acceptance);
  }

  std::vector<double> parVals, momRecVals, rhRecVals;
  MakeBranchesForPars(t, ctx->fullPars, parVals);
  MakeBranchesForMoments(t, ctx->modelsRec, momRecVals);
  MakeBranchesForObservedModel(t, ctx->observed, rhRecVals);

  Chi2FunctionNoGrad fcnNoGrad(ctx);
  Chi2Function fcn(ctx);

  if (cfg.verbose) {
    std::cout << (cfg.useNumericalGradient
                    ? "Using Minuit2 numerical derivatives (analytical gradient disabled)."
                    : "Using analytical gradient.")
              << std::endl;
  }
  if (cfg.verbose && cfg.useMCMCPreScan && cfg.mcmcSteps > 0) {
    std::cout << "Using MCMC pre-scan before each Minuit fit: steps=" << cfg.mcmcSteps
              << ", T=" << cfg.mcmcTemperature
              << ", magSigma=" << cfg.mcmcProposalMagSigma
              << ", phaseSigma=" << cfg.mcmcProposalPhaseSigma << std::endl;
  }

  TBenchmark bench;
  bench.Start("fit");

  const unsigned progressStep = std::max(1u, cfg.nStarts / 10u);
  for (unsigned iStart = 0; iStart < cfg.nStarts; ++iStart) {
    ctx->callCount = 0;
    if (cfg.verbose && (iStart % progressStep == 0 || iStart + 1 == cfg.nStarts)) std::cout << "start " << (iStart + 1) << '/' << cfg.nStarts << std::endl;

    std::vector<double> startVals;
    BuildRandomStartPoint(*ctx, rng, startVals);
    const double startChi2 = EvaluateChi2AtPoint(fcnNoGrad, startVals);
    start_log_val = (startChi2 > 0.0) ? std::log10(startChi2) : -999.0;

    if (cfg.useMCMCPreScan && cfg.mcmcSteps > 0) {
      const auto mcmc = RunMCMCPreScan(*ctx, cfg, fcnNoGrad, rng, startVals);
      startVals = mcmc.xBest;
      prescan_log_val = (mcmc.chi2Best > 0.0) ? std::log10(mcmc.chi2Best) : -999.0;
      mcmc_acceptance = (mcmc.proposed > 0) ? (static_cast<double>(mcmc.accepted) / static_cast<double>(mcmc.proposed)) : 0.0;
    } else {
      prescan_log_val = start_log_val;
      mcmc_acceptance = 0.0;
    }

    std::unique_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    if (!min) throw std::runtime_error("Failed to create Minuit2/Migrad minimizer");
    min->SetMaxFunctionCalls(cfg.maxCalls);
    min->SetMaxIterations(cfg.maxIters);
    min->SetTolerance(cfg.tolerance);
    min->SetStrategy(cfg.strategy);
    min->SetPrintLevel(cfg.printLevel);
    if (cfg.useNumericalGradient) {
      min->SetFunction(fcnNoGrad);
    } else {
      min->SetFunction(fcn);
    }

    for (unsigned i = 0; i < fcn.NDim(); ++i) {
      const int fullIdx = ctx->freeToFull[i];
      const auto& p = ctx->fullPars[fullIdx];
      const double step = (p.step > 0.0) ? p.step : 1e-3;
      min->SetLimitedVariable(i, p.name.c_str(), startVals[i], step, p.low, p.high);
    }

    bool ok = min->Minimize();
    int status = min->Status();
    //if (!ok) continue;
    if (cfg.runHesse) min->Hesse();
    const double chi2 = min->MinValue();
    log_val = std::log10(chi2);
    if (!FillFullFromFree(*ctx, min->X(), parVals)) {
      throw std::runtime_error("Failed to map minimizer coordinates to physical amplitudes");
    }
    EvalAllMoments(*ctx, parVals, momRecVals);
    FillObservedModelValues(*ctx, momRecVals, rhRecVals);
    t->Fill();
  }

  bench.Stop("fit");
  bench.Print("fit");
  fout->Write();
  fout->Close();
  std::cout << "Saved: " << outFile << std::endl;
}

static std::string MakePartFileName(const char* outFile, unsigned workerId) {
  std::string base = outFile ? std::string(outFile) : std::string("resultsGivenMoments_optimized.root");
  const std::string ext = ".root";
  if (base.size() >= ext.size() && base.compare(base.size() - ext.size(), ext.size(), ext) == 0) base.erase(base.size() - ext.size());
  return base + ".part_" + std::to_string(workerId) + ".root";
}

} // namespace

void RunGivenMoments_Chi2Amps_Setup(const char* tableFile="InputFiles/Experiment/experimental_moments_table.root",
                                  const char* treeName = "expMoments",
                                  int bin = 0,
                                  unsigned nStarts = 10000,
                                  const char* outFile = "resultsGivenMoments_chi2_amps.root",
                                  uint32_t seed = 0,
                                  double epsilon = 1.0,
                                  bool useNumericalGradient = true,
                                  bool useMCMCPreScan = false,
                                  unsigned mcmcSteps = 2000,
                                  double mcmcProposalMagSigma = 0.03,
                                  double mcmcProposalPhaseSigma = 0.10,
                                  double mcmcTemperature = 1.0,
                                  double magnitudeStartMean = 0.5,
                                  double magnitudeStartSigma = 0.5,
                                  bool photoProduction = false) {
  chi2_amp_fit_opt::FitConfig cfg;
  cfg.nStarts = nStarts;
  cfg.randomSeed = seed;
  cfg.epsilon = epsilon;
  cfg.useNumericalGradient = useNumericalGradient;
  cfg.useMCMCPreScan = useMCMCPreScan;
  cfg.mcmcSteps = mcmcSteps;
  cfg.mcmcProposalMagSigma = mcmcProposalMagSigma;
  cfg.mcmcProposalPhaseSigma = mcmcProposalPhaseSigma;
  cfg.mcmcTemperature = mcmcTemperature;
  cfg.magnitudeStartMean = magnitudeStartMean;
  cfg.magnitudeStartSigma = magnitudeStartSigma;
  cfg.photoProduction = photoProduction;
  cfg.momentsFile = tableFile ? tableFile : "InputFiles/Experiment/experimental_moments_table.root";
  cfg.momentsTree = treeName ? treeName : "expMoments";
  cfg.bin = bin;
  cfg.verbose = true;

  RunGivenMoments_Chi2Amps_Impl(cfg, outFile);
}

// Macro function, filenames and bin as cml args so that bins can be looped over by  a script
void RunGivenMoments_Chi2Amps(const char* tableFile = "InputFiles/Experiment/e_rho_moments.root", const char* treeName = "expMoments", const int bin = 1, const std::string outFile = "Hermestest.root", const double epsilon = 1, const bool photoProduction = false) {

  // Setup
  // User settings
  const unsigned nStarts = 10000;
  const uint32_t seed = 0;
  //const double epsilon = 0.8; // virtual photon polarisation 0.8 for e rho 0.9 for muon rho and 0.96 for omega
  const unsigned int nCores = 10;

  // Fit / model options
  // photoProduction = true fits only alpha <= 3 and fixes longitudinal amplitudes/phases to 0.
  const bool useNumericalGradient = false;
  const double magnitudeStartMean = 0.5;
  const double magnitudeStartSigma = 0.5;

  // Optional MCMC prescan -- if fit is poor try use this to improve the original posterior distribution
  const bool useMCMCPreScan = false;
  const unsigned mcmcSteps = 2000;
  const double mcmcProposalMagSigma = 0.03;
  const double mcmcProposalPhaseSigma = 0.10;
  const double mcmcTemperature = 1.0;

  // ====================================================================================================================

  // Multicore processing
  // Split up the processes between each core
  if (nStarts == 0) throw std::runtime_error("RunGivenMoments_Chi2Amps: nStarts must be > 0");
  unsigned nWorkers = (nCores == 0) ? std::max(1u, std::thread::hardware_concurrency()) : nCores;
  if (nWorkers > nStarts) nWorkers = nStarts;

  std::string outDir = "./OutputFiles/";
  std::string out = outDir + outFile;

  std::vector<unsigned> workerIds(nWorkers);
  for (unsigned i = 0; i < nWorkers; ++i) workerIds[i] = i;
  ROOT::TProcessExecutor pool(nWorkers);

  // Each core writes to a seperate file
  auto partFiles = pool.Map([=](unsigned workerId) {
    const unsigned baseStarts = nStarts / nWorkers;
    const unsigned extra = nStarts % nWorkers;
    const unsigned myStarts = baseStarts + ((workerId < extra) ? 1u : 0u);
    std::string partFile = MakePartFileName(out.c_str(), workerId);
    const uint32_t workerSeed = (seed == 0) ? (0x9e3779b9u + 100003u * workerId) : (seed + 100003u * workerId);
    RunGivenMoments_Chi2Amps_Setup(tableFile, treeName, bin, myStarts, partFile.c_str(), workerSeed, epsilon,
                                   useNumericalGradient, useMCMCPreScan, mcmcSteps,
                                   mcmcProposalMagSigma, mcmcProposalPhaseSigma, mcmcTemperature,
                                   magnitudeStartMean, magnitudeStartSigma, photoProduction);
    return partFile;
  }, workerIds);

  // Merge the files into final one containing all minimisations
  TFileMerger merger(kTRUE, kFALSE);
  merger.OutputFile(out.c_str(), "RECREATE");
  for (const auto& partFile : partFiles) merger.AddFile(partFile.c_str());
  if (!merger.Merge()) throw std::runtime_error(std::string("Merge failed for output file: ") + outFile);
  for (const auto& partFile : partFiles) gSystem->Unlink(partFile.c_str());
  std::cout << "Merged " << partFiles.size() << " worker files into " << outFile << std::endl;
}

// EXTENSIONS / FUTURE WORK:
// 1) Add if statements for protection against moments that arent included i.e. circ and lin
// For example in GlueX we dont have circ so currently code would set to 0
// Need a flag so that we can skip
// 2) Add setters for all the flags and do a general tidyup on the user input section
// User should just have to specify commands like Set this and Set that with cml for file I/O
// 3) Possibly restructure into different headers for each of the different classes as is standard
// With each class in a different header keep this mean macro for user to only have to set things
// 4) Update verbose to just print once instead of ncores times, simple flag done before