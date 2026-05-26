#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "RunGivenMoments_Chi2Amps.C"

// Need to make function to write up branches for H stuff
namespace fixed_moment_builder {
using namespace chi2_amp_fit_opt;

struct NamedLMValue {
  std::string name;
  int L = 0;
  int M = 0;
  double value = 0.0;
};

static int FindParIndex(const std::shared_ptr<EvalContext>& ctx, const std::string& name) {
  for (size_t i = 0; i < ctx->fullPars.size(); ++i) {
    if (ctx->fullPars[i].name == name) return static_cast<int>(i);
  }
  throw std::runtime_error("Unknown parameter name: " + name);
}

static void SetPar(const std::shared_ptr<EvalContext>& ctx,
                   std::vector<double>& fullVals,
                   const std::string& name,
                   double value) {
  const int idx = FindParIndex(ctx, name);
  const auto& p = ctx->fullPars[static_cast<size_t>(idx)];

  if (p.fixed && std::abs(value - p.init) > 1e-12) {
    throw std::runtime_error(
        "Attempted to change fixed parameter '" + name +
        "'. The fitter fixes it to " + std::to_string(p.init) + ".");
  }
  if (value < p.low - 1e-12 || value > p.high + 1e-12) {
    throw std::runtime_error(
        "Value " + std::to_string(value) + " for parameter '" + name +
        "' is outside allowed range [" + std::to_string(p.low) + ", " +
        std::to_string(p.high) + "].");
  }

  fullVals[static_cast<size_t>(idx)] = value;
}

static void FillUserAmplitudes(const std::shared_ptr<EvalContext>& ctx,
                               std::vector<double>& fullVals) {
  fullVals.assign(ctx->fullPars.size(), 0.0);
  for (size_t i = 0; i < ctx->fullPars.size(); ++i) {
    fullVals[i] = ctx->fullPars[i].init;
  }

  auto set = [&](const char* name, double value) {
    SetPar(ctx, fullVals, name, value);
  };

  // =====================================================================
  // Replace these with any amplitude/phase point you want to test.
  // Keep the fixed parameters fixed: SetPar() will throw if you change one.
  // =====================================================================

  // Positive-reflectivity magnitudes
  set("a_T_0_0", 0.639427);
  set("a_L_0_0", 0.029797);

  set("a_T_1_1", 0.025011);
  set("a_T_1_0", 0.275029);
  set("a_T_1_m1", 0.223211);

  set("a_L_1_1", 0.218638);
  set("a_L_1_0", 0.505355);
  set("a_L_1_m1", 0.026536);

  set("a_T_2_2", 0.736471);
  set("a_T_2_1", 0.676699);
  set("a_T_2_0", 0.892180);
  set("a_T_2_m1", 0.086939);
  set("a_T_2_m2", 0.421922);

  set("a_L_2_2", 0.198838);
  set("a_L_2_1", 0.649884);
  set("a_L_2_0", 0.544941);
  set("a_L_2_m1", 0.220441);
  set("a_L_2_m2", 0.589266);

// Negative-reflectivity magnitudes
  set("b_T_0_0", 0.16947);
  set("b_L_0_0", 0.297);

  set("b_T_1_1", 0.809430);
  set("b_T_1_0", 0.006499);
  set("b_T_1_m1", 0.805819);

  set("b_L_1_1", 0.092746);
  set("b_L_1_0", 0.096716);
  set("b_L_1_m1", 0.847494);

  set("b_T_2_2", 0.698139);
  set("b_T_2_1", 0.340251);
  set("b_T_2_0", 0.155479);
  set("b_T_2_m1", 0.957213);
  set("b_T_2_m2", 0.336595);

  set("b_L_2_2", 0.603726);
  set("b_L_2_1", 0.807128);
  set("b_L_2_0", 0.729732);
  set("b_L_2_m1", 0.536228);
  set("b_L_2_m2", 0.973116);

// Positive-reflectivity phases
  set("aphi_T_0_0", -0.763191);
  set("aphi_L_0_0", -1.323314);

  set("aphi_T_1_1", 0.326981);
  set("aphi_T_1_0", 2.069711);
  set("aphi_T_1_m1", 0.744682);

  set("aphi_L_1_1", -2.640245);
  set("aphi_L_1_0", -1.678924);
  set("aphi_L_1_m1", -2.506982);

  set("aphi_T_2_2", 0.0);
  set("aphi_T_2_1", 0.486018);
  set("aphi_T_2_0", 1.285363);
  set("aphi_T_2_m1", -2.853670);
  set("aphi_T_2_m2", -1.709666);

  set("aphi_L_2_2", -1.395033);
  set("aphi_L_2_1", 0.852531);
  set("aphi_L_2_0", -0.849284);
  set("aphi_L_2_m1", -0.815677);
  set("aphi_L_2_m2", -1.825221);

// Negative-reflectivity phases
  set("bphi_T_0_0", 1.2234);
  set("bphi_L_0_0", -2.3345);

  set("bphi_T_1_1", -1.464122);
  set("bphi_T_1_0", 2.743582);
  set("bphi_T_1_m1", 0.930134);

  set("bphi_L_1_1", 3.075766);
  set("bphi_L_1_0", 0.879644);
  set("bphi_L_1_m1", 0.357826);

  set("bphi_T_2_2", 0.0);
  set("bphi_T_2_1", -2.066297);
  set("bphi_T_2_0", 1.439646);
  set("bphi_T_2_m1", -2.114905);
  set("bphi_T_2_m2", -0.757404);

  set("bphi_L_2_2", 1.159966);
  set("bphi_L_2_1", 2.154202);
  set("bphi_L_2_0", 1.734159);
  set("bphi_L_2_m1", -1.702441);
  set("bphi_L_2_m2", -2.939901);
}

  static void NormalizeAmplitudesToUnitSum(const std::shared_ptr<EvalContext>& ctx, std::vector<double>& fullVals) {
  auto pdefs = ctx->fullPars;
  double sumsq = 0.0;
  std::vector<size_t> magIdx;
  magIdx.reserve(pdefs.size());

  for (size_t i = 0; i < pdefs.size(); ++i) {
    const std::string& n = pdefs[i].name;

    // keep only amplitude magnitudes
    if (n.rfind("aphi_", 0) == 0) continue;  // phase parameter
    if (n.rfind("bphi_", 0) == 0) continue;  // phase parameter
    if (n == "R") continue;                  // not an amplitude
    if (n == "epsR4") continue;              // just in case
    if (n.empty()) continue;  // only amplitude params

    magIdx.push_back(i);
    sumsq += fullVals[i] * fullVals[i];
  }

  if (sumsq <= 0.0) {
    throw std::runtime_error("NormalizeAmplitudesToUnitSum: amplitude norm is zero");
  }

  const double scale = 1.0 / std::sqrt(sumsq);
  for (size_t i : magIdx) fullVals[i] *= scale;

  sumsq = 0;
  for (size_t i = 0; i < pdefs.size(); ++i) {
    const std::string& n = pdefs[i].name;

    // keep only amplitude magnitudes
    if (n.rfind("aphi_", 0) == 0) continue;  // phase parameter
    if (n == "R") continue;                  // not an amplitude
    if (n == "epsR4") continue;              // just in case
    if (n.empty() || n[0] != 'a') continue;  // only amplitude params

    magIdx.push_back(i);
    sumsq += fullVals[i] * fullVals[i];
  }
  std::cout<<sumsq<<std::endl<<std::endl<<std::endl;
}

static int FindModelIndex(const std::shared_ptr<EvalContext>& ctx, int alpha, int L, int M) {
  for (size_t i = 0; i < ctx->modelsRec.size(); ++i) {
    const auto& mm = ctx->modelsRec[i];
    if (mm.alpha == alpha && mm.L == L && mm.M == M) return static_cast<int>(i);
  }
  return -1;
}

static void MakeBranchesForNamedValues(TTree* t,
                                       const std::vector<std::string>& names,
                                       std::vector<double>& storage) {
  storage.assign(names.size(), 0.0);
  for (size_t i = 0; i < names.size(); ++i) {
    t->Branch(names[i].c_str(), &storage[i]);
  }
}

static void PrintNonZeroParameters(const std::shared_ptr<EvalContext>& ctx,
                                   const std::vector<double>& fullVals) {
  std::cout << "\nInput parameters used for synthetic moments\n";
  std::cout << "-----------------------------------------\n";
  for (size_t i = 0; i < ctx->fullPars.size(); ++i) {
    const auto& p = ctx->fullPars[i];
    const double v = fullVals[i];
    if (std::abs(v) < 1e-14) continue;
    std::cout << std::setw(18) << p.name << " = " << std::setw(12) << std::setprecision(8) << v;
    if (p.fixed) std::cout << "   [fixed]";
    std::cout << '\n';
  }
}

static void PrintObservedLikeMoments(const std::shared_ptr<EvalContext>& ctx,
                                     const std::vector<double>& obsVals) {
  std::cout << "\nObserved-like moments from this amplitude point\n";
  std::cout << "----------------------------------------------\n";
  for (size_t i = 0; i < ctx->observed.size(); ++i) {
    std::cout << std::setw(12) << ctx->observed[i].name
              << " = " << std::setw(14) << std::setprecision(10) << obsVals[i] << '\n';
  }
}

} // namespace fixed_moment_builder

void GenerateMomentsFromFixedAmplitudes(std::string outFile = "fixed_input_moments.root",
                                        std::vector<double> Q2vals  = {1.0},
                                        double epsR4 = 1.0,
                                        bool printToScreen = true) {
  using namespace fixed_moment_builder;
  using namespace chi2_amp_fit_opt;

  std::string outDir = "./InputFiles/Generated/";
  auto tmp = outDir + outFile;
  auto out = tmp.c_str();
  std::vector<std::string> Q2name = {"Q2"};
  FitConfig cfg;
  cfg.epsR4 = epsR4;
  cfg.momentsFile = tmp;
  cfg.momentsTree = "genMoments";
  cfg.bin = 0;


  auto ctx = BuildContext(cfg);

  std::vector<double> fullVals;
  FillUserAmplitudes(ctx, fullVals);
  NormalizeAmplitudesToUnitSum(ctx, fullVals);

  // Pure reconstructed moments H^alpha_{LM}
  std::vector<double> Hvals(ctx->modelsRec.size(), 0.0);
  std::vector<double> pairSin, pairCos;

  BuildPhasePairTrigCache(*ctx, fullVals, pairSin, pairCos);

  for (size_t i = 0; i < ctx->modelsRec.size(); ++i) {
    Hvals[i] = EvalMomentOnly(ctx->modelsRec[i], fullVals, pairSin, pairCos);
  }

  // Standard scaled RH^alpha_{LM}
  std::vector<double> RHvals(ctx->modelsRec.size(), 0.0);
  std::vector<double> RHErrs(ctx->modelsRec.size(), 0.0);
  for (size_t i = 0; i < ctx->modelsRec.size(); ++i) {
    RHvals[i] = Ralpha(ctx->modelsRec[i].alpha) * Hvals[i];
    RHErrs[i] = 0.01;

  }

  // Full RH04_{LM} mixed set and the exact observed subset used by the fitter.
  std::vector<NamedLMValue> RH04full;
  std::vector<std::string> RH04names;
  std::vector<double> RH04vals;

  // Skip LM = 00 as we will never have access to it experimentally
  for (int L = 1; L <= 2 * cfg.lmax; ++L) {
    for (int M = 0; M <= L; ++M) {
      const int idx0 = FindModelIndex(ctx, 0, L, M);
      const int idx4 = FindModelIndex(ctx, 4, L, M);
      if (idx0 < 0 || idx4 < 0) continue;

      NamedLMValue x;
      x.L = L;
      x.M = M;
      x.name = std::string("RH04_") + std::to_string(L) + "_" + std::to_string(M);
      x.value = (Hvals[static_cast<size_t>(idx0)] - (cfg.epsR4 * R) * Hvals[static_cast<size_t>(idx4)]) * fac;
      RH04full.push_back(x);
      RH04names.push_back(x.name);
    }
  }


  // Assign RH04
  RH04vals.assign(RH04full.size(), 0.0);
  for (size_t i = 0; i < RH04full.size(); ++i) RH04vals[i] = RH04full[i].value;

  std::vector<std::string> observedNames;
  std::vector<double> observedVals(ctx->observed.size(), 0.0);
  observedNames.reserve(ctx->observed.size());

  std::vector<std::string> errNames;
  std::vector<double> errVals(ctx->observed.size(), 0.01);
  errNames.reserve(ctx->observed.size());
  for (size_t i = 0; i < ctx->observed.size(); ++i) {
    const auto& ob = ctx->observed[i];
    observedNames.push_back(ob.name);
    errNames.push_back(ob.name + "_err");
    if (ob.isMixed04) {
      const int idx0 = (i < ctx->observedModelIdx0.size()) ? ctx->observedModelIdx0[i] : -1;
      const int idx4 = (i < ctx->observedModelIdx4.size()) ? ctx->observedModelIdx4[i] : -1;
      if (idx0 >= 0 && idx4 >= 0) {
        observedVals[i] = (Hvals[static_cast<size_t>(idx0)] - (cfg.epsR4 * R) * Hvals[static_cast<size_t>(idx4)]) * fac;
      }
    } else {
      const int midx = ctx->observedModelIdx[i];
      if (midx >= 0) observedVals[i] = RHvals[static_cast<size_t>(midx)];
    }
  }


  std::unique_ptr<TFile> fout(TFile::Open(out, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error(std::string("Failed to open output file: ") + outFile);
  }

  TTree* t = new TTree("genMoments", "Generated moments built from a fixed amplitude point");
  t->SetDirectory(fout.get());

  double storedR = R;
  t->Branch("R", &storedR);

  std::vector<double> parStorage;
  MakeBranchesForPars(t, ctx->fullPars, parStorage);
  std::vector<double> RH04storage;
  MakeBranchesForNamedValues(t, RH04names, RH04storage);
  std::vector<double> obsStorage;
  MakeBranchesForNamedValues(t, observedNames, obsStorage);
  std::vector<double> errStorage;
  MakeBranchesForNamedValues(t, errNames, errStorage);
  std::vector<double> Q2storage;
  MakeBranchesForNamedValues(t, Q2name, Q2storage);

  parStorage = fullVals;
  RH04storage = RH04vals;
  obsStorage = observedVals;
  Q2storage = Q2vals;
  errStorage = errVals;
  t->Fill();

  fout->Write();
  fout->Close();

  if (printToScreen) {
    PrintNonZeroParameters(ctx, fullVals);
    std::cout << "\nR      = " << std::setprecision(10) << R << '\n';
    PrintObservedLikeMoments(ctx, observedVals);
    std::cout << "\nSaved generated moments to " << outDir<<outFile << std::endl;
  }
}
