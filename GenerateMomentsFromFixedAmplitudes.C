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
  set("a_T_1_1",  0.8);
  set("a_T_1_0",  0.4);
  set("a_T_1_m1", 0.5);
  set("a_L_1_1",  0.2);
  set("a_L_1_m1", 0.25);

  // Negative-reflectivity magnitudes
  set("b_T_1_1",  0.1);
  set("b_T_1_0",  0.6);
  set("b_T_1_m1", 0.18);
  set("b_L_1_1",  0.28);
  set("b_L_1_0",  0.26);
  set("b_L_1_m1", 0.08);

  // Positive-reflectivity phases
  set("aphi_T_1_1",  0.0);
  set("aphi_T_1_0",  -0.3);
  set("aphi_T_1_m1", 2.2);
  set("aphi_L_1_1",  -1.6);
  set("aphi_L_1_m1", -1.95);

  // Negative-reflectivity phases
  set("bphi_T_1_1",  0.0);
  set("bphi_T_1_0",  -0.6);
  set("bphi_T_1_m1", -0.95);
  set("bphi_L_1_1",  1.22);
  set("bphi_L_1_0",  -1.76);
  set("bphi_L_1_m1", 0.4);

  // A very targeted ambiguity test is to duplicate this file and flip only
  // the two suspect phases, e.g.
  //   aphi_T_1_0  -> -aphi_T_1_0
  //   aphi_T_1_m1 -> -aphi_T_1_m1
  // then compare the RH / RH04 branches event-by-event.
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


  // R = -H4_00 / H0_00.
  double H0_00 = 0.0;
  double H4_00 = 0.0;
  if (ctx->idxH0_00 >= 0) H0_00 = Hvals[static_cast<size_t>(ctx->idxH0_00)];
  if (ctx->idxH4_00 >= 0) H4_00 = Hvals[static_cast<size_t>(ctx->idxH4_00)];

  const double R = -H4_00/H0_00;
  const double sqrtR = std::sqrt(R);
  const double fac = 1.0 / (1 + cfg.epsR4 * R);
  auto Ralpha = [&](int alpha) -> double {
    if (alpha >= 0 && alpha <= 3) return 1.0 * fac;
    if (alpha == 4) return cfg.epsR4 * R * fac;
    if (alpha >= 5 && alpha <= 8) return sqrtR * fac;
    return 1.0;
  };

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
