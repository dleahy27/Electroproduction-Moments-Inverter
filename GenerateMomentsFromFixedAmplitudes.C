#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "RunGivenMoments_Chi2Amps.C"

namespace fixed_moment_builder {
using namespace chi2_amp_fit_opt;

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
  // Fixed-amplitude generation is used as an inverter test, so it is allowed
  // to set amplitudes that are fixed or excluded in a particular fit model.
  const int idx = FindParIndex(ctx, name);
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
  // This generator deliberately allows arbitrary amplitude points so that
  // the inverter can be tested independently of physical restrictions.
  // =====================================================================

  // Positive-reflectivity magnitudes
  set("a_T_0_0", 0.278093);
  set("a_L_0_0", 0.012959);

  set("a_T_1_1", 0.010878);
  set("a_T_1_0", 0.119613);
  set("a_T_1_m1", 0.097077);

  set("a_L_1_1", 0.095088);
  set("a_L_1_0", 0.219784);
  set("a_L_1_m1", 0.011541);

  set("a_T_2_2", 0.320298);
  set("a_T_2_1", 0.294303);
  set("a_T_2_0", 0.388018);
  set("a_T_2_m1", 0.037811);
  set("a_T_2_m2", 0.183498);

  set("a_L_2_2", 0.086477);
  set("a_L_2_1", 0.282641);
  set("a_L_2_0", 0.237000);
  set("a_L_2_m1", 0.095872);
  set("a_L_2_m2", 0.256278);

  // Negative-reflectivity magnitudes
  set("b_T_0_0", 0.073704);
  set("b_L_0_0", 0.129168);

  set("b_T_1_1", 0.352029);
  set("b_T_1_0", 0.002826);
  set("b_T_1_m1", 0.350459);

  set("b_L_1_1", 0.040336);
  set("b_L_1_0", 0.042063);
  set("b_L_1_m1", 0.368583);

  set("b_T_2_2", 0.303627);
  set("b_T_2_1", 0.147978);
  set("b_T_2_0", 0.067619);
  set("b_T_2_m1", 0.416301);
  set("b_T_2_m2", 0.146388);

  set("b_L_2_2", 0.262566);
  set("b_L_2_1", 0.351028);
  set("b_L_2_0", 0.317368);
  set("b_L_2_m1", 0.233211);
  set("b_L_2_m2", 0.423218);

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
  std::cout << "Input parameters used for synthetic moments" <<std::endl;
  std::cout << "-----------------------------------------"<<std::endl;
  for (size_t i = 0; i < ctx->fullPars.size(); ++i) {
    const auto& p = ctx->fullPars[i];
    const double v = fullVals[i];
    if (std::abs(v) < 1e-14) continue;
    std::cout << std::setw(18) << p.name << " = " << std::setw(12)<< std::setprecision(8) << v << std::endl;
    if (p.fixed) std::cout << "   [fixed in fit]" << std::endl;
  }
}

static void PrintObservedLikeMoments(const std::shared_ptr<EvalContext>& ctx,
                                     const std::vector<double>& obsVals) {
  std::cout << "Observed Moments written for the inverter";
  std::cout << "--------------------------------"<<std::endl;
  for (size_t i = 0; i < ctx->observed.size(); ++i) {
    std::cout << std::setw(12) << ctx->observed[i].name << " = " << std::setw(14) << std::setprecision(10) << obsVals[i] << std::endl;
  }
}

static void PrintAllMoments(const std::shared_ptr<EvalContext>& ctx,
                                     const std::vector<double>& Hvals) {
  std::cout << "All Moments written for the inverter";
  std::cout << "--------------------------------"<<std::endl;
  for (size_t i = 0; i < ctx->modelsRec.size(); ++i) {
    std::cout << std::setw(12) << ctx->modelsRec[i].name << " = " << std::setw(14) << std::setprecision(10) << Hvals[i] << std::endl;
  }
}

static void CreateSeedInputFile(const std::string& fileName, const std::string& treeName) {
  std::unique_ptr<TFile> f(TFile::Open(fileName.c_str(), "RECREATE"));
  if (!f || f->IsZombie()) throw std::runtime_error("Failed to create seed file: " + fileName);
  TTree t(treeName.c_str(), "Temporary tree used to initialise the moment model");
  double dummy = 0.0;
  t.Branch("dummy", &dummy);
  t.Fill();
  t.Write();
  f->Close();
}

} // namespace fixed_moment_builder

void GenerateMomentsFromFixedAmplitudes(std::string outFile = "fixed_input_moments.root",
                                        std::vector<double> Q2vals = {1.0},
                                        double epsilon = 1.0,
                                        bool printToScreen = true) {
  using namespace fixed_moment_builder;
  using namespace chi2_amp_fit_opt;

  const std::string outDir = "./InputFiles/Generated/";
  gSystem->mkdir(outDir.c_str(), kTRUE);
  const std::string outPath = outDir + outFile;
  const std::string treeName = "genMoments";

  FitConfig cfg;
  cfg.lmax = 2;
  cfg.mmax = 2;
  cfg.negm = true;
  cfg.useNegRef = true;
  cfg.epsilon = epsilon;
  cfg.photoProduction = false;
  cfg.enforceLongitudinalParity = false;
  cfg.momentsFile = outPath;
  cfg.momentsTree = treeName;
  cfg.bin = 0;
  cfg.verbose = false;

  // BuildContext expects an input tree because the fit code maps the observed
  // branch list during context construction.  The file is overwritten below
  // after the synthetic moments have been evaluated.
  CreateSeedInputFile(outPath, treeName);
  auto ctx = BuildContext(cfg);

  std::vector<double> fullVals;
  FillUserAmplitudes(ctx, fullVals);

  std::vector<double> Hvals(ctx->modelsRec.size(), 0.0);
  EvalAllMoments(*ctx, fullVals, Hvals);

  std::vector<std::string> observedNames;
  std::vector<std::string> errNames;
  std::vector<double> observedVals(ctx->observed.size(), 0.0);
  std::vector<double> errVals(ctx->observed.size(), 0.01);
  observedNames.reserve(ctx->observed.size());
  errNames.reserve(ctx->observed.size());

  for (size_t i = 0; i < ctx->observed.size(); ++i) {
    const auto& ob = ctx->observed[i];
    observedNames.push_back(ob.name);
    errNames.push_back(ob.name + "_err");

    if (ob.isMixed04) {
      const int idx0 = ctx->observedModelIdx0[i];
      const int idx4 = ctx->observedModelIdx4[i];
      if (idx0 < 0 || idx4 < 0) throw std::runtime_error("RH04 moment is not mapped to H0/H4");
      observedVals[i] = Hvals[static_cast<size_t>(idx0)]
                      + cfg.epsilon * Hvals[static_cast<size_t>(idx4)];
    } else {
      const int midx = ctx->observedModelIdx[i];
      if (midx < 0) throw std::runtime_error("Observed moment is not mapped to a model moment");
      observedVals[i] = Hvals[static_cast<size_t>(midx)];
    }
  }

  std::unique_ptr<TFile> fout(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error("Failed to open output file: " + outPath);
  }

  TTree* t = new TTree(treeName.c_str(), "Generated moments built from a fixed amplitude point");
  t->SetDirectory(fout.get());

  std::vector<double> parStorage;
  MakeBranchesForPars(t, ctx->fullPars, parStorage);

  std::vector<double> rawMomentStorage;
  MakeBranchesForMoments(t, ctx->modelsRec, rawMomentStorage);

  std::vector<double> obsStorage;
  MakeBranchesForNamedValues(t, observedNames, obsStorage);

  std::vector<double> errStorage;
  MakeBranchesForNamedValues(t, errNames, errStorage);

  std::vector<std::string> q2Name = {"Q2"};
  std::vector<double> q2Storage;
  MakeBranchesForNamedValues(t, q2Name, q2Storage);

  parStorage = fullVals;
  rawMomentStorage = Hvals;
  obsStorage = observedVals;
  errStorage = errVals;
  q2Storage.assign(q2Storage.size(), 0.0);
  for (size_t i = 0; i < q2Storage.size() && i < Q2vals.size(); ++i) q2Storage[i] = Q2vals[i];

  t->Fill();
  fout->Write();
  fout->Close();

  if (printToScreen) {
    PrintNonZeroParameters(ctx, fullVals);
    PrintAllMoments(ctx, Hvals);
    PrintObservedLikeMoments(ctx, observedVals);
    std::cout << "Saved generated moments to " << outPath << std::endl;
  }
}
