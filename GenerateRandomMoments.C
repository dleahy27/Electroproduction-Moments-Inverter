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

struct AmpInfo {
  std::string magName;
  std::string phaseName;
  bool isLongitudinal;
  double mag;
  double phase;
};

namespace random_moment_builder{
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
                                 std::vector<double>& fullVals, double epsilon)
  {
    fullVals.assign(ctx->fullPars.size(), 0.0);
    for (size_t i = 0; i < ctx->fullPars.size(); ++i) {
      fullVals[i] = ctx->fullPars[i].init;
    }

    auto set = [&](std::string name, double value) {
      SetPar(ctx, fullVals, name, value);
    };

    // =====================================================================
    // Replace these with any amplitude/phase point you want to test.
    // This generator deliberately allows arbitrary amplitude points so that
    // the inverter can be tested independently of physical restrictions.
    // =====================================================================

    {
      const double pi = std::acos(-1.0);

      std::random_device rd;
      std::mt19937_64 rng(rd());

      std::uniform_real_distribution<double> magDist(0.0, 1.0);
      std::uniform_real_distribution<double> phaseDist(-pi, pi);

      std::vector<AmpInfo> amps = {
        // Positive-reflectivity amplitudes
        {"a_T_0_0",  "aphi_T_0_0",  false, 0.0, 0.0},
        {"a_L_0_0",  "aphi_L_0_0",  true,  0.0, 0.0},

        {"a_T_1_1",  "aphi_T_1_1",  false, 0.0, 0.0},
        {"a_T_1_0",  "aphi_T_1_0",  false, 0.0, 0.0},
        {"a_T_1_m1", "aphi_T_1_m1", false, 0.0, 0.0},

        {"a_L_1_1",  "aphi_L_1_1",  true, 0.0, 0.0},
        {"a_L_1_0",  "aphi_L_1_0",  true, 0.0, 0.0},

        {"a_T_2_2",  "aphi_T_2_2",  false, 0.0, 0.0},
        {"a_T_2_1",  "aphi_T_2_1",  false, 0.0, 0.0},
        {"a_T_2_0",  "aphi_T_2_0",  false, 0.0, 0.0},
        {"a_T_2_m1", "aphi_T_2_m1", false, 0.0, 0.0},
        {"a_T_2_m2", "aphi_T_2_m2", false, 0.0, 0.0},

        {"a_L_2_2",  "aphi_L_2_2",  true, 0.0, 0.0},
        {"a_L_2_1",  "aphi_L_2_1",  true, 0.0, 0.0},
        {"a_L_2_0",  "aphi_L_2_0",  true, 0.0, 0.0},

        // Negative-reflectivity amplitudes
        {"b_T_0_0",  "bphi_T_0_0",  false, 0.0, 0.0},

        {"b_T_1_1",  "bphi_T_1_1",  false, 0.0, 0.0},
        {"b_T_1_0",  "bphi_T_1_0",  false, 0.0, 0.0},
        {"b_T_1_m1", "bphi_T_1_m1", false, 0.0, 0.0},

        {"b_L_1_1",  "bphi_L_1_1",  true, 0.0, 0.0},

        {"b_T_2_2",  "bphi_T_2_2",  false, 0.0, 0.0},
        {"b_T_2_1",  "bphi_T_2_1",  false, 0.0, 0.0},
        {"b_T_2_0",  "bphi_T_2_0",  false, 0.0, 0.0},
        {"b_T_2_m1", "bphi_T_2_m1", false, 0.0, 0.0},
        {"b_T_2_m2", "bphi_T_2_m2", false, 0.0, 0.0},

        {"b_L_2_2",  "bphi_L_2_2",  true, 0.0, 0.0},
        {"b_L_2_1",  "bphi_L_2_1",  true, 0.0, 0.0},
      };

      double sumT = 0.0;
      double sumL = 0.0;

      // First sample raw magnitudes and phases
      for (auto& amp : amps) {
        amp.mag = magDist(rng);
        amp.phase = phaseDist(rng);

        if (amp.isLongitudinal) {
          sumL += amp.mag * amp.mag;
        } else {
          sumT += amp.mag * amp.mag;
        }
      }

      const double normDenom = sumT + epsilon * sumL;

      const double scale = std::sqrt(1.0 / normDenom);

      //std::cout<<"sumT = "<<sumT<<" and sumL = "<<sumL<<" gives scale = "<<scale<<std::endl;

      // Then set the correctly normalised magnitudes and random phases
      for (const auto& amp : amps) {
        set(amp.magName, amp.mag * scale);
        set(amp.phaseName, amp.phase);
      }

      // Reference Phases
      set("aphi_T_2_2", 0.0);
      set("bphi_T_2_2", 0.0);
    }
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

  void PrintProgress(int i, int nTotal)
  {
    const int barWidth = 50;

    double progress = double(i) / double(nTotal);
    int pos = int(barWidth * progress);

    std::cout << "\r[";
    for (int j = 0; j < barWidth; ++j) {
      if (j < pos) std::cout << "=";
      else if (j == pos) std::cout << ">";
      else std::cout << " ";
    }

    std::cout << "] "
              << std::setw(3) << int(progress * 100.0) << "%"
              << std::flush;
  }

}// namespace fixed_moment_builder

void GenerateRandomMoments(std::string outFile = "random_input_moments.root",
                                        int nEvents = 1E6,
                                        double epsilon = 1.0) {
  using namespace random_moment_builder;
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
  cfg.enforceLongitudinalParity = true;
  cfg.momentsFile = outPath;
  cfg.momentsTree = treeName;
  cfg.bin = 0;
  cfg.verbose = false;

  // BuildContext expects an input tree because the fit code maps the observed
  // branch list during context construction.  The file is overwritten below
  // after the synthetic moments have been evaluated.
  CreateSeedInputFile(outPath, treeName);
  auto ctx = BuildContext(cfg);

  std::unique_ptr<TFile> fout(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error("Failed to open output file: " + outPath);
  }

  TTree* t = new TTree(treeName.c_str(), "Generated moments built from a fixed amplitude point");
  t->SetDirectory(fout.get());
  int i = 0;

  t->Branch("Event", &i, "Event/I");
  std::vector<double> parStorage;
  std::vector<double> rawMomentStorage;
  std::vector<double> obsStorage;
  std::vector<double> errStorage;

  int updateEvery = std::max<int>(1, nEvents / 100);

  for(i; i<nEvents ; i++){
    std::vector<double> fullVals;
    FillUserAmplitudes(ctx, fullVals, cfg.epsilon);

    std::vector<double> Hvals(ctx->modelsRec.size(), 0.0);
    EvalAllMoments(*ctx, fullVals, Hvals);

    std::vector<std::string> observedNames;
    std::vector<std::string> errNames;
    std::vector<double> observedVals(ctx->observed.size(), 0.0);
    observedNames.reserve(ctx->observed.size());

    for (size_t j = 0; j < ctx->observed.size(); ++j) {
      const auto& ob = ctx->observed[j];
      observedNames.push_back(ob.name);

      if (ob.isMixed04) {
        const int idx0 = ctx->observedModelIdx0[j];
        const int idx4 = ctx->observedModelIdx4[j];
        if (idx0 < 0 || idx4 < 0) throw std::runtime_error("RH04 moment is not mapped to H0/H4");
        observedVals[j] = Hvals[static_cast<size_t>(idx0)]
                        + cfg.epsilon * Hvals[static_cast<size_t>(idx4)];
      } else {
        const int midx = ctx->observedModelIdx[j];
        if (midx < 0) throw std::runtime_error("Observed moment is not mapped to a model moment");
        observedVals[j] = Hvals[static_cast<size_t>(midx)];
      }
    }

    // Make the branches and assign the reference
    if (i == 0) {
      MakeBranchesForPars(t, ctx->fullPars, parStorage);
      MakeBranchesForNamedValues(t, observedNames, obsStorage);
    }

    parStorage = fullVals;
    rawMomentStorage = Hvals;
    obsStorage = observedVals;
    t->Fill();
    if (i % updateEvery == 0 || i == nEvents - 1) {
      PrintProgress(i + 1, nEvents);
    }
  }

  fout->Write();
  fout->Close();

  std::cout << "Saved generated moments to " << outPath << std::endl;
}

