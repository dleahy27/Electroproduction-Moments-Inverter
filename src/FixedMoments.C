#include "emi/Runner.h"

#include "Detail.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>


namespace emi {
namespace {
using namespace detail;

static void FillUserAmplitudes(const std::shared_ptr<EvaluationContext>& ctx,
                               std::vector<double>& fullVals) {
  std::vector<double> freeValues(ctx->freeToFull.size(), 0.0);
  for (size_t i = 0; i < ctx->freeToFull.size(); ++i) {
    const auto& parameter = ctx->fullPars[static_cast<size_t>(ctx->freeToFull[i])];
    freeValues[i] = parameter.phase
        ? -2.4 + 4.8 * static_cast<double>((i * 7) % 19) / 18.0
        : 0.025 + 0.005 * static_cast<double>((i * 5) % 7);
  }
  if (!FillFullParameters(*ctx, freeValues.data(), fullVals)) {
    throw std::runtime_error("Could not normalise the generated amplitude point");
  }
}

static void MakeBranchesForNamedValues(TTree* t,
                                       const std::vector<std::string>& names,
                                       std::vector<double>& storage) {
  storage.assign(names.size(), 0.0);
  for (size_t i = 0; i < names.size(); ++i) {
    t->Branch(names[i].c_str(), &storage[i]);
  }
}

static void PrintNonZeroParameters(const std::shared_ptr<EvaluationContext>& ctx,
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

static void PrintObservedLikeMoments(const std::shared_ptr<EvaluationContext>& ctx,
                                     const std::vector<double>& obsVals) {
  std::cout << "Observed Moments written for the inverter";
  std::cout << "--------------------------------"<<std::endl;
  for (size_t i = 0; i < ctx->observed.size(); ++i) {
    std::cout << std::setw(12) << ctx->observed[i].name << " = " << std::setw(14) << std::setprecision(10) << obsVals[i] << std::endl;
  }
}

static void AddSeedMomentBranches(std::vector<std::string>& names,
                                  const std::string& valueName,
                                  const std::string& errName) {
  names.push_back(valueName);
  names.push_back(errName);
}

static std::vector<std::string> RequestedMomentSeedBranches(const InternalConfig& cfg) {
  std::vector<std::string> names;
  int maximumL = 0;
  for (const auto& wave : cfg.waves) maximumL = std::max(maximumL, wave.l);

  auto add = [&](const std::string& valueName) {
    AddSeedMomentBranches(names, valueName, valueName + "_err");
  };
  for (int L = 0; L <= 2 * maximumL; ++L) {
    for (int M = 0; M <= L; ++M) {
      add((cfg.photoproduction ? "RH_0_" : "RH04_") +
          std::to_string(L) + "_" + std::to_string(M));
    }
  }
  const int maximumAlpha = cfg.photoproduction ? 3 : 8;
  for (int alpha = 1; alpha <= maximumAlpha; ++alpha) {
    if (alpha == 4) continue;
    for (int L = 0; L <= 2 * maximumL; ++L) {
      for (int M = 0; M <= L; ++M) {
        if ((alpha == 2 || alpha == 3 || alpha == 6 || alpha == 7) && M == 0) continue;
        add("RH_" + std::to_string(alpha) + "_" +
            std::to_string(L) + "_" + std::to_string(M));
      }
    }
  }
  return names;
}

static void CreateSeedInputFile(const std::string& fileName,
                                const std::string& treeName,
                                const InternalConfig& cfg) {
  std::unique_ptr<TFile> f(TFile::Open(fileName.c_str(), "RECREATE"));
  if (!f || f->IsZombie()) throw std::runtime_error("Failed to create seed file: " + fileName);

  TTree t(treeName.c_str(), "Temporary tree used to initialise the moment model");

  const std::vector<std::string> branchNames = RequestedMomentSeedBranches(cfg);
  if (branchNames.empty()) throw std::runtime_error("No seed moment branch names were requested");

  std::vector<double> storage(branchNames.size(), 0.0);
  for (size_t i = 0; i < branchNames.size(); ++i) {
    if (branchNames[i].size() >= 4 &&
        branchNames[i].compare(branchNames[i].size() - 4, 4, "_err") == 0) {
      storage[i] = 1;
    }
    t.Branch(branchNames[i].c_str(), &storage[i]);
  }

  t.Fill();
  t.Write();
  f->Close();
}

} // namespace

void GenerateFixedMoments(const std::filesystem::path& output,
                          double epsilon,
                          bool printToScreen,
                          const ModelConfig& model) {
  const std::filesystem::path outPath = output.empty()
      ? std::filesystem::path("InputFiles/Generated/fixed_test.root")
      : output;
  if (!outPath.parent_path().empty()) {
    std::filesystem::create_directories(outPath.parent_path());
  }
  const std::string treeName = "genMoments";

  FitConfig fit;
  fit.input = outPath;
  fit.tree = treeName;
  fit.epsilon = epsilon;
  fit.photoproduction = false;
  fit.verbose = false;
  InternalConfig cfg = MakeInternalConfig(fit, model);

  CreateSeedInputFile(outPath.string(), treeName, cfg);
  auto ctx = BuildContext(cfg);

  std::vector<double> fullVals;
  FillUserAmplitudes(ctx, fullVals);

  std::vector<double> Hvals(ctx->modelsRec.size(), 0.0);
  EvaluateAllMoments(*ctx, fullVals, Hvals);

  std::vector<std::string> observedNames;
  std::vector<std::string> errNames;
  std::vector<double> observedVals(ctx->observed.size(), 0.0);
  std::vector<double> errVals(ctx->observed.size(), 1);
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
    throw std::runtime_error("Failed to open output file: " + outPath.string());
  }

  TTree* t = new TTree(treeName.c_str(), "Generated moments built from a fixed amplitude point");
  t->SetDirectory(fout.get());

  std::vector<double> parStorage;
  MakeParameterBranches(t, ctx->fullPars, parStorage);

  std::vector<double> obsStorage;
  MakeBranchesForNamedValues(t, observedNames, obsStorage);

  std::vector<double> errStorage;
  MakeBranchesForNamedValues(t, errNames, errStorage);

  std::vector<std::string> q2Name = {"Q2"};
  std::vector<double> q2Storage;
  MakeBranchesForNamedValues(t, q2Name, q2Storage);

  parStorage = fullVals;
  obsStorage = observedVals;
  errStorage = errVals;
  q2Storage.assign(q2Storage.size(), 1.0);

  t->Fill();
  fout->Write();
  fout->Close();

  if (printToScreen) {
    PrintNonZeroParameters(ctx, fullVals);
    PrintObservedLikeMoments(ctx, observedVals);
    std::cout << "Saved generated moments to " << outPath << std::endl;
  }
}

} // namespace emi
