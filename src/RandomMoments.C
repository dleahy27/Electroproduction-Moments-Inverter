#include "emi/Runner.h"

#include "Detail.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

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
                                 std::vector<double>& fullVals,
                                 TRandom3& random) {
    std::vector<double> freeValues;
    BuildRandomStart(*ctx, random, freeValues);
    if (!FillFullParameters(*ctx, freeValues.data(), fullVals)) {
      throw std::runtime_error("Could not normalise a generated amplitude point");
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

} // namespace

void GenerateRandomMoments(const std::filesystem::path& output,
                           unsigned nEvents,
                           double epsilon,
                           std::uint32_t seed,
                           const ModelConfig& model) {
  const std::filesystem::path outPath = output.empty()
      ? std::filesystem::path("InputFiles/Generated/random_input_moments.root")
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

  GenerateFixedMoments(outPath, epsilon, false, model);
  auto ctx = BuildContext(cfg);

  std::unique_ptr<TFile> fout(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error("Failed to open output file: " + outPath.string());
  }

  TTree* t = new TTree(treeName.c_str(), "Generated moments built from random amplitudes");
  t->SetDirectory(fout.get());
  int i = 0;
  TRandom3 random(seed);
  if (seed == 0) random.SetSeed(0);

  t->Branch("Event", &i, "Event/I");
  std::vector<double> parStorage;
  std::vector<double> obsStorage;
  std::vector<double> errStorage;

  const unsigned updateEvery = std::max(1u, nEvents / 100);

  for (i = 0; i < static_cast<int>(nEvents); ++i) {
    std::vector<double> fullVals;
    FillUserAmplitudes(ctx, fullVals, random);

    std::vector<double> Hvals(ctx->modelsRec.size(), 0.0);
    EvaluateAllMoments(*ctx, fullVals, Hvals);

    std::vector<std::string> observedNames;
    std::vector<std::string> errNames;
    std::vector<double> observedVals(ctx->observed.size(), 0.0);
    observedNames.reserve(ctx->observed.size());

    for (size_t j = 0; j < ctx->observed.size(); ++j) {
      const auto& ob = ctx->observed[j];
      observedNames.push_back(ob.name);
      errNames.push_back(ob.name + "_err");

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

    if (i == 0) {
      MakeParameterBranches(t, ctx->fullPars, parStorage);
      MakeBranchesForNamedValues(t, observedNames, obsStorage);
      MakeBranchesForNamedValues(t, errNames, errStorage);
      std::fill(errStorage.begin(), errStorage.end(), 1.0);
    }

    parStorage = fullVals;
    obsStorage = observedVals;
    t->Fill();
    if (static_cast<unsigned>(i) % updateEvery == 0 ||
        static_cast<unsigned>(i) + 1 == nEvents) {
      PrintProgress(i + 1, nEvents);
    }
  }

  fout->Write();
  fout->Close();

  std::cout << "Saved generated moments to " << outPath << std::endl;
}

} // namespace emi
