#include "emi/Runner.h"

#include "Detail.h"

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace emi {
namespace {
using namespace detail;

void FillRandomAmplitudes(const EvaluationContext& ctx,
                          std::vector<double>& fullVals,
                          TRandom3& random) {
  std::vector<double> freeValues;
  BuildRandomStart(ctx, random, freeValues);
  if (!FillFullParameters(ctx, freeValues.data(), fullVals)) {
    throw std::runtime_error("Could not normalise a generated amplitude point");
  }
}

void PrintProgress(unsigned completed, unsigned total) {
  constexpr int barWidth = 50;

  const double progress = static_cast<double>(completed) / total;
  const int position = static_cast<int>(barWidth * progress);

  std::cout << "\r[";
  for (int j = 0; j < barWidth; ++j) {
    std::cout << (j < position ? '=' : j == position ? '>' : ' ');
  }
  std::cout << "] " << std::setw(3)
            << static_cast<int>(progress * 100.0) << '%' << std::flush;
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

  auto ctx = BuildSyntheticContext(cfg);

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
  std::vector<double> parStorage, obsStorage, errStorage;
  MakeParameterBranches(t, ctx->fullPars, parStorage);

  std::vector<std::string> observedNames, errorNames;
  observedNames.reserve(ctx->observed.size());
  errorNames.reserve(ctx->observed.size());
  for (const auto& observed : ctx->observed) {
    observedNames.push_back(observed.name);
    errorNames.push_back(observed.name + "_err");
  }
  MakeNamedBranches(t, observedNames, obsStorage);
  MakeNamedBranches(t, errorNames, errStorage);
  std::fill(errStorage.begin(), errStorage.end(), 1.0);

  const unsigned updateEvery = std::max(1u, nEvents / 100);
  std::vector<double> fullVals, rawMoments, observedVals;

  for (i = 0; i < static_cast<int>(nEvents); ++i) {
    FillRandomAmplitudes(*ctx, fullVals, random);
    EvaluateAllMoments(*ctx, fullVals, rawMoments);
    EvaluateObservedMoments(*ctx, rawMoments, observedVals);

    parStorage = fullVals;
    obsStorage = observedVals;
    t->Fill();
    if (static_cast<unsigned>(i) % updateEvery == 0 ||
        static_cast<unsigned>(i) + 1 == nEvents) {
      PrintProgress(static_cast<unsigned>(i) + 1, nEvents);
    }
  }

  fout->Write();
  fout->Close();

  if (nEvents > 0) std::cout << '\n';
  std::cout << "Saved generated moments to " << outPath << '\n';
}

} // namespace emi
