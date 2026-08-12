#include "emi/Runner.h"

#include "Detail.h"

#include "ROOT/TProcessExecutor.hxx"
#include "TFile.h"
#include "TFileMerger.h"
#include "TRandom3.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <thread>

namespace emi {
namespace {

struct BestFit {
  bool valid = false;
  int status = -999;
  unsigned validStarts = 0;
  double chi2 = std::numeric_limits<double>::infinity();
  std::vector<double> parameters;
  std::vector<double> outputMoments;
};

void RebuildObservedMap(detail::EvaluationContext& context) {
  context.observedModelIdx.assign(context.observed.size(), -1);
  context.observedModelIdx0.assign(context.observed.size(), -1);
  context.observedModelIdx4.assign(context.observed.size(), -1);
  for (std::size_t i = 0; i < context.observed.size(); ++i) {
    const auto& observed = context.observed[i];
    if (observed.isMixed04) {
      const auto h0 = context.modelIndexByName.find(
          "H_0_" + std::to_string(observed.L) + '_' + std::to_string(observed.M));
      const auto h4 = context.modelIndexByName.find(
          "H_4_" + std::to_string(observed.L) + '_' + std::to_string(observed.M));
      if (h0 == context.modelIndexByName.end() || h4 == context.modelIndexByName.end()) {
        throw std::runtime_error("Could not map bootstrap moment " + observed.name);
      }
      context.observedModelIdx0[i] = static_cast<int>(h0->second);
      context.observedModelIdx4[i] = static_cast<int>(h4->second);
    } else {
      std::string name = observed.name;
      if (name.rfind("RH_", 0) == 0) name = "H_" + name.substr(3);
      const auto model = context.modelIndexByName.find(name);
      if (model == context.modelIndexByName.end()) {
        throw std::runtime_error("Could not map bootstrap moment " + observed.name);
      }
      context.observedModelIdx[i] = static_cast<int>(model->second);
    }
  }
}

std::shared_ptr<detail::EvaluationContext> SampleContext(
    const detail::EvaluationContext& nominal, TRandom3& random) {
  auto sampled = std::make_shared<detail::EvaluationContext>(nominal);
  for (auto& moment : sampled->observed) {
    if (std::isfinite(moment.value) && std::isfinite(moment.sigma) && moment.sigma > 0.0) {
      moment.value = random.Gaus(moment.value, std::abs(moment.sigma));
    }
  }
  RebuildObservedMap(*sampled);
  sampled->callCount = 0;
  return sampled;
}

BestFit FitToy(const detail::InternalConfig& config,
               const std::shared_ptr<detail::EvaluationContext>& context,
               TRandom3& random) {
  BestFit best;
  best.parameters.assign(context->fullPars.size(),
                         std::numeric_limits<double>::quiet_NaN());
  best.outputMoments.assign(context->outputMoments.size(),
                            std::numeric_limits<double>::quiet_NaN());

  for (unsigned startIndex = 0; startIndex < config.starts; ++startIndex) {
    std::vector<double> start;
    detail::BuildRandomStart(*context, random, start);
    auto result = detail::Minimize(context, start, config.runHesse);
    if (!result.valid) continue;
    ++best.validStarts;
    if (result.chi2 >= best.chi2) continue;

    std::vector<double> full;
    if (!detail::FillFullParameters(*context, result.freeValues.data(), full)) continue;
    best.valid = true;
    best.status = result.status;
    best.chi2 = result.chi2;
    best.parameters = std::move(full);
    std::vector<double> rawMoments;
    detail::EvaluateAllMoments(*context, best.parameters, rawMoments);
    detail::EvaluateOutputMoments(*context, rawMoments, best.outputMoments);
  }
  return best;
}

std::filesystem::path PartFile(const std::filesystem::path& output, unsigned worker) {
  return output.parent_path() /
         (output.stem().string() + ".part_" + std::to_string(worker) +
          output.extension().string());
}

void RunBootstrapWorker(const detail::InternalConfig& config, unsigned toys,
                        const std::filesystem::path& output) {
  auto nominal = detail::BuildContext(config);
  TRandom3 random(config.seed);
  if (config.seed == 0) random.SetSeed(0);
  std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "RECREATE"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Could not create output file " + output.string());
  }

  TTree tree("PartialWaves", "Best amplitude fit for each bootstrap sample");
  int toy = 0;
  int valid = 0;
  int status = -999;
  int ndf = static_cast<int>(nominal->observed.size()) -
            static_cast<int>(nominal->freeToFull.size());
  unsigned validStarts = 0;
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  double value = std::numeric_limits<double>::quiet_NaN();
  double logValue = std::numeric_limits<double>::quiet_NaN();
  double log10Chi2 = std::numeric_limits<double>::quiet_NaN();
  tree.Branch("toy", &toy);
  tree.Branch("valid", &valid);
  tree.Branch("status", &status);
  tree.Branch("ndf", &ndf);
  tree.Branch("valid_starts", &validStarts);
  tree.Branch("chi2", &chi2);
  tree.Branch("val", &value);
  tree.Branch("log_val", &logValue);
  tree.Branch("log10_chi2", &log10Chi2);

  std::vector<double> parameters, outputMoments;
  detail::MakeParameterBranches(&tree, nominal->fullPars, parameters);
  detail::MakeOutputMomentBranches(&tree, nominal->outputMoments, outputMoments);

  const unsigned progressStep = std::max(1u, toys / 10u);
  for (toy = 0; toy < static_cast<int>(toys); ++toy) {
    if (config.verbose &&
        (static_cast<unsigned>(toy) % progressStep == 0 ||
         static_cast<unsigned>(toy) + 1 == toys)) {
      std::cout << "toy " << toy + 1 << '/' << toys << '\n';
    }
    auto sampled = SampleContext(*nominal, random);
    const BestFit best = FitToy(config, sampled, random);
    valid = best.valid;
    status = best.status;
    validStarts = best.validStarts;
    if (best.valid) {
      chi2 = best.chi2;
      value = ndf > 0 ? chi2 / static_cast<double>(ndf)
                      : std::numeric_limits<double>::quiet_NaN();
      logValue = value > 0.0 ? std::log(value) : -999.0;
      log10Chi2 = chi2 > 0.0 ? std::log10(chi2) : -999.0;
      parameters = best.parameters;
      outputMoments = best.outputMoments;
    } else {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      chi2 = value = logValue = log10Chi2 = nan;
      std::fill(parameters.begin(), parameters.end(), nan);
      std::fill(outputMoments.begin(), outputMoments.end(), nan);
    }
    tree.Fill();
  }
  tree.Write();
  file->Close();
}

} // namespace

void RunBootstrap(const BootstrapConfig& config, const ModelConfig& model) {
  if (config.fit.input.empty()) throw std::invalid_argument("Bootstrap input is required");
  if (config.fit.output.empty()) throw std::invalid_argument("Bootstrap output is required");
  if (config.toys == 0) throw std::invalid_argument("Bootstrap toys must be positive");
  if (config.startsPerToy == 0) {
    throw std::invalid_argument("Bootstrap starts per toy must be positive");
  }
  std::filesystem::create_directories(config.fit.output.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : config.fit.output.parent_path());

  unsigned workers = config.fit.workers == 0
                         ? std::max(1u, std::thread::hardware_concurrency())
                         : config.fit.workers;
  workers = std::min(workers, config.toys);
  FitConfig fit = config.fit;
  fit.starts = config.startsPerToy;
  if (workers == 1) {
    RunBootstrapWorker(detail::MakeInternalConfig(fit, model), config.toys,
                       fit.output);
    return;
  }

  std::vector<unsigned> ids(workers);
  for (unsigned i = 0; i < workers; ++i) ids[i] = i;
  ROOT::TProcessExecutor pool(workers);
  const auto parts = pool.Map([=](unsigned worker) {
    FitConfig workerFit = fit;
    workerFit.seed = fit.seed == 0 ? 0x9e3779b9u + 100003u * worker
                                   : fit.seed + 100003u * worker;
    workerFit.verbose = fit.verbose && worker == 0;
    const unsigned workerToys = config.toys / workers +
                                (worker < config.toys % workers);
    const auto part = PartFile(fit.output, worker);
    RunBootstrapWorker(detail::MakeInternalConfig(workerFit, model), workerToys, part);
    return part.string();
  }, ids);

  TFileMerger merger(true, false);
  merger.OutputFile(fit.output.c_str(), "RECREATE");
  for (const auto& part : parts) merger.AddFile(part.c_str());
  if (!merger.Merge()) throw std::runtime_error("Could not merge bootstrap output files");
  for (const auto& part : parts) std::filesystem::remove(part);
  if (fit.verbose) std::cout << "Saved " << fit.output << '\n';
}

} // namespace emi
