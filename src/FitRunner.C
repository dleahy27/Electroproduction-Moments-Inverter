#include "emi/Runner.h"

#include "Detail.h"

#include "ROOT/TProcessExecutor.hxx"
#include "TBenchmark.h"
#include "TFile.h"
#include "TFileMerger.h"
#include "TRandom3.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>

namespace emi {
namespace {

template <typename Definition>
void MakeErrorBranches(TTree* tree, const std::vector<Definition>& definitions,
                       std::vector<double>& storage) {
  storage.assign(definitions.size(), std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i = 0; i < definitions.size(); ++i) {
    tree->Branch(("err__" + definitions[i].name).c_str(), &storage[i]);
  }
}

std::filesystem::path PartFile(const std::filesystem::path& output, unsigned worker) {
  const auto stem = output.stem().string() + ".part_" + std::to_string(worker);
  return output.parent_path() / (stem + output.extension().string());
}

void RunWorker(const detail::InternalConfig& config,
               const std::filesystem::path& output) {
  auto context = detail::BuildContext(config);
  TRandom3 random(config.seed);
  if (config.seed == 0) random.SetSeed(0);

  std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "RECREATE"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Could not create output file " + output.string());
  }

  TTree tree("fitResults", "Amplitude fits from independent random starts");
  bool fitOk = false;
  bool hesseOk = false;
  int status = -999;
  int covarianceStatus = -1;
  int ndf = static_cast<int>(context->observed.size()) -
            static_cast<int>(context->freeToFull.size());
  double edm = std::numeric_limits<double>::quiet_NaN();
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  double value = std::numeric_limits<double>::quiet_NaN();
  double logValue = std::numeric_limits<double>::quiet_NaN();
  tree.Branch("fit_ok", &fitOk);
  tree.Branch("hesse_ok", &hesseOk);
  tree.Branch("status", &status);
  tree.Branch("cov_status", &covarianceStatus);
  tree.Branch("ndf", &ndf);
  tree.Branch("edm", &edm);
  tree.Branch("chi2", &chi2);
  tree.Branch("val", &value);
  tree.Branch("log_val", &logValue);

  std::vector<double> parameters, rawMoments, outputMoments;
  std::vector<double> parameterErrors, outputMomentErrors;
  double ratio = std::numeric_limits<double>::quiet_NaN();
  double ratioError = std::numeric_limits<double>::quiet_NaN();
  detail::MakeParameterBranches(&tree, context->fullPars, parameters);
  detail::MakeOutputMomentBranches(&tree, context->outputMoments, outputMoments);
  if (config.runHesse) {
    MakeErrorBranches(&tree, context->fullPars, parameterErrors);
    MakeErrorBranches(&tree, context->outputMoments, outputMomentErrors);
  }
  tree.Branch("R", &ratio);
  if (config.runHesse) tree.Branch("err__R", &ratioError);

  TBenchmark benchmark;
  benchmark.Start("fit");
  const unsigned progressStep = std::max(1u, config.starts / 10u);
  for (unsigned startIndex = 0; startIndex < config.starts; ++startIndex) {
    if (config.verbose &&
        (startIndex % progressStep == 0 || startIndex + 1 == config.starts)) {
      std::cout << "start " << startIndex + 1 << '/' << config.starts << '\n';
    }

    std::vector<double> start;
    detail::BuildRandomStart(*context, random, start);
    auto result = detail::Minimize(context, start, config.runHesse);
    fitOk = result.fitOk;
    hesseOk = result.hesseOk;
    status = result.status;
    covarianceStatus = result.covarianceStatus;
    edm = result.edm;
    chi2 = result.chi2;
    value = ndf > 0 ? chi2 / static_cast<double>(ndf)
                    : std::numeric_limits<double>::quiet_NaN();
    logValue = value > 0.0 ? std::log(value)
                           : -std::numeric_limits<double>::infinity();

    if (!detail::FillFullParameters(*context, result.freeValues.data(), parameters)) {
      throw std::runtime_error("Minuit result violates amplitude normalisation");
    }
    detail::EvaluateAllMoments(*context, parameters, rawMoments);
    detail::EvaluateOutputMoments(*context, rawMoments, outputMoments);
    ratio = std::numeric_limits<double>::quiet_NaN();
    if (context->idxH0_00 >= 0 && context->idxH4_00 >= 0) {
      const double h0 = rawMoments[static_cast<size_t>(context->idxH0_00)];
      if (std::abs(h0) > 1e-15) {
        ratio = rawMoments[static_cast<size_t>(context->idxH4_00)] / h0;
      }
    }
    if (config.runHesse) {
      detail::FillHessianProducts(
          *context, *result.minimizer, parameters,
          hesseOk && covarianceStatus >= 1, parameterErrors,
          outputMomentErrors, ratio, ratioError);
    }
    tree.Fill();
  }

  benchmark.Stop("fit");
  if (config.verbose) benchmark.Print("fit");
  tree.Write();
  file->Close();
}

} // namespace

void RunFit(const FitConfig& config, const ModelConfig& model) {
  if (config.input.empty()) throw std::invalid_argument("Fit input is required");
  if (config.output.empty()) throw std::invalid_argument("Fit output is required");
  if (config.starts == 0) throw std::invalid_argument("Fit starts must be positive");
  std::filesystem::create_directories(config.output.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : config.output.parent_path());

  unsigned workers = config.workers == 0
                         ? std::max(1u, std::thread::hardware_concurrency())
                         : config.workers;
  workers = std::min(workers, config.starts);
  if (workers == 1) {
    RunWorker(detail::MakeInternalConfig(config, model), config.output);
    return;
  }

  std::vector<unsigned> ids(workers);
  for (unsigned i = 0; i < workers; ++i) ids[i] = i;
  ROOT::TProcessExecutor pool(workers);
  const auto parts = pool.Map([=](unsigned worker) {
    FitConfig workerConfig = config;
    workerConfig.starts = config.starts / workers + (worker < config.starts % workers);
    workerConfig.seed = config.seed == 0 ? 0x9e3779b9u + 100003u * worker
                                         : config.seed + 100003u * worker;
    workerConfig.verbose = config.verbose && worker == 0;
    const auto part = PartFile(config.output, worker);
    RunWorker(detail::MakeInternalConfig(workerConfig, model), part);
    return part.string();
  }, ids);

  TFileMerger merger(true, false);
  merger.OutputFile(config.output.c_str(), "RECREATE");
  for (const auto& part : parts) merger.AddFile(part.c_str());
  if (!merger.Merge()) throw std::runtime_error("Could not merge worker output files");
  for (const auto& part : parts) std::filesystem::remove(part);
  if (config.verbose) {
    std::cout << "Saved " << config.output << " (" << workers << " workers)\n";
  }
}

} // namespace emi
