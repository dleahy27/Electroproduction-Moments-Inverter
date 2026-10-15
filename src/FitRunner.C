// Coordinate multi-start amplitude fits.  Starts can run in separate ROOT
// processes, after which all starts and their diagnostics are merged into
// the final ROOT tree in deterministic worker order.
// Multiple starts are essential because the amplitude-to-moment map is
// nonlinear and can contain discrete ambiguities as well as local minima.
#include "emi/Runner.h"

#include "Detail.h"
#include "Progress.h"
#include "WorkerOutput.h"

#include "ROOT/TProcessExecutor.hxx"
#include "TBenchmark.h"
#include "TFile.h"
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
  // ROOT branches keep pointers to these vector elements until each Fill().
  // Resize once before creating branches so later reallocation cannot
  // invalidate those addresses.
  storage.assign(definitions.size(), std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i = 0; i < definitions.size(); ++i) {
    tree->Branch(("err__" + definitions[i].name).c_str(), &storage[i]);
  }
}

template <typename Definition>
void MakeCorrelationBranches(TTree* tree, const std::vector<Definition>& definitions,
                             std::vector<double>& storage) {
  const std::size_t count = definitions.size();
  storage.assign(count * (count - 1) / 2, std::numeric_limits<double>::quiet_NaN());
  std::size_t pairIndex = 0;
  for (std::size_t first = 0; first < count; ++first) {
    for (std::size_t second = first + 1; second < count; ++second) {
      const std::string name =
          "corr__" + definitions[first].name + "__" + definitions[second].name;
      tree->Branch(name.c_str(), &storage[pairIndex++]);
    }
  }
}

void RunWorker(const detail::InternalConfig& config,
               const std::filesystem::path& output,
               detail::ProgressCounter progress = {}) {
  auto context = detail::BuildContext(config);
  TRandom3 random(config.seed);
  if (config.seed == 0)
    random.SetSeed(0);

  std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "RECREATE"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Could not create output file " + output.string());
  }

  TTree tree("fitResults", "Amplitude fits from independent random starts");
  // One row represents one independent start, including failed starts. This
  // makes convergence diagnostics possible and lets analyses choose the best
  // accepted minimum rather than trusting the final attempted start.
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
  std::vector<double> parameterErrors, parameterCorrelations;
  std::vector<double> outputMomentErrors;
  double ratio = std::numeric_limits<double>::quiet_NaN();
  double ratioError = std::numeric_limits<double>::quiet_NaN();
  detail::MakeParameterBranches(&tree, context->fullPars, parameters);
  detail::MakeOutputMomentBranches(&tree, context->outputMoments, outputMoments);
  if (config.runHesse) {
    MakeErrorBranches(&tree, context->fullPars, parameterErrors);
    MakeCorrelationBranches(&tree, context->fullPars, parameterCorrelations);
    MakeErrorBranches(&tree, context->outputMoments, outputMomentErrors);
  }
  tree.Branch("R", &ratio);
  if (config.runHesse)
    tree.Branch("err__R", &ratioError);

  TBenchmark benchmark;
  benchmark.Start("fit");
  for (unsigned startIndex = 0; startIndex < config.starts; ++startIndex) {

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
    logValue = value > 0.0 ? std::log(value) : -std::numeric_limits<double>::infinity();

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
      // HESSE is meaningful only at the returned minimum. Failed covariance
      // estimates remain NaN rather than being mistaken for zero uncertainty.
      detail::FillHessianProducts(*context, *result.minimizer, parameters,
                                  hesseOk && covarianceStatus >= 1, parameterErrors,
                                  parameterCorrelations, outputMomentErrors, ratio,
                                  ratioError);
    }
    tree.Fill();
    progress.Advance();
  }

  benchmark.Stop("fit");
  if (config.verbose)
    benchmark.Print("fit");
  tree.Write();
  file->Close();
}

} // namespace

void RunFit(const FitConfig& config, const ModelConfig& model) {
  // Each worker builds its own equivalent context and explores disjoint
  // random starts, avoiding shared mutable ROOT and Minuit objects.
  if (config.input.empty())
    throw std::invalid_argument("Fit input is required");
  if (config.output.empty())
    throw std::invalid_argument("Fit output is required");
  if (config.starts == 0)
    throw std::invalid_argument("Fit starts must be positive");
  std::filesystem::create_directories(config.output.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : config.output.parent_path());

  unsigned workers = config.workers == 0
                         ? std::max(1u, std::thread::hardware_concurrency())
                         : config.workers;
  workers = std::min(workers, config.starts);
  if (workers == 1) {
    detail::ProgressMonitor progress(config.verbose, config.starts, workers, "fit");
    RunWorker(detail::MakeInternalConfig(config, model), config.output,
              progress.Counter());
    return;
  }

  std::vector<unsigned> ids(workers);
  for (unsigned i = 0; i < workers; ++i)
    ids[i] = i;
  ROOT::TProcessExecutor pool(workers);
  detail::ProgressMonitor progress(config.verbose, config.starts, workers, "fit");
  const auto counter = progress.Counter();
  // Large, separated seed offsets prevent workers from replaying the same
  // random starts. Seed zero still remains nondeterministic from the user's
  // perspective, while each worker receives a distinct concrete seed.
  const auto parts = pool.Map(
      [=](unsigned worker) {
        FitConfig workerConfig = config;
        workerConfig.starts =
            config.starts / workers + (worker < config.starts % workers);
        workerConfig.seed = config.seed == 0 ? 0x9e3779b9u + 100003u * worker
                                             : config.seed + 100003u * worker;
        workerConfig.verbose = false;
        const auto part = detail::WorkerPartFile(config.output, worker);
        RunWorker(detail::MakeInternalConfig(workerConfig, model), part, counter);
        return part.string();
      },
      ids);

  detail::MergeWorkerFiles(config.output, parts);
  for (const auto& part : parts)
    std::filesystem::remove(part);
  if (config.verbose) {
    std::cout << "Saved " << config.output << " (" << workers << " workers)\n";
  }
}

} // namespace emi
