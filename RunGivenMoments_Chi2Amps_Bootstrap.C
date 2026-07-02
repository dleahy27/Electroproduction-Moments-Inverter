#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// Load the minimiser implementation.  The first include is the name I would use
// in the repository; the second allows this uploaded copy to run without being
// renamed first.
#if __has_include("RunGivenMoments_Chi2Amps_normalised_analyticgrad.C")
#include "RunGivenMoments_Chi2Amps_normalised_analyticgrad.C"
#elif __has_include("RunGivenMoments_Chi2Amps_normalised_analyticgrad(1).C")
#include "RunGivenMoments_Chi2Amps_normalised_analyticgrad(1).C"
#else
#error "Could not find RunGivenMoments_Chi2Amps_normalised_analyticgrad.C"
#endif

namespace chi2_amp_fit_toys {
using namespace chi2_amp_fit_opt;

static constexpr double kFitValDenominator = 7.0; // Keep the same val/log_val convention as the new main fitter.

static std::string SafeBranchName(const std::string& prefix, const std::string& name) {
  std::string out = prefix + name;
  for (char& c : out) {
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_')) c = '_';
  }
  return out;
}

static void MakePrefixedBranches(TTree* t,
                                 const std::string& prefix,
                                 const std::vector<ObservedMoment>& observed,
                                 std::vector<double>& storage) {
  storage.assign(observed.size(), 0.0);
  for (size_t i = 0; i < observed.size(); ++i) {
    const std::string branchName = SafeBranchName(prefix, observed[i].name);
    t->Branch(branchName.c_str(), &storage[i]);
  }
}

static void ThrowObservedMoment(ObservedMoment& moment, TRandom3& rng) {
  const double sigma = std::abs(moment.sigma);
  if (!std::isfinite(moment.value) || !(sigma > 0.0) || !std::isfinite(sigma)) return;
  moment.value = rng.Gaus(moment.value, sigma);
}

static std::vector<ObservedMoment> SampleObservedMoments(const std::vector<ObservedMoment>& nominal,
                                                         TRandom3& rng) {
  std::vector<ObservedMoment> sampled = nominal;
  for (auto& moment : sampled) ThrowObservedMoment(moment, rng);
  return sampled;
}

static void FillObservedValues(const std::vector<ObservedMoment>& observed,
                               std::vector<double>& storage) {
  storage.assign(observed.size(), 0.0);
  for (size_t i = 0; i < observed.size(); ++i) storage[i] = observed[i].value;
}

static void RebuildObservedIndexMaps(EvalContext& ctx) {
  ctx.observedModelIdx.assign(ctx.observed.size(), -1);
  ctx.observedModelIdx0.assign(ctx.observed.size(), -1);
  ctx.observedModelIdx4.assign(ctx.observed.size(), -1);

  for (size_t i = 0; i < ctx.observed.size(); ++i) {
    const auto& ob = ctx.observed[i];

    if (ob.isMixed04) {
      std::ostringstream os0, os4;
      os0 << "H_0_" << ob.L << '_' << ob.M;
      os4 << "H_4_" << ob.L << '_' << ob.M;

      auto it0 = ctx.modelIndexByName.find(os0.str());
      auto it4 = ctx.modelIndexByName.find(os4.str());
      if (it0 != ctx.modelIndexByName.end()) ctx.observedModelIdx0[i] = static_cast<int>(it0->second);
      if (it4 != ctx.modelIndexByName.end()) ctx.observedModelIdx4[i] = static_cast<int>(it4->second);
      if (ctx.observedModelIdx0[i] < 0 || ctx.observedModelIdx4[i] < 0) {
        throw std::runtime_error("Failed to map " + ob.name + " to H0/H4 model moments");
      }
      continue;
    }

    std::string key = ob.name;
    if (key.rfind("RH_", 0) == 0) key = "H_" + key.substr(3);
    auto it = ctx.modelIndexByName.find(key);
    if (it != ctx.modelIndexByName.end()) ctx.observedModelIdx[i] = static_cast<int>(it->second);
    if (ctx.observedModelIdx[i] < 0) {
      throw std::runtime_error("Failed to map observed moment " + ob.name);
    }
  }
}

static std::shared_ptr<EvalContext> BuildContextWithObserved(const std::shared_ptr<const EvalContext>& templateCtx,
                                                             const std::vector<ObservedMoment>& observed) {
  if (!templateCtx) throw std::runtime_error("BuildContextWithObserved: null template context");
  auto ctx = std::make_shared<EvalContext>(*templateCtx);
  ctx->observed = observed;
  RebuildObservedIndexMaps(*ctx);

  // These members are per-fit scratch/output state and should not be shared
  // between toys copied from the template context.
  ctx->iterTree = nullptr;
  ctx->iter_parVals.clear();
  ctx->iter_Hrec.clear();
  ctx->iter_log_val = 0.0;
  ctx->callCount = 0;
  return ctx;
}

struct BestFitResult {
  bool valid = false;
  int status = -999;
  unsigned validStarts = 0;
  double chi2 = std::numeric_limits<double>::infinity();
  double val = std::numeric_limits<double>::quiet_NaN();
  double log_val = std::numeric_limits<double>::quiet_NaN();
  double log10_chi2 = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> parVals;
  std::vector<double> momVals;
  std::vector<double> fittedObservedVals;
};

static BestFitResult FitBestForObserved(const FitConfig& cfg,
                                        const std::shared_ptr<const EvalContext>& templateCtx,
                                        const std::vector<ObservedMoment>& observed,
                                        TRandom3& rng) {
  auto ctx = BuildContextWithObserved(templateCtx, observed);
  Chi2FunctionNoGrad fcnNoGrad(ctx);
  Chi2Function fcn(ctx);

  const unsigned nDim = fcnNoGrad.NDim();

  BestFitResult best;
  best.parVals.assign(ctx->fullPars.size(), std::numeric_limits<double>::quiet_NaN());
  best.momVals.assign(ctx->modelsRec.size(), std::numeric_limits<double>::quiet_NaN());
  best.fittedObservedVals.assign(ctx->observed.size(), std::numeric_limits<double>::quiet_NaN());

  for (unsigned iStart = 0; iStart < cfg.nStarts; ++iStart) {
    ctx->callCount = 0;

    std::vector<double> startVals;
    BuildRandomStartPoint(*ctx, rng, startVals);
    if (cfg.useMCMCPreScan && cfg.mcmcSteps > 0) {
      const auto mcmc = RunMCMCPreScan(*ctx, cfg, fcnNoGrad, rng, startVals);
      startVals = mcmc.xBest;
    }

    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    if (!min) throw std::runtime_error("Failed to create Minuit2/Migrad minimizer");

    min->SetMaxFunctionCalls(cfg.maxCalls);
    min->SetMaxIterations(cfg.maxIters);
    min->SetTolerance(cfg.tolerance);
    min->SetStrategy(cfg.strategy);
    min->SetPrintLevel(cfg.printLevel);

    if (cfg.useNumericalGradient) {
      min->SetFunction(fcnNoGrad);
    } else {
      min->SetFunction(fcn);
    }

    for (unsigned i = 0; i < nDim; ++i) {
      const int fullIdx = ctx->freeToFull[i];
      const auto& p = ctx->fullPars[static_cast<size_t>(fullIdx)];
      const double step = (p.step > 0.0) ? p.step : 1e-3;
      min->SetLimitedVariable(i, p.name.c_str(), startVals[i], step, p.low, p.high);
    }

    const bool ok = min->Minimize();
    const int status = min->Status();
    if (cfg.runHesse && ok) min->Hesse();

    std::vector<double> xBest(nDim, 0.0);
    for (unsigned i = 0; i < nDim; ++i) xBest[i] = min->X()[i];

    const double chi2 = EvaluateChi2AtPoint(fcnNoGrad, xBest);
    if (!std::isfinite(chi2) || !(chi2 < 1e299)) continue;
    ++best.validStarts;
    if (chi2 >= best.chi2) continue;

    std::vector<double> fullVals;
    if (!FillFullFromFree(*ctx, xBest.data(), fullVals)) continue;

    best.valid = true;
    best.status = status;
    best.chi2 = chi2;
    best.val = chi2 / kFitValDenominator;
    best.log_val = (best.val > 0.0) ? std::log(best.val) : -999.0;
    best.log10_chi2 = (chi2 > 0.0) ? std::log10(chi2) : -999.0;
    best.parVals = std::move(fullVals);

    EvalAllMoments(*ctx, best.parVals, best.momVals);
    FillObservedModelValues(*ctx, best.momVals, best.fittedObservedVals);
  }

  return best;
}

static void FillWithNaNs(std::vector<double>& values) {
  std::fill(values.begin(), values.end(), std::numeric_limits<double>::quiet_NaN());
}

static void ApplyChi2AmpsDefaults(FitConfig& cfg) {
  // Match the current normalised analytic-gradient minimiser defaults.
  cfg.useNumericalGradient = false;
  cfg.useMCMCPreScan = false;
  cfg.mcmcSteps = 2000;
  cfg.mcmcProposalMagSigma = 0.03;
  cfg.mcmcProposalPhaseSigma = 0.10;
  cfg.mcmcTemperature = 1.0;
  cfg.magnitudeStartMean = 0.5;
  cfg.magnitudeStartSigma = 0.5;

  // HESSE is unnecessary for bootstrap toys because the toy spread gives the
  // uncertainty.  Keeping this false makes the bootstrap much faster.  It can
  // still be re-enabled through the setup wrapper below.
  cfg.runHesse = false;
}

static void RunGivenMoments_Chi2Amps_Toys_Impl(const FitConfig& cfg,
                                               unsigned nToys,
                                               const char* outFile) {
  if (cfg.momentsFile.empty()) throw std::runtime_error("Toy fit requires cfg.momentsFile");
  if (cfg.momentsTree.empty()) throw std::runtime_error("Toy fit requires cfg.momentsTree");
  if (nToys == 0) throw std::runtime_error("Toy fit requires nToys > 0");

  TRandom3 rng(cfg.randomSeed);
  if (cfg.randomSeed == 0) rng.SetSeed(0);

  auto templateCtx = BuildContext(cfg);
  if (!templateCtx) throw std::runtime_error("Failed to build toy template context");

  std::unique_ptr<TFile> fout(TFile::Open(outFile, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error(std::string("Failed to open output file: ") + outFile);
  }

  TTree* t = new TTree("PartialWaves", "Best fit per Gaussian-thrown moment in the bootstrap");
  int toy = 0;
  int valid = 0;
  int status = -999;
  unsigned valid_starts = 0;
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  double val = std::numeric_limits<double>::quiet_NaN();
  double log_val = std::numeric_limits<double>::quiet_NaN();
  double log10_chi2 = std::numeric_limits<double>::quiet_NaN();

  t->Branch("toy", &toy);
  t->Branch("valid", &valid);
  t->Branch("status", &status);
  t->Branch("valid_starts", &valid_starts);
  t->Branch("chi2", &chi2);
  t->Branch("val", &val);
  t->Branch("log_val", &log_val);
  t->Branch("log10_chi2", &log10_chi2);

  std::vector<double> thrownObservedVals;
  MakePrefixedBranches(t, "throw_", templateCtx->observed, thrownObservedVals);

  std::vector<double> fittedObservedVals;
  MakePrefixedBranches(t, "fit_", templateCtx->observed, fittedObservedVals);

  std::vector<double> parVals;
  MakeBranchesForPars(t, templateCtx->fullPars, parVals);

  std::vector<double> momRecVals;
  MakeBranchesForMoments(t, templateCtx->modelsRec, momRecVals);

  TBenchmark bench;
  bench.Start("toyfit");

  const unsigned progressStep = std::max(1u, nToys / 10u);
  for (toy = 0; toy < static_cast<int>(nToys); ++toy) {
    if (cfg.verbose && (static_cast<unsigned>(toy) % progressStep == 0 || static_cast<unsigned>(toy) + 1 == nToys)) {
      std::cout << "toy " << (toy + 1) << '/' << nToys << std::endl;
    }

    const std::vector<ObservedMoment> observed = SampleObservedMoments(templateCtx->observed, rng);
    FillObservedValues(observed, thrownObservedVals);

    const BestFitResult best = FitBestForObserved(cfg, templateCtx, observed, rng);

    valid = best.valid ? 1 : 0;
    status = best.status;
    valid_starts = best.validStarts;

    if (best.valid) {
      chi2 = best.chi2;
      val = best.val;
      log_val = best.log_val;
      log10_chi2 = best.log10_chi2;
      parVals = best.parVals;
      momRecVals = best.momVals;
      fittedObservedVals = best.fittedObservedVals;
    } else {
      chi2 = std::numeric_limits<double>::quiet_NaN();
      val = std::numeric_limits<double>::quiet_NaN();
      log_val = std::numeric_limits<double>::quiet_NaN();
      log10_chi2 = std::numeric_limits<double>::quiet_NaN();
      FillWithNaNs(parVals);
      FillWithNaNs(momRecVals);
      FillWithNaNs(fittedObservedVals);
    }

    t->Fill();
  }

  bench.Stop("toyfit");
  if (cfg.verbose) bench.Print("toyfit");

  fout->Write();
  fout->Close();

  if (cfg.verbose) {
    std::cout << "Saved: " << outFile << std::endl;
  }
}

} // namespace chi2_amp_fit_toys

void RunGivenMoments_Chi2Amps_Toys_Setup(const char* tableFile = "InputFiles/Experiment/e_rho_moments.root",
                                         const char* treeName = "expMoments",
                                         int bin = 1,
                                         unsigned nToys = 1000,
                                         unsigned nStartsPerToy = 1,
                                         const char* outFile = "resultsGivenMoments_chi2_amps_toys.root",
                                         uint32_t seed = 0,
                                         double epsilon = 1.0,
                                         bool verbose = false,
                                         bool photoProduction = false,
                                         bool useNumericalGradient = false,
                                         bool useMCMCPreScan = false,
                                         unsigned mcmcSteps = 2000,
                                         double mcmcProposalMagSigma = 0.03,
                                         double mcmcProposalPhaseSigma = 0.10,
                                         double mcmcTemperature = 1.0,
                                         double magnitudeStartMean = 0.5,
                                         double magnitudeStartSigma = 0.5,
                                         bool runHesse = false) {
  chi2_amp_fit_opt::FitConfig cfg;
  chi2_amp_fit_toys::ApplyChi2AmpsDefaults(cfg);
  cfg.nStarts = nStartsPerToy;
  cfg.randomSeed = seed;
  cfg.epsilon = epsilon;
  cfg.verbose = verbose;
  cfg.photoProduction = photoProduction;
  cfg.useNumericalGradient = useNumericalGradient;
  cfg.useMCMCPreScan = useMCMCPreScan;
  cfg.mcmcSteps = mcmcSteps;
  cfg.mcmcProposalMagSigma = mcmcProposalMagSigma;
  cfg.mcmcProposalPhaseSigma = mcmcProposalPhaseSigma;
  cfg.mcmcTemperature = mcmcTemperature;
  cfg.magnitudeStartMean = magnitudeStartMean;
  cfg.magnitudeStartSigma = magnitudeStartSigma;
  cfg.runHesse = runHesse;
  cfg.momentsFile = tableFile ? tableFile : "InputFiles/Experiment/e_rho_moments.root";
  cfg.momentsTree = treeName ? treeName : "expMoments";
  cfg.bin = bin;

  chi2_amp_fit_toys::RunGivenMoments_Chi2Amps_Toys_Impl(cfg, nToys, outFile);
}

namespace {
std::string MakeToyPartFileName(const char* outFile, unsigned workerId) {
  std::string base = outFile ? std::string(outFile) : std::string("toyFits.root");
  const std::string ext = ".root";
  if (base.size() >= ext.size() &&
      base.compare(base.size() - ext.size(), ext.size(), ext) == 0) {
    base.erase(base.size() - ext.size());
  }

  std::ostringstream os;
  os << base << ".part_" << workerId << ".root";
  return os.str();
}
}

void RunGivenMoments_Chi2Amps_Bootstrap(const char* tableFile = "InputFiles/Experiment/e_rho_moments.root",
                                        const char* treeName = "expMoments",
                                        int bin = 1,
                                        const char* outFile = "toyFits.root",
                                        double epsilon = 1.0,
                                        bool photoProduction = false,
                                        unsigned nToys = 1000,
                                        unsigned nStartsPerToy = 1000,
                                        unsigned nCores = 14,
                                        uint32_t seed = 0,
                                        bool useNumericalGradient = false,
                                        bool useMCMCPreScan = false,
                                        unsigned mcmcSteps = 2000,
                                        double mcmcProposalMagSigma = 0.03,
                                        double mcmcProposalPhaseSigma = 0.10,
                                        double mcmcTemperature = 1.0,
                                        double magnitudeStartMean = 0.5,
                                        double magnitudeStartSigma = 0.5,
                                        bool runHesse = false) {
  if (nToys == 0) throw std::runtime_error("RunGivenMoments_Chi2Amps_Bootstrap: nToys must be > 0");
  if (nStartsPerToy == 0) throw std::runtime_error("RunGivenMoments_Chi2Amps_Bootstrap: nStartsPerToy must be > 0");

  unsigned nWorkers = nCores ? nCores : std::max(1u, std::thread::hardware_concurrency());
  if (nWorkers > nToys) nWorkers = nToys;

  std::vector<unsigned> workerIds(nWorkers);
  for (unsigned i = 0; i < nWorkers; ++i) workerIds[i] = i;

  ROOT::TProcessExecutor pool(nWorkers);

  auto partFiles = pool.Map([=](unsigned workerId) {
    const unsigned base = nToys / nWorkers;
    const unsigned extra = nToys % nWorkers;
    const unsigned myNToys = base + (workerId < extra ? 1u : 0u);

    const uint32_t workerSeed =
      (seed == 0) ? (0x9e3779b9u + 100003u * workerId)
                  : (seed + 100003u * workerId);

    const std::string partFile = MakeToyPartFileName(outFile, workerId);
    const bool verbose = (workerId == 0);

    RunGivenMoments_Chi2Amps_Toys_Setup(tableFile,
                                        treeName,
                                        bin,
                                        myNToys,
                                        nStartsPerToy,
                                        partFile.c_str(),
                                        workerSeed,
                                        epsilon,
                                        verbose,
                                        photoProduction,
                                        useNumericalGradient,
                                        useMCMCPreScan,
                                        mcmcSteps,
                                        mcmcProposalMagSigma,
                                        mcmcProposalPhaseSigma,
                                        mcmcTemperature,
                                        magnitudeStartMean,
                                        magnitudeStartSigma,
                                        runHesse);

    return partFile;
  }, workerIds);

  TFileMerger merger(kTRUE, kFALSE);
  merger.OutputFile(outFile, "RECREATE");
  for (const auto& f : partFiles) merger.AddFile(f.c_str());

  if (!merger.Merge()) {
    throw std::runtime_error("RunGivenMoments_Chi2Amps_Bootstrap: merge failed");
  }

  for (const auto& f : partFiles) gSystem->Unlink(f.c_str());

  std::cout << "Merged " << partFiles.size() << " files into " << outFile << std::endl;
}
