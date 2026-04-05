#include <limits>

#include "RunGivenMoments_Chi2Amps.C"

namespace chi2_amp_fit_toys {
using namespace chi2_amp_fit_opt;

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
      continue;
    }

    std::string key = ob.name;
    if (key.rfind("RH_", 0) == 0) key = "H_" + key.substr(3);
    auto it = ctx.modelIndexByName.find(key);
    if (it != ctx.modelIndexByName.end()) ctx.observedModelIdx[i] = static_cast<int>(it->second);
  }
}

static std::shared_ptr<EvalContext> BuildContextWithObserved(const std::shared_ptr<const EvalContext>& templateCtx,
                                                             const std::vector<ObservedMoment>& observed) {
  if (!templateCtx) throw std::runtime_error("BuildContextWithObserved: null template context");
  auto ctx = std::make_shared<EvalContext>(*templateCtx);
  ctx->observed = observed;
  RebuildObservedIndexMaps(*ctx);
  ctx->iterTree = nullptr;
  ctx->callCount = 0;
  return ctx;
}

struct BestFitResult {
  bool valid = false;
  double chi2 = std::numeric_limits<double>::infinity();
  double log_val = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> parVals;
  std::vector<double> momVals;
};

static BestFitResult FitBestForObserved(const FitConfig& cfg,
                                        const std::shared_ptr<const EvalContext>& templateCtx,
                                        const std::vector<ObservedMoment>& observed,
                                        TRandom3& rng) {
  auto ctx = BuildContextWithObserved(templateCtx, observed);
  Chi2FunctionNoGrad fcnNoGrad(ctx);
  Chi2Function fcn(ctx);

  BestFitResult best;
  best.parVals.assign(ctx->fullPars.size(), std::numeric_limits<double>::quiet_NaN());
  best.momVals.assign(ctx->modelsRec.size(), std::numeric_limits<double>::quiet_NaN());

  for (unsigned iStart = 0; iStart < cfg.nStarts; ++iStart) {
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

    std::vector<double> startVals;
    BuildRandomStartPoint(*ctx, rng, startVals);
    if (cfg.useMCMCPreScan && cfg.mcmcSteps > 0) {
      const auto mcmc = RunMCMCPreScan(*ctx, cfg, fcnNoGrad, rng, startVals);
      startVals = mcmc.xBest;
    }

    for (unsigned i = 0; i < fcn.NDim(); ++i) {
      const int fullIdx = ctx->freeToFull[i];
      const auto& p = ctx->fullPars[static_cast<size_t>(fullIdx)];
      if (p.isPhase) {
        const double step = (p.step > 0.0) ? p.step : 1e-3;
        min->SetLimitedVariable(i, p.name.c_str(), startVals[i], step, p.low, p.high);
      } else {
        const std::string coordName = std::string("logit_") + p.name;
        const double step = (cfg.simplexLogitStep > 0.0) ? cfg.simplexLogitStep : 0.2;
        min->SetVariable(i, coordName.c_str(), startVals[i], step);
      }
    }

    const bool ok = min->Minimize();
    if (cfg.runHesse && ok) min->Hesse();

    const double chi2 = min->MinValue();
    if (!std::isfinite(chi2) || chi2 >= best.chi2) continue;

    best.valid = true;
    best.chi2 = chi2;
    best.log_val = (chi2 > 0.0) ? std::log10(chi2) : -999.0;

    if (!FillFullFromFree(*ctx, min->X(), best.parVals)) {
      throw std::runtime_error("Failed to map simplex coordinates to physical amplitudes in toy fit");
    }
    EvalAllMoments(*ctx, best.parVals, best.momVals);
  }

  return best;
}

static void FillWithNaNs(std::vector<double>& values) {
  std::fill(values.begin(), values.end(), std::numeric_limits<double>::quiet_NaN());
}

static void ApplyChi2AmpsDefaults(FitConfig& cfg) {
  cfg.depNormMagName = "a_T_1_1";
  cfg.useNumericalGradient = false;
  cfg.dirichletMagnitudeAlpha = 1.0;
  cfg.simplexLogitStep = 0.01;
}

static void RunGivenMoments_Chi2Amps_Toys_Impl(const FitConfig& cfg,
                                               unsigned nToys,
                                               const char* outFile) {
  if (cfg.momentsFile.empty()) throw std::runtime_error("Toy fit requires cfg.momentsFile");
  if (cfg.momentsTree.empty()) throw std::runtime_error("Toy fit requires cfg.momentsTree");

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
  double log_val = std::numeric_limits<double>::quiet_NaN();
  t->Branch("toy", &toy);
  t->Branch("log_val", &log_val);

  std::vector<double> toyObservedVals;
  MakeBranchesForObservedRH(t, templateCtx->observed, toyObservedVals);
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
    FillObservedValues(observed, toyObservedVals);
    const BestFitResult best = FitBestForObserved(cfg, templateCtx, observed, rng);

    if (best.valid) {
      log_val = best.log_val;
      parVals = best.parVals;
      momRecVals = best.momVals;
    } else {
      log_val = std::numeric_limits<double>::quiet_NaN();
      FillWithNaNs(parVals);
      FillWithNaNs(momRecVals);
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
                                         double epsR4 = 1.0,
                                         bool verbose = false,
                                         bool photoProduction = false) {
  chi2_amp_fit_opt::FitConfig cfg;
  chi2_amp_fit_toys::ApplyChi2AmpsDefaults(cfg);
  cfg.nStarts = nStartsPerToy;
  cfg.randomSeed = seed;
  cfg.epsR4 = epsR4;
  cfg.verbose = verbose;
  cfg.photoProduction = photoProduction;
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
                                   double epsR4 = 1.0,
                                   bool photoProduction = false){

  // SETUP -- majority should already be set inside RunGivenMoments_Chi2A

  unsigned nToys = 1000;
  unsigned nStartsPerToy = 1000;
  unsigned nCores = 10;
  const uint32_t seed = 0;

  // STUFF STARTS BELOW

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
                                        epsR4,
                                        verbose,
                                        photoProduction);

    return partFile;
  }, workerIds);

  TFileMerger merger(kTRUE, kFALSE);
  merger.OutputFile(outFile, "RECREATE");
  for (const auto& f : partFiles) merger.AddFile(f.c_str());

  if (!merger.Merge()) {
    throw std::runtime_error("RunGivenMoments_Chi2Amps_Toys: merge failed");
  }

  for (const auto& f : partFiles) gSystem->Unlink(f.c_str());

  std::cout << "Merged " << partFiles.size() << " files into " << outFile << std::endl;
}