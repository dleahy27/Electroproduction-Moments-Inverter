#include <limits>

#include "RunGivenMoments_Chi2Amps.C"

namespace chi2_amp_fit_toys {
using namespace chi2_amp_fit;

static void ThrowValErr2(ValErr2& x, TRandom3& rng) {
  x.v = rng.Gaus(x.v, Comb(x));
}

static ProtonSDMEsTable SampleProtonSDMEs_TableQ2Bin(int q2bin, TRandom3& rng) {
  ProtonSDMEsTable p = GetProtonSDMEs_TableQ2Bin(q2bin);

  ThrowValErr2(p.r00_04, rng);
  ThrowValErr2(p.re_r10_04, rng);
  ThrowValErr2(p.r1m1_04, rng);

  ThrowValErr2(p.r1m1_1, rng);
  ThrowValErr2(p.re_r10_1, rng);
  ThrowValErr2(p.im_r10_2, rng);
  ThrowValErr2(p.r00_1, rng);
  ThrowValErr2(p.im_r10_3, rng);
  ThrowValErr2(p.r00_8, rng);
  ThrowValErr2(p.r11_5, rng);
  ThrowValErr2(p.r1m1_5, rng);
  ThrowValErr2(p.im_r1m1_6, rng);
  ThrowValErr2(p.im_r1m1_7, rng);
  ThrowValErr2(p.r11_8, rng);
  ThrowValErr2(p.r1m1_8, rng);
  ThrowValErr2(p.r11_1, rng);
  ThrowValErr2(p.im_r1m1_3, rng);
  ThrowValErr2(p.im_r1m1_2, rng);
  ThrowValErr2(p.re_r10_5, rng);
  ThrowValErr2(p.im_r10_6, rng);
  ThrowValErr2(p.im_r10_7, rng);
  ThrowValErr2(p.re_r10_8, rng);
  ThrowValErr2(p.r00_5, rng);

  return p;
}

static std::vector<ObservedMoment> BuildObservedProtonMomentsFromTable(const ProtonSDMEsTable& p) {
  std::vector<ObservedMoment> obs;
  obs.reserve(24);

  auto add = [&](int alpha, int L, int M, double val, double sig) {
    ObservedMoment m;
    m.alpha = alpha;
    m.L = L;
    m.M = M;
    m.value = val;
    m.sigma = sig;
    std::ostringstream os;
    os << "RH_" << alpha << '_' << L << '_' << M;
    m.name = os.str();
    obs.push_back(m);
  };

  auto add04 = [&](int L, int M, double val, double sig) {
    ObservedMoment m;
    m.alpha = 0;
    m.L = L;
    m.M = M;
    m.value = val;
    m.sigma = sig;
    m.isMixed04 = true;
    std::ostringstream os;
    os << "RH04_" << L << '_' << M;
    m.name = os.str();
    obs.push_back(m);
  };

  const double s_r00_04    = Comb(p.r00_04);
  const double s_re_r10_04 = Comb(p.re_r10_04);
  const double s_r1m1_04   = Comb(p.r1m1_04);
  const double s_r1m1_1    = Comb(p.r1m1_1);
  const double s_re_r10_1  = Comb(p.re_r10_1);
  const double s_r00_1     = Comb(p.r00_1);
  const double s_r11_1     = Comb(p.r11_1);
  const double s_im_r10_2  = Comb(p.im_r10_2);
  const double s_im_r1m1_2 = Comb(p.im_r1m1_2);
  const double s_im_r10_3  = Comb(p.im_r10_3);
  const double s_im_r1m1_3 = Comb(p.im_r1m1_3);
  const double s_r00_5     = Comb(p.r00_5);
  const double s_r11_5     = Comb(p.r11_5);
  const double s_re_r10_5  = Comb(p.re_r10_5);
  const double s_r1m1_5    = Comb(p.r1m1_5);
  const double s_im_r10_6  = Comb(p.im_r10_6);
  const double s_im_r1m1_6 = Comb(p.im_r1m1_6);
  const double s_im_r10_7  = Comb(p.im_r10_7);
  const double s_im_r1m1_7 = Comb(p.im_r1m1_7);
  const double s_r00_8     = Comb(p.r00_8);
  const double s_r11_8     = Comb(p.r11_8);
  const double s_re_r10_8  = Comb(p.re_r10_8);
  const double s_r1m1_8    = Comb(p.r1m1_8);

  const double f21 = TMath::Sqrt(12.0) / 5.0;
  const double f22 = TMath::Sqrt(6.0)  / 5.0;

  add04(2, 0, 0.2 * (3.0 * p.r00_04.v - 1.0),  TMath::Abs(0.6) * s_r00_04);
  add04(2, 1, f21 * p.re_r10_04.v,              TMath::Abs(f21) * s_re_r10_04);
  add04(2, 2, -f22 * p.r1m1_04.v,               TMath::Abs(f22) * s_r1m1_04);

  add(1, 0, 0, -(2.0 * p.r11_1.v + p.r00_1.v),
      TMath::Sqrt((2 * s_r11_1) * (2 * s_r11_1) + s_r00_1 * s_r00_1));
  add(1, 2, 0, 0.4 * (p.r11_1.v - p.r00_1.v),
      0.4 * TMath::Sqrt(s_r11_1 * s_r11_1 + s_r00_1 * s_r00_1));
  add(1, 2, 1, -f21 * p.re_r10_1.v,             TMath::Abs(f21) * s_re_r10_1);
  add(1, 2, 2,  f22 * p.r1m1_1.v,               TMath::Abs(f22) * s_r1m1_1);

  add(2, 2, 1, -f21 * p.im_r10_2.v,             TMath::Abs(f21) * s_im_r10_2);
  add(2, 2, 2,  f22 * p.im_r1m1_2.v,            TMath::Abs(f22) * s_im_r1m1_2);

  add(3, 2, 1, -f21 * p.im_r10_3.v,             TMath::Abs(f21) * s_im_r10_3);
  add(3, 2, 2,  f22 * p.im_r1m1_3.v,            TMath::Abs(f22) * s_im_r1m1_3);

  add(5, 0, 0, -(2.0 * p.r11_5.v + p.r00_5.v),
      TMath::Sqrt((2 * s_r11_5) * (2 * s_r11_5) + s_r00_5 * s_r00_5));
  add(5, 2, 0, 0.4 * (p.r11_5.v - p.r00_5.v),
      0.4 * TMath::Sqrt(s_r11_5 * s_r11_5 + s_r00_5 * s_r00_5));
  add(5, 2, 1, -f21 * p.re_r10_5.v,             TMath::Abs(f21) * s_re_r10_5);
  add(5, 2, 2,  f22 * p.r1m1_5.v,               TMath::Abs(f22) * s_r1m1_5);

  add(6, 2, 1, -f21 * p.im_r10_6.v,             TMath::Abs(f21) * s_im_r10_6);
  add(6, 2, 2,  f22 * p.im_r1m1_6.v,            TMath::Abs(f22) * s_im_r1m1_6);

  add(7, 2, 1, -f21 * p.im_r10_7.v,             TMath::Abs(f21) * s_im_r10_7);
  add(7, 2, 2,  f22 * p.im_r1m1_7.v,            TMath::Abs(f22) * s_im_r1m1_7);

  add(8, 0, 0, -(2.0 * p.r11_8.v + p.r00_8.v),
      TMath::Sqrt((2 * s_r11_8) * (2 * s_r11_8) + s_r00_8 * s_r00_8));
  add(8, 2, 0, 0.4 * (p.r11_8.v - p.r00_8.v),
      0.4 * TMath::Sqrt(s_r11_8 * s_r11_8 + s_r00_8 * s_r00_8));
  add(8, 2, 1, -f21 * p.re_r10_8.v,             TMath::Abs(f21) * s_re_r10_8);
  add(8, 2, 2,  f22 * p.r1m1_8.v,               TMath::Abs(f22) * s_r1m1_8);

  return obs;
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

static std::shared_ptr<EvalContext> BuildContextWithObserved(const FitConfig& cfg,
                                                             const std::vector<ObservedMoment>& observed) {
  auto ctx = BuildContext(cfg);
  ctx->observed = observed;
  RebuildObservedIndexMaps(*ctx);
  ctx->iterTree = nullptr;
  return ctx;
}

struct BestFitResult {
  bool valid = false;
  double chi2 = std::numeric_limits<double>::infinity();
  double log_val = -999.0;
  std::vector<double> parVals;
  std::vector<double> momVals;
};

static BestFitResult FitBestForObserved(const FitConfig& cfg,
                                        const std::vector<ObservedMoment>& observed,
                                        TRandom3& rngStarts) {
  auto ctx = BuildContextWithObserved(cfg, observed);
  Chi2Function fcn(ctx);

  BestFitResult best;
  best.parVals.assign(ctx->fullPars.size(), std::numeric_limits<double>::quiet_NaN());
  best.momVals.assign(ctx->modelsRec.size(), std::numeric_limits<double>::quiet_NaN());

  for (unsigned iStart = 0; iStart < cfg.nStarts; ++iStart) {
    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

    min->SetMaxFunctionCalls(cfg.maxCalls);
    min->SetMaxIterations(cfg.maxIters);
    min->SetTolerance(cfg.tolerance);
    min->SetStrategy(cfg.strategy);
    min->SetPrintLevel(cfg.printLevel);
    min->SetFunction(fcn);

    for (unsigned i = 0; i < fcn.NDim(); ++i) {
      const int fullIdx = ctx->freeToFull[i];
      const auto& p = ctx->fullPars[static_cast<size_t>(fullIdx)];

      double start = p.isPhase ? rngStarts.Gaus(0.0, 0.1) : rngStarts.Gaus(0.5, 0.1);
      if (start < p.low) start = p.low;
      if (start > p.high) start = p.high;

      min->SetLimitedVariable(i, p.name.c_str(), start,
                              (p.step > 0.0 ? p.step : 1e-3), p.low, p.high);
    }

    const bool ok = min->Minimize();
    if (cfg.runHesse && ok) min->Hesse();

    const double chi2 = min->MinValue();
    if (!std::isfinite(chi2) || chi2 >= best.chi2) continue;

    best.valid = true;
    best.chi2 = chi2;
    best.log_val = (chi2 > 0.0) ? std::log10(chi2) : -999.0;

    FillFullFromFree(*ctx, min->X(), best.parVals);
    for (size_t j = 0; j < ctx->fullPars.size(); ++j) {
      if (ctx->fullPars[j].isPhase) best.parVals[j] = best.parVals[j];
    }

    EvalAllMoments(ctx->modelsRec, best.parVals, best.momVals);
  }

  return best;
}

static void FillWithNaNs(std::vector<double>& values) {
  std::fill(values.begin(), values.end(), std::numeric_limits<double>::quiet_NaN());
}

static void RunGivenMoments_Chi2Amps_Toys_Impl(const FitConfig& cfg,
                                               unsigned nToys,
                                               const char* outFile) {
  TRandom3 rng(cfg.randomSeed);
  if (cfg.randomSeed == 0) rng.SetSeed(0);

  auto templateCtx = BuildContext(cfg);

  std::unique_ptr<TFile> fout(TFile::Open(outFile, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error(std::string("Failed to open output file: ") + outFile);
  }

  TTree* t = new TTree("toyFitResults", "Best fit per Gaussian-thrown SDME toy");
  int toy = 0;
  double log_val = -999.0;
  t->Branch("toy", &toy);
  t->Branch("log_val", &log_val);

  std::vector<double> parVals;
  MakeBranchesForPars(t, templateCtx->fullPars, parVals);
  std::vector<double> momRecVals;
  MakeBranchesForMoments(t, templateCtx->modelsRec, momRecVals);

  TBenchmark bench;
  bench.Start("toyfit");

  for (toy = 0; toy < static_cast<int>(nToys); ++toy) {
    if (toy % 10 == 0 && cfg.verbose) std::cout << "toy " << toy*cfg.cores << '/' << nToys*cfg.cores << std::endl;

    const ProtonSDMEsTable sampled = SampleProtonSDMEs_TableQ2Bin(cfg.observedInput.q2bin, rng);
    const std::vector<ObservedMoment> observed = BuildObservedProtonMomentsFromTable(sampled);
    const BestFitResult best = FitBestForObserved(cfg, observed, rng);

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

  if (cfg.verbose){
    std::cout << "Saving intermediate files."<< std::endl;
    std::cout << "Now compiling to: " << outFile << std::endl;
  }
}

} // namespace chi2_amp_fit_toys

void RunGivenMoments_Chi2Amps_Setup(double q2Value,
                                   unsigned nToys = 1000,
                                   unsigned nStartsPerToy = 1,
                                   const char* outFile = "resultsGivenMoments_chi2_amps_toys.root",
                                   uint32_t seed = 0,
                                   double epsR4 = 1.0, bool verbose = false, unsigned int cores = 1) {
  chi2_amp_fit::FitConfig cfg;
  cfg.nStarts = nStartsPerToy;
  cfg.randomSeed = seed;
  cfg.epsR4 = epsR4;
  cfg.recordEvery = 0;
  cfg.observedInput.source = chi2_amp_fit::ObservedMomentsSource::kQ2Table;
  cfg.observedInput.q2bin = chi2_amp_fit::GetClosestQ2Bin(q2Value);
  cfg.verbose = verbose;
  cfg.cores = cores;

  chi2_amp_fit_toys::RunGivenMoments_Chi2Amps_Toys_Impl(cfg, nToys, outFile);
}

// void RunGivenMoments_Chi2Amps_Setup(int q2bin,
//                                          unsigned nToys = 1000,
//                                          unsigned nStartsPerToy = 1,
//                                          const char* outFile = "resultsGivenMoments_chi2_amps_toys.root",
//                                          uint32_t seed = 0,
//                                          double epsR4 = 1.0) {
//   chi2_amp_fit::FitConfig cfg;
//   cfg.nStarts = nStartsPerToy;
//   cfg.randomSeed = seed;
//   cfg.epsR4 = epsR4;
//   cfg.recordEvery = 0;
//   cfg.observedInput.source = chi2_amp_fit::ObservedMomentsSource::kQ2Table;
//   cfg.observedInput.q2bin = q2bin;
//
//   chi2_amp_fit_toys::RunGivenMoments_Chi2Amps_Toys_Impl(cfg, nToys, outFile);
// }

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

// Parallel wrapper around your existing serial toy macro
void RunGivenMoments_Chi2Amps_Toys(double q2Value,
                                      unsigned nToys = 1000,
                                      unsigned nStartsPerToy = 1,
                                      const char* outFile = "toyFits.root",
                                      uint32_t seed = 0,
                                      double epsR4 = 1.0,
                                      unsigned nCores = 0)
{
  if (nToys == 0) {
    throw std::runtime_error("RunGivenMoments_Chi2Amps_Toys_MP: nToys must be > 0");
  }

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

    bool verbose = false;
    if (workerId == 0) verbose = true;
    // call your EXISTING serial toy function here
    RunGivenMoments_Chi2Amps_Setup(q2Value,
                                  myNToys,
                                  nStartsPerToy,
                                  partFile.c_str(),
                                  workerSeed,
                                  epsR4, verbose, nCores);

    return partFile;
  }, workerIds);

  TFileMerger merger(kTRUE, kFALSE);
  merger.OutputFile(outFile, "RECREATE");
  for (const auto& f : partFiles) merger.AddFile(f.c_str());

  if (!merger.Merge()) {
    throw std::runtime_error("RunGivenMoments_Chi2Amps_Toys_MP: merge failed");
  }

  for (const auto& f : partFiles) gSystem->Unlink(f.c_str());

  std::cout << "Merged " << partFiles.size() << " files into " << outFile << std::endl;
}
