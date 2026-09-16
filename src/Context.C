// Build the immutable context used by every objective-function evaluation.
// Expensive index maps and covariance factors are prepared once here rather
// than recomputed at each Minuit step.
#include "Detail.h"

#include "Math/Minimizer.h"
#include "TRandom3.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace emi::detail {

InternalConfig MakeInternalConfig(const FitConfig& fit, const ModelConfig& model) {
  // Public configuration is copied into a compact internal form so worker
  // processes do not depend on the CLI or on Cling-owned objects.
  InternalConfig config;
  config.waves = model.waves;
  config.usePositiveReflectivity = model.usePositiveReflectivity;
  config.useNegativeReflectivity = model.useNegativeReflectivity;
  config.enforceLongitudinalParity = model.enforceLongitudinalParity;
  config.nucleonPolarization = model.nucleonPolarization;
  config.kMixingGauge = model.kMixingGauge;
  if (config.kMixingGauge && config.enforceLongitudinalParity &&
      config.kMixingGauge->orientation == 'L') {
    config.kMixingGauge->m = std::abs(config.kMixingGauge->m);
  }
  config.normalisationMomentTarget = model.normalisationMoment;
  config.epsilon = fit.epsilon;
  config.photoproduction = fit.photoproduction;
  config.starts = fit.starts;
  config.seed = fit.seed;
  config.runHesse = fit.hesse;
  config.useNumericalGradients = fit.useNumericalGradients;
  config.verbose = fit.verbose;
  config.momentsFile = fit.input.string();
  config.momentsTree = fit.tree;
  config.bin = fit.bin;
  return config;
}

std::shared_ptr<EvaluationContext> BuildContext(const InternalConfig& cfg) {
  auto ctx = std::make_shared<EvaluationContext>();
  ctx->cfg = cfg;
  ctx->fullPars = BuildParameters(cfg);
  ApplyPolarizationGauge(*ctx);
  MarkAmplitudeNormalisationParameter(*ctx);

  std::unordered_map<long long, int> paramIndex;
  // Moment terms refer to parameters by integer index. Encoding the quantum
  // labels once here avoids repeated string searches in the objective.
  paramIndex.reserve(ctx->fullPars.size() * 2);
  for (int i = 0; i < static_cast<int>(ctx->fullPars.size()); ++i) {
    const auto& p = ctx->fullPars[static_cast<size_t>(i)];
    const auto label = ParseParameterLabel(p.name);
    if (!label.valid) throw std::runtime_error("Could not parse parameter label: " + p.name);
    paramIndex.emplace(MakeParameterKey(label.reflectivity, label.orientation,
                                        label.k, label.l, label.m, p.phase), i);
  }

  if (cfg.verbose) {
    std::cout << "Production mode: "
              << (cfg.photoproduction
                      ? "photoproduction (alpha <= 3, no longitudinal amplitudes)"
                      : "electroproduction/full")
              << std::endl;
    std::cout << "Longitudinal negative-m waves are folded onto m >= 0 using parity." << std::endl;
    if (cfg.nucleonPolarization != NucleonPolarization::None) {
      const char* mode =
          cfg.nucleonPolarization == NucleonPolarization::Initial
              ? "initial"
              : cfg.nucleonPolarization == NucleonPolarization::Recoil
                    ? "recoil" : "initial and recoil";
      std::cout << "Nucleon polarization: " << mode
                << " (explicit k=+/-1 amplitudes)." << std::endl;
      const auto& phaseReference =
          ctx->fullPars[static_cast<size_t>(ctx->phaseReferenceMagFullIdx)];
      std::cout << "Overall phase fixed by " << phaseReference.name
                << " real and nonnegative." << std::endl;
      if (ctx->rotationGaugeMagFullIdx >= 0) {
        const auto& rotationReference =
            ctx->fullPars[static_cast<size_t>(ctx->rotationGaugeMagFullIdx)];
        std::cout << "The k-basis mixing angle is fixed by Im("
                  << rotationReference.name << ") = 0; its stored magnitude "
                  << "is a signed real coordinate." << std::endl;
      }
    }
  }

  ctx->fullToFree.assign(ctx->fullPars.size(), -1);
  // `fullPars` describes the physical output convention. Minuit sees only the
  // unfixed coordinates, so these two maps translate in both directions.
  for (int i = 0; i < static_cast<int>(ctx->fullPars.size()); ++i) {
    if (ctx->fullPars[static_cast<size_t>(i)].fixed) continue;
    ctx->fullToFree[static_cast<size_t>(i)] = static_cast<int>(ctx->freeToFull.size());
    ctx->freeToFull.push_back(i);
  }

  const auto inputMoments = ReadObservedMoments(cfg);
  if (inputMoments.empty()) {
    throw std::runtime_error("No requested moment branches were found in the input ROOT file");
  }

  ctx->modelsRec = BuildMomentModels(cfg, paramIndex, {}, ctx->phasePairs);
  ctx->modelIndexByName.reserve(ctx->modelsRec.size() * 2);
  for (size_t i = 0; i < ctx->modelsRec.size(); ++i) ctx->modelIndexByName.emplace(ctx->modelsRec[i].name, i);

  const std::string h0Name = MakeMomentName(cfg, 0, 0, 0, 0, 0);
  const std::string h4Name = MakeMomentName(cfg, 4, 0, 0, 0, 0);
  auto it0 = ctx->modelIndexByName.find(h0Name);
  auto it4 = ctx->modelIndexByName.find(h4Name);
  if (it0 != ctx->modelIndexByName.end()) ctx->idxH0_00 = static_cast<int>(it0->second);
  if (it4 != ctx->modelIndexByName.end()) ctx->idxH4_00 = static_cast<int>(it4->second);
  SetAmplitudeNormalisationWeights(*ctx);

  for (const auto& ob : inputMoments) {
    // H04 is the experimentally inseparable H0 + epsilon H4 response. Store
    // both component indices so it can be evaluated without inventing a model
    // parameter for the combined quantity.
    int modelIndex = -1;
    int modelIndex0 = -1;
    int modelIndex4 = -1;
    if (ob.isMixed04) {
      const std::string key0 = MakeMomentName(cfg, 0, ob.beta, ob.delta, ob.L, ob.M);
      const std::string key4 = MakeMomentName(cfg, 4, ob.beta, ob.delta, ob.L, ob.M);
      auto itObs0 = ctx->modelIndexByName.find(key0);
      auto itObs4 = ctx->modelIndexByName.find(key4);
      if (itObs0 != ctx->modelIndexByName.end()) modelIndex0 = static_cast<int>(itObs0->second);
      if (itObs4 != ctx->modelIndexByName.end()) modelIndex4 = static_cast<int>(itObs4->second);
    } else {
      const std::string key = MakeMomentName(cfg, ob.alpha, ob.beta, ob.delta,
                                             ob.L, ob.M);
      auto it = ctx->modelIndexByName.find(key);
      if (it != ctx->modelIndexByName.end()) modelIndex = static_cast<int>(it->second);
    }

    const bool constructible = ob.isMixed04
        ? modelIndex0 >= 0 && modelIndex4 >= 0
        : modelIndex >= 0;
    if (!constructible) {
      if (cfg.verbose) std::cout << "Skipping unconstructible moment " << ob.name << '\n';
      continue;
    }
    ctx->observed.push_back(ob);
    ctx->observedModelIdx.push_back(modelIndex);
    ctx->observedModelIdx0.push_back(modelIndex0);
    ctx->observedModelIdx4.push_back(modelIndex4);
  }
  if (ctx->observed.empty()) {
    throw std::runtime_error("None of the input moments can be constructed from the selected waves");
  }

  for (const auto& moment : ctx->modelsRec) {
    // Output branches are a reader-friendly view of the complete internal
    // model. Electroproduction writes H0, H4, and H04; photoproduction has no
    // longitudinal H4 term and therefore writes only H04.
    if (moment.alpha == 0) {
      std::string suffix;
      if (cfg.nucleonPolarization != NucleonPolarization::None) {
        suffix = std::to_string(moment.beta) + "_" + std::to_string(moment.delta) + "_";
      }
      suffix += std::to_string(moment.L) + "_" + std::to_string(moment.M);
      int second = -1;
      if (!cfg.photoproduction) {
        const auto h4 = ctx->modelIndexByName.find(
            MakeMomentName(cfg, 4, moment.beta, moment.delta, moment.L, moment.M));
        if (h4 == ctx->modelIndexByName.end()) continue;
        second = static_cast<int>(h4->second);
      }
      ctx->outputMoments.push_back({"H04_" + suffix,
                                    static_cast<int>(ctx->modelIndexByName.at(moment.name)),
                                    second, cfg.epsilon});
      if (!cfg.photoproduction) {
        ctx->outputMoments.push_back({
            moment.name,
            static_cast<int>(ctx->modelIndexByName.at(moment.name)),
            -1, 0.0});
        const auto& h4Moment = ctx->modelsRec[static_cast<size_t>(second)];
        ctx->outputMoments.push_back({h4Moment.name, second, -1, 0.0});
      }
    }
  }
  for (const auto& moment : ctx->modelsRec) {
    if (moment.alpha == 0 || moment.alpha == 4) continue;
    ctx->outputMoments.push_back({moment.name,
                                  static_cast<int>(ctx->modelIndexByName.at(moment.name)),
                                  -1, 0.0});
  }

  return ctx;
}

void MakeParameterBranches(TTree* t, const std::vector<Parameter>& pars, std::vector<double>& storage) {
  storage.assign(pars.size(), 0.0);
  for (size_t i = 0; i < pars.size(); ++i) t->Branch(pars[i].name.c_str(), &storage[i]);
}

void MakeNamedBranches(TTree* t, const std::vector<std::string>& names,
                       std::vector<double>& storage) {
  storage.assign(names.size(), 0.0);
  for (size_t i = 0; i < names.size(); ++i) {
    t->Branch(names[i].c_str(), &storage[i]);
  }
}

void MakeOutputMomentBranches(TTree* t, const std::vector<OutputMoment>& models,
                              std::vector<double>& storage) {
  storage.assign(models.size(), 0.0);
  for (size_t i = 0; i < models.size(); ++i) t->Branch(models[i].name.c_str(), &storage[i]);
}

void EvaluateOutputMoments(const EvaluationContext& ctx,
                           const std::vector<double>& rawMoments,
                           std::vector<double>& values) {
  values.assign(ctx.outputMoments.size(), 0.0);
  for (size_t i = 0; i < ctx.outputMoments.size(); ++i) {
    const auto& output = ctx.outputMoments[i];
    values[i] = rawMoments[static_cast<size_t>(output.firstModel)];
    if (output.secondModel >= 0) {
      values[i] += output.secondScale *
                   rawMoments[static_cast<size_t>(output.secondModel)];
    }
  }
}

static double VarianceFromGradient(const std::vector<double>& gradient,
                                   const std::vector<double>& covariance) {
  // First-order error propagation: Var(f) = grad(f)^T Cov grad(f).
  const size_t n = gradient.size();
  if (covariance.size() != n * n) return std::numeric_limits<double>::quiet_NaN();
  double variance = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      variance += gradient[i] * covariance[i * n + j] * gradient[j];
    }
  }
  if (variance < 0.0 && variance > -1e-12) variance = 0.0;
  return (variance >= 0.0 && std::isfinite(variance))
           ? variance
           : std::numeric_limits<double>::quiet_NaN();
}

void FillHessianProducts(const EvaluationContext& ctx,
                                const ROOT::Math::Minimizer& min,
                                const std::vector<double>& fullVals,
                                bool covarianceIsUsable,
                                std::vector<double>& parErrors,
                                std::vector<double>& momentErrors,
                                double& ratioR,
                                double& ratioRError) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::fill(parErrors.begin(), parErrors.end(), nan);
  std::fill(momentErrors.begin(), momentErrors.end(), nan);
  ratioRError = nan;

  const size_t nFree = ctx.freeToFull.size();
  if (!covarianceIsUsable || nFree == 0) return;

  std::vector<double> freeCov(nFree * nFree, 0.0);
  for (size_t i = 0; i < nFree; ++i) {
    for (size_t j = 0; j < nFree; ++j) {
      freeCov[i * nFree + j] = min.CovMatrix(static_cast<unsigned>(i), static_cast<unsigned>(j));
    }
  }

  // Jacobian from the independent Minuit coordinates to the complete physical
  // amplitude basis.  This includes the amplitude fixed by normalisation.
  std::vector<double> fullJac(ctx.fullPars.size() * nFree, 0.0);
  for (size_t p = 0; p < ctx.fullPars.size(); ++p) {
    const int freeIndex = ctx.fullToFree[p];
    if (freeIndex >= 0) fullJac[p * nFree + static_cast<size_t>(freeIndex)] = 1.0;
    if (static_cast<int>(p) == ctx.normalisedMagFullIdx) {
      for (size_t i = 0; i < nFree; ++i) {
        fullJac[p * nFree + i] = NormalisedMagnitudeDerivative(
            ctx, fullVals, ctx.freeToFull[i]);
      }
    }
  }

  std::vector<double> fullCov(ctx.fullPars.size() * ctx.fullPars.size(), 0.0);
  // Transform the covariance with J C J^T. This restores uncertainties for
  // derived coordinates, especially the magnitude eliminated by H00
  // normalisation, before the result is written to the ROOT tree.
  for (size_t p = 0; p < ctx.fullPars.size(); ++p) {
    for (size_t q = 0; q < ctx.fullPars.size(); ++q) {
      double value = 0.0;
      for (size_t i = 0; i < nFree; ++i) {
        for (size_t j = 0; j < nFree; ++j) {
          value += fullJac[p * nFree + i] * freeCov[i * nFree + j]
                 * fullJac[q * nFree + j];
        }
      }
      fullCov[p * ctx.fullPars.size() + q] = value;
    }
  }

  for (size_t p = 0; p < ctx.fullPars.size(); ++p) {
    double variance = fullCov[p * ctx.fullPars.size() + p];
    if (variance < 0.0 && variance > -1e-12) variance = 0.0;
    parErrors[p] = (variance >= 0.0 && std::isfinite(variance)) ? std::sqrt(variance) : nan;
  }
  std::vector<double> pairSin, pairCos;
  BuildPhasePairTrigCache(ctx, fullVals, pairSin, pairCos);
  std::vector<double> dFull, dFull2, gradient(nFree, 0.0), gradient2(nFree, 0.0);

  auto momentGradient = [&](int modelIndex, std::vector<double>& out) {
    // Moment derivatives are initially with respect to the full amplitudes;
    // the chain rule folds the derived normalisation magnitude back onto the
    // independent Minuit coordinates.
    EvaluateMomentAndDerivative(ctx, ctx.modelsRec[static_cast<size_t>(modelIndex)],
                          fullVals, pairSin, pairCos, dFull);
    out.assign(nFree, 0.0);
    for (size_t i = 0; i < nFree; ++i) {
      out[i] = ApplyNormalisationChainRule(ctx, fullVals, dFull, ctx.freeToFull[i]);
    }
  };

  for (size_t m = 0; m < ctx.outputMoments.size(); ++m) {
    const auto& output = ctx.outputMoments[m];
    momentGradient(output.firstModel, gradient);
    if (output.secondModel >= 0) {
      momentGradient(output.secondModel, gradient2);
      for (size_t i = 0; i < nFree; ++i) {
        gradient[i] += output.secondScale * gradient2[i];
      }
    }
    const double variance = VarianceFromGradient(gradient, freeCov);
    momentErrors[m] = std::isfinite(variance) ? std::sqrt(variance) : nan;
  }

  const auto h0It = ctx.modelIndexByName.find(
      MakeMomentName(ctx.cfg, 0, 0, 0, 0, 0));
  const auto h4It = ctx.modelIndexByName.find(
      MakeMomentName(ctx.cfg, 4, 0, 0, 0, 0));
  if (h0It != ctx.modelIndexByName.end() && h4It != ctx.modelIndexByName.end()) {
    // R = H4/H0, so its gradient uses the ordinary quotient rule.
    const double h0 = EvaluateMomentAndDerivative(ctx, ctx.modelsRec[h0It->second],
                                             fullVals, pairSin, pairCos, dFull);
    const double h4 = EvaluateMomentAndDerivative(ctx, ctx.modelsRec[h4It->second],
                                             fullVals, pairSin, pairCos, dFull2);
    if (std::abs(h0) > 1e-15) {
      ratioR = h4 / h0;
      for (size_t i = 0; i < nFree; ++i) {
        const int fullIndex = ctx.freeToFull[i];
        const double dh0 = ApplyNormalisationChainRule(ctx, fullVals, dFull, fullIndex);
        const double dh4 = ApplyNormalisationChainRule(ctx, fullVals, dFull2, fullIndex);
        gradient[i] = (dh4 * h0 - h4 * dh0) / (h0 * h0);
      }
      const double variance = VarianceFromGradient(gradient, freeCov);
      ratioRError = std::isfinite(variance) ? std::sqrt(variance) : nan;
    }
  }
}

void BuildRandomStart(const EvaluationContext& ctx,
                                  TRandom3& rng,
                                  std::vector<double>& xStart) {
  xStart.assign(ctx.freeToFull.size(), 0.0);

  for (unsigned i = 0; i < ctx.freeToFull.size(); ++i) {
    const int fullIdx = ctx.freeToFull[i];
    const auto& p = ctx.fullPars[static_cast<size_t>(fullIdx)];
    if (p.phase) {
      // Phases have no preferred direction, hence a uniform draw across their
      // full allowed interval. Magnitudes use the configured truncated normal.
      xStart[i] = rng.Uniform(p.low, p.high);
    } else {
      double value = rng.Gaus(ctx.cfg.magnitudeStartMean, ctx.cfg.magnitudeStartSigma);
      while (value < p.low || value > p.high) {
        value = rng.Gaus(ctx.cfg.magnitudeStartMean, ctx.cfg.magnitudeStartSigma);
      }
      xStart[i] = value;
    }
  }

  if (ctx.normalisedMagFullIdx >= 0) {
    std::vector<double> fullVals;
    if (!FillFullParameters(ctx, xStart.data(), fullVals)) {
      // If the random point gives an invalid normalisation square-root, shrink
      // all free magnitudes together and then calculate the derived a_T wave.
      // This enforces sum |T|^2 + epsilon sum |L|^2 < target/2.
      FillFullParametersRaw(ctx, xStart.data(), fullVals);
      const double base = 0.5 * ctx.cfg.normalisationMomentTarget;
      const double sum = EvalNormalisationSum(ctx, fullVals);
      if (sum > 0.0) {
        const double scale = std::sqrt(0.95 * base / sum);
        for (unsigned i = 0; i < ctx.freeToFull.size(); ++i) {
          const int fullIdx = ctx.freeToFull[i];
          const auto& p = ctx.fullPars[static_cast<size_t>(fullIdx)];
          if (!p.phase) xStart[i] *= scale;
        }
      }
    }
  }
}

} // namespace emi::detail
