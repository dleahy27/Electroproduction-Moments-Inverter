// Propagate Minuit's covariance into physical amplitudes, moments, and R.
// Minuit fits independent coordinates; the normalization constraint supplies
// one derived magnitude. Its Jacobian restores the full output covariance.
#include "Detail.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace emi::detail {

static double VarianceFromGradient(const std::vector<double>& gradient,
                                   const std::vector<double>& covariance) {
  // First-order error propagation: Var(f) = grad(f)^T Cov grad(f).
  const size_t n = gradient.size();
  if (covariance.size() != n * n)
    return std::numeric_limits<double>::quiet_NaN();
  double variance = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      variance += gradient[i] * covariance[i * n + j] * gradient[j];
    }
  }
  if (variance < 0.0 && variance > -1e-12)
    variance = 0.0;
  return (variance >= 0.0 && std::isfinite(variance))
             ? variance
             : std::numeric_limits<double>::quiet_NaN();
}

void FillHessianProducts(const EvaluationContext& ctx, const ROOT::Math::Minimizer& min,
                         const std::vector<double>& fullVals, bool covarianceIsUsable,
                         std::vector<double>& parErrors,
                         std::vector<double>& parCorrelations,
                         std::vector<double>& momentErrors, double& ratioR,
                         double& ratioRError) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::fill(parErrors.begin(), parErrors.end(), nan);
  std::fill(parCorrelations.begin(), parCorrelations.end(), nan);
  std::fill(momentErrors.begin(), momentErrors.end(), nan);
  ratioRError = nan;

  const size_t nFree = ctx.freeToFull.size();
  if (!covarianceIsUsable || nFree == 0)
    return;

  std::vector<double> freeCov(nFree * nFree, 0.0);
  for (size_t i = 0; i < nFree; ++i) {
    for (size_t j = 0; j < nFree; ++j) {
      freeCov[i * nFree + j] =
          min.CovMatrix(static_cast<unsigned>(i), static_cast<unsigned>(j));
    }
  }

  // Jacobian from the independent Minuit coordinates to the complete physical
  // amplitude basis.  This includes the amplitude fixed by normalisation.
  std::vector<double> fullJac(ctx.fullPars.size() * nFree, 0.0);
  for (size_t p = 0; p < ctx.fullPars.size(); ++p) {
    const int freeIndex = ctx.fullToFree[p];
    if (freeIndex >= 0)
      fullJac[p * nFree + static_cast<size_t>(freeIndex)] = 1.0;
    if (static_cast<int>(p) == ctx.normalisedMagFullIdx) {
      for (size_t i = 0; i < nFree; ++i) {
        fullJac[p * nFree + i] =
            NormalisedMagnitudeDerivative(ctx, fullVals, ctx.freeToFull[i]);
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
          value +=
              fullJac[p * nFree + i] * freeCov[i * nFree + j] * fullJac[q * nFree + j];
        }
      }
      fullCov[p * ctx.fullPars.size() + q] = value;
    }
  }

  for (size_t p = 0; p < ctx.fullPars.size(); ++p) {
    double variance = fullCov[p * ctx.fullPars.size() + p];
    if (variance < 0.0 && variance > -1e-12)
      variance = 0.0;
    parErrors[p] =
        (variance >= 0.0 && std::isfinite(variance)) ? std::sqrt(variance) : nan;
  }

  // Store only the upper triangle.  The analysis reconstructs the symmetric
  // matrix, while omitting the diagonal avoids hundreds of redundant branches.
  std::size_t pairIndex = 0;
  for (size_t first = 0; first < ctx.fullPars.size(); ++first) {
    for (size_t second = first + 1; second < ctx.fullPars.size(); ++second) {
      const double denominator = parErrors[first] * parErrors[second];
      if (denominator > 0.0 && std::isfinite(denominator)) {
        double correlation =
            fullCov[first * ctx.fullPars.size() + second] / denominator;
        if (std::isfinite(correlation) && std::abs(correlation) <= 1.0 + 1.0e-10) {
          parCorrelations[pairIndex] = std::clamp(correlation, -1.0, 1.0);
        }
      }
      ++pairIndex;
    }
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

  const auto h0It = ctx.modelIndexByName.find(MakeMomentName(ctx.cfg, 0, 0, 0, 0, 0));
  const auto h4It = ctx.modelIndexByName.find(MakeMomentName(ctx.cfg, 4, 0, 0, 0, 0));
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
        const double dh4 =
            ApplyNormalisationChainRule(ctx, fullVals, dFull2, fullIndex);
        gradient[i] = (dh4 * h0 - h4 * dh0) / (h0 * h0);
      }
      const double variance = VarianceFromGradient(gradient, freeCov);
      ratioRError = std::isfinite(variance) ? std::sqrt(variance) : nan;
    }
  }
}

} // namespace emi::detail
