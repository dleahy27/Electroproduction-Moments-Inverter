#include "Detail.h"

#include "TRandom3.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace emi::detail {

namespace {

void FastSinCos(double angle, double& sine, double& cosine) {
#if defined(__GLIBC__) || defined(__APPLE__)
  ::sincos(angle, &sine, &cosine);
#else
  sine = std::sin(angle);
  cosine = std::cos(angle);
#endif
}

int SelectNormalisationIndex(const EvaluationContext& ctx) {
  const auto& pars = ctx.fullPars;
  const bool hasPositive = std::any_of(pars.begin(), pars.end(), [](const Parameter& p) {
    const auto label = ParseParameterLabel(p.name);
    return !p.phase && !p.fixed && label.valid &&
           label.reflectivity == 'a' && label.orientation == 'T';
  });
  const char selectedReflectivity = hasPositive ? 'a' : 'b';
  const int selectedK =
      ctx.cfg.nucleonPolarization == NucleonPolarization::None ? 0 : 1;
  int bestIdx = -1;
  int bestL = -1;
  int bestM = -999999;

  for (int i = 0; i < static_cast<int>(pars.size()); ++i) {
    const auto& p = pars[static_cast<size_t>(i)];
    if (p.phase || p.fixed) continue;

    const auto label = ParseParameterLabel(p.name);
    if (!label.valid) continue;
    if (label.reflectivity != selectedReflectivity ||
        label.orientation != 'T' || label.k != selectedK) continue;

    if (label.l > bestL || (label.l == bestL && label.m > bestM)) {
      bestIdx = i;
      bestL = label.l;
      bestM = label.m;
    }
  }
  return bestIdx;
}

} // namespace

void ApplyPolarizationGauge(EvaluationContext& ctx) {
  if (ctx.cfg.nucleonPolarization == NucleonPolarization::None) return;

  auto select = [&](char reflectivity, char orientation, int l, int m,
                    bool chooseHighest) {
    int magnitude = -1;
    int phase = -1;
    int bestL = -1;
    int bestM = -999999;
    for (int i = 0; i < static_cast<int>(ctx.fullPars.size()); ++i) {
      const auto label = ParseParameterLabel(ctx.fullPars[static_cast<size_t>(i)].name);
      if (!label.valid || label.k != 1 || label.reflectivity != reflectivity ||
          label.orientation != orientation ||
          (!chooseHighest && (label.l != l || label.m != m))) continue;
      if (chooseHighest &&
          (label.l < bestL || (label.l == bestL && label.m < bestM))) continue;
      if (chooseHighest && (label.l > bestL || label.m > bestM)) {
        magnitude = -1;
        phase = -1;
        bestL = label.l;
        bestM = label.m;
      }
      if (ctx.fullPars[static_cast<size_t>(i)].phase) phase = i;
      else magnitude = i;
    }
    return std::pair{magnitude, phase};
  };

  auto fixPhase = [&](int phaseIndex) {
    auto& phase = ctx.fullPars[static_cast<size_t>(phaseIndex)];
    phase.init = 0.0;
    phase.fixed = true;
    phase.low = 0.0;
    phase.high = 0.0;
    phase.step = 0.0;
  };

  const auto natural = select('a', 'T', 0, 0, true);
  if (natural.first < 0 || natural.second < 0) {
    throw std::runtime_error(
        "The polarized gauge requires a natural transverse non-flip wave");
  }
  fixPhase(natural.second);
  ctx.phaseReferenceMagFullIdx = natural.first;

  if (ctx.cfg.nucleonPolarization == NucleonPolarization::Both) return;

  std::pair<int, int> rotation;
  if (ctx.cfg.kMixingGauge) {
    const auto& gauge = *ctx.cfg.kMixingGauge;
    if ((gauge.reflectivity != 'a' && gauge.reflectivity != 'b') ||
        (gauge.orientation != 'T' && gauge.orientation != 'L') ||
        gauge.l < 0 || std::abs(gauge.m) > gauge.l) {
      throw std::invalid_argument("Invalid k-mixing gauge reference");
    }
    rotation = select(gauge.reflectivity, gauge.orientation,
                      gauge.l, gauge.m, false);
  } else {
    rotation = select('b', 'T', 0, 0, true);
  }
  if (rotation.first < 0 || rotation.second < 0) {
    throw std::runtime_error(
        ctx.cfg.kMixingGauge
            ? "The requested non-flip k-mixing gauge wave is not in the model"
            : "The polarized gauge requires an unnatural transverse non-flip wave");
  }
  if (rotation.first == natural.first) {
    throw std::runtime_error(
        "The phase and k-mixing gauges must use different reference waves");
  }

  fixPhase(rotation.second);
  auto& realCoordinate = ctx.fullPars[static_cast<size_t>(rotation.first)];
  realCoordinate.low = -ctx.cfg.magnitudeMax;
  realCoordinate.high = ctx.cfg.magnitudeMax;
  ctx.rotationGaugeMagFullIdx = rotation.first;
}

void MarkAmplitudeNormalisationParameter(EvaluationContext& ctx) {
  const int idx = SelectNormalisationIndex(ctx);
  if (idx < 0) {
    throw std::runtime_error("Could not find a transverse magnitude for amplitude normalisation");
  }

  auto& p = ctx.fullPars[static_cast<size_t>(idx)];
  p.init = 0.0;
  p.fixed = true;
  p.low = 0.0;
  p.high = ctx.cfg.magnitudeMax;
  p.step = 0.0;
  ctx.normalisedMagFullIdx = idx;
}

void SetAmplitudeNormalisationWeights(EvaluationContext& ctx) {
  ctx.normalisationWeights.assign(ctx.fullPars.size(), 0.0);
  auto addMoment = [&](int momentIndex, double scale) {
    if (momentIndex < 0) return;
    for (const auto& term :
         ctx.modelsRec[static_cast<size_t>(momentIndex)].terms) {
      if (term.idxMag1 == term.idxMag2 && term.trig == TrigKind::kCos) {
        ctx.normalisationWeights[static_cast<size_t>(term.idxMag1)] +=
            scale * term.coeff;
      }
    }
  };

  // H04_00 = H0_00 + epsilon H4_00 is diagonal in the amplitudes. Dividing
  // its coefficient by two retains the conventional target/2 equation while
  // also handling explicit wave lists that do not contain every +/-m partner.
  addMoment(ctx.idxH0_00, 0.5);
  addMoment(ctx.idxH4_00, 0.5 * ctx.cfg.epsilon);
  if (ctx.normalisedMagFullIdx < 0 ||
      !(ctx.normalisationWeights[
          static_cast<size_t>(ctx.normalisedMagFullIdx)] > 0.0)) {
    throw std::runtime_error(
        "The normalization amplitude does not contribute to H04_00");
  }
}

static double RawNormalisationWeight(const EvaluationContext& ctx, int fullIdx) {
  return ctx.normalisationWeights[static_cast<size_t>(fullIdx)];
}

double EvalNormalisationSum(const EvaluationContext& ctx, const std::vector<double>& fullVals) {
  double sum = 0.0;
  for (int i = 0; i < static_cast<int>(ctx.fullPars.size()); ++i) {
    if (i == ctx.normalisedMagFullIdx) continue;
    const double w = RawNormalisationWeight(ctx, i);
    if (w == 0.0) continue;
    const double v = fullVals[static_cast<size_t>(i)];
    sum += w * v * v;
  }
  return sum;
}

bool ApplyAmplitudeNormalisation(const EvaluationContext& ctx, std::vector<double>& fullVals) {
  if (ctx.normalisedMagFullIdx < 0) return true;

  const double base = 0.5 * ctx.cfg.normalisationMomentTarget;
  const double rest = EvalNormalisationSum(ctx, fullVals);
  const double normWeight = RawNormalisationWeight(ctx, ctx.normalisedMagFullIdx);
  if (!(normWeight > 0.0) || !std::isfinite(normWeight)) return false;
  const double norm2 = (base - rest) / normWeight;
  if (!(norm2 >= 0.0) || !std::isfinite(norm2)) return false;

  const double magnitude = std::sqrt(norm2);
  fullVals[static_cast<size_t>(ctx.normalisedMagFullIdx)] = magnitude;
  return true;
}

double NormalisedMagnitudeDerivative(const EvaluationContext& ctx,
                                            const std::vector<double>& fullVals,
                                            int wrtFullIdx) {
  if (ctx.normalisedMagFullIdx < 0) return 0.0;
  if (wrtFullIdx == ctx.normalisedMagFullIdx) return 0.0;

  const double w = RawNormalisationWeight(ctx, wrtFullIdx);
  if (w == 0.0) return 0.0;

  const double normMag = fullVals[static_cast<size_t>(ctx.normalisedMagFullIdx)];
  if (!(normMag > 0.0) || !std::isfinite(normMag)) return 0.0;

  const double normWeight = RawNormalisationWeight(ctx, ctx.normalisedMagFullIdx);
  if (!(normWeight > 0.0) || !std::isfinite(normWeight)) return 0.0;
  return -w * fullVals[static_cast<size_t>(wrtFullIdx)] /
         (normWeight * normMag);
}

double ApplyNormalisationChainRule(const EvaluationContext& ctx,
                                          const std::vector<double>& fullVals,
                                          const std::vector<double>& dHdFull,
                                          int wrtFullIdx) {
  double dH = dHdFull[static_cast<size_t>(wrtFullIdx)];
  if (ctx.normalisedMagFullIdx >= 0) {
    const int normIdx = ctx.normalisedMagFullIdx;
    const double dHdNorm = dHdFull[static_cast<size_t>(normIdx)];
    dH += dHdNorm * NormalisedMagnitudeDerivative(ctx, fullVals, wrtFullIdx);
  }
  return dH;
}

void EvaluateObservedMoments(const EvaluationContext& ctx,
                             const std::vector<double>& rawMoments,
                             std::vector<double>& values) {
  values.resize(ctx.observed.size());
  for (size_t i = 0; i < ctx.observed.size(); ++i) {
    if (ctx.observed[i].isMixed04) {
      const int idx0 = ctx.observedModelIdx0[i];
      const int idx4 = ctx.observedModelIdx4[i];
      if (idx0 < 0 || idx4 < 0) {
        throw std::runtime_error("RH04 moment is not mapped to H0/H4");
      }
      values[i] = rawMoments[static_cast<size_t>(idx0)] +
                  ctx.cfg.epsilon * rawMoments[static_cast<size_t>(idx4)];
    } else {
      const int index = ctx.observedModelIdx[i];
      if (index < 0) {
        throw std::runtime_error("Observed moment is not mapped to a model moment");
      }
      values[i] = rawMoments[static_cast<size_t>(index)];
    }
  }
}


bool FillFullParametersRaw(const EvaluationContext& ctx,
                                const double* x,
                                std::vector<double>& fullVals) {
  if (fullVals.size() != ctx.fullPars.size()) fullVals.resize(ctx.fullPars.size());
  for (size_t i = 0; i < ctx.fullPars.size(); ++i) fullVals[i] = ctx.fullPars[i].init;

  for (size_t i = 0; i < ctx.freeToFull.size(); ++i) {
    const int fullIdx = ctx.freeToFull[i];
    const auto& p = ctx.fullPars[static_cast<size_t>(fullIdx)];
    const double value = x[i];
    if (!std::isfinite(value)) return false;
    if (!p.phase && (value < p.low || value > p.high)) return false;
    fullVals[static_cast<size_t>(fullIdx)] = value;
  }
  return true;
}

bool FillFullParameters(const EvaluationContext& ctx,
                             const double* x,
                             std::vector<double>& fullVals) {
  if (!FillFullParametersRaw(ctx, x, fullVals)) return false;
  return ApplyAmplitudeNormalisation(ctx, fullVals);
}

static inline void EnsureSize(std::vector<double>& v, size_t n, double fill = 0.0) {
  if (v.size() != n) v.assign(n, fill);
}

void BuildPhasePairTrigCache(const EvaluationContext& ctx,
                                    const std::vector<double>& fullVals,
                                    std::vector<double>& pairSin,
                                    std::vector<double>& pairCos) {
  EnsureSize(pairSin, ctx.phasePairs.size());
  EnsureSize(pairCos, ctx.phasePairs.size());
  for (size_t i = 0; i < ctx.phasePairs.size(); ++i) {
    const auto& pp = ctx.phasePairs[i];
    FastSinCos(fullVals[pp.idxPhi1] - fullVals[pp.idxPhi2], pairSin[i], pairCos[i]);
  }
}

static double EvalMomentOnly(const MomentModel& mm,
                             const std::vector<double>& fullVals,
                             const std::vector<double>& pairSin,
                             const std::vector<double>& pairCos) {
  double H = 0.0;
  for (const auto& t : mm.terms) {
    const double trig = t.ignorePhase
                            ? (t.trig == TrigKind::kCos ? 1.0 : 0.0)
                            : (t.trig == TrigKind::kCos
                                   ? pairCos[t.phasePairIdx]
                                   : pairSin[t.phasePairIdx]);
    H += t.coeff * fullVals[t.idxMag1] * fullVals[t.idxMag2] * trig;
  }
  return H;
}

double EvaluateMomentAndDerivative(const EvaluationContext& ctx,
                                     const MomentModel& mm,
                                     const std::vector<double>& fullVals,
                                     const std::vector<double>& pairSin,
                                     const std::vector<double>& pairCos,
                                     std::vector<double>& dHdFull) {
  EnsureSize(dHdFull, ctx.fullPars.size());
  std::fill(dHdFull.begin(), dHdFull.end(), 0.0);
  double H = 0.0;
  for (const auto& t : mm.terms) {
    const auto& pp = ctx.phasePairs[t.phasePairIdx];
    const double m1 = fullVals[t.idxMag1]; // Magnitudes
    const double m2 = fullVals[t.idxMag2];

    // const double dtrig = (t.trig == TrigKind::kCos) ? -pairSin[t.phasePairIdx] : pairCos[t.phasePairIdx];

    double trig;
    double dtrig;
    if (t.ignorePhase) {
      trig = t.trig == TrigKind::kCos ? 1.0 : 0.0;
      dtrig = 0.0;
    } else {
      trig = t.trig == TrigKind::kCos ? pairCos[t.phasePairIdx]
                                      : pairSin[t.phasePairIdx];
      dtrig = t.trig == TrigKind::kCos ? -pairSin[t.phasePairIdx]
                                       : pairCos[t.phasePairIdx];
    }
    const double val = t.coeff * m1 * m2 * trig;
    H += val;
    dHdFull[t.idxMag1] += t.coeff * m2 * trig;
    dHdFull[t.idxMag2] += t.coeff * m1 * trig;
    const double common = t.coeff * m1 * m2 * dtrig;
    dHdFull[pp.idxPhi1] += common;
    dHdFull[pp.idxPhi2] -= common;
  }
  return H;
}

void EvaluateAllMoments(const EvaluationContext& ctx,
                           const std::vector<double>& fullVals,
                           std::vector<double>& values) {
  EnsureSize(values, ctx.modelsRec.size());
  std::vector<double> pairSin, pairCos;
  BuildPhasePairTrigCache(ctx, fullVals, pairSin, pairCos);
  for (size_t i = 0; i < ctx.modelsRec.size(); ++i) values[i] = EvalMomentOnly(ctx.modelsRec[i], fullVals, pairSin, pairCos);
}

class Chi2Function final : public ROOT::Math::IMultiGradFunction {
public:
  explicit Chi2Function(std::shared_ptr<EvaluationContext> ctx) : ctx_(std::move(ctx)) {
    if (!ctx_) throw std::runtime_error("Chi2Function: null context");
  }

  unsigned int NDim() const override { return static_cast<unsigned>(ctx_->freeToFull.size()); }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override { return new Chi2Function(ctx_); }

  double DoEval(const double* x) const override {
    if (!FillFullParameters(*ctx_, x, fullVals_)) return 1e300;
    BuildPhasePairTrigCache(*ctx_, fullVals_, pairSin_, pairCos_);
    ++ctx_->callCount;

    EnsureSize(buf_momRec_, ctx_->modelsRec.size());
    if (buf_momSeen_.size() != ctx_->modelsRec.size()) buf_momSeen_.assign(ctx_->modelsRec.size(), 0);
    std::fill(buf_momSeen_.begin(), buf_momSeen_.end(), static_cast<unsigned char>(0));

    auto getMomentRaw = [&](int idx) -> double {
      if (idx < 0) throw std::runtime_error("Observed moment is not mapped to a model moment");
      const size_t uidx = static_cast<size_t>(idx);
      if (!buf_momSeen_[uidx]) {
        buf_momRec_[uidx] = EvalMomentOnly(ctx_->modelsRec[uidx], fullVals_, pairSin_, pairCos_);
        buf_momSeen_[uidx] = 1;
      }
      return buf_momRec_[uidx];
    };

    double chi2 = 0.0;
    for (size_t i = 0; i < ctx_->observed.size(); ++i) {
      const auto& ob = ctx_->observed[i];
      double model = 0.0;
      if (ob.isMixed04) {
        const double H0 = getMomentRaw(ctx_->observedModelIdx0[i]);
        const double H4 = getMomentRaw(ctx_->observedModelIdx4[i]);
        model = H0 + ctx_->cfg.epsilon * H4;
      } else {
        model = getMomentRaw(ctx_->observedModelIdx[i]);
      }

      const double sigma = ob.sigma;
      const double r = (ob.value - model) / sigma;
      chi2 += r * r;
    }

    return chi2;
  }

  void Gradient(const double* x, double* grad) const override {
    std::fill(grad, grad + NDim(), 0.0);
    if (!FillFullParameters(*ctx_, x, fullVals_)) return;
    BuildPhasePairTrigCache(*ctx_, fullVals_, pairSin_, pairCos_);

    const size_t nFull = ctx_->fullPars.size();
    EnsureSize(buf_fullA_, nFull);
    EnsureSize(buf_fullB_, nFull);

    for (size_t iobs = 0; iobs < ctx_->observed.size(); ++iobs) {
      const auto& ob = ctx_->observed[iobs];
      const double sigma = ob.sigma;
      double model = 0.0;

      if (ob.isMixed04) {
        const int idx0 = ctx_->observedModelIdx0[iobs];
        const int idx4 = ctx_->observedModelIdx4[iobs];
        if (idx0 < 0 || idx4 < 0) throw std::runtime_error("RH04 moment is not mapped to H0/H4 model moments");

        const double H0 = EvaluateMomentAndDerivative(*ctx_, ctx_->modelsRec[idx0], fullVals_, pairSin_, pairCos_, buf_fullA_);
        const double H4 = EvaluateMomentAndDerivative(*ctx_, ctx_->modelsRec[idx4], fullVals_, pairSin_, pairCos_, buf_fullB_);
        model = H0 + ctx_->cfg.epsilon * H4;
        const double pull = (ob.value - model) / sigma;

        for (unsigned i = 0; i < NDim(); ++i) {
          const int fullIdx = ctx_->freeToFull[i];
          const double dH0 = ApplyNormalisationChainRule(*ctx_, fullVals_, buf_fullA_, fullIdx);
          const double dH4 = ApplyNormalisationChainRule(*ctx_, fullVals_, buf_fullB_, fullIdx);
          const double dModel = dH0 + ctx_->cfg.epsilon * dH4;
          grad[i] += -2.0 * pull * dModel / sigma;
        }
      } else {
        const int midx = ctx_->observedModelIdx[iobs];
        if (midx < 0) throw std::runtime_error("Observed moment is not mapped to a model moment");

        model = EvaluateMomentAndDerivative(*ctx_, ctx_->modelsRec[midx], fullVals_, pairSin_, pairCos_, buf_fullA_);
        const double pull = (ob.value - model) / sigma;

        for (unsigned i = 0; i < NDim(); ++i) {
          const int fullIdx = ctx_->freeToFull[i];
          const double dModel = ApplyNormalisationChainRule(*ctx_, fullVals_, buf_fullA_, fullIdx);
          grad[i] += -2.0 * pull * dModel / sigma;
        }
      }
    }
  }

  double DoDerivative(const double* x, unsigned int icoord) const override {
    EnsureSize(buf_gradTmp_, NDim());
    Gradient(x, buf_gradTmp_.data());
    return (icoord < NDim()) ? buf_gradTmp_[icoord] : 0.0;
  }

private:
  std::shared_ptr<EvaluationContext> ctx_;
  mutable std::vector<double> fullVals_;
  mutable std::vector<double> pairSin_;
  mutable std::vector<double> pairCos_;
  mutable std::vector<double> buf_momRec_;
  mutable std::vector<unsigned char> buf_momSeen_;
  mutable std::vector<double> buf_fullA_, buf_fullB_;
  mutable std::vector<double> buf_gradTmp_;
};

class Chi2FunctionNoGrad final : public ROOT::Math::IBaseFunctionMultiDim {
public:
  explicit Chi2FunctionNoGrad(std::shared_ptr<EvaluationContext> ctx) : fcn_(std::move(ctx)) {}

  unsigned int NDim() const override { return fcn_.NDim(); }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override { return new Chi2FunctionNoGrad(*this); }
  double DoEval(const double* x) const override { return fcn_.DoEval(x); }

private:
  Chi2Function fcn_;
};

std::unique_ptr<ROOT::Math::IBaseFunctionMultiDim>
MakeNumericalChi2(std::shared_ptr<EvaluationContext> context) {
  return std::make_unique<Chi2FunctionNoGrad>(std::move(context));
}

std::unique_ptr<ROOT::Math::IMultiGradFunction>
MakeAnalyticChi2(std::shared_ptr<EvaluationContext> context) {
  return std::make_unique<Chi2Function>(std::move(context));
}

double EvaluateChi2(const ROOT::Math::IBaseFunctionMultiDim& function,
                    const std::vector<double>& point) {
  if (point.empty()) return 1e300;
  const double chi2 = function(point.data());
  return std::isfinite(chi2) ? chi2 : 1e300;
}

} // namespace emi::detail
