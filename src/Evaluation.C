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

int SelectNormalisationIndex(const std::vector<Parameter>& pars) {
  const bool hasPositive = std::any_of(pars.begin(), pars.end(), [](const Parameter& p) {
    const auto label = ParseParameterLabel(p.name);
    return !p.phase && !p.fixed && label.valid &&
           label.reflectivity == 'a' && label.orientation == 'T';
  });
  const char selectedReflectivity = hasPositive ? 'a' : 'b';
  int bestIdx = -1;
  int bestL = -1;
  int bestM = -999999;

  for (int i = 0; i < static_cast<int>(pars.size()); ++i) {
    const auto& p = pars[static_cast<size_t>(i)];
    if (p.phase || p.fixed) continue;

    const auto label = ParseParameterLabel(p.name);
    if (!label.valid) continue;
    if (label.reflectivity != selectedReflectivity || label.orientation != 'T') continue;

    if (label.l > bestL || (label.l == bestL && label.m > bestM)) {
      bestIdx = i;
      bestL = label.l;
      bestM = label.m;
    }
  }
  return bestIdx;
}

} // namespace

void MarkAmplitudeNormalisationParameter(EvaluationContext& ctx) {
  if (ctx.cfg.photoproduction) return;

  const int idx = SelectNormalisationIndex(ctx.fullPars);
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

double NormalisationWeightForMagnitude(const EvaluationContext& ctx, int fullIdx) {
  if (fullIdx == ctx.normalisedMagFullIdx) return 0.0;

  const auto& p = ctx.fullPars[static_cast<size_t>(fullIdx)];
  if (p.phase) return 0.0;

  const auto label = ParseParameterLabel(p.name);
  if (!label.valid) return 0.0;

  // The normalisation condition is
  //   |a_T_norm|^2 = target/2 - (sum |T|^2 + epsilon sum |L|^2).
  if (label.orientation == 'T') return 1.0;
  if (label.orientation == 'L') {
    double w = ctx.cfg.epsilon;
    if (ctx.cfg.enforceLongitudinalParity && label.m > 0) w *= 2.0;
    return w;
  }
  return 0.0;
}

double EvalNormalisationSum(const EvaluationContext& ctx, const std::vector<double>& fullVals) {
  double sum = 0.0;
  for (int i = 0; i < static_cast<int>(ctx.fullPars.size()); ++i) {
    const double w = NormalisationWeightForMagnitude(ctx, i);
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
  const double norm2 = base - rest;
  if (!(norm2 >= 0.0) || !std::isfinite(norm2)) return false;

  fullVals[static_cast<size_t>(ctx.normalisedMagFullIdx)] = std::sqrt(norm2);
  return true;
}

double NormalisedMagnitudeDerivative(const EvaluationContext& ctx,
                                            const std::vector<double>& fullVals,
                                            int wrtFullIdx) {
  if (ctx.normalisedMagFullIdx < 0) return 0.0;

  const double w = NormalisationWeightForMagnitude(ctx, wrtFullIdx);
  if (w == 0.0) return 0.0;

  const double normMag = fullVals[static_cast<size_t>(ctx.normalisedMagFullIdx)];
  if (!(normMag > 0.0) || !std::isfinite(normMag)) return 0.0;

  // a_norm = sqrt(target/2 - sum_i w_i a_i^2)
  // d a_norm / d a_i = -w_i a_i / a_norm.
  return -w * fullVals[static_cast<size_t>(wrtFullIdx)] / normMag;
}

double ApplyNormalisationChainRule(const EvaluationContext& ctx,
                                          const std::vector<double>& fullVals,
                                          const std::vector<double>& dHdFull,
                                          int wrtFullIdx) {
  double dH = dHdFull[static_cast<size_t>(wrtFullIdx)];
  if (ctx.normalisedMagFullIdx >= 0) {
    const int normIdx = ctx.normalisedMagFullIdx;
    dH += dHdFull[static_cast<size_t>(normIdx)]
        * NormalisedMagnitudeDerivative(ctx, fullVals, wrtFullIdx);
  }
  return dH;
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
    if (t.ignorePhase && (mm.alpha==3 || mm.alpha==7 || mm.alpha==8)) continue; // No imaginary parts for these
    const double trig = t.ignorePhase ? 1.0 : ((t.trig == TrigKind::kCos) ? pairCos[t.phasePairIdx] : pairSin[t.phasePairIdx]);
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
    if (t.ignorePhase && (mm.alpha==3 || mm.alpha==7 || mm.alpha==8))
    {
       trig = 0.0;
       dtrig = 1.0;
    }else
    {
      trig = t.ignorePhase ? 1.0 : ((t.trig == TrigKind::kCos) ? pairCos[t.phasePairIdx] : pairSin[t.phasePairIdx]);
      dtrig = t.ignorePhase ? 0.0 : ((t.trig == TrigKind::kCos) ? -pairSin[t.phasePairIdx] : pairCos[t.phasePairIdx]);
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
