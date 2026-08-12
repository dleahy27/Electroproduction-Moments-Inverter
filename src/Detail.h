#pragma once

#include "emi/Config.h"

#include "Math/IFunction.h"
#include "Math/Minimizer.h"

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

class TTree;
class TRandom3;

namespace emi::detail {

inline constexpr double kPi = 3.1415926535897932384626433832795;
inline constexpr double kTwoPi = 2.0 * kPi;
inline constexpr double kSqrt2 = 1.4142135623730950488016887242097;
inline constexpr double kInvSqrt2 = 1.0 / kSqrt2;
inline constexpr double kSqrt6Over5 = 0.48989794855663561963945681494118;
inline constexpr double kSqrt12Over5 = 0.69282032302755092063339055356909;
inline constexpr int kMaximumInputBins = 18;

struct ObservedMoment {
  int alpha = 0;
  int L = 0;
  int M = 0;
  double value = 0.0;
  double sigma = 1.0;
  bool isMixed04 = false;
  std::string name;
};

struct Parameter {
  std::string name;
  double init = 0.0;
  double step = 1e-3;
  double low = -1.0;
  double high = 1.0;
  bool fixed = false;
  bool phase = false;
};

struct ParameterLabel {
  char reflectivity = '\0';
  char orientation = '\0';
  int l = -1;
  int m = 0;
  bool phase = false;
  bool valid = false;
};

struct InternalConfig {
  std::vector<Wave> waves;
  bool usePositiveReflectivity = true;
  bool useNegativeReflectivity = true;
  bool enforceLongitudinalParity = true;
  double epsilon = 1.0;
  bool photoproduction = false;
  unsigned starts = 10000;
  unsigned maxCalls = 1000000;
  unsigned maxIterations = 1000000;
  double tolerance = 1e-8;
  int strategy = 2;
  int printLevel = 0;
  bool runHesse = true;
  bool useNumericalGradients = false;
  std::uint32_t seed = 0;
  bool verbose = true;
  double magnitudeStartMean = 0.5;
  double magnitudeStartSigma = 0.5;
  double magnitudeMin = 0.0;
  double magnitudeMax = 5.0;
  double normalisationMomentTarget = 2.0;
  std::string momentsFile;
  std::string momentsTree;
  int bin = 0;
};

enum class TrigKind : std::uint8_t { kCos, kSin };

struct PhasePair {
  int idxPhi1 = -1;
  int idxPhi2 = -1;
};

struct Term {
  double coeff = 0.0;
  int idxMag1 = -1;
  int idxMag2 = -1;
  int phasePairIdx = -1;
  TrigKind trig = TrigKind::kCos;
  bool ignorePhase = false;
};

struct MomentModel {
  int alpha = 0;
  int L = 0;
  int M = 0;
  std::string name;
  std::vector<Term> terms;
};

struct OutputMoment {
  std::string name;
  int firstModel = -1;
  int secondModel = -1;
  double secondScale = 0.0;
};

struct EvaluationContext {
  InternalConfig cfg;
  std::vector<Parameter> fullPars;
  std::vector<int> freeToFull;
  std::vector<int> fullToFree;
  std::vector<ObservedMoment> observed;
  std::vector<MomentModel> modelsRec;
  std::vector<OutputMoment> outputMoments;
  std::unordered_map<std::string, std::size_t> modelIndexByName;
  std::vector<int> observedModelIdx;
  std::vector<int> observedModelIdx0;
  std::vector<int> observedModelIdx4;
  std::vector<PhasePair> phasePairs;
  int idxH0_00 = -1;
  int idxH4_00 = -1;
  int normalisedMagFullIdx = -1;
  mutable unsigned callCount = 0;
};

struct MinimizerResult {
  bool valid = false;
  bool fitOk = false;
  bool hesseOk = false;
  int status = -999;
  int covarianceStatus = -1;
  double edm = 0.0;
  double chi2 = 1e300;
  std::vector<double> freeValues;
  std::unique_ptr<ROOT::Math::Minimizer> minimizer;
};

InternalConfig MakeInternalConfig(const FitConfig& fit, const ModelConfig& model);
std::vector<ObservedMoment> ReadObservedMoments(const InternalConfig& config);
ParameterLabel ParseParameterLabel(const std::string& name);
long long MakeParameterKey(char reflectivity, char orientation, int l, int m, bool phase);
std::vector<Parameter> BuildParameters(const InternalConfig& config);
std::vector<MomentModel> BuildMomentModels(
    const InternalConfig& config,
    const std::unordered_map<long long, int>& parameterIndex,
    const std::unordered_set<std::string>& neededMoments,
    std::vector<PhasePair>& phasePairs);
std::shared_ptr<EvaluationContext> BuildContext(const InternalConfig& config);

void MarkAmplitudeNormalisationParameter(EvaluationContext& context);
double EvalNormalisationSum(const EvaluationContext& context,
                            const std::vector<double>& fullValues);
bool FillFullParametersRaw(const EvaluationContext& context, const double* freeValues,
                           std::vector<double>& fullValues);
bool FillFullParameters(const EvaluationContext& context, const double* freeValues,
                        std::vector<double>& fullValues);
double NormalisedMagnitudeDerivative(const EvaluationContext& context,
                                     const std::vector<double>& fullValues,
                                     int parameterIndex);
double ApplyNormalisationChainRule(const EvaluationContext& context,
                                   const std::vector<double>& fullValues,
                                   const std::vector<double>& fullDerivative,
                                   int parameterIndex);
double EvaluateMomentAndDerivative(const EvaluationContext& context,
                                   const MomentModel& moment,
                                   const std::vector<double>& fullValues,
                                   const std::vector<double>& pairSin,
                                   const std::vector<double>& pairCos,
                                   std::vector<double>& derivative);
void EvaluateAllMoments(const EvaluationContext& context,
                        const std::vector<double>& fullValues,
                        std::vector<double>& values);
void EvaluateOutputMoments(const EvaluationContext& context,
                           const std::vector<double>& rawMoments,
                           std::vector<double>& values);
void BuildPhasePairTrigCache(const EvaluationContext& context,
                             const std::vector<double>& fullValues,
                             std::vector<double>& pairSin,
                             std::vector<double>& pairCos);
void BuildRandomStart(const EvaluationContext& context, TRandom3& random,
                      std::vector<double>& values);

std::unique_ptr<ROOT::Math::IBaseFunctionMultiDim>
MakeNumericalChi2(std::shared_ptr<EvaluationContext> context);
std::unique_ptr<ROOT::Math::IMultiGradFunction>
MakeAnalyticChi2(std::shared_ptr<EvaluationContext> context);
double EvaluateChi2(const ROOT::Math::IBaseFunctionMultiDim& function,
                    const std::vector<double>& point);
MinimizerResult Minimize(const std::shared_ptr<EvaluationContext>& context,
                         const std::vector<double>& start,
                         bool runHesse);

void MakeParameterBranches(TTree* tree, const std::vector<Parameter>& parameters,
                           std::vector<double>& storage);
void MakeOutputMomentBranches(TTree* tree,
                              const std::vector<OutputMoment>& moments,
                              std::vector<double>& storage);
void FillHessianProducts(
    const EvaluationContext& context, const ROOT::Math::Minimizer& minimizer,
    const std::vector<double>& fullValues, bool covarianceIsUsable,
    std::vector<double>& parameterErrors, std::vector<double>& momentErrors,
    double& ratio, double& ratioError);

} // namespace emi::detail
