// Simultaneous mass-dependent moment fit.  The angular algebra is shared with
// the ordinary inversion; this layer replaces its free bin-by-bin complex
// amplitudes with coherent Breit-Wigner, Flatte, and background functions.
#include "emi/Runner.h"

#include "Detail.h"
#include "MassModels/MassModel.h"
#include "MassModels/OrderedMasses.h"
#include "Progress.h"
#include "WorkerOutput.h"

#include "Math/Factory.h"
#include "ROOT/TProcessExecutor.hxx"
#include "TFile.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace emi {
namespace {

using Complex = std::complex<double>;

// One mass point: values/errors follow angular->observed branch order.
struct DataPoint {
  double mass = 0.0;
  std::vector<double> values, errors;
};
// One shared pole, coupling, or background coordinate supplied to Minuit.
using Coordinate = MassModels::FitCoordinate;
// Indices into the wave and coordinate tables, built once before fitting.
struct ResonanceTerm {
  int wave = -1, mass = -1, width = -1, magnitude = -1, phase = -1;
  MassLineShape lineShape = MassLineShape::ConstantWidthBreitWigner;
  std::vector<MassModels::FlatteChannel> channels;
  std::vector<int> channelCoordinates;
};
struct BackgroundTerm {
  int wave = -1, re = -1, im = -1;
};

// A mass fit reuses the angular bilinear algebra, but fits common line-shape
// coordinates across all mass points instead of independent bin amplitudes.
struct Problem {
  std::shared_ptr<detail::EvaluationContext> angular;
  std::vector<DataPoint> data;
  std::vector<int> amplitudeForMagnitude;
  std::vector<Coordinate> coordinates; // Independent MINUIT coordinates.
  MassModels::OrderedMasses massMapping;
  std::vector<ResonanceTerm> resonances;
  std::vector<BackgroundTerm> backgrounds;
  std::unordered_map<long long, int> waveIndex;
};

long long WaveKey(const MassWave& wave) {
  return detail::MakeParameterKey(wave.reflectivity, wave.orientation, wave.k, wave.l,
                                  wave.m, false);
}

void ValidateWave(const MassWave& wave) {
  if ((wave.reflectivity != 'a' && wave.reflectivity != 'b') ||
      (wave.orientation != 'T' && wave.orientation != 'L') || wave.l < 0 ||
      std::abs(wave.m) > wave.l) {
    throw std::invalid_argument("Invalid mass-model wave label");
  }
}

double LeafValue(TTree& tree, const std::string& name) {
  auto* leaf = tree.GetLeaf(name.c_str());
  if (!leaf)
    throw std::runtime_error("Missing scalar branch '" + name + "'");
  if (leaf->GetNdata() != 1) {
    throw std::runtime_error("Mass-dependent fits require scalar branch '" + name +
                             "'");
  }
  return leaf->GetValue(0);
}

// Resolve selection before discovering moments so excluded rows cannot affect
// the fitted branch set. Preserve tree order even for an unordered index list.
std::vector<Long64_t> SelectMassRows(const FitMassModelConfig& config) {
  if (config.yieldBranch &&
      (config.yieldBranch->empty() || !std::isfinite(config.normalisation) ||
       config.normalisation <= 0.0))
    throw std::invalid_argument("Yield scaling requires a branch and positive finite normalisation");
  if ((config.minMassGeV && !std::isfinite(*config.minMassGeV)) ||
      (config.maxMassGeV && !std::isfinite(*config.maxMassGeV)) ||
      (config.minMassGeV && config.maxMassGeV &&
       *config.minMassGeV > *config.maxMassGeV))
    throw std::invalid_argument("Mass selection requires finite, ordered bounds");
  std::unique_ptr<TFile> file(TFile::Open(config.fit.input.c_str(), "READ"));
  if (!file || file->IsZombie())
    throw std::runtime_error("Failed to open " + config.fit.input.string());
  auto* tree = dynamic_cast<TTree*>(file->Get(config.fit.tree.c_str()));
  if (!tree || tree->GetEntries() == 0)
    throw std::runtime_error("Mass input tree is missing or empty");
  std::unordered_set<Long64_t> indices;
  for (auto bin : config.bins) {
    if (bin < 0 || bin >= tree->GetEntries())
      throw std::invalid_argument("Mass bin index is outside the input tree: " +
                                  std::to_string(bin));
    indices.insert(bin);
  }
  std::vector<Long64_t> rows;
  for (Long64_t row = 0; row < tree->GetEntries(); ++row) {
    if (!indices.empty() && indices.count(row) == 0)
      continue;
    tree->GetEntry(row);
    const double mass = LeafValue(*tree, config.massBranch);
    if (!std::isfinite(mass) || mass <= 0.0)
      throw std::runtime_error("Mass branch contains an invalid value at entry " +
                               std::to_string(row));
    if ((config.minMassGeV && mass < *config.minMassGeV) ||
        (config.maxMassGeV && mass > *config.maxMassGeV))
      continue;
    rows.push_back(row);
  }
  if (rows.empty())
    throw std::invalid_argument("Mass selection contains no input entries");
  return rows;
}

std::vector<DataPoint> ReadMassData(const FitMassModelConfig& config,
                                    const detail::EvaluationContext& angular,
                                    const std::vector<Long64_t>& rows) {
  std::unique_ptr<TFile> file(TFile::Open(config.fit.input.c_str(), "READ"));
  if (!file || file->IsZombie())
    throw std::runtime_error("Failed to open " + config.fit.input.string());
  auto* tree = dynamic_cast<TTree*>(file->Get(config.fit.tree.c_str()));
  if (!tree || tree->GetEntries() == 0)
    throw std::runtime_error("Mass input tree is missing or empty");
  std::vector<DataPoint> result;
  result.reserve(rows.size());
  for (Long64_t row : rows) {
    // Read each ROOT entry once. All leaves below refer to this loaded row.
    tree->GetEntry(row);
    DataPoint point;
    point.mass = LeafValue(*tree, config.massBranch);
    if (!std::isfinite(point.mass) || point.mass <= 0.0)
      throw std::runtime_error("Mass branch contains an invalid value");
    double scale = 1.0;
    if (config.yieldBranch) {
      const double yield = LeafValue(*tree, *config.yieldBranch);
      scale = yield / config.normalisation;
      if (!std::isfinite(yield) || yield <= 0.0 ||
          !std::isfinite(scale) || scale <= 0.0)
        throw std::runtime_error("Yield must be positive and finite at entry " +
                                 std::to_string(row));
    }
    for (const auto& moment : angular.observed) {
      const double value = scale * LeafValue(*tree, moment.name);
      const double error = scale * LeafValue(*tree, moment.name + "_err");
      if (!std::isfinite(value) || !std::isfinite(error) || error <= 0.0) {
        throw std::runtime_error("Moment '" + moment.name +
                                 "' has an invalid value or uncertainty");
      }
      point.values.push_back(value);
      point.errors.push_back(error);
    }
    result.push_back(std::move(point));
  }
  return result;
}

int AddCoordinate(Problem& problem, const std::string& name,
                  const MassParameter& value) {
  if (!std::isfinite(value.initial) || (value.bounded &&
      (!std::isfinite(value.lower) || !std::isfinite(value.upper) || value.lower > value.upper ||
       value.initial < value.lower || value.initial > value.upper))) {
    throw std::invalid_argument("Invalid limits for mass-model coordinate " + name);
  }
  for (const auto& existing : problem.coordinates)
    if (existing.name == name)
      throw std::invalid_argument("Duplicate mass coordinate " + name);
  problem.coordinates.push_back({name, value});
  return static_cast<int>(problem.coordinates.size() - 1);
}

// Generation and fitting must assign exactly the same coordinate names and
// wave indices. Keep that mapping here so neither path can drift independently.
Problem BuildAmplitudeProblem(const MassModelConfig& massModel,
                              std::shared_ptr<detail::EvaluationContext> angular) {
  if (massModel.resonances.empty() && massModel.backgrounds.empty())
    throw std::invalid_argument("Mass model has no amplitudes");
  Problem problem;
  problem.angular = std::move(angular);
  problem.amplitudeForMagnitude.assign(problem.angular->fullPars.size(), -1);
  for (size_t i = 0; i < problem.angular->fullPars.size(); ++i) {
    const auto label = detail::ParseParameterLabel(problem.angular->fullPars[i].name);
    if (!label.valid || problem.angular->fullPars[i].phase)
      continue;
    const long long key = detail::MakeParameterKey(
        label.reflectivity, label.orientation, label.k, label.l, label.m, false);
    const int index = static_cast<int>(problem.waveIndex.size());
    auto inserted = problem.waveIndex.emplace(key, index);
    problem.amplitudeForMagnitude[i] = inserted.first->second;
  }
  auto addWave = [&](const MassWave& wave) {
    ValidateWave(wave);
    const auto found = problem.waveIndex.find(WaveKey(wave));
    if (found == problem.waveIndex.end())
      throw std::invalid_argument("Mass-model wave is absent from ModelConfig");
    return found->second;
  };
  std::unordered_map<std::string, int> sharedPoleCoordinates;
  auto addShared = [&](const std::string& name, const MassParameter& control) {
    const auto found = sharedPoleCoordinates.find(name);
    if (found != sharedPoleCoordinates.end()) {
      const auto& previous = problem.coordinates[found->second].control;
      if (previous.initial != control.initial || previous.bounded != control.bounded ||
          (control.bounded && (previous.lower != control.lower || previous.upper != control.upper)) ||
          previous.fixed != control.fixed)
        throw std::invalid_argument("Inconsistent shared pole control " + name);
      return found->second;
    }
    const int index = AddCoordinate(problem, name, control);
    sharedPoleCoordinates.emplace(name, index);
    return index;
  };
  std::unordered_map<std::string, const MassResonanceConfig*> sharedShapes;
  for (const auto& resonance : massModel.resonances) {
    if (resonance.name.empty())
      throw std::invalid_argument("Every mass resonance needs a name");
    if (resonance.poleMass.initial <= 0.0 ||
        (resonance.poleMass.bounded && resonance.poleMass.lower <= 0.0) ||
        (resonance.lineShape == MassLineShape::ConstantWidthBreitWigner &&
         (resonance.width.initial <= 0.0 || (resonance.width.bounded && resonance.width.lower <= 0.0))) ||
        (resonance.couplingMagnitude.bounded && resonance.couplingMagnitude.lower < 0.0))
      throw std::invalid_argument(
          "Mass/width bounds must be positive and magnitudes nonnegative");
    const std::string prefix = resonance.name;
    const std::string couplingPrefix = prefix + "_" + resonance.wave.reflectivity +
                                       (resonance.wave.orientation == 'L' ? "_L" : "") +
                                       "_" + std::to_string(resonance.wave.l) + "_" +
                                       detail::IndexString(resonance.wave.m) +
                                       (resonance.wave.k == 0 ? "" :
                                        "_k" + detail::IndexString(resonance.wave.k));
    const auto shared = sharedShapes.emplace(prefix, &resonance);
    if (!shared.second) {
      const auto& previous = *shared.first->second;
      if (previous.lineShape != resonance.lineShape ||
          previous.channels.size() != resonance.channels.size())
        throw std::invalid_argument("Inconsistent shared line shape " + prefix);
      for (size_t i = 0; i < resonance.channels.size(); ++i) {
        const auto& first = previous.channels[i];
        const auto& second = resonance.channels[i];
        if (first.name != second.name || first.firstMassGeV != second.firstMassGeV ||
            first.secondMassGeV != second.secondMassGeV)
          throw std::invalid_argument("Inconsistent shared channel definition " + prefix);
      }
    }
    ResonanceTerm term;
    term.wave = addWave(resonance.wave);
    term.mass = addShared(prefix + "_mass", resonance.poleMass);
    term.lineShape = resonance.lineShape;
    switch (resonance.lineShape) {
    case MassLineShape::ConstantWidthBreitWigner:
      term.width = addShared(prefix + "_width", resonance.width);
      break;
    case MassLineShape::Flatte: {
      if (resonance.channels.size() < 2)
        throw std::invalid_argument("Flatte requires at least two decay channels");
      std::unordered_set<std::string> channelNames;
      for (const auto& channel : resonance.channels) {
        if (channel.name.empty() || !channelNames.insert(channel.name).second ||
            !std::isfinite(channel.firstMassGeV) || channel.firstMassGeV < 0.0 ||
            !std::isfinite(channel.secondMassGeV) || channel.secondMassGeV < 0.0 ||
            channel.couplingSquared.lower < 0.0)
          throw std::invalid_argument("Invalid Flatte channel definition for " + prefix);
        term.channels.push_back({channel.firstMassGeV, channel.secondMassGeV,
                                 channel.couplingSquared.initial});
        term.channelCoordinates.push_back(
            addShared(prefix + "_" + channel.name + "_g2", channel.couplingSquared));
      }
      break;
    }
    default:
      throw std::invalid_argument("Unknown mass line shape");
    }
    term.magnitude = AddCoordinate(problem, couplingPrefix + "_mag", resonance.couplingMagnitude);
    term.phase = AddCoordinate(problem, couplingPrefix + "_phase", resonance.couplingPhase);
    problem.resonances.push_back(std::move(term));
  }
  for (const auto& background : massModel.backgrounds) {
    if (background.name.empty())
      throw std::invalid_argument("Every mass background needs a name");
    const std::string prefix = "bg_" + background.name;
    problem.backgrounds.push_back(
        {addWave(background.wave),
         AddCoordinate(problem, prefix + "_re", background.real),
         AddCoordinate(problem, prefix + "_im", background.imag)});
  }
  problem.massMapping = MassModels::OrderedMasses(massModel, problem.coordinates);
  problem.coordinates = problem.massMapping.fitCoordinates;
  return problem;
}

Problem BuildProblem(const FitMassModelConfig& config, const ModelConfig& model) {
  if (config.fit.input.empty() || config.fit.output.empty())
    throw std::invalid_argument("Mass fit requires input and output paths");
  if (config.fit.starts == 0)
    throw std::invalid_argument("Mass fit starts must be positive");
  FitConfig angularFit = config.fit;
  angularFit.bin = 0; // Mass-fit inputs have scalar leaves.
  angularFit.verbose = false;
  const auto rows = SelectMassRows(config);
  auto internal = detail::MakeInternalConfig(angularFit, model);
  internal.momentEntry = rows.front();
  auto angular = detail::BuildContext(internal);
  auto problem = BuildAmplitudeProblem(config.model, std::move(angular));
  problem.data = ReadMassData(config, *problem.angular, rows);
  return problem;
}

void AmplitudesAndDerivatives(const Problem& problem, double mass, const double* x,
                              std::vector<Complex>& amplitude,
                              std::vector<std::vector<Complex>>& derivative,
                              bool computeDerivatives = true,
                              const MassModels::PhysicalPoint* mapped = nullptr) {
  const auto storage = mapped ? MassModels::PhysicalPoint{} : problem.massMapping.Evaluate(x);
  const auto& point = mapped ? *mapped : storage;
  const auto* physical = point.values.data();
  amplitude.assign(problem.waveIndex.size(), {});
  if (computeDerivatives)
    derivative.assign(point.values.size(),
                      std::vector<Complex>(problem.waveIndex.size(), Complex{}));
  for (const auto& term : problem.resonances) {
    const double pole = physical[term.mass];
    Complex shape, dPole, dWidth;
    std::vector<Complex> dChannels;
    if (term.lineShape == MassLineShape::Flatte) {
      auto channels = term.channels;
      for (size_t i = 0; i < channels.size(); ++i)
        channels[i].couplingSquaredGeV2 = physical[term.channelCoordinates[i]];
      const auto evaluated = MassModels::Flatte(pole, mass, channels);
      shape = evaluated.value;
      dPole = evaluated.poleMassDerivative;
      dChannels = evaluated.couplingSquaredDerivatives;
    } else {
      const double width = physical[term.width];
      const double numerator = pole * width;
      const Complex denominator{pole * pole - mass * mass, -pole * width};
      shape = numerator / denominator;
      if (computeDerivatives) {
        dPole = (width * denominator - numerator * Complex{2.0 * pole, -width}) /
                (denominator * denominator);
        dWidth = (pole * denominator - numerator * Complex{0.0, -pole}) /
                 (denominator * denominator);
      }
    }
    const Complex coupling = std::polar(physical[term.magnitude], physical[term.phase]);
    amplitude[term.wave] += coupling * shape;
    if (!computeDerivatives) continue;
    derivative[term.mass][term.wave] += coupling * dPole;
    if (term.width >= 0)
      derivative[term.width][term.wave] += coupling * dWidth;
    for (size_t i = 0; i < term.channelCoordinates.size(); ++i)
      derivative[term.channelCoordinates[i]][term.wave] += coupling * dChannels[i];
    derivative[term.magnitude][term.wave] += std::polar(1.0, physical[term.phase]) * shape;
    derivative[term.phase][term.wave] += Complex{0.0, 1.0} * coupling * shape;
  }
  for (const auto& term : problem.backgrounds) {
    amplitude[term.wave] += Complex{physical[term.re], physical[term.im]};
    if (!computeDerivatives) continue;
    derivative[term.re][term.wave] += 1.0;
    derivative[term.im][term.wave] += Complex{0.0, 1.0};
  }
  if (!computeDerivatives) return;
  std::vector<std::vector<Complex>> fitDerivative(
      problem.coordinates.size(), std::vector<Complex>(amplitude.size(), Complex{}));
  for (size_t p = 0; p < point.jacobian.size(); ++p)
    for (const auto& entry : point.jacobian[p])
      for (size_t wave = 0; wave < amplitude.size(); ++wave)
        fitDerivative[entry.first][wave] += entry.second * derivative[p][wave];
  derivative.swap(fitDerivative);
}

double MomentValue(const Problem& problem, const detail::MomentModel& moment,
                   const std::vector<Complex>& amplitudes) {
  double value = 0.0;
  for (const auto& term : moment.terms) {
    const int first = problem.amplitudeForMagnitude[term.idxMag1];
    const int second = problem.amplitudeForMagnitude[term.idxMag2];
    if (first < 0 || second < 0)
      continue;
    const Complex product = amplitudes[first] * std::conj(amplitudes[second]);
    value += term.coeff *
             (term.trig == detail::TrigKind::kCos ? product.real() : product.imag());
  }
  return value;
}

// Accumulate derivatives with respect to complex wave amplitudes first,
// then apply the chain rule to the BW coordinates. For an imaginary bilinear
// the two wave contributions have opposite signs (the second is conjugated).
void MomentAndGradient(const Problem& problem, const detail::MomentModel& moment,
                       const std::vector<Complex>& amplitudes,
                       const std::vector<std::vector<Complex>>& dAmplitude,
                       double& value, std::vector<double>& gradient) {
  value = 0.0;
  std::vector<Complex> waveDerivative(amplitudes.size(), Complex{});
  for (const auto& term : moment.terms) {
    const int first = problem.amplitudeForMagnitude[term.idxMag1];
    const int second = problem.amplitudeForMagnitude[term.idxMag2];
    if (first < 0 || second < 0)
      continue;
    const Complex firstAmplitude = amplitudes[first],
                  secondAmplitude = amplitudes[second];
    const Complex product = firstAmplitude * std::conj(secondAmplitude);
    if (term.trig == detail::TrigKind::kCos) {
      value += term.coeff * product.real();
      waveDerivative[first] += term.coeff * std::conj(secondAmplitude);
      waveDerivative[second] += term.coeff * std::conj(firstAmplitude);
    } else {
      value += term.coeff * product.imag();
      waveDerivative[first] += Complex{0.0, -term.coeff} * std::conj(secondAmplitude);
      waveDerivative[second] += Complex{0.0, term.coeff} * std::conj(firstAmplitude);
    }
  }
  gradient.assign(problem.coordinates.size(), 0.0);
  for (size_t parameter = 0; parameter < gradient.size(); ++parameter) {
    double derivative = 0.0;
    for (size_t wave = 0; wave < waveDerivative.size(); ++wave)
      derivative += (dAmplitude[parameter][wave] * waveDerivative[wave]).real();
    gradient[parameter] = derivative;
  }
}

// Measured H04 is H0 + epsilon H4; ordinary observations use one raw moment.
// Use the same contraction for fitting, generated truth, and saved predictions.
void ObservedMomentAndGradient(const Problem& problem, size_t index,
                               const std::vector<Complex>& amplitudes,
                               const std::vector<std::vector<Complex>>& derivatives,
                               double& value, std::vector<double>& gradient,
                               std::vector<double>& secondaryGradient) {
  const auto& angular = *problem.angular;
  const bool mixed = angular.observed[index].isMixed04;
  const int primary =
      mixed ? angular.observedModelIdx0[index] : angular.observedModelIdx[index];
  MomentAndGradient(problem, angular.modelsRec[primary], amplitudes, derivatives, value,
                    gradient);
  if (mixed) {
    double h4 = 0.0;
    MomentAndGradient(problem, angular.modelsRec[angular.observedModelIdx4[index]],
                      amplitudes, derivatives, h4, secondaryGradient);
    value += angular.cfg.epsilon * h4;
    for (size_t p = 0; p < gradient.size(); ++p)
      gradient[p] += angular.cfg.epsilon * secondaryGradient[p];
  }
}

// ROOT requires value, full-gradient, single-derivative, and clone entry points.
// Keep reusable scratch buffers local to each objective/worker.
class MassChi2 final : public ROOT::Math::IMultiGradFunction {
public:
  explicit MassChi2(std::shared_ptr<Problem> problem) : problem_(std::move(problem)) {}
  unsigned int NDim() const override {
    return static_cast<unsigned>(problem_->coordinates.size());
  }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override {
    return new MassChi2(problem_);
  }
  double DoEval(const double* x) const {
    const auto mapped = problem_->massMapping.Evaluate(x);
    if (!problem_->massMapping.Valid(mapped)) return 1.e300;
    double chi2 = 0.0;
    for (const auto& point : problem_->data) {
      AmplitudesAndDerivatives(*problem_, point.mass, x, amplitude_, derivative_, false, &mapped);
      for (size_t i = 0; i < problem_->angular->observed.size(); ++i) {
        const auto& observed = problem_->angular->observed[i];
        double model = MomentValue(
            *problem_,
            problem_->angular->modelsRec[observed.isMixed04
                                             ? problem_->angular->observedModelIdx0[i]
                                             : problem_->angular->observedModelIdx[i]],
            amplitude_);
        if (observed.isMixed04)
          model +=
              problem_->angular->cfg.epsilon *
              MomentValue(
                  *problem_,
                  problem_->angular->modelsRec[problem_->angular->observedModelIdx4[i]],
                  amplitude_);
        const double pull = (point.values[i] - model) / point.errors[i];
        chi2 += pull * pull;
      }
    }
    return std::isfinite(chi2) ? chi2 : 1.e300;
  }
  void Gradient(const double* x, double* output) const override {
    std::fill(output, output + NDim(), 0.0);
    const auto mapped = problem_->massMapping.Evaluate(x);
    if (!problem_->massMapping.Valid(mapped)) return;
    const auto& angular = *problem_->angular;
    for (const auto& point : problem_->data) {
      AmplitudesAndDerivatives(*problem_, point.mass, x, amplitude_, derivative_, true, &mapped);
      std::vector<Complex> waveGradient(amplitude_.size(), Complex{});
      // Sum the chi-square derivative in wave space before applying the
      // line-shape chain rule. This avoids repeating the full coordinate-by-
      // wave contraction independently for every observed moment.
      auto accumulate = [&](const detail::MomentModel& moment, double weight) {
        for (const auto& term : moment.terms) {
          const int first = problem_->amplitudeForMagnitude[term.idxMag1];
          const int second = problem_->amplitudeForMagnitude[term.idxMag2];
          if (first < 0 || second < 0) continue;
          const double coefficient = weight * term.coeff;
          if (term.trig == detail::TrigKind::kCos) {
            waveGradient[first] += coefficient * std::conj(amplitude_[second]);
            waveGradient[second] += coefficient * std::conj(amplitude_[first]);
          } else {
            waveGradient[first] += Complex{0.0, -coefficient} * std::conj(amplitude_[second]);
            waveGradient[second] += Complex{0.0, coefficient} * std::conj(amplitude_[first]);
          }
        }
      };
      for (size_t i = 0; i < angular.observed.size(); ++i) {
        const bool mixed = angular.observed[i].isMixed04;
        const auto& primary = angular.modelsRec[mixed ? angular.observedModelIdx0[i]
                                                     : angular.observedModelIdx[i]];
        double model = MomentValue(*problem_, primary, amplitude_);
        if (mixed)
          model += angular.cfg.epsilon * MomentValue(*problem_, angular.modelsRec[angular.observedModelIdx4[i]], amplitude_);
        const double factor = -2.0 * (point.values[i] - model) / (point.errors[i] * point.errors[i]);
        accumulate(primary, factor);
        if (mixed) accumulate(angular.modelsRec[angular.observedModelIdx4[i]], angular.cfg.epsilon * factor);
      }
      for (size_t p = 0; p < derivative_.size(); ++p)
        for (size_t wave = 0; wave < waveGradient.size(); ++wave)
          output[p] += (derivative_[p][wave] * waveGradient[wave]).real();
    }
  }
  double DoDerivative(const double* x, unsigned int coordinate) const override {
    if (coordinate >= NDim()) return 0.0;
    // Some ROOT adapters request components separately. Reuse the complete
    // gradient for an unchanged point instead of rebuilding all wave and
    // line-shape derivatives for each coordinate.
    if (derivativePoint_.size() != NDim() ||
        !std::equal(derivativePoint_.begin(), derivativePoint_.end(), x)) {
      derivativePoint_.assign(x, x + NDim());
      derivativeCache_.resize(NDim());
      Gradient(x, derivativeCache_.data());
    }
    return derivativeCache_[coordinate];
  }

private:
  std::shared_ptr<Problem> problem_;
  mutable std::vector<Complex> amplitude_;
  mutable std::vector<std::vector<Complex>> derivative_;
  mutable std::vector<double> derivativePoint_, derivativeCache_;
};

class MassChi2NoGrad final : public ROOT::Math::IBaseFunctionMultiDim {
public:
  explicit MassChi2NoGrad(std::shared_ptr<Problem> problem)
      : objective_(std::move(problem)) {}
  unsigned int NDim() const override {
    return objective_.NDim();
  }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override {
    return new MassChi2NoGrad(*this);
  }
  double DoEval(const double* x) const override {
    return objective_.DoEval(x);
  }

private:
  MassChi2 objective_;
};

void WritePredictions(const Problem& problem, const std::vector<double>& x,
                      const std::string& massBranch, const char* treeName,
                      const std::vector<double>& masses) {
  TTree predictions(treeName, "Best-fit model evaluated from resonance parameters");
  double mass = 0.0;
  predictions.Branch(massBranch.c_str(), &mass);
  std::vector<double> values(problem.angular->observed.size());
  for (size_t i = 0; i < values.size(); ++i)
    predictions.Branch(("model__" + problem.angular->observed[i].name).c_str(), &values[i]);
  std::vector<Complex> amplitudes;
  std::vector<std::vector<Complex>> derivatives;
  for (double pointMass : masses) {
    mass = pointMass;
    AmplitudesAndDerivatives(problem, mass, x.data(), amplitudes, derivatives, false);
    for (size_t i = 0; i < values.size(); ++i) {
      const auto& angular = *problem.angular;
      const auto& observed = angular.observed[i];
      values[i] = MomentValue(problem, angular.modelsRec[observed.isMixed04
          ? angular.observedModelIdx0[i] : angular.observedModelIdx[i]], amplitudes);
      if (observed.isMixed04)
        values[i] += angular.cfg.epsilon * MomentValue(
            problem, angular.modelsRec[angular.observedModelIdx4[i]], amplitudes);
    }
    predictions.Fill();
  }
  predictions.Write();
}

// Single-k reflectivities do not interfere, so each populated sector needs
// its own reference. Explicit-k polarized models retain one overall phase.
std::vector<int> FindPhaseGauges(const MassModelConfig& massModel,
                                  const Problem& problem) {
  std::vector<MassPhaseGauge> overrides = massModel.phaseGauges;
  if (massModel.phaseGauge &&
      std::none_of(overrides.begin(), overrides.end(), [&](const auto& gauge) {
        return gauge.wave.reflectivity == massModel.phaseGauge->wave.reflectivity;
      }))
    overrides.push_back(*massModel.phaseGauge);
  std::unordered_set<char> overriddenSectors;
  for (const auto& gauge : overrides) {
    ValidateWave(gauge.wave);
    if (!overriddenSectors.insert(gauge.wave.reflectivity).second)
      throw std::invalid_argument("Duplicate mass phase gauge for a reflectivity");
    const auto found = std::find_if(
        massModel.resonances.begin(), massModel.resonances.end(),
        [&](const auto& resonance) {
          return resonance.name == gauge.resonance &&
                 WaveKey(resonance.wave) == WaveKey(gauge.wave);
        });
    if (found == massModel.resonances.end())
      throw std::invalid_argument("Requested mass phase gauge resonance is absent");
  }
  const bool singleK =
      problem.angular->cfg.nucleonPolarization == NucleonPolarization::None;
  if (!singleK && overrides.size() > 1)
    throw std::invalid_argument("Explicit-k mass fits accept one overall phase gauge");

  auto select = [&](char reflectivity, const MassPhaseGauge* gauge) {
    int selected = -1, bestL = -1, bestM = -999999;
    for (size_t i = 0; i < massModel.resonances.size(); ++i) {
      const auto& resonance = massModel.resonances[i];
      const bool matches = gauge
          ? resonance.name == gauge->resonance &&
                WaveKey(resonance.wave) == WaveKey(gauge->wave)
          : resonance.wave.reflectivity == reflectivity &&
                resonance.wave.orientation == 'T';
      if (matches && (gauge || resonance.wave.l > bestL ||
                      (resonance.wave.l == bestL && resonance.wave.m > bestM))) {
        selected = problem.massMapping.DirectFitIndex(problem.resonances[i].phase);
        bestL = resonance.wave.l;
        bestM = resonance.wave.m;
      }
    }
    if (selected < 0)
      throw std::invalid_argument(
          std::string("Mass phase gauge needs a transverse resonance in sector ") +
          reflectivity + " or an explicit override");
    return selected;
  };
  if (!singleK)
    return {select('a', overrides.empty() ? nullptr : &overrides.front())};

  std::vector<int> result;
  for (char reflectivity : {'a', 'b'}) {
    const bool populated = std::any_of(
        massModel.resonances.begin(), massModel.resonances.end(),
        [&](const auto& r) { return r.wave.reflectivity == reflectivity; }) ||
        std::any_of(massModel.backgrounds.begin(), massModel.backgrounds.end(),
        [&](const auto& b) { return b.wave.reflectivity == reflectivity; });
    if (!populated)
      continue;
    const auto override = std::find_if(overrides.begin(), overrides.end(),
        [&](const auto& gauge) { return gauge.wave.reflectivity == reflectivity; });
    result.push_back(select(reflectivity,
                           override == overrides.end() ? nullptr : &*override));
  }
  return result;
}

void WriteMassMetadata(const Problem& problem) {
  std::string physicalNames, fitNames, ordering;
  for (const auto& coordinate : problem.massMapping.physicalCoordinates)
    physicalNames += coordinate.name + "\n";
  for (const auto& coordinate : problem.coordinates)
    fitNames += coordinate.name + "\n";
  for (const auto& node : problem.massMapping.nodes)
    ordering += node.group + " " + node.resonance + " " +
                std::to_string(node.physical) + " " + std::to_string(node.previous) +
                " " + std::to_string(node.fit) + " " +
                std::to_string(static_cast<int>(node.kind)) + "\n";
  TObjString(physicalNames.c_str()).Write("mass_fit_physical_coordinate_names");
  TObjString(fitNames.c_str()).Write("mass_fit_coordinate_names");
  TObjString(ordering.c_str()).Write("mass_fit_ordering_map");
  TParameter<double>("mass_fit_minimum_gap_GeV", problem.massMapping.minimumGap, 'f').Write();
}

void RunMassWorker(const FitMassModelConfig& config, const ModelConfig& model,
                   const std::filesystem::path& output, unsigned worker,
                   unsigned firstStart, detail::ProgressCounter progress = {}) {
  auto problem = std::make_shared<Problem>(BuildProblem(config, model));
  MassChi2 objective(problem);
  const auto phaseGauges = FindPhaseGauges(config.model, *problem);
  auto isPhaseGauge = [&](size_t i) {
    return std::find(phaseGauges.begin(), phaseGauges.end(), static_cast<int>(i)) !=
           phaseGauges.end();
  };
  unsigned freeCoordinates = 0;
  for (size_t i = 0; i < problem->coordinates.size(); ++i)
    if (!problem->coordinates[i].control.fixed && !isPhaseGauge(i))
      ++freeCoordinates;
  if (freeCoordinates == 0)
    throw std::invalid_argument("Mass model has no free fit coordinates");

  TRandom3 random(config.fit.seed);
  if (config.fit.seed == 0)
    random.SetSeed(0);
  std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "RECREATE"));
  if (!file || file->IsZombie())
    throw std::runtime_error("Could not create " + output.string());
  TTree summary("massFitResults", "Simultaneous mass-dependent moment fit");
  bool ok = false, hesse = false;
  int status = -999;
  int startIndex = 0, workerIndex = static_cast<int>(worker);
  unsigned seed = config.fit.seed;
  int ndf = static_cast<int>(problem->data.size() * problem->angular->observed.size()) -
            static_cast<int>(freeCoordinates);
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  summary.Branch("fit_ok", &ok);
  summary.Branch("hesse_ok", &hesse);
  summary.Branch("status", &status);
  summary.Branch("worker", &workerIndex);
  summary.Branch("seed", &seed);
  summary.Branch("start", &startIndex);
  double initialChi2 = 0.0;
  summary.Branch("initial_chi2", &initialChi2);
  summary.Branch("chi2", &chi2);
  summary.Branch("ndf", &ndf);
  const size_t physicalSize = problem->massMapping.physicalCoordinates.size();
  const size_t fitSize = problem->coordinates.size();
  std::vector<double> errors(physicalSize), result(physicalSize);
  std::vector<double> fitResult(fitSize), fitErrors(fitSize);
  std::vector<double> fitCovariance, physicalCovariance, massJacobian;
  bool covarianceOk = false;
  summary.Branch("covariance_ok", &covarianceOk);
  summary.Branch("fit_covariance", &fitCovariance);
  summary.Branch("physical_covariance", &physicalCovariance);
  summary.Branch("mass_jacobian", &massJacobian);
  for (size_t i = 0; i < physicalSize; ++i) {
    const auto& name = problem->massMapping.physicalCoordinates[i].name;
    summary.Branch(name.c_str(), &result[i]);
    summary.Branch(("err__" + name).c_str(), &errors[i]);
  }
  // Non-mass fit coordinates already have their physical branches above.
  for (size_t i = 0; i < fitSize; ++i) {
    const auto& name = problem->coordinates[i].name;
    if (name.rfind("ordered__", 0) != 0 && name.rfind("positive__", 0) != 0)
      continue;
    summary.Branch(name.c_str(), &fitResult[i]);
    summary.Branch(("err__" + name).c_str(), &fitErrors[i]);
  }
  std::vector<double> orderedMasses(problem->massMapping.nodes.size());
  std::vector<double> orderedErrors(orderedMasses.size());
  for (size_t i = 0; i < orderedMasses.size(); ++i) {
    const auto& node = problem->massMapping.nodes[i];
    const std::string name = "ordered__" + node.group + "__" +
        (node.previous < 0 ? "base_GeV" : "gap_" + node.resonance + "_GeV");
    summary.Branch(name.c_str(), &orderedMasses[i]);
    summary.Branch(("err__" + name).c_str(), &orderedErrors[i]);
  }
  for (unsigned localStart = 0; localStart < config.fit.starts; ++localStart) {
    // Each trial owns fresh minimizer and covariance state.
    std::unique_ptr<ROOT::Math::IBaseFunctionMultiDim> numericalObjective;
    std::unique_ptr<ROOT::Math::Minimizer> minimizer(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    if (!minimizer)
      throw std::runtime_error("Could not create Minuit2/Migrad");
    if (config.fit.useNumericalGradients) {
      numericalObjective = std::make_unique<MassChi2NoGrad>(problem);
      minimizer->SetFunction(*numericalObjective);
    } else {
      minimizer->SetFunction(objective);
    }
    // Delta(chi2)=1 defines one-sigma errors. A tight tolerance matters for
    // noiseless closure: the default can stop while weak couplings are biased.
    minimizer->SetErrorDef(1.0);
    minimizer->SetStrategy(1);
    minimizer->SetTolerance(1.e-8);
    minimizer->SetMaxFunctionCalls(config.maxFunctionCalls);
    minimizer->SetPrintLevel(-1); // The parent progress monitor owns terminal output.
    for (size_t i = 0; i < problem->coordinates.size(); ++i) {
      const auto& coordinate = problem->coordinates[i];
      if (coordinate.control.fixed || isPhaseGauge(i)) {
        minimizer->SetFixedVariable(i, coordinate.name, coordinate.control.initial);
      } else if (!coordinate.control.bounded) {
        minimizer->SetVariable(i, coordinate.name, coordinate.control.initial, 1.e-3);
      } else {
        const double step = std::min(1.e-3,
            0.01 * (coordinate.control.upper - coordinate.control.lower));
        minimizer->SetLimitedVariable(i, coordinate.name, coordinate.control.initial,
                                      step, coordinate.control.lower, coordinate.control.upper);
      }
    }
    std::vector<double> startValues(problem->coordinates.size());
    for (size_t i = 0; i < problem->coordinates.size(); ++i) {
      const auto& coordinate = problem->coordinates[i];
      const bool coupling =
          ((coordinate.name.size() >= 4 &&
            coordinate.name.compare(coordinate.name.size() - 4, 4, "_mag") == 0) ||
           (coordinate.name.size() >= 6 &&
            coordinate.name.compare(coordinate.name.size() - 6, 6, "_phase") == 0));
      double value = coordinate.control.initial;
      // Start zero uses the steering card exactly. Only subsequent global
      // starts explore random couplings; pole starts remain configured.
      if (firstStart + localStart > 0 && coupling && !coordinate.control.fixed &&
          !isPhaseGauge(i))
        value = coordinate.control.bounded
            ? random.Uniform(coordinate.control.lower, coordinate.control.upper)
            : random.Gaus(coordinate.control.initial, 1.0);
      // Fix the reference at its supplied phase. Setting it to zero alone
      // would alter relative phases and spoil an otherwise exact truth start.
      startValues[i] = value;
      minimizer->SetVariableValue(i, value);
    }
    initialChi2 = objective.DoEval(startValues.data());
    ok = minimizer->Minimize();
    hesse = config.fit.hesse && ok ? minimizer->Hesse() : false;
    status = minimizer->Status();
    chi2 = minimizer->MinValue();
    startIndex = static_cast<int>(firstStart + localStart);
    const double* values = minimizer->X();
    ok = ok && problem->massMapping.Valid(values) && std::isfinite(chi2) && chi2 < 1.e299;
    std::copy(values, values + fitSize, fitResult.begin());
    const auto physical = problem->massMapping.Evaluate(values);
    std::copy(physical.values.begin(), physical.values.end(), result.begin());
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::fill(errors.begin(), errors.end(), missing);
    std::fill(fitErrors.begin(), fitErrors.end(), missing);
    fitCovariance.clear();
    physicalCovariance.clear();
    massJacobian.clear();
    covarianceOk = false;
    if (hesse) {
      fitCovariance.resize(fitSize * fitSize);
      covarianceOk = minimizer->GetCovMatrix(fitCovariance.data());
      if (covarianceOk) {
        physicalCovariance = MassModels::OrderedMasses::Propagate(physical, fitCovariance, fitSize);
        massJacobian.assign(physicalSize * fitSize, 0.0);
        for (size_t i = 0; i < physicalSize; ++i) {
          errors[i] = std::sqrt(std::max(0.0, physicalCovariance[i * physicalSize + i]));
          for (const auto& entry : physical.jacobian[i])
            massJacobian[i * fitSize + entry.first] = entry.second;
        }
        for (size_t i = 0; i < fitSize; ++i)
          fitErrors[i] = std::sqrt(std::max(0.0, fitCovariance[i * fitSize + i]));
      } else {
        fitCovariance.clear();
      }
    }
    for (size_t i = 0; i < orderedMasses.size(); ++i) {
      const auto& node = problem->massMapping.nodes[i];
      const int p = node.physical, q = node.previous;
      orderedMasses[i] = result[p] - (q < 0 ? 0.0 : result[q]);
      orderedErrors[i] = missing;
      if (covarianceOk) {
        double variance = physicalCovariance[p * physicalSize + p];
        if (q >= 0)
          variance += physicalCovariance[q * physicalSize + q] -
                      2.0 * physicalCovariance[p * physicalSize + q];
        orderedErrors[i] = std::sqrt(std::max(0.0, variance));
      }
    }
    summary.Fill();
    progress.Advance();
  }
  summary.Write();
  WriteMassMetadata(*problem);
  TParameter<int>("mass_fit_analytic_gradient",
                  config.fit.useNumericalGradients ? 0 : 1, 'f')
      .Write();
  // Gauge metadata is identical across workers; preserve it when merging.
  // Keep the original metadata key for consumers expecting the first gauge.
  TParameter<int>("mass_fit_phase_gauge_coordinate", phaseGauges.front(), 'f').Write();
  TParameter<int>("mass_fit_phase_gauge_count", static_cast<int>(phaseGauges.size()), 'f')
      .Write();
  for (size_t i = 0; i < phaseGauges.size(); ++i)
    TParameter<int>(("mass_fit_phase_gauge_coordinate_" + std::to_string(i)).c_str(),
                    phaseGauges[i], 'f').Write();
  TParameter<int>("mass_fit_points", static_cast<int>(problem->data.size()), 'f').Write();
  file->Close();
}

void WriteGlobalBestPrediction(const FitMassModelConfig& config,
                               const ModelConfig& model) {
  auto problem = BuildProblem(config, model);
  std::unique_ptr<TFile> file(TFile::Open(config.fit.output.c_str(), "UPDATE"));
  if (!file || file->IsZombie())
    throw std::runtime_error("Could not reopen merged mass-fit output " +
                             config.fit.output.string());
  auto* summary = dynamic_cast<TTree*>(file->Get("massFitResults"));
  if (!summary || summary->GetEntries() == 0)
    throw std::runtime_error("Merged mass-fit output has no start records");
  bool ok = false;
  double chi2 = std::numeric_limits<double>::infinity();
  summary->SetBranchAddress("fit_ok", &ok);
  summary->SetBranchAddress("chi2", &chi2);
  std::vector<double> values(problem.coordinates.size());
  for (size_t i = 0; i < values.size(); ++i)
    summary->SetBranchAddress(problem.coordinates[i].name.c_str(), &values[i]);
  double bestChi2 = std::numeric_limits<double>::infinity();
  std::vector<double> best;
  for (Long64_t entry = 0; entry < summary->GetEntries(); ++entry) {
    summary->GetEntry(entry);
    if (ok && std::isfinite(chi2) && chi2 < bestChi2) {
      bestChi2 = chi2;
      best = values;
    }
  }
  if (best.empty())
    throw std::runtime_error("No successful finite mass-fit minimum was recorded");
  file->cd();
  std::vector<double> masses;
  for (const auto& point : problem.data) masses.push_back(point.mass);
  WritePredictions(problem, best, config.massBranch, "massFitPredictions", masses);
  if (config.curvePoints > 0) {
    const auto bounds = std::minmax_element(masses.begin(), masses.end());
    std::vector<double> curveMasses(config.curvePoints);
    for (size_t i = 0; i < curveMasses.size(); ++i)
      curveMasses[i] = *bounds.first + (*bounds.second - *bounds.first) *
          (curveMasses.size() > 1 ? static_cast<double>(i) / (curveMasses.size() - 1) : 0.0);
    WritePredictions(problem, best, config.massBranch, "massFitCurves", curveMasses);
  }
  TParameter<double>("mass_fit_best_chi2", bestChi2).Write();
  TParameter<Long64_t>("mass_fit_recorded_starts", summary->GetEntries()).Write();
  file->Close();
}

} // namespace

void FitMassModel(const FitMassModelConfig& config, const ModelConfig& model) {
  if (config.fit.input.empty() || config.fit.output.empty())
    throw std::invalid_argument("Mass fit requires input and output paths");
  if (config.fit.starts == 0)
    throw std::invalid_argument("Mass fit starts must be positive");
  if (config.maxFunctionCalls == 0)
    throw std::invalid_argument("Mass fit maxFunctionCalls must be positive");
  std::filesystem::create_directories(config.fit.output.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : config.fit.output.parent_path());
  unsigned workers = config.fit.workers == 0
                         ? std::max(1u, std::thread::hardware_concurrency())
                         : config.fit.workers;
  workers = std::min(workers, config.fit.starts);
  if (workers == 1) {
    detail::ProgressMonitor progress(config.fit.verbose, config.fit.starts, workers,
                                     "mass-fit");
    RunMassWorker(config, model, config.fit.output, 0, 0, progress.Counter());
    WriteGlobalBestPrediction(config, model);
    return;
  }
  std::vector<unsigned> ids(workers), offsets(workers);
  unsigned offset = 0;
  for (unsigned worker = 0; worker < workers; ++worker) {
    ids[worker] = worker;
    offsets[worker] = offset;
    offset += config.fit.starts / workers + (worker < config.fit.starts % workers);
  }
  ROOT::TProcessExecutor pool(workers);
  detail::ProgressMonitor progress(config.fit.verbose, config.fit.starts, workers,
                                   "mass-fit");
  const auto counter = progress.Counter();
  const auto parts = pool.Map(
      [=](unsigned worker) {
        FitMassModelConfig child = config;
        child.fit.workers = 1;
        child.fit.starts =
            config.fit.starts / workers + (worker < config.fit.starts % workers);
        child.fit.seed = config.fit.seed == 0 ? 0x9e3779b9u + worker * 100003u
                                              : config.fit.seed + worker * 100003u;
        child.fit.verbose = false;
        const auto part = detail::WorkerPartFile(config.fit.output, worker);
        RunMassWorker(child, model, part, worker, offsets[worker], counter);
        return part.string();
      },
      ids);
  detail::MergeWorkerFiles(config.fit.output, parts, false);
  WriteGlobalBestPrediction(config, model);
  for (const auto& part : parts)
    std::filesystem::remove(part);
  if (config.fit.verbose)
    std::cout << "Saved " << config.fit.output << " (" << workers << " workers)\n";
}

void GenerateMassModelMoments(const MassMomentGenerationConfig& config) {
  if (config.massesGeV.empty())
    throw std::invalid_argument("Mass generator needs at least one mass");
  if (!std::isfinite(config.relativeUncertainty) ||
      !std::isfinite(config.absoluteUncertainty) || config.relativeUncertainty < 0.0 ||
      config.absoluteUncertainty <= 0.0)
    throw std::invalid_argument("Invalid mass-generator uncertainties");
  std::filesystem::create_directories(config.output.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : config.output.parent_path());
  FitConfig fit;
  fit.input = config.output;
  fit.tree = config.tree;
  fit.epsilon = config.epsilon;
  fit.photoproduction = config.photoproduction;
  fit.verbose = false;
  const auto internal = detail::MakeInternalConfig(fit, config.model);
  auto angular = detail::BuildSyntheticContext(internal);
  auto problem = BuildAmplitudeProblem(config.massModel, angular);
  std::vector<double> truth;
  for (const auto& c : problem.coordinates)
    truth.push_back(c.control.initial);
  std::unique_ptr<TFile> output(TFile::Open(config.output.c_str(), "RECREATE"));
  if (!output || output->IsZombie())
    throw std::runtime_error("Could not create mass moments output");
  TTree tree(config.tree.c_str(), "Synthetic mass-dependent moments");
  double mass = 0.0;
  tree.Branch(config.massBranch.c_str(), &mass);
  std::vector<double> value(angular->observed.size()), error(angular->observed.size());
  for (size_t i = 0; i < value.size(); ++i) {
    tree.Branch(angular->observed[i].name.c_str(), &value[i]);
    tree.Branch((angular->observed[i].name + "_err").c_str(), &error[i]);
  }
  std::vector<Complex> amplitudes;
  std::vector<std::vector<Complex>> derivatives;
  std::vector<double> gradient, secondaryGradient;
  for (double pointMass : config.massesGeV) {
    if (!std::isfinite(pointMass) || pointMass <= 0.0)
      throw std::invalid_argument("Generated masses must be finite and positive");
    mass = pointMass;
    AmplitudesAndDerivatives(problem, mass, truth.data(), amplitudes, derivatives);
    for (size_t i = 0; i < value.size(); ++i) {
      ObservedMomentAndGradient(problem, i, amplitudes, derivatives, value[i], gradient,
                                secondaryGradient);
      error[i] = std::max(config.absoluteUncertainty,
                          std::abs(value[i]) * config.relativeUncertainty);
    }
    tree.Fill();
  }
  tree.Write();
  const auto physical = problem.massMapping.Evaluate(truth.data());
  for (size_t i = 0; i < physical.values.size(); ++i)
    TParameter<double>(problem.massMapping.physicalCoordinates[i].name.c_str(),
                       physical.values[i]).Write();
  for (size_t i = 0; i < truth.size(); ++i)
    if (problem.coordinates[i].name.rfind("ordered__", 0) == 0 ||
        problem.coordinates[i].name.rfind("positive__", 0) == 0)
      TParameter<double>(problem.coordinates[i].name.c_str(), truth[i]).Write();
  WriteMassMetadata(problem);
  output->Close();
}

} // namespace emi
