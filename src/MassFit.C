// Simultaneous mass-dependent moment fit.  The angular algebra is shared with
// the ordinary inversion; this layer replaces its free bin-by-bin complex
// amplitudes with coherent Breit-Wigner and constant-background functions.
#include "emi/Runner.h"

#include "Detail.h"
#include "Progress.h"
#include "WorkerOutput.h"

#include "Math/Factory.h"
#include "ROOT/TProcessExecutor.hxx"
#include "TFile.h"
#include "TLeaf.h"
#include "TParameter.h"
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
struct Coordinate {
  std::string name;
  MassParameter control;
};
// Indices into the wave and coordinate tables, built once before fitting.
struct ResonanceTerm {
  int wave = -1, mass = -1, width = -1, magnitude = -1, phase = -1;
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
  std::vector<Coordinate> coordinates;
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

std::vector<DataPoint> ReadMassData(const FitMassModelConfig& config,
                                    const detail::EvaluationContext& angular) {
  std::unique_ptr<TFile> file(TFile::Open(config.fit.input.c_str(), "READ"));
  if (!file || file->IsZombie())
    throw std::runtime_error("Failed to open " + config.fit.input.string());
  auto* tree = dynamic_cast<TTree*>(file->Get(config.fit.tree.c_str()));
  if (!tree || tree->GetEntries() == 0)
    throw std::runtime_error("Mass input tree is missing or empty");
  std::vector<DataPoint> result;
  result.reserve(static_cast<size_t>(tree->GetEntries()));
  for (Long64_t row = 0; row < tree->GetEntries(); ++row) {
    // Read each ROOT entry once. All leaves below refer to this loaded row.
    tree->GetEntry(row);
    DataPoint point;
    point.mass = LeafValue(*tree, config.massBranch);
    if (!std::isfinite(point.mass) || point.mass <= 0.0)
      throw std::runtime_error("Mass branch contains an invalid value");
    for (const auto& moment : angular.observed) {
      const double value = LeafValue(*tree, moment.name);
      const double error = LeafValue(*tree, moment.name + "_err");
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
  if (!std::isfinite(value.initial) || !std::isfinite(value.lower) ||
      !std::isfinite(value.upper) || value.lower > value.upper ||
      value.initial < value.lower || value.initial > value.upper) {
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
      if (previous.initial != control.initial || previous.lower != control.lower ||
          previous.upper != control.upper || previous.fixed != control.fixed)
        throw std::invalid_argument("Inconsistent shared pole control " + name);
      return found->second;
    }
    const int index = AddCoordinate(problem, name, control);
    sharedPoleCoordinates.emplace(name, index);
    return index;
  };
  for (const auto& resonance : massModel.resonances) {
    if (resonance.name.empty())
      throw std::invalid_argument("Every mass resonance needs a name");
    if (resonance.poleMass.lower <= 0.0 || resonance.width.lower <= 0.0 ||
        resonance.couplingMagnitude.lower < 0.0)
      throw std::invalid_argument(
          "Mass/width bounds must be positive and magnitudes nonnegative");
    const std::string prefix = "bw_" + resonance.name;
    const std::string couplingPrefix = prefix + "_" + resonance.wave.reflectivity +
                                       (resonance.wave.orientation == 'L' ? "_L" : "") +
                                       "_" + std::to_string(resonance.wave.l) + "_" +
                                       std::to_string(resonance.wave.m) + "_k" +
                                       std::to_string(resonance.wave.k);
    problem.resonances.push_back(
        {addWave(resonance.wave), addShared(prefix + "_mass", resonance.poleMass),
         addShared(prefix + "_width", resonance.width),
         AddCoordinate(problem, couplingPrefix + "_mag", resonance.couplingMagnitude),
         AddCoordinate(problem, couplingPrefix + "_phase", resonance.couplingPhase)});
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
  return problem;
}

Problem BuildProblem(const FitMassModelConfig& config, const ModelConfig& model) {
  if (config.fit.input.empty() || config.fit.output.empty())
    throw std::invalid_argument("Mass fit requires input and output paths");
  if (config.fit.starts == 0)
    throw std::invalid_argument("Mass fit starts must be positive");
  FitConfig angularFit = config.fit;
  angularFit.bin = 0; // BuildContext discovers the schema from entry zero.
  angularFit.verbose = false;
  auto angular = detail::BuildContext(detail::MakeInternalConfig(angularFit, model));
  auto problem = BuildAmplitudeProblem(config.model, std::move(angular));
  problem.data = ReadMassData(config, *problem.angular);
  return problem;
}

void AmplitudesAndDerivatives(const Problem& problem, double mass, const double* x,
                              std::vector<Complex>& amplitude,
                              std::vector<std::vector<Complex>>& derivative) {
  amplitude.assign(problem.waveIndex.size(), {});
  derivative.assign(problem.coordinates.size(),
                    std::vector<Complex>(problem.waveIndex.size(), Complex{}));
  for (const auto& term : problem.resonances) {
    const double pole = x[term.mass], width = x[term.width];
    const double numerator = pole * width;
    const Complex denominator{pole * pole - mass * mass, -pole * width};
    const Complex bw = numerator / denominator;
    const Complex coupling = std::polar(x[term.magnitude], x[term.phase]);
    amplitude[term.wave] += coupling * bw;
    const Complex dPole =
        (width * denominator - numerator * Complex{2.0 * pole, -width}) /
        (denominator * denominator);
    const Complex dWidth = (pole * denominator - numerator * Complex{0.0, -pole}) /
                           (denominator * denominator);
    derivative[term.mass][term.wave] += coupling * dPole;
    derivative[term.width][term.wave] += coupling * dWidth;
    derivative[term.magnitude][term.wave] += std::polar(1.0, x[term.phase]) * bw;
    derivative[term.phase][term.wave] += Complex{0.0, 1.0} * coupling * bw;
  }
  for (const auto& term : problem.backgrounds) {
    amplitude[term.wave] += Complex{x[term.re], x[term.im]};
    derivative[term.re][term.wave] += 1.0;
    derivative[term.im][term.wave] += Complex{0.0, 1.0};
  }
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

double MomentDerivative(const Problem& problem, const detail::MomentModel& moment,
                        const std::vector<Complex>& amplitudes,
                        const std::vector<std::vector<Complex>>& dAmplitude,
                        size_t coordinate) {
  double derivative = 0.0;
  for (const auto& term : moment.terms) {
    const int first = problem.amplitudeForMagnitude[term.idxMag1];
    const int second = problem.amplitudeForMagnitude[term.idxMag2];
    if (first < 0 || second < 0)
      continue;
    const Complex dp = dAmplitude[coordinate][first] * std::conj(amplitudes[second]) +
                       amplitudes[first] * std::conj(dAmplitude[coordinate][second]);
    derivative +=
        term.coeff * (term.trig == detail::TrigKind::kCos ? dp.real() : dp.imag());
  }
  return derivative;
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
    double chi2 = 0.0;
    for (const auto& point : problem_->data) {
      AmplitudesAndDerivatives(*problem_, point.mass, x, amplitude_, derivative_);
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
    for (const auto& point : problem_->data) {
      AmplitudesAndDerivatives(*problem_, point.mass, x, amplitude_, derivative_);
      for (size_t i = 0; i < problem_->angular->observed.size(); ++i) {
        double model = 0.0;
        ObservedMomentAndGradient(*problem_, i, amplitude_, derivative_, model,
                                  momentGradient_, momentGradient2_);
        const double factor =
            -2.0 * (point.values[i] - model) / (point.errors[i] * point.errors[i]);
        for (size_t p = 0; p < momentGradient_.size(); ++p)
          output[p] += factor * momentGradient_[p];
      }
    }
  }
  double DoDerivative(const double* x, unsigned int coordinate) const override {
    if (coordinate >= NDim())
      return 0.0;
    double output = 0.0;
    for (const auto& point : problem_->data) {
      AmplitudesAndDerivatives(*problem_, point.mass, x, amplitude_, derivative_);
      for (size_t i = 0; i < problem_->angular->observed.size(); ++i) {
        const auto& observed = problem_->angular->observed[i];
        const auto& primary =
            problem_->angular->modelsRec[observed.isMixed04
                                             ? problem_->angular->observedModelIdx0[i]
                                             : problem_->angular->observedModelIdx[i]];
        double model = MomentValue(*problem_, primary, amplitude_);
        double derivative =
            MomentDerivative(*problem_, primary, amplitude_, derivative_, coordinate);
        if (observed.isMixed04) {
          const auto& secondary =
              problem_->angular->modelsRec[problem_->angular->observedModelIdx4[i]];
          model += problem_->angular->cfg.epsilon *
                   MomentValue(*problem_, secondary, amplitude_);
          derivative += problem_->angular->cfg.epsilon *
                        MomentDerivative(*problem_, secondary, amplitude_, derivative_,
                                         coordinate);
        }
        output += -2.0 * (point.values[i] - model) * derivative /
                  (point.errors[i] * point.errors[i]);
      }
    }
    return output;
  }

private:
  std::shared_ptr<Problem> problem_;
  mutable std::vector<Complex> amplitude_;
  mutable std::vector<std::vector<Complex>> derivative_;
  mutable std::vector<double> momentGradient_, momentGradient2_;
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
                      const std::string& massBranch) {
  TTree predictions("massFitPredictions",
                    "Moment-model predictions at each mass point");
  double mass = 0.0;
  predictions.Branch(massBranch.c_str(), &mass);
  std::vector<double> values(problem.angular->observed.size());
  for (size_t i = 0; i < values.size(); ++i)
    predictions.Branch(("model__" + problem.angular->observed[i].name).c_str(),
                       &values[i]);
  std::vector<Complex> amplitudes;
  std::vector<std::vector<Complex>> derivatives;
  std::vector<double> gradient, secondaryGradient;
  for (const auto& point : problem.data) {
    mass = point.mass;
    AmplitudesAndDerivatives(problem, mass, x.data(), amplitudes, derivatives);
    for (size_t i = 0; i < values.size(); ++i) {
      ObservedMomentAndGradient(problem, i, amplitudes, derivatives, values[i],
                                gradient, secondaryGradient);
    }
    predictions.Fill();
  }
  predictions.Write();
}

// An explicit gauge wins. Otherwise use the positive-reflectivity transverse
// resonance with the largest (l,m), preserving its configured phase value.
int FindPhaseGauge(const MassModelConfig& massModel, const Problem& problem) {
  int phaseGauge = -1;
  auto matchesGauge = [&](const MassResonanceConfig& resonance) {
    if (massModel.phaseGauge) {
      const auto& gauge = *massModel.phaseGauge;
      return resonance.name == gauge.resonance &&
             WaveKey(resonance.wave) == WaveKey(gauge.wave);
    }
    return resonance.wave.reflectivity == 'a' && resonance.wave.orientation == 'T';
  };
  int bestL = -1, bestM = -999999;
  for (size_t i = 0; i < massModel.resonances.size(); ++i) {
    const auto& resonance = massModel.resonances[i];
    if (!matchesGauge(resonance))
      continue;
    if (massModel.phaseGauge || resonance.wave.l > bestL ||
        (resonance.wave.l == bestL && resonance.wave.m > bestM)) {
      phaseGauge = problem.resonances[i].phase;
      bestL = resonance.wave.l;
      bestM = resonance.wave.m;
    }
  }
  if (phaseGauge < 0)
    throw std::invalid_argument("Mass phase gauge needs a positive-reflectivity "
                                "transverse resonance or an explicit override");
  return phaseGauge;
}

void RunMassWorker(const FitMassModelConfig& config, const ModelConfig& model,
                   const std::filesystem::path& output, unsigned worker,
                   unsigned firstStart, detail::ProgressCounter progress = {}) {
  auto problem = std::make_shared<Problem>(BuildProblem(config, model));
  MassChi2 objective(problem);
  const int phaseGauge = FindPhaseGauge(config.model, *problem);
  unsigned freeCoordinates = 0;
  for (size_t i = 0; i < problem->coordinates.size(); ++i)
    if (!problem->coordinates[i].control.fixed && static_cast<int>(i) != phaseGauge)
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
  std::vector<double> errors(problem->coordinates.size());
  std::vector<double> result(problem->coordinates.size());
  for (size_t i = 0; i < result.size(); ++i) {
    summary.Branch(problem->coordinates[i].name.c_str(), &result[i]);
    summary.Branch(("err__" + problem->coordinates[i].name).c_str(), &errors[i]);
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
      minimizer->SetLimitedVariable(i, coordinate.name, coordinate.control.initial,
                                    1.e-3, coordinate.control.lower,
                                    coordinate.control.upper);
      if (coordinate.control.fixed || static_cast<int>(i) == phaseGauge)
        minimizer->FixVariable(i);
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
          static_cast<int>(i) != phaseGauge)
        value = random.Uniform(coordinate.control.lower, coordinate.control.upper);
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
    result.assign(values, values + problem->coordinates.size());
    errors.assign(problem->coordinates.size(),
                  std::numeric_limits<double>::quiet_NaN());
    if (hesse) {
      const double* fitErrors = minimizer->Errors();
      errors.assign(fitErrors, fitErrors + problem->coordinates.size());
    }
    summary.Fill();
    progress.Advance();
  }
  summary.Write();
  TParameter<int>("mass_fit_analytic_gradient",
                  config.fit.useNumericalGradients ? 0 : 1)
      .Write();
  TParameter<int>("mass_fit_phase_gauge_coordinate", phaseGauge).Write();
  TParameter<int>("mass_fit_points", static_cast<int>(problem->data.size())).Write();
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
  WritePredictions(problem, best, config.massBranch);
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
  for (size_t i = 0; i < truth.size(); ++i)
    TParameter<double>(problem.coordinates[i].name.c_str(), truth[i]).Write();
  output->Close();
}

} // namespace emi
