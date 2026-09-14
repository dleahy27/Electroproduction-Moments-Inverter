#include "emi/Runner.h"

#include "Detail.h"
#include "MassModels/PhotoTest.h"

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TObjString.h"
#include "TParameter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace emi {
namespace {
using namespace detail;

static void FillExampleAmplitudes(const EvaluationContext& ctx,
                                  std::vector<double>& fullVals) {
  std::vector<double> freeValues(ctx.freeToFull.size(), 0.0);
  for (size_t i = 0; i < ctx.freeToFull.size(); ++i) {
    const auto& parameter = ctx.fullPars[static_cast<size_t>(ctx.freeToFull[i])];
    if (parameter.phase) {
      const double candidate =
          -2.4 + 4.8 * static_cast<double>((i * 7) % 19) / 18.0;
      freeValues[i] = std::clamp(candidate, parameter.low, parameter.high);
    } else {
      freeValues[i] = 0.025 + 0.005 * static_cast<double>((i * 5) % 7);
    }
  }
  if (!FillFullParameters(ctx, freeValues.data(), fullVals)) {
    throw std::runtime_error("Could not normalise the generated amplitude point");
  }
}

static std::string AmplitudeName(char reflectivity, char orientation,
                                 int l, int m, int k, bool phase) {
  auto index = [](int value) {
    return value < 0 ? "m" + std::to_string(-value)
                     : std::to_string(value);
  };
  std::string name(1, reflectivity);
  if (phase) name += "phi";
  name += '_';
  name += orientation;
  name += '_';
  name += std::to_string(l) + '_' + index(m);
  if (k != 0) name += '_' + index(k);
  return name;
}

static bool IsDominantKMode(FixedGenerationMode mode) {
  return mode == FixedGenerationMode::KPositiveDominant ||
         mode == FixedGenerationMode::KReflectivitySplit;
}

struct PhotoTestFillInfo {
  double referencePhaseBeforeRotation = 0.0;
  double rawH000 = 0.0;
  double amplitudeScale = 1.0;
};

static void ValidateRandomRange(double minimum, double maximum,
                                double magnitudeLimit) {
  if (!std::isfinite(minimum) || !std::isfinite(maximum) ||
      minimum < 0.0 || maximum < minimum || maximum > magnitudeLimit) {
    throw std::invalid_argument("Invalid random-amplitude range");
  }
}

static void FillCustomAmplitudes(const EvaluationContext& ctx,
                                 const FixedMomentsConfig& generation,
                                 TRandom3& random,
                                 std::vector<double>& fullVals) {
  fullVals.resize(ctx.fullPars.size());
  std::unordered_map<std::string, int> parameterIndices;
  for (int i = 0; i < static_cast<int>(ctx.fullPars.size()); ++i) {
    fullVals[static_cast<size_t>(i)] = ctx.fullPars[static_cast<size_t>(i)].init;
    parameterIndices.emplace(ctx.fullPars[static_cast<size_t>(i)].name, i);
  }

  std::unordered_set<std::string> assigned;
  const bool usesK =
      ctx.cfg.nucleonPolarization != NucleonPolarization::None;

  auto assign = [&](char reflectivity, char orientation, int l, int m, int k,
                    double magnitude, double phase, bool randomPhase) {
    const std::string inputName =
        AmplitudeName(reflectivity, orientation, l, m, k, false);
    if ((reflectivity != 'a' && reflectivity != 'b') ||
        (orientation != 'T' && orientation != 'L')) {
      throw std::invalid_argument(
          "Invalid reflectivity or orientation in " + inputName);
    }
    if ((usesK && k != -1 && k != 1) || (!usesK && k != 0)) {
      throw std::invalid_argument(
          inputName + (usesK ? " must use k=+1 or k=-1" : " must use k=0"));
    }
    if (!std::isfinite(magnitude) || magnitude < 0.0 ||
        magnitude > ctx.cfg.magnitudeMax || !std::isfinite(phase) ||
        phase < -kPi || phase > kPi) {
      throw std::invalid_argument("Invalid magnitude or phase for " + inputName);
    }
    if (ctx.cfg.photoproduction && orientation == 'L') {
      throw std::invalid_argument(
          "Photoproduction has no longitudinal amplitude " + inputName);
    }
    if (ctx.cfg.enforceLongitudinalParity && orientation == 'L') {
      m = std::abs(m);
    }

    const std::string magnitudeName =
        AmplitudeName(reflectivity, orientation, l, m, k, false);
    const std::string phaseName =
        AmplitudeName(reflectivity, orientation, l, m, k, true);
    if (!assigned.insert(magnitudeName).second) {
      throw std::invalid_argument(
          "Amplitude " + inputName +
          " cannot be both fixed and random or listed twice");
    }
    const auto magnitudeIt = parameterIndices.find(magnitudeName);
    const auto phaseIt = parameterIndices.find(phaseName);
    if (magnitudeIt == parameterIndices.end() ||
        phaseIt == parameterIndices.end()) {
      throw std::invalid_argument(
          "Amplitude " + inputName + " is not in the selected model");
    }

    const int magnitudeIndex = magnitudeIt->second;
    const int phaseIndex = phaseIt->second;
    if (magnitudeIndex == ctx.normalisedMagFullIdx) {
      throw std::invalid_argument(
          "Do not list normalization amplitude " + inputName +
          "; EMI derives its magnitude");
    }
    const auto& magnitudeParameter =
        ctx.fullPars[static_cast<size_t>(magnitudeIndex)];
    if (magnitudeParameter.fixed && magnitude != magnitudeParameter.init) {
      throw std::invalid_argument(
          "Amplitude " + inputName + " is fixed to zero by the model");
    }
    if (ctx.fullPars[static_cast<size_t>(phaseIndex)].fixed) {
      if (!randomPhase && std::abs(phase) > 1e-12) {
        throw std::invalid_argument(
            "The reference phase of " + inputName + " must be zero");
      }
      phase = 0.0;
    }

    fullVals[static_cast<size_t>(magnitudeIndex)] = magnitude;
    fullVals[static_cast<size_t>(phaseIndex)] = phase;
  };

  for (const auto& amplitude : generation.fixedAmplitudes) {
    assign(amplitude.reflectivity, amplitude.orientation,
           amplitude.l, amplitude.m, amplitude.k,
           amplitude.magnitude, amplitude.phase, false);
  }
  for (const auto& amplitude : generation.randomAmplitudes) {
    ValidateRandomRange(amplitude.minimum, amplitude.maximum,
                        ctx.cfg.magnitudeMax);
    assign(amplitude.reflectivity, amplitude.orientation,
           amplitude.l, amplitude.m, amplitude.k,
           random.Uniform(amplitude.minimum, amplitude.maximum),
           random.Uniform(-kPi, kPi), true);
  }

  if (!ApplyAmplitudeNormalisation(ctx, fullVals)) {
    throw std::invalid_argument(
        "Specified amplitudes exceed the requested normalization");
  }
}

static double ModeScale(FixedGenerationMode mode,
                        char reflectivity, int k, double suppression) {
  if (mode == FixedGenerationMode::KPositiveDominant) {
    return k == 1 ? 1.0 : suppression;
  }
  if (mode == FixedGenerationMode::KReflectivitySplit) {
    const bool dominant = (reflectivity == 'a' && k == 1) ||
                          (reflectivity == 'b' && k == -1);
    return dominant ? 1.0 : suppression;
  }
  return 1.0;
}

static void FillRandomModeAmplitudes(const EvaluationContext& ctx,
                                     const FixedMomentsConfig& generation,
                                     TRandom3& random,
                                     std::vector<double>& fullVals) {
  if (!generation.fixedAmplitudes.empty() ||
      !generation.randomAmplitudes.empty()) {
    throw std::invalid_argument(
        "Amplitude lists must be empty in a non-custom generation mode");
  }
  if (IsDominantKMode(generation.mode) &&
      ctx.cfg.nucleonPolarization != NucleonPolarization::Both) {
    throw std::invalid_argument(
        "K generation modes require both nucleon polarizations");
  }
  ValidateRandomRange(generation.randomMinimum, generation.randomMaximum,
                      ctx.cfg.magnitudeMax);
  if (IsDominantKMode(generation.mode) &&
      (!std::isfinite(generation.suppression) ||
       generation.suppression < 0.0 || generation.suppression > 1.0)) {
    throw std::invalid_argument(
        "Suppression must be between zero and one");
  }

  fullVals.assign(ctx.fullPars.size(), 0.0);
  std::unordered_map<long long, double> pairedMagnitudes;
  for (size_t i = 0; i < ctx.fullPars.size(); ++i) {
    const auto& parameter = ctx.fullPars[i];
    const auto label = ParseParameterLabel(parameter.name);
    if (parameter.phase) {
      fullVals[i] = parameter.fixed ? parameter.init
                                    : random.Uniform(-kPi, kPi);
      continue;
    }
    if (parameter.fixed && static_cast<int>(i) != ctx.normalisedMagFullIdx) {
      fullVals[i] = parameter.init;
      continue;
    }

    double magnitude;
    if (generation.mode == FixedGenerationMode::KPositiveDominant ||
        generation.mode == FixedGenerationMode::KReflectivitySplit) {
      const long long key = MakeParameterKey(
          label.reflectivity, label.orientation, 0, label.l, label.m, false);
      const auto position = pairedMagnitudes.emplace(
          key, random.Uniform(generation.randomMinimum,
                              generation.randomMaximum)).first;
      magnitude = position->second;
    } else {
      magnitude = random.Uniform(generation.randomMinimum,
                                 generation.randomMaximum);
    }
    fullVals[i] = magnitude * ModeScale(
        generation.mode, label.reflectivity, label.k,
        generation.suppression);
  }

  std::vector<double> rawMoments;
  EvaluateAllMoments(ctx, fullVals, rawMoments);
  double rawNorm = rawMoments[static_cast<size_t>(ctx.idxH0_00)];
  if (ctx.idxH4_00 >= 0) {
    rawNorm += ctx.cfg.epsilon *
               rawMoments[static_cast<size_t>(ctx.idxH4_00)];
  }
  if (!(rawNorm > 0.0)) {
    throw std::runtime_error("Random amplitude point has zero normalization");
  }
  const double scale =
      std::sqrt(ctx.cfg.normalisationMomentTarget / rawNorm);
  for (size_t i = 0; i < ctx.fullPars.size(); ++i) {
    if (!ctx.fullPars[i].phase) fullVals[i] *= scale;
  }
}

static PhotoTestFillInfo FillPhotoTestAmplitudes(
    const EvaluationContext& ctx,
    const FixedMomentsConfig& generation,
    std::vector<double>& fullVals) {
  if (!generation.photoproduction || !ctx.cfg.photoproduction) {
    throw std::invalid_argument(
        "PhotoTest is a photoproduction mass model");
  }
  if (ctx.cfg.nucleonPolarization != NucleonPolarization::Both) {
    throw std::invalid_argument(
        "PhotoTest requires NucleonPolarization::Both for generated truth");
  }
  if (!ctx.cfg.usePositiveReflectivity ||
      !ctx.cfg.useNegativeReflectivity) {
    throw std::invalid_argument(
        "PhotoTest requires both reflectivities");
  }
  if (!generation.fixedAmplitudes.empty() ||
      !generation.randomAmplitudes.empty()) {
    throw std::invalid_argument(
        "Amplitude lists must be empty in PhotoTest mode");
  }
  if (!generation.massModel.massGeV) {
    throw std::invalid_argument(
        "PhotoTest requires an invariant mass in GeV");
  }

  constexpr std::array<Wave, 9> requiredWaves{{
      {0, 0},
      {1, -1}, {1, 0}, {1, +1},
      {2, -2}, {2, -1}, {2, 0}, {2, +1}, {2, +2},
  }};
  for (const auto& required : requiredWaves) {
    const bool present = std::any_of(
        ctx.cfg.waves.begin(), ctx.cfg.waves.end(),
        [&](const Wave& wave) {
          return wave.l == required.l && wave.m == required.m;
        });
    if (!present) {
      throw std::invalid_argument(
          "PhotoTest requires wave (l,m)=(" +
          std::to_string(required.l) + "," +
          std::to_string(required.m) + ")");
    }
  }

  MassModels::PhotoTest model({
      generation.massModel.backgroundEnabled,
      generation.massModel.kMinusScale,
  });
  auto amplitudes = model.Evaluate(*generation.massModel.massGeV);

  const auto reference = std::find_if(
      amplitudes.begin(), amplitudes.end(),
      [](const MassModels::ComplexAmplitude& amplitude) {
        const auto& key = amplitude.key;
        return key.reflectivity == 'a' && key.orientation == 'T' &&
               key.l == 2 && key.m == 2 && key.k == 1;
      });
  if (reference == amplitudes.end() ||
      !(std::abs(reference->value) > 1e-14)) {
    throw std::runtime_error(
        "PhotoTest phase reference a_T_2_2_1 is zero");
  }

  PhotoTestFillInfo info;
  info.referencePhaseBeforeRotation = std::arg(reference->value);
  const std::complex<double> phaseRotation =
      std::polar(1.0, -info.referencePhaseBeforeRotation);
  for (auto& amplitude : amplitudes) {
    amplitude.value *= phaseRotation;
  }

  std::unordered_map<long long, int> parameterIndices;
  for (int i = 0; i < static_cast<int>(ctx.fullPars.size()); ++i) {
    const auto& parameter = ctx.fullPars[static_cast<size_t>(i)];
    const auto label = ParseParameterLabel(parameter.name);
    if (!label.valid) {
      throw std::runtime_error(
          "Could not parse amplitude parameter " + parameter.name);
    }
    parameterIndices.emplace(
        MakeParameterKey(label.reflectivity, label.orientation, label.k,
                         label.l, label.m, parameter.phase),
        i);
  }

  fullVals.assign(ctx.fullPars.size(), 0.0);
  std::unordered_set<long long> assignedMagnitudes;
  for (const auto& amplitude : amplitudes) {
    const auto& key = amplitude.key;
    const long long magnitudeKey = MakeParameterKey(
        key.reflectivity, key.orientation, key.k, key.l, key.m, false);
    const long long phaseKey = MakeParameterKey(
        key.reflectivity, key.orientation, key.k, key.l, key.m, true);
    const auto magnitude = parameterIndices.find(magnitudeKey);
    const auto phase = parameterIndices.find(phaseKey);
    if (magnitude == parameterIndices.end() ||
        phase == parameterIndices.end()) {
      throw std::invalid_argument(
          "PhotoTest amplitude " +
          AmplitudeName(key.reflectivity, key.orientation, key.l, key.m,
                        key.k, false) +
          " is not in the selected generation model");
    }
    if (!assignedMagnitudes.insert(magnitudeKey).second) {
      throw std::runtime_error("PhotoTest produced a duplicate amplitude key");
    }

    const double absoluteValue = std::abs(amplitude.value);
    fullVals[static_cast<size_t>(magnitude->second)] = absoluteValue;
    const bool isReference =
        key.reflectivity == 'a' && key.orientation == 'T' &&
        key.l == 2 && key.m == 2 && key.k == 1;
    fullVals[static_cast<size_t>(phase->second)] =
        absoluteValue == 0.0 || isReference
            ? 0.0 : std::arg(amplitude.value);
  }
  if (assignedMagnitudes.size() != amplitudes.size() ||
      amplitudes.size() != 36) {
    throw std::runtime_error(
        "PhotoTest did not produce its complete 36-amplitude basis");
  }

  std::vector<double> rawMoments;
  EvaluateAllMoments(ctx, fullVals, rawMoments);
  if (ctx.idxH0_00 < 0) {
    throw std::runtime_error("PhotoTest could not find H_0_0_0_0_0");
  }
  info.rawH000 = rawMoments[static_cast<size_t>(ctx.idxH0_00)];
  if (!(info.rawH000 > 0.0) || !std::isfinite(info.rawH000)) {
    throw std::runtime_error(
        "PhotoTest has a non-positive raw normalization moment");
  }
  info.amplitudeScale = std::sqrt(
      ctx.cfg.normalisationMomentTarget / info.rawH000);
  for (size_t i = 0; i < ctx.fullPars.size(); ++i) {
    if (!ctx.fullPars[i].phase) fullVals[i] *= info.amplitudeScale;
  }
  return info;
}

static std::optional<PhotoTestFillInfo> FillConfiguredAmplitudes(
    const EvaluationContext& ctx,
    const FixedMomentsConfig& generation,
    std::vector<double>& fullVals) {
  TRandom3 random(generation.seed);
  if (generation.seed == 0) random.SetSeed(0);
  if (generation.mode == FixedGenerationMode::Custom) {
    FillCustomAmplitudes(ctx, generation, random, fullVals);
    return std::nullopt;
  }
  if (generation.mode == FixedGenerationMode::PhotoTest) {
    return FillPhotoTestAmplitudes(ctx, generation, fullVals);
  }
  FillRandomModeAmplitudes(ctx, generation, random, fullVals);
  return std::nullopt;
}

static void PrintAmplitudes(const EvaluationContext& ctx,
                            const std::vector<double>& fullVals) {
  std::unordered_map<std::string, double> values;
  for (size_t i = 0; i < ctx.fullPars.size(); ++i) {
    values.emplace(ctx.fullPars[i].name, fullVals[i]);
  }

  std::cout << "Amplitudes used for synthetic moments\n"
            << "-------------------------------------\n";
  for (size_t i = 0; i < ctx.fullPars.size(); ++i) {
    const auto& parameter = ctx.fullPars[i];
    if (parameter.phase) continue;
    const auto label = ParseParameterLabel(parameter.name);
    const double magnitude = fullVals[i];
    std::cout << std::setw(18) << parameter.name << " = " << std::setw(12)
              << std::setprecision(8) << magnitude;
    if (std::abs(magnitude) > 1e-14) {
      const auto phase = values.find(AmplitudeName(
          label.reflectivity, label.orientation, label.l, label.m,
          label.k, true));
      if (phase != values.end()) {
        std::cout << "   phase = " << std::setw(12) << phase->second;
      }
    }
    if (static_cast<int>(i) == ctx.normalisedMagFullIdx) {
      std::cout << "   [normalization]";
    } else if (parameter.fixed) {
      std::cout << "   [fixed by model]";
    }
    std::cout << '\n';
  }
}

static std::vector<std::string> RequestedMomentSeedBranches(const InternalConfig& cfg) {
  std::vector<std::string> names;
  int maximumL = 0;
  for (const auto& wave : cfg.waves) maximumL = std::max(maximumL, wave.l);
  const int betaMax = (cfg.nucleonPolarization == NucleonPolarization::Initial ||
                       cfg.nucleonPolarization == NucleonPolarization::Both) ? 3 : 0;
  const int deltaMax = (cfg.nucleonPolarization == NucleonPolarization::Recoil ||
                        cfg.nucleonPolarization == NucleonPolarization::Both) ? 3 : 0;

  auto add = [&](const std::string& valueName) {
    names.push_back(valueName);
    names.push_back(valueName + "_err");
  };
  auto suffix = [&](int beta, int delta, int L, int M) {
    std::string value;
    if (cfg.nucleonPolarization != NucleonPolarization::None) {
      value = std::to_string(beta) + '_' + std::to_string(delta) + '_';
    }
    return value + std::to_string(L) + '_' + std::to_string(M);
  };
  for (int beta = 0; beta <= betaMax; ++beta) {
    for (int delta = 0; delta <= deltaMax; ++delta) {
      for (int L = 0; L <= 2 * maximumL; ++L) {
        for (int M = 0; M <= L; ++M) {
          add((cfg.photoproduction ? "RH_0_" : "RH04_") +
              suffix(beta, delta, L, M));
        }
      }
    }
  }
  const int maximumAlpha = cfg.photoproduction ? 3 : 8;
  for (int alpha = 1; alpha <= maximumAlpha; ++alpha) {
    if (alpha == 4) continue;
    for (int beta = 0; beta <= betaMax; ++beta) {
      for (int delta = 0; delta <= deltaMax; ++delta) {
        const int parity = ((alpha == 2 || alpha == 3 || alpha == 6 || alpha == 7) ? -1 : 1)
                         * (beta >= 2 ? -1 : 1) * (delta >= 2 ? -1 : 1);
        for (int L = 0; L <= 2 * maximumL; ++L) {
          for (int M = 0; M <= L; ++M) {
            if (parity < 0 && M == 0) continue;
            add("RH_" + std::to_string(alpha) + "_" +
                suffix(beta, delta, L, M));
          }
        }
      }
    }
  }
  return names;
}

static void CreateSeedInputFile(const std::string& fileName,
                                const std::string& treeName,
                                const InternalConfig& cfg) {
  std::unique_ptr<TFile> f(TFile::Open(fileName.c_str(), "RECREATE"));
  if (!f || f->IsZombie()) throw std::runtime_error("Failed to create seed file: " + fileName);

  TTree t(treeName.c_str(), "Temporary tree used to initialise the moment model");

  const std::vector<std::string> branchNames = RequestedMomentSeedBranches(cfg);
  if (branchNames.empty()) throw std::runtime_error("No seed moment branch names were requested");

  std::vector<double> storage(branchNames.size(), 0.0);
  for (size_t i = 0; i < branchNames.size(); ++i) {
    if (branchNames[i].size() >= 4 &&
        branchNames[i].compare(branchNames[i].size() - 4, 4, "_err") == 0) {
      storage[i] = 1;
    }
    t.Branch(branchNames[i].c_str(), &storage[i]);
  }

  t.Fill();
  t.Write();
  f->Close();
}

} // namespace

namespace detail {

std::shared_ptr<EvaluationContext> BuildSyntheticContext(
    const InternalConfig& config) {
  CreateSeedInputFile(config.momentsFile, config.momentsTree, config);
  return BuildContext(config);
}

} // namespace detail

static void GenerateFixedPoint(const FixedMomentsConfig& generation,
                               bool useExampleAmplitudes) {
  const std::filesystem::path outPath = generation.output.empty()
      ? std::filesystem::path("InputFiles/Generated/fixed_test.root")
      : generation.output;
  if (!outPath.parent_path().empty()) {
    std::filesystem::create_directories(outPath.parent_path());
  }
  const std::string treeName = "genMoments";

  FitConfig fit;
  fit.input = outPath;
  fit.tree = treeName;
  fit.epsilon = generation.epsilon;
  fit.photoproduction = generation.photoproduction;
  fit.verbose = false;
  InternalConfig cfg = MakeInternalConfig(fit, generation.model);

  auto ctx = BuildSyntheticContext(cfg);

  std::vector<double> fullVals;
  std::optional<PhotoTestFillInfo> photoTestInfo;
  if (useExampleAmplitudes) {
    FillExampleAmplitudes(*ctx, fullVals);
  } else {
    photoTestInfo = FillConfiguredAmplitudes(*ctx, generation, fullVals);
  }

  std::vector<double> Hvals(ctx->modelsRec.size(), 0.0);
  EvaluateAllMoments(*ctx, fullVals, Hvals);

  std::vector<std::string> observedNames;
  std::vector<std::string> errNames;
  std::vector<double> observedVals;
  std::vector<double> errVals(ctx->observed.size(), 1);
  observedNames.reserve(ctx->observed.size());
  errNames.reserve(ctx->observed.size());

  EvaluateObservedMoments(*ctx, Hvals, observedVals);
  for (size_t i = 0; i < ctx->observed.size(); ++i) {
    const auto& ob = ctx->observed[i];
    observedNames.push_back(ob.name);
    errNames.push_back(ob.name + "_err");
  }

  // H04 is the observable supplied to an electroproduction fit. Keep its
  // separate H0 and H4 values as diagnostics.
  std::vector<std::string> raw04Names;
  std::vector<double> raw04Vals;
  if (!cfg.photoproduction) {
    for (size_t i = 0; i < ctx->modelsRec.size(); ++i) {
      const auto& moment = ctx->modelsRec[i];
      if (moment.alpha != 0 && moment.alpha != 4) continue;
      raw04Names.push_back("R" + moment.name);
      raw04Vals.push_back(Hvals[i]);
    }
  }

  std::unique_ptr<TFile> fout(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error("Failed to open output file: " + outPath.string());
  }

  TTree* t = new TTree(treeName.c_str(), "Generated moments built from a fixed amplitude point");
  t->SetDirectory(fout.get());

  if (photoTestInfo) {
    TObjString("PhotoTest").Write("mass_model");
    TParameter<double>("invariant_mass_GeV",
                       *generation.massModel.massGeV).Write();
    TParameter<double>("k_minus_scale",
                       generation.massModel.kMinusScale).Write();
    TParameter<int>("background_enabled",
                    generation.massModel.backgroundEnabled ? 1 : 0).Write();
    TParameter<double>("reference_phase_before_rotation",
                       photoTestInfo->referencePhaseBeforeRotation).Write();
    TParameter<double>("raw_H_0_0_0",
                       photoTestInfo->rawH000).Write();
    TParameter<double>("amplitude_normalisation_scale",
                       photoTestInfo->amplitudeScale).Write();
  }

  std::vector<double> parStorage;
  MakeParameterBranches(t, ctx->fullPars, parStorage);

  std::vector<double> obsStorage;
  MakeNamedBranches(t, observedNames, obsStorage);

  std::vector<double> errStorage;
  MakeNamedBranches(t, errNames, errStorage);

  std::vector<double> raw04Storage;
  MakeNamedBranches(t, raw04Names, raw04Storage);

  std::vector<std::string> q2Name = {"Q2"};
  std::vector<double> q2Storage;
  MakeNamedBranches(t, q2Name, q2Storage);

  parStorage = fullVals;
  obsStorage = observedVals;
  errStorage = errVals;
  raw04Storage = raw04Vals;
  q2Storage.assign(q2Storage.size(), cfg.photoproduction ? 0.0 : 1.0);

  t->Fill();
  fout->Write();
  fout->Close();

  if (generation.printAmplitudes) {
    PrintAmplitudes(*ctx, fullVals);
    std::cout << "Saved generated moments to " << outPath << std::endl;
  }
}

void GenerateFixedMoments(const FixedMomentsConfig& config) {
  GenerateFixedPoint(config, false);
}

void GenerateFixedMoments(const std::filesystem::path& output,
                          double epsilon,
                          bool printAmplitudes,
                          const ModelConfig& model,
                          bool photoproduction) {
  FixedMomentsConfig config;
  config.output = output;
  config.epsilon = epsilon;
  config.printAmplitudes = printAmplitudes;
  config.model = model;
  config.photoproduction = photoproduction;
  GenerateFixedPoint(config, true);
}

} // namespace emi
