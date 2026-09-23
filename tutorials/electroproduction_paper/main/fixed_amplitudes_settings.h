#pragma once

// Deterministic full S/P/D electroproduction closure point. The active values
// are the normalized "large longitudinal" test from the legacy
// GenerateMomentsFromFixedAmplitudes.C. Other legacy blocks (random closure,
// large a_T/b_L, large b_T/a_L, large natural/unnatural, and large transverse)
// remain inactive in the legacy source rather than being mixed into this case.

#include "emi/Config.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <vector>

namespace emi::user {

inline ModelConfig FixedModel() {
  ModelConfig model;
  model.SetWaves({
      {0, 0},
      {1, -1}, {1, 0}, {1, 1},
      {2, -2}, {2, -1}, {2, 0}, {2, 1}, {2, 2},
  }).UseReflectivities(true, true)
    .EnforceLongitudinalParity(true)
    .UseNucleonPolarization(NucleonPolarization::None);
  model.normalisationMoment = 2.0;
  return model;
}

// The legacy file stored convenient raw hypotheses and normalized the complete
// vector afterwards. This helper preserves that format for the inactive blocks
// below, then removes a_T_2_2 because EMI derives it during generation.
inline std::vector<FixedAmplitude> NormalizeLegacy(
    std::initializer_list<FixedAmplitude> raw, double epsilon) {
  double norm = 0.0;
  for (const auto& amplitude : raw) {
    const double parityWeight =
        amplitude.orientation == 'L' && amplitude.m != 0 ? 2.0 : 1.0;
    const double polarizationWeight =
        amplitude.orientation == 'L' ? epsilon : 1.0;
    norm += parityWeight * polarizationWeight *
            amplitude.magnitude * amplitude.magnitude;
  }
  const double scale = 1.0 / std::sqrt(norm);
  std::vector<FixedAmplitude> result;
  result.reserve(raw.size());
  for (auto amplitude : raw) {
    const bool normalizationAmplitude =
        amplitude.reflectivity == 'a' && amplitude.orientation == 'T' &&
        amplitude.l == 2 && amplitude.m == 2 && amplitude.k == 0;
    if (normalizationAmplitude) continue;
    amplitude.magnitude *= scale;
    result.push_back(amplitude);
  }
  return result;
}

inline FixedMomentsConfig FixedMoments() {
  FixedMomentsConfig generation;
  generation.photoproduction = false;
  generation.epsilon = 0.8;
  generation.seed = 12345;
  generation.printAmplitudes = true;
  generation.model = FixedModel();
  generation.mode = FixedGenerationMode::Custom;

  // Magnitudes include the legacy common scale 0.3711113985402077. The
  // normalization amplitude a_T_2_2 (truth magnitude 0.007422228) is omitted
  // because EMI derives it from H04_00.
  generation.fixedAmplitudes = {
      {'b', 'T', 0,  0, 0, 0.031544469, -0.763191},
      {'b', 'T', 1,  1, 0, 0.033400026,  0.000000},
      {'b', 'T', 1,  0, 0, 0.028946689,  2.069711},
      {'b', 'T', 1, -1, 0, 0.030431135,  0.744682},
      {'b', 'T', 2,  2, 0, 0.035255583,  0.000000},
      {'b', 'T', 2,  1, 0, 0.032657803,  0.486018},
      {'b', 'T', 2,  0, 0, 0.008164451,  1.285363},
      {'b', 'T', 2, -1, 0, 0.176277914, -2.853670},
      {'b', 'T', 2, -2, 0, 0.018555570, -1.709666},
      {'b', 'L', 1,  1, 0, 0.304311347, -2.640245},
      {'b', 'L', 2,  2, 0, 0.192977927, -1.395033},
      {'b', 'L', 2,  1, 0, 0.241222409,  0.852531},

      {'a', 'T', 0,  0, 0, 0.005566671,  1.223400},
      {'a', 'T', 1,  1, 0, 0.007422228,  0.000000},
      {'a', 'T', 1,  0, 0, 0.003711114,  2.743582},
      {'a', 'T', 1, -1, 0, 0.009277785,  0.930134},
      {'a', 'T', 2,  1, 0, 0.005566671, -2.066297},
      {'a', 'T', 2,  0, 0, 0.003711114,  1.439646},
      {'a', 'T', 2, -1, 0, 0.009277785, -2.114905},
      {'a', 'T', 2, -2, 0, 0.005566671, -0.757404},
      {'a', 'L', 0,  0, 0, 0.207822383, -3.098700},
      {'a', 'L', 1,  1, 0, 0.315444689,  3.075766},
      {'a', 'L', 1,  0, 0, 0.155866787, -0.345900},
      {'a', 'L', 2,  2, 0, 0.352555829,  1.159966},
      {'a', 'L', 2,  1, 0, 0.334000259,  2.154202},
      {'a', 'L', 2,  0, 0, 0.296889119, -1.789620},
  };

  // -----------------------------------------------------------------------
  // Inactive legacy alternatives. To select one, comment out the active
  // assignment above and uncomment exactly one complete assignment below.
  // Values here are the original unscaled values; NormalizeLegacy applies the
  // old common normalization and drops the derived a_T_2_2 magnitude.
  // -----------------------------------------------------------------------

  // Random closure point.
  // generation.fixedAmplitudes = NormalizeLegacy({
  //     {'a', 'T', 0,  0, 0, 0.278093, -0.763191},
  //     {'a', 'L', 0,  0, 0, 0.001290, -1.323314},
  //     {'a', 'T', 1,  1, 0, 0.187800,  0.326981},
  //     {'a', 'T', 1,  0, 0, 0.119613,  2.069711},
  //     {'a', 'T', 1, -1, 0, 0.097077,  0.744682},
  //     {'a', 'L', 1,  1, 0, 0.095088, -2.640245},
  //     {'a', 'L', 1,  0, 0, 0.219784, -1.678924},
  //     {'a', 'T', 2,  2, 0, 0.320298,  0.000000},
  //     {'a', 'T', 2,  1, 0, 0.294303,  0.486018},
  //     {'a', 'T', 2,  0, 0, 0.388018,  1.285363},
  //     {'a', 'T', 2, -1, 0, 0.037811, -2.853670},
  //     {'a', 'T', 2, -2, 0, 0.183498, -1.709666},
  //     {'a', 'L', 2,  2, 0, 0.086477, -1.395033},
  //     {'a', 'L', 2,  1, 0, 0.282641,  0.852531},
  //     {'a', 'L', 2,  0, 0, 0.237000, -0.849284},
  //     {'b', 'T', 0,  0, 0, 0.073704,  1.223400},
  //     {'b', 'T', 1,  1, 0, 0.352029, -1.464122},
  //     {'b', 'T', 1,  0, 0, 0.002826,  2.743582},
  //     {'b', 'T', 1, -1, 0, 0.350459,  0.930134},
  //     {'b', 'L', 1,  1, 0, 0.040336,  1.075766},
  //     {'b', 'T', 2,  2, 0, 0.303627,  0.000000},
  //     {'b', 'T', 2,  1, 0, 0.147978, -2.066297},
  //     {'b', 'T', 2,  0, 0, 0.067619,  1.439646},
  //     {'b', 'T', 2, -1, 0, 0.416301, -2.114905},
  //     {'b', 'T', 2, -2, 0, 0.146388, -0.757404},
  //     {'b', 'L', 2,  2, 0, 0.262566,  1.159966},
  //     {'b', 'L', 2,  1, 0, 0.351028,  2.154202},
  // }, generation.epsilon);

  // Natural transverse / unnatural longitudinal dominant (large a_T, b_L).
  // generation.fixedAmplitudes = NormalizeLegacy({
  //     {'a', 'T', 0,  0, 0, 0.85,  -0.763191},
  //     {'a', 'L', 0,  0, 0, 0.01,  -1.323314},
  //     {'a', 'T', 1,  1, 0, 0.90,   0.000000},
  //     {'a', 'T', 1,  0, 0, 0.78,   2.069711},
  //     {'a', 'T', 1, -1, 0, 0.82,   0.744682},
  //     {'a', 'L', 1,  1, 0, 0.02,  -2.640245},
  //     {'a', 'L', 1,  0, 0, 0.015, -1.678924},
  //     {'a', 'T', 2,  2, 0, 0.95,   0.000000},
  //     {'a', 'T', 2,  1, 0, 0.88,   0.486018},
  //     {'a', 'T', 2,  0, 0, 0.92,   1.285363},
  //     {'a', 'T', 2, -1, 0, 0.75,  -2.853670},
  //     {'a', 'T', 2, -2, 0, 0.80,  -1.709666},
  //     {'a', 'L', 2,  2, 0, 0.02,  -1.395033},
  //     {'a', 'L', 2,  1, 0, 0.025,  0.852531},
  //     {'a', 'L', 2,  0, 0, 0.015, -0.849284},
  //     {'b', 'T', 0,  0, 0, 0.015,  1.223400},
  //     {'b', 'T', 1,  1, 0, 0.02,   0.000000},
  //     {'b', 'T', 1,  0, 0, 0.01,   2.743582},
  //     {'b', 'T', 1, -1, 0, 0.025,  0.930134},
  //     {'b', 'L', 1,  1, 0, 0.85,   3.075766},
  //     {'b', 'T', 2,  2, 0, 0.02,   0.000000},
  //     {'b', 'T', 2,  1, 0, 0.015, -2.066297},
  //     {'b', 'T', 2,  0, 0, 0.01,   1.439646},
  //     {'b', 'T', 2, -1, 0, 0.025, -2.114905},
  //     {'b', 'T', 2, -2, 0, 0.015, -0.757404},
  //     {'b', 'L', 2,  2, 0, 0.95,   1.159966},
  //     {'b', 'L', 2,  1, 0, 0.90,   2.154202},
  // }, generation.epsilon);

  // Unnatural transverse / natural longitudinal dominant (large b_T, a_L).
  // generation.fixedAmplitudes = NormalizeLegacy({
  //     {'b', 'T', 0,  0, 0, 0.85,  -0.763191},
  //     {'b', 'T', 1,  1, 0, 0.90,   0.000000},
  //     {'b', 'T', 1,  0, 0, 0.78,   2.069711},
  //     {'b', 'T', 1, -1, 0, 0.82,   0.744682},
  //     {'b', 'L', 1,  1, 0, 0.02,  -2.640245},
  //     {'b', 'T', 2,  2, 0, 0.95,   0.000000},
  //     {'b', 'T', 2,  1, 0, 0.88,   0.486018},
  //     {'b', 'T', 2,  0, 0, 0.92,   1.285363},
  //     {'b', 'T', 2, -1, 0, 0.75,  -2.853670},
  //     {'b', 'T', 2, -2, 0, 0.80,  -1.709666},
  //     {'b', 'L', 2,  2, 0, 0.02,  -1.395033},
  //     {'b', 'L', 2,  1, 0, 0.025,  0.852531},
  //     {'a', 'T', 0,  0, 0, 0.015,  1.223400},
  //     {'a', 'L', 0,  0, 0, 0.56,  -3.098700},
  //     {'a', 'T', 1,  1, 0, 0.02,   0.000000},
  //     {'a', 'T', 1,  0, 0, 0.01,   2.743582},
  //     {'a', 'T', 1, -1, 0, 0.025,  0.930134},
  //     {'a', 'L', 1,  1, 0, 0.85,   3.075766},
  //     {'a', 'L', 1,  0, 0, 0.42,  -0.345900},
  //     {'a', 'T', 2,  2, 0, 0.02,   0.000000},
  //     {'a', 'T', 2,  1, 0, 0.015, -2.066297},
  //     {'a', 'T', 2,  0, 0, 0.01,   1.439646},
  //     {'a', 'T', 2, -1, 0, 0.025, -2.114905},
  //     {'a', 'T', 2, -2, 0, 0.015, -0.757404},
  //     {'a', 'L', 2,  2, 0, 0.95,   1.159966},
  //     {'a', 'L', 2,  1, 0, 0.90,   2.154202},
  //     {'a', 'L', 2,  0, 0, 0.80,  -1.789620},
  // }, generation.epsilon);

  // Unnatural-reflectivity dominant (large b).
  // generation.fixedAmplitudes = NormalizeLegacy({
  //     {'b', 'T', 0,  0, 0, 0.85,  -0.763191},
  //     {'b', 'T', 1,  1, 0, 0.90,   0.000000},
  //     {'b', 'T', 1,  0, 0, 0.78,   2.069711},
  //     {'b', 'T', 1, -1, 0, 0.82,   0.744682},
  //     {'b', 'L', 1,  1, 0, 0.12,  -2.640245},
  //     {'b', 'T', 2,  2, 0, 0.95,   0.000000},
  //     {'b', 'T', 2,  1, 0, 0.88,   0.486018},
  //     {'b', 'T', 2,  0, 0, 0.92,   1.285363},
  //     {'b', 'T', 2, -1, 0, 0.75,  -2.853670},
  //     {'b', 'T', 2, -2, 0, 0.80,  -1.709666},
  //     {'b', 'L', 2,  2, 0, 0.72,  -1.395033},
  //     {'b', 'L', 2,  1, 0, 0.425,  0.852531},
  //     {'a', 'T', 0,  0, 0, 0.015,  1.223400},
  //     {'a', 'L', 0,  0, 0, 0.056, -3.098700},
  //     {'a', 'T', 1,  1, 0, 0.02,   0.000000},
  //     {'a', 'T', 1,  0, 0, 0.01,   2.743582},
  //     {'a', 'T', 1, -1, 0, 0.025,  0.930134},
  //     {'a', 'L', 1,  1, 0, 0.085,  3.075766},
  //     {'a', 'L', 1,  0, 0, 0.042, -0.345900},
  //     {'a', 'T', 2,  2, 0, 0.02,   0.000000},
  //     {'a', 'T', 2,  1, 0, 0.015, -2.066297},
  //     {'a', 'T', 2,  0, 0, 0.01,   1.439646},
  //     {'a', 'T', 2, -1, 0, 0.025, -2.114905},
  //     {'a', 'T', 2, -2, 0, 0.015, -0.757404},
  //     {'a', 'L', 2,  2, 0, 0.095,  1.159966},
  //     {'a', 'L', 2,  1, 0, 0.090,  2.154202},
  //     {'a', 'L', 2,  0, 0, 0.080, -1.789620},
  // }, generation.epsilon);

  // Natural-reflectivity dominant (large a).
  // generation.fixedAmplitudes = NormalizeLegacy({
  //     {'b', 'T', 0,  0, 0, 0.085, -0.763191},
  //     {'b', 'T', 1,  1, 0, 0.090,  0.000000},
  //     {'b', 'T', 1,  0, 0, 0.078,  2.069711},
  //     {'b', 'T', 1, -1, 0, 0.082,  0.744682},
  //     {'b', 'L', 1,  1, 0, 0.020, -2.640245},
  //     {'b', 'T', 2,  2, 0, 0.095,  0.000000},
  //     {'b', 'T', 2,  1, 0, 0.188,  0.486018},
  //     {'b', 'T', 2,  0, 0, 0.092,  1.285363},
  //     {'b', 'T', 2, -1, 0, 0.075, -2.853670},
  //     {'b', 'T', 2, -2, 0, 0.080, -1.709666},
  //     {'b', 'L', 2,  2, 0, 0.020, -1.395033},
  //     {'b', 'L', 2,  1, 0, 0.025,  0.852531},
  //     {'a', 'T', 0,  0, 0, 0.15,   1.223400},
  //     {'a', 'L', 0,  0, 0, 0.56,  -3.098700},
  //     {'a', 'T', 1,  1, 0, 0.20,   0.000000},
  //     {'a', 'T', 1,  0, 0, 0.31,   2.743582},
  //     {'a', 'T', 1, -1, 0, 0.25,   0.930134},
  //     {'a', 'L', 1,  1, 0, 0.85,   3.075766},
  //     {'a', 'L', 1,  0, 0, 0.42,  -0.345900},
  //     {'a', 'T', 2,  2, 0, 0.52,   0.000000},
  //     {'a', 'T', 2,  1, 0, 0.15,  -2.066297},
  //     {'a', 'T', 2,  0, 0, 0.10,   1.439646},
  //     {'a', 'T', 2, -1, 0, 0.29,  -2.114905},
  //     {'a', 'T', 2, -2, 0, 0.11,  -0.757404},
  //     {'a', 'L', 2,  2, 0, 0.95,   1.159966},
  //     {'a', 'L', 2,  1, 0, 0.90,   2.154202},
  //     {'a', 'L', 2,  0, 0, 0.80,  -1.789620},
  // }, generation.epsilon);
  //
  // Transverse dominant (large T).
  // generation.fixedAmplitudes = NormalizeLegacy({
  //     {'b', 'T', 0,  0, 0, 0.85,  -0.763191},
  //     {'b', 'T', 1,  1, 0, 0.90,   0.000000},
  //     {'b', 'T', 1,  0, 0, 0.78,   2.069711},
  //     {'b', 'T', 1, -1, 0, 0.82,   0.744682},
  //     {'b', 'L', 1,  1, 0, 0.02,  -2.640245},
  //     {'b', 'T', 2,  2, 0, 0.95,   0.000000},
  //     {'b', 'T', 2,  1, 0, 0.88,   0.486018},
  //     {'b', 'T', 2,  0, 0, 0.92,   1.285363},
  //     {'b', 'T', 2, -1, 0, 0.75,  -2.853670},
  //     {'b', 'T', 2, -2, 0, 0.80,  -1.709666},
  //     {'b', 'L', 2,  2, 0, 0.02,  -1.395033},
  //     {'b', 'L', 2,  1, 0, 0.025,  0.852531},
  //     {'a', 'T', 0,  0, 0, 0.455,  1.223400},
  //     {'a', 'L', 0,  0, 0, 0.056, -3.098700},
  //     {'a', 'T', 1,  1, 0, 0.234,  0.000000},
  //     {'a', 'T', 1,  0, 0, 0.12,   2.743582},
  //     {'a', 'T', 1, -1, 0, 0.45,   0.930134},
  //     {'a', 'L', 1,  1, 0, 0.085,  3.075766},
  //     {'a', 'L', 1,  0, 0, 0.042, -0.345900},
  //     {'a', 'T', 2,  2, 0, 0.20,   0.000000},
  //     {'a', 'T', 2,  1, 0, 0.17,  -2.066297},
  //     {'a', 'T', 2,  0, 0, 0.199,  1.439646},
  //     {'a', 'T', 2, -1, 0, 0.67,  -2.114905},
  //     {'a', 'T', 2, -2, 0, 0.89,  -0.757404},
  //     {'a', 'L', 2,  2, 0, 0.095,  1.159966},
  //     {'a', 'L', 2,  1, 0, 0.090,  2.154202},
  //     {'a', 'L', 2,  0, 0, 0.080, -1.789620},
  // }, generation.epsilon);

  return generation;
}

inline std::vector<InitialAmplitude> NearTruthStarts(
    const std::vector<FixedAmplitude>& truth, double epsilon) {
  struct Choice {
    std::size_t index;
    double plusChange;
    double minusChange;
  };

  constexpr double offset = 0.01;
  constexpr double sigma = 0.01;
  std::vector<double> means(truth.size());
  std::vector<Choice> choices;
  double normalizationChange = 0.0;
  for (std::size_t index = 0; index < truth.size(); ++index) {
    const auto& amplitude = truth[index];
    const double weight = amplitude.orientation == 'L'
        ? epsilon * (amplitude.m == 0 ? 1.0 : 2.0)
        : 1.0;
    const double plus = amplitude.magnitude + offset;
    const double plusChange =
        weight * (plus * plus - amplitude.magnitude * amplitude.magnitude);
    means[index] = plus;
    if (amplitude.magnitude < offset) {
      normalizationChange += plusChange;
      continue;
    }
    const double minus = amplitude.magnitude - offset;
    const double minusChange =
        weight * (minus * minus - amplitude.magnitude * amplitude.magnitude);
    choices.push_back({index, plusChange, minusChange});
  }

  // Assign the largest changes first and greedily balance the +/- offsets.
  // This keeps the omitted normalization amplitude close to its truth value.
  std::sort(choices.begin(), choices.end(), [](const Choice& first,
                                                const Choice& second) {
    return std::max(std::abs(first.plusChange), std::abs(first.minusChange)) >
           std::max(std::abs(second.plusChange), std::abs(second.minusChange));
  });
  for (const auto& choice : choices) {
    const bool useMinus =
        std::abs(normalizationChange + choice.minusChange) <
        std::abs(normalizationChange + choice.plusChange);
    normalizationChange += useMinus ? choice.minusChange : choice.plusChange;
    means[choice.index] = truth[choice.index].magnitude +
                          (useMinus ? -offset : offset);
  }

  std::vector<InitialAmplitude> starts;
  starts.reserve(truth.size());
  for (std::size_t index = 0; index < truth.size(); ++index) {
    const auto& amplitude = truth[index];
    starts.push_back({
        amplitude.reflectivity, amplitude.orientation,
        amplitude.l, amplitude.m, amplitude.k, means[index], sigma,
    });
  }
  return starts;
}

inline FitConfig Fit(const std::vector<FixedAmplitude>& truth,
                     double epsilon) {
  FitConfig fit;
  fit.photoproduction = false;
  fit.epsilon = epsilon;
  fit.SetStarts(50000).SetWorkers(10).SetSeed(12345)
      .UseHesse(false).UseNumericalGradients(false);

  // Every independent magnitude is displaced from truth by exactly +/-0.03.
  // Signs are balanced against the normalization sum, and all sigmas are
  // 0.01. Selecting another generation block updates these starts with it.
  fit.initialAmplitudes = NearTruthStarts(truth, epsilon);
  return fit;
}

inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.fixedMoments = FixedMoments();
  settings.model = FixedModel();
  settings.fit = Fit(settings.fixedMoments.fixedAmplitudes,
                     settings.fixedMoments.epsilon);
  settings.bootstrap.fit = settings.fit;
  return settings;
}

} // namespace emi::user
