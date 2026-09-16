#pragma once

// Public, dependency-light configuration objects.  Fluent Set... methods let
// C++, the CLI settings header, and ROOT's interpreter express the same model
// without exposing implementation classes from src/.

#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace emi {

struct FitConfig {
  // One TTree entry is one kinematic bin; all random starts fit that same row.
  std::filesystem::path input;
  std::string tree = "expMoments";
  int bin = 0;
  std::filesystem::path output;
  double epsilon = 1.0; // Virtual-photon L/T polarization.
  // Photoproduction uses alpha <= 3 and fixes longitudinal amplitudes to zero.
  bool photoproduction = false;

  // A worker count of zero lets ROOT choose. Seed zero delegates seed creation
  // to ROOT; use a non-zero value for reproducible comparisons.
  unsigned starts = 10000;
  unsigned workers = 0;
  std::uint32_t seed = 0;
  bool hesse = true;
  bool useNumericalGradients = false;
  bool verbose = true;

  FitConfig& SetStarts(unsigned value) { starts = value; return *this; }
  FitConfig& SetWorkers(unsigned value) { workers = value; return *this; }
  FitConfig& SetSeed(std::uint32_t value) { seed = value; return *this; }
  FitConfig& UseHesse(bool value = true) { hesse = value; return *this; }
  FitConfig& UseNumericalGradients(bool value = true) {
    useNumericalGradients = value;
    return *this;
  }
};

struct BootstrapConfig {
  // `fit.starts` is irrelevant here: every toy instead receives
  // `startsPerToy` random minimizer starts.
  FitConfig fit;
  unsigned toys = 1000;
  unsigned startsPerToy = 1000;
};

struct Wave {
  // Orbital angular momentum and its projection, with |m| <= l.
  int l;
  int m;
};

enum class NucleonPolarization {
  // Initial and Recoil resolve one nucleon spin index; Both resolves the full
  // pair and therefore removes the single-polarization SO(2) ambiguity.
  None,
  Initial,
  Recoil,
  Both,
};

// Reference used to fix the unobserved k-basis rotation for single polarization.
struct KMixingGauge {
  char reflectivity;
  char orientation;
  int l;
  int m;
};

struct ModelConfig {
  std::vector<Wave> waves = {{1, -1}, {1, 0}, {1, 1}};
  bool usePositiveReflectivity = true;
  bool useNegativeReflectivity = true;
  // When true, negative-m longitudinal waves are derived from parity rather
  // than introduced as independent fit parameters.
  bool enforceLongitudinalParity = true;
  NucleonPolarization nucleonPolarization = NucleonPolarization::None;
  std::optional<KMixingGauge> kMixingGauge;
  double normalisationMoment = 2.0;

  ModelConfig& SetWaves(std::initializer_list<Wave> values) {
    waves.assign(values);
    return *this;
  }
  ModelConfig& SetWaves(std::vector<Wave> values) {
    waves = std::move(values);
    return *this;
  }
  ModelConfig& AddWave(int l, int m) { waves.push_back({l, m}); return *this; }
  ModelConfig& UseReflectivities(bool positive, bool negative) {
    usePositiveReflectivity = positive;
    useNegativeReflectivity = negative;
    return *this;
  }
  ModelConfig& EnforceLongitudinalParity(bool value = true) {
    enforceLongitudinalParity = value;
    return *this;
  }
  ModelConfig& UseNucleonPolarization(NucleonPolarization value) {
    nucleonPolarization = value;
    return *this;
  }
  ModelConfig& SetKMixingGauge(char reflectivity, char orientation,
                               int l, int m) {
    kMixingGauge = KMixingGauge{reflectivity, orientation, l, m};
    return *this;
  }
};

// One complex partial-wave amplitude used to make an exact synthetic point.
// Use k=0 without nucleon polarization and k=+/-1 with polarization.
struct FixedAmplitude {
  char reflectivity;
  char orientation;
  int l;
  int m;
  int k;
  double magnitude;
  double phase;
};

// A magnitude sampled uniformly in [minimum, maximum]. Its phase is sampled
// uniformly in [-pi, pi], unless fixed by a gauge convention.
struct RandomAmplitude {
  char reflectivity;
  char orientation;
  int l;
  int m;
  int k;
  double minimum;
  double maximum;
};

enum class FixedGenerationMode {
  // Custom uses the two explicit lists below. The other modes populate the
  // complete selected model according to their documented rule.
  Custom,
  AllRandom,
  KPositiveDominant,
  KReflectivitySplit,
  PhotoTest,
};

// Shared controls for deterministic mass-dependent generation modes. A mode
// may reject an option that has no meaning for that particular model.
struct MassModelGenerationConfig {
  // No arbitrary default: a mass-dependent mode must be given a mass.
  std::optional<double> massGeV;
  bool backgroundEnabled = true;
  double kMinusScale = 1.0;
};

// Settings used by `emi generate-fixed`. In Custom mode, amplitudes absent
// from both lists are zero and the normalization amplitude must be absent.
struct FixedMomentsConfig {
  // Destination for the generated one-entry ROOT tree.
  std::filesystem::path output = "fixed_test.root";

  // Process and reproducibility settings shared by every generation mode.
  bool photoproduction = false;
  double epsilon = 0.8;
  std::uint32_t seed = 12345;
  bool printAmplitudes = true;

  // The selected wave basis and the rule used to populate its amplitudes.
  ModelConfig model;
  FixedGenerationMode mode = FixedGenerationMode::Custom;

  // Uniform magnitude interval and relative weak-sector scale used by random
  // modes. Custom and PhotoTest generation ignore these values.
  double randomMinimum = 0.2;
  double randomMaximum = 1.0;
  double suppression = 0.1;

  // Used by deterministic mass-dependent modes such as PhotoTest.
  MassModelGenerationConfig massModel;

  std::vector<FixedAmplitude> fixedAmplitudes;
  std::vector<RandomAmplitude> randomAmplitudes;
};

// Everything read from app/UserSettings.h when EMI starts. Command-line
// values, when supplied, override the corresponding fit defaults.
struct UserSettingsConfig {
  ModelConfig model;
  FitConfig fit;
  BootstrapConfig bootstrap;
  FixedMomentsConfig fixedMoments;
};

} // namespace emi
