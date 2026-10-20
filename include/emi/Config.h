#pragma once

// Public, dependency-light configuration objects.  Fluent Set... methods let
// C++, the CLI settings header, and ROOT's interpreter express the same model
// without exposing implementation classes from src/.
//
// All structures are value types: constructing or copying one has no I/O side
// effects.  Paths, tree names, bin choices, and numerical controls are checked
// only when a Runner operation begins, where useful context can be included in
// any error message.

#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace emi {

// Gaussian used to initialize one independent amplitude magnitude. Use k=0
// without nucleon polarization and k=+/-1 with polarization. The amplitude
// eliminated by normalization must not be listed.
struct InitialAmplitude {
  char reflectivity;
  char orientation;
  int l;
  int m;
  int k;
  double mean;
  double standardDeviation;
};

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

  // These distributions are used only when every independent magnitude other
  // than the normalization amplitude is present. An empty or partial list
  // retains the default broad random initialization for all magnitudes.
  std::vector<InitialAmplitude> initialAmplitudes;

  FitConfig& SetStarts(unsigned value) {
    starts = value;
    return *this;
  }
  FitConfig& SetWorkers(unsigned value) {
    workers = value;
    return *this;
  }
  FitConfig& SetSeed(std::uint32_t value) {
    seed = value;
    return *this;
  }
  FitConfig& UseHesse(bool value = true) {
    hesse = value;
    return *this;
  }
  FitConfig& UseNumericalGradients(bool value = true) {
    useNumericalGradients = value;
    return *this;
  }
};

// A real minimizer coordinate used by the mass-dependent model.  Bounds are
// optional. Unbounded positive shape parameters use free log coordinates,
// retaining physical positivity without MINUIT limits.
struct MassParameter {
  double initial = 0.0;
  double lower = -10.0;
  double upper = 10.0;
  bool fixed = false;
  bool bounded = true;

  // Keep the supplied initial value and fixed flag; ignore lower/upper limits.
  MassParameter& Unbounded() {
    bounded = false;
    return *this;
  }
};

// Identifies one complex partial wave in a mass model.  The key has the same
// convention as FixedAmplitude; use k=0 for an unpolarized analysis.
struct MassWave {
  char reflectivity = 'a';
  char orientation = 'T';
  int l = 0;
  int m = 0;
  int k = 0;
};

enum class MassLineShape {
  ConstantWidthBreitWigner,
  Flatte,
};

// One S-wave two-body Flatte channel. Daughter masses are fixed in GeV;
// couplingSquared is the coefficient g_i^2 in GeV^2 multiplying rho_i=2q_i/m.
// Channel names identify shared fitted coordinates across projections/sectors.
struct MassDecayChannel {
  std::string name;
  double firstMassGeV = 0.0;
  double secondMassGeV = 0.0;
  MassParameter couplingSquared{0.1, 1.e-6, 4.0, false};
};

// A selectable line shape times a complex production coupling. The default
// retains the constant-width relativistic BW used by existing steering cards.
struct MassResonanceConfig {
  std::string name;
  MassWave wave;
  MassParameter poleMass{1.3, 0.1, 4.0, false};
  MassParameter width{0.15, 1.e-4, 2.0, false};
  // Same polar convention as ordinary EMI amplitudes.
  MassParameter couplingMagnitude{1.0, 0.0, 10.0, false};
  MassParameter couplingPhase{0.0, -3.141592653589793, 3.141592653589793, false};
  MassLineShape lineShape = MassLineShape::ConstantWidthBreitWigner;
  // Flatte uses poleMass and these channels; width is unused and not fitted.
  // F(m) = (1 GeV^2)/(m0^2-m^2-i sum_i g_i^2 rho_i(m)).
  std::vector<MassDecayChannel> channels;
  // Empty uses the angular-momentum sector l. Set a common explicit group for
  // resonances sharing further quantum numbers absent from this wave schema.
  // A shared resonance name must use the same l and group in every amplitude.
  std::string massOrderingGroup;

  MassResonanceConfig& UseFlatte(std::initializer_list<MassDecayChannel> values) {
    lineShape = MassLineShape::Flatte;
    channels.assign(values);
    return *this;
  }
};

// A coherent complex constant term. It is useful for non-resonant production
// and is evaluated in precisely the same phase convention as the resonances.
struct MassBackgroundConfig {
  std::string name;
  MassWave wave;
  MassParameter real{0.0, -10.0, 10.0, false};
  MassParameter imag{0.0, -10.0, 10.0, false};
};

// Phase-reference coupling held at its configured initial value. Single-k
// fits default to the highest-(l,m) transverse resonance in each populated
// reflectivity sector; explicit-k fits use one overall reference (a by default).
struct MassPhaseGauge {
  std::string resonance;
  MassWave wave;
};

struct MassModelConfig {
  std::vector<MassResonanceConfig> resonances;
  std::vector<MassBackgroundConfig> backgrounds;
  // Legacy single-reference override; SetPhaseGauge overrides it for the same
  // reflectivity. Retained for existing steering cards.
  std::optional<MassPhaseGauge> phaseGauge;
  std::vector<MassPhaseGauge> phaseGauges;
  // Strict minimum adjacent mass separation in each ordering group (GeV).
  // This removes equality, without assigning a physical resonance separation.
  double minimumMassGapGeV = 1.e-6;

  // Select a reference resonance independently in a and b for k=0. The phase
  // value comes from that resonance's couplingPhase.initial. Explicit-k fits
  // have only one overall phase freedom and accept only one override.
  MassModelConfig& SetPhaseGauge(std::string resonance, MassWave wave) {
    if (wave.reflectivity != 'a' && wave.reflectivity != 'b')
      throw std::invalid_argument("Mass phase gauge reflectivity must be a or b");
    for (auto& gauge : phaseGauges) {
      if (gauge.wave.reflectivity == wave.reflectivity) {
        gauge = {std::move(resonance), wave};
        return *this;
      }
    }
    phaseGauges.push_back({std::move(resonance), wave});
    return *this;
  }
};

// Fits selected entries of one moment tree together. `fit.bin` is deliberately
// not used: every entry is a mass point and `massBranch` supplies its invariant
// mass in GeV. Empty selections use all entries. The usual FitConfig controls
// Minuit starts and diagnostics.
struct FitMassModelConfig {
  FitConfig fit;
  std::string massBranch = "mass_GeV";
  // Optional scalar branch containing the unnormalised zeroth moment (yield).
  // When set, multiply every measured moment and error by yield / normalisation.
  // The yield is treated as exact; its uncertainty/correlations are not propagated.
  std::optional<std::string> yieldBranch;
  double normalisation = 2.0;
  // Zero-based original TTree entry indices, independent of mass ordering.
  // Entries are fitted once in tree order; duplicates do not repeat observations.
  std::vector<std::int64_t> bins;
  // Inclusive mass bounds in GeV; either bound may be omitted. When combined
  // with `bins`, an entry must satisfy both the index and mass selections.
  std::optional<double> minMassGeV;
  std::optional<double> maxMassGeV;
  // A finite budget makes high-dimensional multistart scans predictable. A
  // start that reaches this limit is retained with its Minuit status.
  unsigned maxFunctionCalls = 100000;
  // Dense best-model curves, evaluated from resonance parameters, not by
  // interpolating the fitted data points. Zero disables this additional tree.
  unsigned curvePoints = 1001;
  MassModelConfig model;
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
  ModelConfig& AddWave(int l, int m) {
    waves.push_back({l, m});
    return *this;
  }
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
  ModelConfig& SetKMixingGauge(char reflectivity, char orientation, int l, int m) {
    kMixingGauge = KMixingGauge{reflectivity, orientation, l, m};
    return *this;
  }
};

// Closure-test generator for a complete mass scan. It writes scalar moment
// branches per entry, ready to be consumed directly by FitMassModel.
struct MassMomentGenerationConfig {
  std::filesystem::path output = "mass_model_moments.root";
  std::string tree = "massMoments";
  std::vector<double> massesGeV;
  double epsilon = 1.0;
  bool photoproduction = true;
  ModelConfig model;
  MassModelConfig massModel;
  std::string massBranch = "mass_GeV";
  double relativeUncertainty = 0.02;
  double absoluteUncertainty = 1.e-5;
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
  // Signed k-sector asymmetry for PhotoTest. s=+1 selects k=+1, s=-1
  // selects k=-1, and s=0 gives both sectors equal amplitude factors.
  double kAsymmetry = 0.0;
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

  // When present, smear every non-normalization generated moment with a
  // Gaussian of relative width `noise`. An absent value keeps exact moments.
  std::optional<double> noise;

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
  FitMassModelConfig massFit;
  MassMomentGenerationConfig massGeneration;
  BootstrapConfig bootstrap;
  FixedMomentsConfig fixedMoments;
};

} // namespace emi
