// Concrete gamma p -> eta pi0 p toy model used by the mass-scan tutorial.
// A signed asymmetry splits each wave between k=+1 and k=-1, making the effect
// of omitting a helicity sector directly measurable in closure fits.
// This is a pedagogical toy model, not a data-driven resonance fit.  Its fixed
// couplings create S-, P-, and D-wave interference with enough structure to
// test mass-dependent inversion and polarization assumptions.
#include "MassModels/PhotoTest.h"
#include "emi/PhotoTest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>

namespace emi::MassModels {
namespace {

using Complex = std::complex<double>;

enum class ResonanceId : std::size_t {
  A0_980,
  Pi1_1600,
  A2_1320,
  A2_1700,
};

struct Sector {
  char reflectivity;
  int k;
};

struct Wave {
  int l;
  int m;
};

struct Coupling {
  // A real production coefficient multiplies a complex line shape.  Coherent
  // addition of contributions creates mass-dependent phases and interference.
  ResonanceId resonance;
  int m;
  // k=+1 couplings for positive and negative reflectivity.
  std::array<double, 2> coefficient;
};

constexpr std::array<Sector, 4> kSectors{{
    {'a', +1},
    {'a', -1},
    {'b', +1},
    {'b', -1},
}};

constexpr std::array<Wave, 9> kWaves{{
    {0, 0},
    {1, -1},
    {1, 0},
    {1, +1},
    {2, -2},
    {2, -1},
    {2, 0},
    {2, +1},
    {2, +2},
}};

constexpr std::array<Resonance, 4> kResonances{{
    {"a0_980", 0, 0.980, 0.075},
    {"pi1_1600", 1, 1.564, 0.492},
    {"a2_1320", 2, 1.306, 0.114},
    {"a2_1700", 2, 1.722, 0.247},
}};

constexpr std::array<Coupling, 14> kCouplings{{
    {ResonanceId::A0_980, 0, {+0.32270, +0.08000}},

    {ResonanceId::Pi1_1600, -1, {+0.00600, +0.00300}},
    {ResonanceId::Pi1_1600, 0, {+0.03033, +0.01200}},
    {ResonanceId::Pi1_1600, +1, {-0.03000, -0.00900}},

    {ResonanceId::A2_1320, -2, {+0.00300, +0.00150}},
    {ResonanceId::A2_1320, -1, {-0.01000, +0.00400}},
    {ResonanceId::A2_1320, 0, {+0.05279, -0.01800}},
    {ResonanceId::A2_1320, +1, {-0.10900, -0.03000}},
    {ResonanceId::A2_1320, +2, {+0.05279, +0.01500}},

    {ResonanceId::A2_1700, -2, {-0.00150, +0.00080}},
    {ResonanceId::A2_1700, -1, {+0.00400, -0.00200}},
    {ResonanceId::A2_1700, 0, {+0.01322, +0.00600}},
    {ResonanceId::A2_1700, +1, {-0.03600, +0.01200}},
    {ResonanceId::A2_1700, +2, {+0.01322, -0.00500}},
}};

const std::array<Complex, 2> kBackground{{
    {+0.020, +0.020},
    {+0.008, -0.006},
}};

// Fixed k=-1 phase offsets, in the order of kWaves. They are deliberately
// independent of reflectivity, mass, and k-sector asymmetry so each partial
// wave keeps the same relative k-sector phase throughout a scan.
constexpr std::array<double, kWaves.size()> kMinusPhaseOffsets{{
    +0.73, // S(0, 0)
    -0.41, // P(1,-1)
    +1.19, // P(1, 0)
    -2.06, // P(1,+1)
    +0.37, // D(2,-2)
    -1.46, // D(2,-1)
    +2.28, // D(2, 0)
    -0.92, // D(2,+1)
    +1.64, // D(2,+2)
}};

constexpr std::size_t ToIndex(ResonanceId resonance) {
  return static_cast<std::size_t>(resonance);
}

} // namespace

PhotoTest::PhotoTest() = default;

PhotoTest::PhotoTest(Options options) : options_(options) {
  if (!std::isfinite(options_.kAsymmetry) || options_.kAsymmetry < -1.0 ||
      options_.kAsymmetry > 1.0) {
    throw std::invalid_argument(
        "PhotoTest k asymmetry must be finite and between -1 and 1");
  }
}

const std::array<Resonance, 4>& PhotoTest::Resonances() {
  // Expose the immutable table for validation and metadata without duplicating
  // pole masses or widths outside this model.
  return kResonances;
}

AmplitudeSet PhotoTest::Evaluate(double massGeV) const {
  if (!std::isfinite(massGeV) || massGeV <= 0.0) {
    throw std::invalid_argument("PhotoTest invariant mass must be finite and positive");
  }

  std::array<Complex, kResonances.size()> lineShapes;
  // Evaluate each resonance once per mass. Its several m projections reuse
  // the same line shape with different real production couplings.
  for (std::size_t i = 0; i < kResonances.size(); ++i) {
    lineShapes[i] = ConstantWidthBreitWigner(kResonances[i], massGeV);
  }

  // Enumerate the complete rectangular basis so waves without a coupling are
  // represented explicitly by zero. Stable ordering simplifies branch creation
  // and point-by-point comparisons across a mass scan.
  AmplitudeSet amplitudes;
  amplitudes.reserve(kWaves.size() * kSectors.size());
  for (std::size_t waveIndex = 0; waveIndex < kWaves.size(); ++waveIndex) {
    const auto& wave = kWaves[waveIndex];
    for (std::size_t sectorIndex = 0; sectorIndex < kSectors.size(); ++sectorIndex) {
      const std::size_t reflectivityIndex =
          kSectors[sectorIndex].reflectivity == 'a' ? 0 : 1;
      Complex value{0.0, 0.0};
      for (const auto& coupling : kCouplings) {
        const auto resonanceIndex = ToIndex(coupling.resonance);
        if (kResonances[resonanceIndex].spin == wave.l && coupling.m == wave.m) {
          value += coupling.coefficient[reflectivityIndex] * lineShapes[resonanceIndex];
        }
      }

      if (options_.backgroundEnabled && wave.l == 0 && wave.m == 0) {
        value += kBackground[reflectivityIndex];
      }
      const bool negativeK = kSectors[sectorIndex].k == -1;
      const double sectorFactor = negativeK ? (1.0 - options_.kAsymmetry) / 2.0
                                            : (1.0 + options_.kAsymmetry) / 2.0;
      // Apply the split after resonances and background are combined so each
      // k sector has the same wave shape. The negative sector receives its
      // fixed partial-wave phase shift.
      value *= std::sqrt(sectorFactor);
      if (negativeK) {
        value *= std::polar(1.0, kMinusPhaseOffsets[waveIndex]);
      }

      amplitudes.push_back({
          {kSectors[sectorIndex].reflectivity, 'T', wave.l, wave.m,
           kSectors[sectorIndex].k},
          value,
      });
    }
  }
  return amplitudes;
}

} // namespace emi::MassModels

namespace emi::MassModels {

// Convert the truth table to fitted BW coordinates. Pole coordinates with the
// same resonance name are shared across all wave projections and k sectors.
// Negative real coefficients become positive magnitudes with a pi phase.
// The effective k=0 model uses unsplit couplings as an initial guess; those
// coordinates are not predicted to equal either fitted two-k sector.
MassModelConfig MakePhotoTestMassModelImpl(double kAsymmetry, bool bothKSectors) {
  if (!std::isfinite(kAsymmetry) || kAsymmetry < -1.0 || kAsymmetry > 1.0) {
    throw std::invalid_argument("PhotoTest k asymmetry must be in [-1, 1]");
  }
  MassModelConfig model;
  const int firstK = bothKSectors ? -1 : 0;
  const int lastK = bothKSectors ? 1 : 0;
  for (const auto& coupling : kCouplings) {
    const auto& resonance = kResonances[ToIndex(coupling.resonance)];
    for (int reflectivityIndex = 0; reflectivityIndex < 2; ++reflectivityIndex) {
      for (int k = firstK; k <= lastK; k += bothKSectors ? 2 : 1) {
        const bool minus = k == -1;
        const double factor = bothKSectors ? std::sqrt(minus ? (1.0 - kAsymmetry) / 2.0
                                                             : (1.0 + kAsymmetry) / 2.0)
                                           : 1.0;
        const auto waveIt =
            std::find_if(kWaves.begin(), kWaves.end(), [&](const Wave& wave) {
              return wave.l == resonance.spin && wave.m == coupling.m;
            });
        const auto waveIndex = static_cast<size_t>(waveIt - kWaves.begin());
        const double phase = minus ? kMinusPhaseOffsets[waveIndex] : 0.0;
        MassResonanceConfig entry;
        entry.name = resonance.name;
        entry.wave = {reflectivityIndex == 0 ? 'a' : 'b', 'T', resonance.spin,
                      coupling.m, k};
        entry.poleMass = {resonance.poleMassGeV, 0.5, 2.5, false};
        entry.width = {resonance.widthGeV, 1.e-3, 1.0, false};
        entry.couplingMagnitude = {std::abs(coupling.coefficient[reflectivityIndex]) *
                                       factor,
                                   0.0, 2.0, false};
        const double couplingPhase =
            std::atan2(std::sin(phase + (coupling.coefficient[reflectivityIndex] < 0.0
                                             ? 3.141592653589793
                                             : 0.0)),
                       std::cos(phase + (coupling.coefficient[reflectivityIndex] < 0.0
                                             ? 3.141592653589793
                                             : 0.0)));
        entry.couplingPhase = {couplingPhase, -3.141592653589793, 3.141592653589793,
                               false};
        model.resonances.push_back(std::move(entry));
      }
    }
  }
  return model;
}

} // namespace emi::MassModels

namespace emi {
MassModelConfig MakePhotoTestMassModel(double kAsymmetry, bool bothKSectors) {
  return MassModels::MakePhotoTestMassModelImpl(kAsymmetry, bothKSectors);
}
} // namespace emi
