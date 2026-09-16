// Concrete gamma p -> eta pi0 p toy model used by the mass-scan tutorial.
// Resonance couplings define k=+1; k=-1 is a controlled scaled copy, making
// the effect of omitting a helicity sector directly measurable in closure fits.
#include "MassModels/PhotoTest.h"

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
  ResonanceId resonance;
  int m;
  // k=+1 couplings for positive and negative reflectivity.
  std::array<double, 2> coefficient;
};

constexpr std::array<Sector, 4> kSectors{{
    // Couplings below are indexed only by reflectivity; k=-1 is derived from
    // the matching k=+1 complex amplitude after all contributions are summed.
    {'a', +1},
    {'a', -1},
    {'b', +1},
    {'b', -1},
}};

constexpr std::array<Wave, 9> kWaves{{
    {0, 0},
    {1, -1}, {1, 0}, {1, +1},
    {2, -2}, {2, -1}, {2, 0}, {2, +1}, {2, +2},
}};

constexpr std::array<Resonance, 4> kResonances{{
    {"a0_980", 0, 0.980, 0.075},
    {"pi1_1600", 1, 1.564, 0.492},
    {"a2_1320", 2, 1.306, 0.114},
    {"a2_1700", 2, 1.722, 0.247},
}};

constexpr std::array<Coupling, 14> kCouplings{{
    {ResonanceId::A0_980, 0,
     {+0.32270, +0.08000}},

    {ResonanceId::Pi1_1600, -1,
     {+0.00600, +0.00300}},
    {ResonanceId::Pi1_1600, 0,
     {+0.03033, +0.01200}},
    {ResonanceId::Pi1_1600, +1,
     {-0.03000, -0.00900}},

    {ResonanceId::A2_1320, -2,
     {+0.00300, +0.00150}},
    {ResonanceId::A2_1320, -1,
     {-0.01000, +0.00400}},
    {ResonanceId::A2_1320, 0,
     {+0.05279, -0.01800}},
    {ResonanceId::A2_1320, +1,
     {-0.10900, -0.03000}},
    {ResonanceId::A2_1320, +2,
     {+0.05279, +0.01500}},

    {ResonanceId::A2_1700, -2,
     {-0.00150, +0.00080}},
    {ResonanceId::A2_1700, -1,
     {+0.00400, -0.00200}},
    {ResonanceId::A2_1700, 0,
     {+0.01322, +0.00600}},
    {ResonanceId::A2_1700, +1,
     {-0.03600, +0.01200}},
    {ResonanceId::A2_1700, +2,
     {+0.01322, -0.00500}},
}};

const std::array<Complex, 2> kBackground{{
    {+0.020, +0.020},
    {+0.008, -0.006},
}};

constexpr std::size_t ToIndex(ResonanceId resonance) {
  return static_cast<std::size_t>(resonance);
}

} // namespace

PhotoTest::PhotoTest() = default;

PhotoTest::PhotoTest(Options options) : options_(options) {
  if (!std::isfinite(options_.kMinusScale) ||
      options_.kMinusScale < 0.0) {
    throw std::invalid_argument(
        "PhotoTest k-minus scale must be finite and nonnegative");
  }
}

const std::array<Resonance, 4>& PhotoTest::Resonances() {
  // Expose the immutable table for validation and metadata without duplicating
  // pole masses or widths outside this model.
  return kResonances;
}

AmplitudeSet PhotoTest::Evaluate(double massGeV) const {
  if (!std::isfinite(massGeV) || massGeV <= 0.0) {
    throw std::invalid_argument(
        "PhotoTest invariant mass must be finite and positive");
  }

  std::array<Complex, kResonances.size()> lineShapes;
  // Evaluate each resonance once per mass. Its several m projections reuse
  // the same line shape with different real production couplings.
  for (std::size_t i = 0; i < kResonances.size(); ++i) {
    lineShapes[i] = ConstantWidthBreitWigner(kResonances[i], massGeV);
  }

  AmplitudeSet amplitudes;
  amplitudes.reserve(kWaves.size() * kSectors.size());
  for (const auto& wave : kWaves) {
    for (std::size_t sectorIndex = 0;
         sectorIndex < kSectors.size(); ++sectorIndex) {
      const std::size_t reflectivityIndex =
          kSectors[sectorIndex].reflectivity == 'a' ? 0 : 1;
      Complex value{0.0, 0.0};
      for (const auto& coupling : kCouplings) {
        const auto resonanceIndex = ToIndex(coupling.resonance);
        if (kResonances[resonanceIndex].spin == wave.l &&
            coupling.m == wave.m) {
          value += coupling.coefficient[reflectivityIndex] *
                   lineShapes[resonanceIndex];
        }
      }

      if (options_.backgroundEnabled && wave.l == 0 && wave.m == 0) {
        value += kBackground[reflectivityIndex];
      }
      if (kSectors[sectorIndex].k == -1) {
        // Scale after resonances and background are combined, making the whole
        // negative-k complex amplitude an exact partner of positive k.
        value *= options_.kMinusScale;
      }

      amplitudes.push_back({
          {kSectors[sectorIndex].reflectivity, 'T',
           wave.l, wave.m, kSectors[sectorIndex].k},
          value,
      });
    }
  }
  return amplitudes;
}

} // namespace emi::MassModels
