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
  // epsilon=+,k=+; epsilon=+,k=-; epsilon=-,k=+; epsilon=-,k=-
  std::array<double, 4> coefficient;
};

constexpr std::array<Sector, 4> kSectors{{
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
     {+0.32270, +0.12000, +0.08000, -0.04500}},

    {ResonanceId::Pi1_1600, -1,
     {+0.00600, -0.00400, +0.00300, +0.00200}},
    {ResonanceId::Pi1_1600, 0,
     {+0.03033, -0.01800, +0.01200, +0.00700}},
    {ResonanceId::Pi1_1600, +1,
     {-0.03000, -0.01200, -0.00900, +0.00600}},

    {ResonanceId::A2_1320, -2,
     {+0.00300, -0.00200, +0.00150, +0.00100}},
    {ResonanceId::A2_1320, -1,
     {-0.01000, +0.00700, +0.00400, -0.00300}},
    {ResonanceId::A2_1320, 0,
     {+0.05279, +0.02500, -0.01800, +0.01200}},
    {ResonanceId::A2_1320, +1,
     {-0.10900, +0.04500, -0.03000, -0.02000}},
    {ResonanceId::A2_1320, +2,
     {+0.05279, -0.02000, +0.01500, +0.00900}},

    {ResonanceId::A2_1700, -2,
     {-0.00150, +0.00100, +0.00080, -0.00050}},
    {ResonanceId::A2_1700, -1,
     {+0.00400, +0.00250, -0.00200, +0.00150}},
    {ResonanceId::A2_1700, 0,
     {+0.01322, -0.00900, +0.00600, +0.00400}},
    {ResonanceId::A2_1700, +1,
     {-0.03600, -0.01800, +0.01200, -0.00800}},
    {ResonanceId::A2_1700, +2,
     {+0.01322, +0.00700, -0.00500, +0.00300}},
}};

const std::array<Complex, 4> kBackground{{
    {+0.020, +0.020},
    {-0.012, +0.016},
    {+0.008, -0.006},
    {-0.005, -0.008},
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
  return kResonances;
}

AmplitudeSet PhotoTest::Evaluate(double massGeV) const {
  if (!std::isfinite(massGeV) || massGeV <= 0.0) {
    throw std::invalid_argument(
        "PhotoTest invariant mass must be finite and positive");
  }

  std::array<Complex, kResonances.size()> lineShapes;
  for (std::size_t i = 0; i < kResonances.size(); ++i) {
    lineShapes[i] = ConstantWidthBreitWigner(kResonances[i], massGeV);
  }

  AmplitudeSet amplitudes;
  amplitudes.reserve(kWaves.size() * kSectors.size());
  for (const auto& wave : kWaves) {
    for (std::size_t sectorIndex = 0;
         sectorIndex < kSectors.size(); ++sectorIndex) {
      Complex value{0.0, 0.0};
      for (const auto& coupling : kCouplings) {
        const auto resonanceIndex = ToIndex(coupling.resonance);
        if (kResonances[resonanceIndex].spin == wave.l &&
            coupling.m == wave.m) {
          value += coupling.coefficient[sectorIndex] *
                   lineShapes[resonanceIndex];
        }
      }

      if (options_.backgroundEnabled && wave.l == 0 && wave.m == 0) {
        value += kBackground[sectorIndex];
      }
      if (kSectors[sectorIndex].k == -1) {
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
