#pragma once

#include <complex>
#include <string_view>
#include <vector>

namespace emi::MassModels {

// Internal value types shared by deterministic mass-dependent generators.
// Keeping this header under src/ makes the model layer private to emi_core.
struct Resonance {
  std::string_view name;
  int spin = 0;
  double poleMassGeV = 0.0;
  double widthGeV = 0.0;
};

struct AmplitudeKey {
  char reflectivity = 'a';
  char orientation = 'T';
  int l = 0;
  int m = 0;
  int k = 1;
};

struct ComplexAmplitude {
  AmplitudeKey key;
  std::complex<double> value{};
};

using AmplitudeSet = std::vector<ComplexAmplitude>;

[[nodiscard]] std::complex<double> ConstantWidthBreitWigner(
    const Resonance& resonance, double massGeV);

} // namespace emi::MassModels
