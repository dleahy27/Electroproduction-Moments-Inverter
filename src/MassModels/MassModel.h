#pragma once

// Private line-shape types shared by deterministic truth generators.
// Models return complex partial waves before fit-gauge rotation and moment
// normalization. Public steering types belong in include/emi/Config.h.

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

// Physical-sheet two-body phase space, continued below threshold.
[[nodiscard]] std::complex<double> TwoBodyPhaseSpace(double massGeV,
                                                     double firstMassGeV,
                                                     double secondMassGeV);

struct FlatteChannel {
  double firstMassGeV;
  double secondMassGeV;
  double couplingSquaredGeV2;
};

struct FlatteValue {
  std::complex<double> value;
  std::complex<double> poleMassDerivative;
  std::vector<std::complex<double>> couplingSquaredDerivatives;
};

// Fixed 1 GeV^2 numerator; production magnitude supplies the amplitude scale.
[[nodiscard]] FlatteValue Flatte(double poleMassGeV, double massGeV,
                                const std::vector<FlatteChannel>& channels);

[[nodiscard]] std::complex<double> ConstantWidthBreitWigner(const Resonance& resonance,
                                                            double massGeV);

} // namespace emi::MassModels
