// Reusable resonance line shapes for deterministic mass-dependent generators.
// Complex phases arise from the Breit-Wigner denominator and are retained so
// interference between resonances is represented before moments are formed.
// The line shape is dimensionless; production couplings supply the amplitude
// scale, while pole masses, widths, and the scan mass are all expressed in GeV.
#include "MassModels/MassModel.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

namespace emi::MassModels {

std::complex<double> ConstantWidthBreitWigner(const Resonance& resonance,
                                              double massGeV) {
  // The numerator fixes the value to +i at the pole. Extra magnitude and phase
  // therefore live in the production coupling or explicit complex background.
  if (!std::isfinite(massGeV) || massGeV <= 0.0) {
    throw std::invalid_argument("Invariant mass must be finite and positive");
  }
  if (!std::isfinite(resonance.poleMassGeV) || !std::isfinite(resonance.widthGeV) ||
      resonance.poleMassGeV <= 0.0 || resonance.widthGeV <= 0.0) {
    throw std::invalid_argument("Invalid resonance parameters for " +
                                std::string(resonance.name));
  }

  const double numerator = resonance.poleMassGeV * resonance.widthGeV;
  const std::complex<double> denominator{
      resonance.poleMassGeV * resonance.poleMassGeV - massGeV * massGeV,
      -resonance.poleMassGeV * resonance.widthGeV,
  };
  return numerator / denominator;
}

std::complex<double> TwoBodyPhaseSpace(double massGeV, double firstMassGeV,
                                        double secondMassGeV) {
  if (!std::isfinite(massGeV) || massGeV <= 0.0 ||
      !std::isfinite(firstMassGeV) || firstMassGeV < 0.0 ||
      !std::isfinite(secondMassGeV) || secondMassGeV < 0.0)
    throw std::invalid_argument("Invalid two-body phase-space masses");
  const double sum = firstMassGeV + secondMassGeV;
  const double difference = firstMassGeV - secondMassGeV;
  const double s = massGeV * massGeV;
  // Separate roots select positive real rho above threshold and positive
  // imaginary rho between the threshold and pseudothreshold. Closed channels
  // shift the real denominator; they must not be clamped to zero.
  return std::sqrt(std::complex<double>{1.0 - sum * sum / s, 0.0}) *
         std::sqrt(std::complex<double>{1.0 - difference * difference / s, 0.0});
}

FlatteValue Flatte(double poleMassGeV, double massGeV,
                   const std::vector<FlatteChannel>& channels) {
  if (!std::isfinite(poleMassGeV) || poleMassGeV <= 0.0 ||
      !std::isfinite(massGeV) || massGeV <= 0.0 || channels.size() < 2)
    throw std::invalid_argument("Flatte requires positive masses and at least two channels");
  const std::complex<double> imaginary{0.0, 1.0};
  std::complex<double> denominator{poleMassGeV * poleMassGeV - massGeV * massGeV, 0.0};
  std::vector<std::complex<double>> rho;
  rho.reserve(channels.size());
  for (const auto& channel : channels) {
    if (!std::isfinite(channel.couplingSquaredGeV2) || channel.couplingSquaredGeV2 < 0.0)
      throw std::invalid_argument("Flatte channel couplings squared must be nonnegative");
    rho.push_back(TwoBodyPhaseSpace(massGeV, channel.firstMassGeV,
                                   channel.secondMassGeV));
    denominator -= imaginary * channel.couplingSquaredGeV2 * rho.back();
  }
  FlatteValue result;
  result.value = 1.0 / denominator;
  const auto inverseSquared = result.value * result.value;
  result.poleMassDerivative = -2.0 * poleMassGeV * inverseSquared;
  for (const auto& phaseSpace : rho)
    result.couplingSquaredDerivatives.push_back(imaginary * phaseSpace * inverseSquared);
  return result;
}

} // namespace emi::MassModels
