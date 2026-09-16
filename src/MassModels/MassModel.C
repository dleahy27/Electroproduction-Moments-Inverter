// Reusable resonance line shapes for deterministic mass-dependent generators.
// Complex phases arise from the Breit-Wigner denominator and are retained so
// interference between resonances is represented before moments are formed.
#include "MassModels/MassModel.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

namespace emi::MassModels {

std::complex<double> ConstantWidthBreitWigner(
    const Resonance& resonance, double massGeV) {
  // The numerator fixes the value to +i at the pole. Extra magnitude and phase
  // therefore live in the production coupling or explicit complex background.
  if (!std::isfinite(massGeV) || massGeV <= 0.0) {
    throw std::invalid_argument(
        "Invariant mass must be finite and positive");
  }
  if (!std::isfinite(resonance.poleMassGeV) ||
      !std::isfinite(resonance.widthGeV) ||
      resonance.poleMassGeV <= 0.0 || resonance.widthGeV <= 0.0) {
    throw std::invalid_argument(
        "Invalid resonance parameters for " + std::string(resonance.name));
  }

  const double numerator = resonance.poleMassGeV * resonance.widthGeV;
  const std::complex<double> denominator{
      resonance.poleMassGeV * resonance.poleMassGeV - massGeV * massGeV,
      -resonance.poleMassGeV * resonance.widthGeV,
  };
  return numerator / denominator;
}

} // namespace emi::MassModels
