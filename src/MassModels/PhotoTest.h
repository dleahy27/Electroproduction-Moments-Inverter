#pragma once

// Interface to the deterministic photoproduction closure-test model.  Keeping
// its options explicit makes every generated scan point reproducible.

#include "MassModels/MassModel.h"

#include <array>

namespace emi::MassModels {

// Hard-coded gamma p -> eta pi0 p truth model. The k=+1 amplitudes define the
// model and each k=-1 amplitude is a scaled copy of its k=+1 partner.
class PhotoTest final {
public:
  struct Options {
    bool backgroundEnabled = true;
    double kMinusScale = 1.0;
  };

  PhotoTest();
  explicit PhotoTest(Options options);

  [[nodiscard]] AmplitudeSet Evaluate(double massGeV) const;

  [[nodiscard]] static const std::array<Resonance, 4>& Resonances();

private:
  Options options_;
};

} // namespace emi::MassModels
