#pragma once

// Interface to the deterministic photoproduction closure-test model.  Keeping
// its options explicit makes every generated scan point reproducible.

#include "MassModels/MassModel.h"

#include <array>

namespace emi::MassModels {

// Hard-coded gamma p -> eta pi0 p truth model. The shared wave model is split
// between k=+1 and k=-1 by a signed asymmetry, with a fixed relative phase.
class PhotoTest final {
public:
  struct Options {
    bool backgroundEnabled = true;
    double kAsymmetry = 0.0;
  };

  PhotoTest();
  explicit PhotoTest(Options options);

  [[nodiscard]] AmplitudeSet Evaluate(double massGeV) const;

  [[nodiscard]] static const std::array<Resonance, 4>& Resonances();

private:
  Options options_;
};

} // namespace emi::MassModels
