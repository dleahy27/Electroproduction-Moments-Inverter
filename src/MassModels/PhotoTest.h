#pragma once

#include "MassModels/MassModel.h"

#include <array>

namespace emi::MassModels {

// Hard-coded gamma p -> eta pi0 p truth model. It returns raw complex
// amplitudes; generation owns phase conventions, normalisation, and output.
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
