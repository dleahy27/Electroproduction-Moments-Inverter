#pragma once

#include "emi/Config.h"

namespace emi {

// Build the mass-fit representation of the deterministic PhotoTest truth.
// `bothKSectors` retains explicit k=+/-1 amplitudes; false creates the
// single-k effective model used to study an unpolarized inversion.
MassModelConfig MakePhotoTestMassModel(double kAsymmetry, bool bothKSectors = true);

} // namespace emi
