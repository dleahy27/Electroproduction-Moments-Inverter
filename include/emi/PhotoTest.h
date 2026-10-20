#pragma once

#include "emi/Config.h"

namespace emi {

// Build the mass-fit representation of the deterministic PhotoTest truth.
// `bothKSectors` retains explicit k=+/-1 amplitudes; false creates the
// unsplit k=0 model for either mass-moment generation or effective single-k
// fitting. kAsymmetry affects only the explicit-k model.
MassModelConfig MakePhotoTestMassModel(double kAsymmetry, bool bothKSectors = true);

} // namespace emi
