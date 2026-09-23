#pragma once

#include "emi/Config.h"

namespace emi::user {
inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.model.SetWaves({
      {0, 0}, {1, -1}, {1, 0}, {1, 1},
      {2, -2}, {2, -1}, {2, 0}, {2, 1}, {2, 2},
  }).UseReflectivities(true, true)
    .EnforceLongitudinalParity(true)
    .UseNucleonPolarization(NucleonPolarization::Both);
  settings.model.normalisationMoment = 2.0;
  settings.fit.photoproduction = false;
  settings.fit.epsilon = 0.8;
  settings.fit.SetStarts(20000).SetWorkers(10).SetSeed(0)
      .UseHesse(false);
  return settings;
}
} // namespace emi::user
