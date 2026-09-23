#pragma once

#include "emi/Config.h"

namespace emi::user {
inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.model.SetWaves({{0, 0}, {1, -1}, {1, 0}, {1, 1}})
      .UseReflectivities(true, true)
      .UseNucleonPolarization(NucleonPolarization::None);
  settings.model.normalisationMoment = 2.0;
  settings.fit.photoproduction = true;
  settings.fit.SetStarts(10000).SetWorkers(8).SetSeed(0).UseHesse(true);
  settings.bootstrap.fit = settings.fit;
  settings.bootstrap.fit.UseHesse(false);
  settings.bootstrap.toys = 1000;
  settings.bootstrap.startsPerToy = 1000;
  return settings;
}
} // namespace emi::user
