#pragma once

// P-wave model used for the published electroproduction data. Run-specific
// paths, bins, epsilon values, fit counts, and Hessian/bootstrap choices are
// supplied by the Python runners as command-line overrides.

#include "emi/Config.h"

namespace emi::user {

inline ModelConfig ExperimentalModel() {
  ModelConfig model;
  model.SetWaves({{1, -1}, {1, 0}, {1, 1}})
      .UseReflectivities(true, true)
      .EnforceLongitudinalParity(true)
      .UseNucleonPolarization(NucleonPolarization::None);
  model.normalisationMoment = 2.0;
  return model;
}

inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.model = ExperimentalModel();
  settings.fit.photoproduction = false;
  settings.fit.epsilon = 0.8;
  settings.fit.SetStarts(10000).SetWorkers(10).SetSeed(0)
      .UseHesse(true).UseNumericalGradients(false);
  settings.bootstrap.fit = settings.fit;
  settings.bootstrap.fit.UseHesse(false);
  settings.bootstrap.toys = 1000;
  settings.bootstrap.startsPerToy = 1000;
  return settings;
}

} // namespace emi::user
