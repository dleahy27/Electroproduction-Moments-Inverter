#pragma once

#include "emi/Config.h"

namespace emi::user {

inline ModelConfig Model() {
  ModelConfig model;
  model
      .SetWaves({
          {1, -1},
          {1,  0},
          {1,  1},
      })
      .UseReflectivities(true, true)
      .EnforceLongitudinalParity(true);
  model.normalisationMoment = 2.0;
  return model;
}

inline ModelConfig GenerationModel() {
  ModelConfig model;
  model
      .SetWaves({
          {0,  0},
          {1, -1}, {1, 0}, {1, 1},
          {2, -2}, {2, -1}, {2, 0}, {2, 1}, {2, 2},
      })
      .UseReflectivities(true, true)
      .EnforceLongitudinalParity(true);
  return model;
}

inline void ConfigureFit(FitConfig& fit) {
  fit
      .SetStarts(10000)
      .SetWorkers(0)
      .SetSeed(0)
      .UseHesse(true)
      .UseNumericalGradients(false);
}

inline void ConfigureBootstrap(BootstrapConfig& bootstrap) {
  bootstrap.toys = 1000;
  bootstrap.startsPerToy = 1000;
  bootstrap.fit.UseHesse(false);
}

} // namespace emi::user
