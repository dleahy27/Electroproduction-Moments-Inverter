#pragma once
#include "emi/Config.h"
#include "emi/PhotoTest.h"

namespace emi::user {
inline ModelConfig Model() {
  ModelConfig model;
  model
      .SetWaves(
          {{0, 0}, {1, -1}, {1, 0}, {1, 1}, {2, -2}, {2, -1}, {2, 0}, {2, 1}, {2, 2}})
      .UseReflectivities(true, true)
      .UseNucleonPolarization(NucleonPolarization::None);
  return model;
}

inline MassMomentGenerationConfig MassGeneration() {
  MassMomentGenerationConfig config;
  config.photoproduction = true;
  config.tree = "massMoments";
  config.model = Model();
  config.model.UseNucleonPolarization(NucleonPolarization::Both);
  config.massModel = MakePhotoTestMassModel(0.0, true);
  // Exact moments, with absolute errors defining the chi-square weights.
  config.relativeUncertainty = 0.0;
  config.absoluteUncertainty = 1;
  for (int i = 0; i < 20; ++i)
    config.massesGeV.push_back(0.8 + i / 19.0);
  return config;
}

inline FitMassModelConfig MassFit() {
  FitMassModelConfig config;
  config.fit.photoproduction = true;
  config.fit.tree = "massMoments";
  config.fit.starts = 50000;
  config.fit.workers = 14;
  config.fit.seed = 1441;
  config.fit.hesse = false;
  config.model = MakePhotoTestMassModel(0.0, false);
  return config;
}

inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.model = Model();
  settings.massGeneration = MassGeneration();
  settings.massFit = MassFit();
  settings.fit = settings.massFit.fit;
  return settings;
}
} // namespace emi::user
