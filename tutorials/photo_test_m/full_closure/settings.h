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
      .UseNucleonPolarization(NucleonPolarization::Both);
  return model;
}

inline MassMomentGenerationConfig MassGeneration() {
  MassMomentGenerationConfig config;
  config.output = "tutorials/photo_test_m/full_closure/output/moments.root";
  config.tree = "massMoments";
  config.photoproduction = true;
  config.model = Model();
  config.massModel = MakePhotoTestMassModel(0.0, true);
  config.relativeUncertainty = 0.0;
  config.absoluteUncertainty = 0.001;
  for (int i = 0; i < 21; ++i)
    config.massesGeV.push_back(0.8 + 0.05 * i);
  return config;
}

inline FitMassModelConfig MassFit() {
  FitMassModelConfig config;
  config.fit.input = "tutorials/photo_test_m/full_closure/output/moments.root";
  config.fit.tree = "massMoments";
  config.fit.output = "tutorials/photo_test_m/full_closure/output/fit.root";
  config.fit.photoproduction = true;
  // One supplied truth start checks the full model. Increase starts for a
  // random-start recovery study; that is substantially more expensive.
  config.fit.starts = 1;
  config.fit.workers = 1;
  config.fit.seed = 1441;
  config.fit.hesse = false;
  config.fit.verbose = false;
  config.maxFunctionCalls = 100000;
  config.model = MakePhotoTestMassModel(0.0, true);
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
