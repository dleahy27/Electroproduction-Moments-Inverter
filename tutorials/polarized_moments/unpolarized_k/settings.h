#pragma once

#include "emi/Config.h"

namespace emi::user {
inline ModelConfig StudyModel() {
  ModelConfig model;
  model.SetWaves({
      {0, 0}, {1, -1}, {1, 0}, {1, 1},
      {2, -2}, {2, -1}, {2, 0}, {2, 1}, {2, 2},
  }).UseReflectivities(true, true).EnforceLongitudinalParity(true);
  model.normalisationMoment = 2.0;
  return model;
}

inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.model = StudyModel();
  settings.fit.photoproduction = false;
  settings.fit.epsilon = 0.8;
  settings.fit.SetStarts(1000).SetWorkers(8).SetSeed(12345)
      .UseHesse(false);
  settings.fixedMoments.model = StudyModel();
  settings.fixedMoments.model.UseNucleonPolarization(
      NucleonPolarization::Both);
  settings.fixedMoments.photoproduction = false;
  settings.fixedMoments.epsilon = 0.8;
  settings.fixedMoments.mode = FixedGenerationMode::KPositiveDominant;
  settings.fixedMoments.randomMinimum = 0.2;
  settings.fixedMoments.randomMaximum = 1.0;
  settings.fixedMoments.suppression = 0.1;
  return settings;
}
} // namespace emi::user
