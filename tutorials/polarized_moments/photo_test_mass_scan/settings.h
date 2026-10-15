#pragma once

#include "emi/Config.h"

namespace emi::user {
inline ModelConfig ScanModel() {
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
  settings.model = ScanModel();
  settings.fit.photoproduction = true;
  settings.fit.SetStarts(500000).SetWorkers(10).SetSeed(1441)
      .UseHesse(false);
  settings.fixedMoments.model = ScanModel();
  settings.fixedMoments.model.UseNucleonPolarization(
      NucleonPolarization::Both);
  settings.fixedMoments.photoproduction = true;
  settings.fixedMoments.mode = FixedGenerationMode::PhotoTest;
  // settings.fixedMoments.noise = 0.05;
  settings.fixedMoments.massModel.massGeV = 1.3;
  settings.fixedMoments.massModel.backgroundEnabled = false;
  settings.fixedMoments.massModel.kAsymmetry = 0.0;
  return settings;
}
} // namespace emi::user
