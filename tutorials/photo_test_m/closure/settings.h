#pragma once
#include "emi/Config.h"

namespace emi::user {
inline ModelConfig Model() {
  ModelConfig model;
  model.SetWaves({{2, 1}}).UseReflectivities(true, false);
  return model;
}

inline MassModelConfig TruthModel() {
  MassModelConfig model;
  // Coordinates are {initial, lower, upper, fixed}; masses/widths are GeV.
  // A single wave measures the coupling magnitude, but not its overall phase.
  model.resonances = {{"D2",
                       {'a', 'T', 2, 1, 0},
                       {1.320, 1.0, 1.6, true},
                       {0.120, 0.02, 0.4, true},
                       {0.220, 0.0, 2.0, true},
                       {-0.100, -3.14, 3.14, true}}};
  return model;
}

inline MassMomentGenerationConfig MassGeneration() {
  MassMomentGenerationConfig config;
  config.output = "tutorials/photo_test_m/closure/output/moments.root";
  config.tree = "massMoments";
  config.model = Model();
  config.massModel = TruthModel();
  config.photoproduction = true;
  config.relativeUncertainty = 0.0;
  config.absoluteUncertainty = 0.001; // No smearing; errors define fit weights.
  for (int i = 0; i < 41; ++i)
    config.massesGeV.push_back(1.0 + 0.015 * i);
  return config;
}

inline FitMassModelConfig MassFit() {
  FitMassModelConfig config;
  config.fit.input = "tutorials/photo_test_m/closure/output/moments.root";
  config.fit.tree = "massMoments";
  config.fit.output = "tutorials/photo_test_m/closure/output/fit.root";
  config.fit.photoproduction = true;
  config.fit.starts = 8;
  config.fit.workers = 2;
  config.fit.seed = 1441;
  config.fit.hesse = true;
  config.fit.verbose = false;
  // Optional experimental input controls (leave unset for this closure test):
  // config.yieldBranch = "yield";
  // config.normalisation = 2.0; // Scale every moment/error by yield / 2.
  // config.bins = {0, 2, 5};    // Original zero-based TTree entry indices.
  // config.minMassGeV = 1.1;   // Inclusive; combines with bins by intersection.
  // config.maxMassGeV = 1.5;
  config.model = TruthModel();
  auto &resonance = config.model.resonances.front();
  // Deliberately displaced starts check recovery, not just truth evaluation.
  resonance.poleMass.initial = 1.27;
  resonance.poleMass.fixed = false;
  resonance.width.initial = 0.18;
  resonance.width.fixed = false;
  resonance.couplingMagnitude.initial = 0.16;
  resonance.couplingMagnitude.fixed = false;
  resonance.couplingPhase.initial = 0.0;
  resonance.couplingPhase.fixed = false; // The fitter holds this gauge phase.
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
