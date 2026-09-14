#pragma once

#include "emi/Config.h"

namespace emi::user {

// Model used for fits and bootstraps.
inline ModelConfig Model() {
  ModelConfig model;
  model
      .SetWaves({
          {0, 0},
          {1, -1},
          {1,  0},
          {1,  1},
          {2, -2},
          {2, -1},
          {2,  0},
          {2,  1},
          {2,  2},
      })
      .UseReflectivities(true, true)
      .EnforceLongitudinalParity(true)
      .UseNucleonPolarization(NucleonPolarization::None);
  model.normalisationMoment = 2.0;
  return model;
}

// Editable configuration used by `emi generate-fixed`.
inline FixedMomentsConfig FixedMoments() {
  FixedMomentsConfig generation;
  generation.output = "InputFiles/Generated/fixed_photomoments.root";
  generation.photoproduction = true; // false selects electroproduction.
  generation.epsilon = 0.8;          // Ignored for photoproduction.
  generation.seed = 12345;
  generation.printAmplitudes = true;

  // Generate with the same waves as the fit, but retain both k sectors in
  // the truth sample. Other choices are None, Initial, and Recoil.
  generation.model = Model();
  generation.model.UseNucleonPolarization(NucleonPolarization::Both);

  // PhotoTest evaluates a deterministic S/P/D Breit-Wigner model at one mass.
  generation.mode = FixedGenerationMode::PhotoTest;
  generation.massModel.massGeV = 1.300;
  generation.massModel.backgroundEnabled = true;
  generation.massModel.kMinusScale = 1.0;

  // Used only by the non-custom random modes.
  generation.randomMinimum = 0.20;
  generation.randomMaximum = 1.00;
  generation.suppression = 0.10;

  // Fixed/random lists are unused by deterministic mass-model modes.
  generation.fixedAmplitudes = {};
  generation.randomAmplitudes = {};

  return generation;
}

// Defaults for `emi fit`. Input/output paths and physics-mode command-line
// arguments may still override these values for individual runs.
inline FitConfig Fit() {
  FitConfig fit;
  fit.input = "InputFiles/Generated/fixed_photomoments.root";
  fit.tree = "genMoments";
  fit.output = "OutputFiles/fixed_photomoments_unpolarized_fit.root";
  fit.photoproduction = true;
  fit.epsilon = 0.8; // Ignored for photoproduction.
  fit
      .SetStarts(10000)
      .SetWorkers(1)
      .SetSeed(12345)
      .UseHesse(false)
      .UseNumericalGradients(false);
  return fit;
}

inline BootstrapConfig Bootstrap() {
  BootstrapConfig bootstrap;
  bootstrap.fit = Fit();
  bootstrap.toys = 1000;
  bootstrap.startsPerToy = 1000;
  bootstrap.fit.UseHesse(false);
  return bootstrap;
}

// EMI loads this aggregate from this header at run time. Editing this file
// therefore does not require rebuilding the executable.
inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.model = Model();
  settings.fit = Fit();
  settings.bootstrap = Bootstrap();
  settings.fixedMoments = FixedMoments();
  return settings;
}

} // namespace emi::user
