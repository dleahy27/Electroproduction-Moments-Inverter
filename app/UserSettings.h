#pragma once

// Editable example configuration consumed by `emi --settings`.  Functions are
// inline because ROOT's interpreter includes this header as a small C++ unit.
//
// Treat this file as the human-readable steering card for an analysis.  Wave
// content, fit cost, generation choices, and bootstrap settings are collected
// here so tutorial scripts can share one reproducible physics definition.

#include "emi/Config.h"

namespace emi::user {

// Model used for fits and bootstraps.
inline ModelConfig Model() {
  ModelConfig model;
  model
      .SetWaves({
          {0, 0},
          {1, -1},
          {1, 0},
          {1, 1},
          {2, -2},
          {2, -1},
          {2, 0},
          {2, 1},
          {2, 2},
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
  generation.output = "fixed_recoil_moments.root";
  generation.photoproduction = true; // false selects electroproduction.
  generation.epsilon = 0.8;          // Ignored for photoproduction.
  generation.seed = 44335;
  generation.printAmplitudes = true;
  // generation.noise = 0.05; // Relative Gaussian width; absent keeps truth exact.

  // Generate with the same waves as the fit, but retain both k sectors in
  // the truth sample. Other choices are None, Initial, and Recoil.
  generation.model = Model();
  generation.model.UseNucleonPolarization(NucleonPolarization::Both);

  // PhotoTest evaluates a deterministic S/P/D Breit-Wigner model at one mass.
  generation.mode = FixedGenerationMode::PhotoTest;
  generation.massModel.massGeV = 1.300;
  generation.massModel.backgroundEnabled = false;
  generation.massModel.kAsymmetry = 0.0;
  //
  //   // Used only by the non-custom random modes.
  //   generation.randomMinimum = 0.20;
  //   generation.randomMaximum = 1.00;
  //   generation.suppression = 0.10;

  // Fixed/random lists are unused by deterministic mass-model modes.
  generation.fixedAmplitudes = {};
  generation.randomAmplitudes = {};

  return generation;
}

// Defaults for `emi fit`. Input/output paths and physics-mode command-line
// arguments may still override these values for individual runs.
inline FitConfig Fit() {
  FitConfig fit;
  fit.input = "fixed_recoil_moments.root";
  fit.tree = "genMoments";
  fit.output = "fixed_photomoments_unpolarized_fit.root";
  fit.photoproduction = true;
  fit.epsilon = 0.8; // Ignored for photoproduction.
  fit.SetStarts(10000)
      .SetWorkers(1)
      .SetSeed(12345)
      .UseHesse(true)
      .UseNumericalGradients(false);

  // Optionally give every independent amplitude magnitude (except the one
  // derived from normalization) a Gaussian start distribution. Use k=0 for
  // unpolarized fits. If this list is incomplete, EMI uses its default broad
  // random initialization for the complete start instead.
  // fit.initialAmplitudes = {
  //     {'a', 'T', 0, 0, 0, 0.30, 0.05},
  //     // ...one row for every other independent magnitude...
  // };
  return fit;
}

// A coherent Breit-Wigner-plus-background model for `emi mass-fit`.  The
// selected waves must also appear in Model().  Each entry supplies the start
// value and physical bounds for the fitted pole and complex production term.
inline FitMassModelConfig MassFit() {
  FitMassModelConfig massFit;
  massFit.fit = Fit();
  massFit.fit.input = "mass_model_moments.root";
  massFit.fit.tree = "massMoments";
  massFit.fit.output = "mass_model_fit.root";
  massFit.fit.photoproduction = true;
  massFit.fit.SetStarts(1).SetWorkers(1).UseHesse(true);
  massFit.massBranch = "mass_GeV";
  massFit.model.resonances = {
      {"D2",
       {'a', 'T', 2, 1, 0},
       {1.28, 1.0, 1.6, false},
       {0.14, 0.02, 0.4, false},
       {0.15, -2.0, 2.0, false},
       {0.0, -2.0, 2.0, false}},
  };
  massFit.model.backgrounds = {
      {"D2_nr", {'a', 'T', 2, 1, 0}, {0.02, -1.0, 1.0, false}, {0.0, -1.0, 1.0, false}},
  };
  return massFit;
}

inline MassMomentGenerationConfig MassGeneration() {
  MassMomentGenerationConfig generation;
  generation.output = "mass_model_moments.root";
  generation.tree = "massMoments";
  generation.photoproduction = true;
  generation.epsilon = 1.0;
  generation.model = Model();
  generation.model.UseReflectivities(true, false);
  generation.relativeUncertainty = 0.02;
  for (int bin = 0; bin < 31; ++bin) {
    generation.massesGeV.push_back(1.0 + 0.02 * bin);
  }
  generation.massModel.resonances = {
      {"D2",
       {'a', 'T', 2, 1, 0},
       {1.32, 1.0, 1.6, true},
       {0.12, 0.02, 0.4, true},
       {0.22, -2.0, 2.0, true},
       {-0.10, -2.0, 2.0, true}},
  };
  generation.massModel.backgrounds = {
      {"D2_nr",
       {'a', 'T', 2, 1, 0},
       {0.035, -1.0, 1.0, true},
       {0.015, -1.0, 1.0, true}},
  };
  return generation;
}

inline BootstrapConfig Bootstrap() {
  // A bootstrap refits each fluctuated sample from several random starts. Its
  // uncertainty comes from the ensemble, so a per-toy Hessian is normally
  // unnecessary and comparatively expensive.
  BootstrapConfig bootstrap;
  bootstrap.fit = Fit();
  bootstrap.toys = 1000;
  bootstrap.startsPerToy = 1000;
  bootstrap.fit.UseHesse(true);
  return bootstrap;
}

// EMI loads this aggregate from this header at run time. Editing this file
// therefore does not require rebuilding the executable.
inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.model = Model();
  settings.fit = Fit();
  settings.massFit = MassFit();
  settings.massGeneration = MassGeneration();
  settings.bootstrap = Bootstrap();
  settings.fixedMoments = FixedMoments();
  return settings;
}

} // namespace emi::user
