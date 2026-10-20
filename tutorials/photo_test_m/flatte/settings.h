#pragma once
#include "emi/Config.h"
#include <cstdlib>

namespace emi::user {
// PDG 2025 f2(1270) BW averages; f0 Flatte parameters from the BES
// ABLIKIM 05 measurement cited by PDG. See README.md for conventions.
inline ModelConfig Model() {
  ModelConfig model;
  model.SetWaves({{0, 0}, {2, -2}, {2, -1}, {2, 0}, {2, 1}, {2, 2}})
      .UseReflectivities(true, true)
      .UseNucleonPolarization(NucleonPolarization::None);
  return model;
}

inline MassModelConfig TruthModel() {
  MassModelConfig model;
  for (char reflectivity : {'a', 'b'}) {
    MassResonanceConfig scalar;
    scalar.name = "f0_980";
    scalar.wave = {reflectivity, 'T', 0, 0, 0};
    scalar.poleMass.initial = 0.965;
    scalar.poleMass.Unbounded();
    scalar.couplingMagnitude = {reflectivity == 'a' ? 0.035 : 0.025, 0.0, 0.15, false};
    scalar.couplingPhase = {0.0, -3.141592653589793, 3.141592653589793, false};
    // D = m0^2 - s - i (g_pi^2 rho_pi + g_K^2 rho_K).
    // BES's coefficients g1,g2 enter directly here, not their squares.
    // Simple charged-pion / charged-kaon two-channel approximation.
    scalar.UseFlatte({
        {"pipi", 0.13957039, 0.13957039, {0.165, 0.0, 1.0, false}},
        {"KK", 0.493677, 0.493677, {0.165 * 4.21, 0.0, 2.0, false}},
    });
    for (auto& channel : scalar.channels) channel.couplingSquared.Unbounded();
    model.resonances.push_back(scalar);
    for (int m = -2; m <= 2; ++m) {
      MassResonanceConfig tensor;
      tensor.name = "f2_1270";
      tensor.wave = {reflectivity, 'T', 2, m, 0};
      tensor.poleMass.initial = 1.2754;
      tensor.width.initial = 0.1866;
      tensor.poleMass.Unbounded();
      tensor.width.Unbounded();
      tensor.couplingMagnitude = {(reflectivity == 'a' ? 0.22 : 0.16) / (std::abs(m) + 1), 0.0, 0.6, false};
      tensor.couplingPhase = {0.35 + 0.4 * m + (reflectivity == 'b' ? 0.2 : 0.0),
                             -3.141592653589793, 3.141592653589793, false};
      model.resonances.push_back(tensor);
    }
  }
  model.SetPhaseGauge("f0_980", {'a', 'T', 0, 0, 0})
       .SetPhaseGauge("f0_980", {'b', 'T', 0, 0, 0});
  return model;
}

inline UserSettingsConfig Settings() {
  UserSettingsConfig settings;
  settings.model = Model();
  auto& generation = settings.massGeneration;
  generation.model = Model();
  generation.massModel = TruthModel();
  generation.photoproduction = true;
  // Exact deterministic points; these errors only define chi-square weights.
  generation.relativeUncertainty = 0.0;
  generation.absoluteUncertainty = 0.001;
  for (int i = 0; i <= 80; ++i) generation.massesGeV.push_back(0.8 + i * 0.01);
  generation.massesGeV.push_back(2 * 0.493677);

  auto& fit = settings.massFit;
  fit.model = TruthModel();
  fit.fit.tree = "massMoments";
  fit.fit.photoproduction = true;
  fit.fit.starts = 20;
  fit.fit.workers = 10;
  fit.fit.seed = 1441;
  fit.fit.hesse = false;
  fit.curvePoints = 2001;
  for (auto& resonance : fit.model.resonances) {
    if (resonance.name == "f0_980") {
      resonance.poleMass.initial = 0.950;
      resonance.channels[0].couplingSquared.initial = 0.140;
      resonance.channels[1].couplingSquared.initial = 0.580;
    } else {
      resonance.poleMass.initial = 1.290;
      resonance.width.initial = 0.170;
    }
    resonance.couplingMagnitude.initial *= 0.85;
  }
  settings.fit = fit.fit;
  return settings;
}
} // namespace emi::user
