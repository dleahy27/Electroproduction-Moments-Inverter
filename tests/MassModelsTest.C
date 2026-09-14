#include "emi/Runner.h"

#include "MassModels/MassModel.h"
#include "MassModels/PhotoTest.h"

#include "TFile.h"
#include "TLeaf.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TTree.h"

#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using emi::MassModels::ComplexAmplitude;
using emi::MassModels::PhotoTest;

bool Near(double actual, double expected, double tolerance = 1e-12) {
  return std::abs(actual - expected) <= tolerance;
}

bool Near(std::complex<double> actual, std::complex<double> expected,
          double tolerance = 1e-12) {
  return std::abs(actual - expected) <= tolerance;
}

void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

const ComplexAmplitude& Find(const emi::MassModels::AmplitudeSet& amplitudes,
                             char reflectivity, int l, int m, int k) {
  for (const auto& amplitude : amplitudes) {
    const auto& key = amplitude.key;
    if (key.reflectivity == reflectivity && key.orientation == 'T' &&
        key.l == l && key.m == m && key.k == k) {
      return amplitude;
    }
  }
  throw std::runtime_error("Could not find requested PhotoTest amplitude");
}

std::string Index(int value) {
  return value < 0 ? "m" + std::to_string(-value) : std::to_string(value);
}

std::string AmplitudeName(const ComplexAmplitude& amplitude, bool phase) {
  const auto& key = amplitude.key;
  std::string name(1, key.reflectivity);
  if (phase) name += "phi";
  name += "_T_" + std::to_string(key.l) + '_' + Index(key.m) + '_' +
          Index(key.k);
  return name;
}

double LeafValue(TTree& tree, const std::string& name) {
  TLeaf* leaf = tree.GetLeaf(name.c_str());
  if (!leaf) throw std::runtime_error("Missing ROOT leaf " + name);
  return leaf->GetValue();
}

void CheckMassModel() {
  const auto& resonances = PhotoTest::Resonances();
  for (const auto& resonance : resonances) {
    const auto value = emi::MassModels::ConstantWidthBreitWigner(
        resonance, resonance.poleMassGeV);
    Require(Near(value, {0.0, 1.0}),
            "Breit-Wigner is not +i at the pole for " +
                std::string(resonance.name));
  }

  constexpr double mass = 1.410;
  const PhotoTest withoutBackground({false, 1.0});
  const auto amplitudes = withoutBackground.Evaluate(mass);
  Require(amplitudes.size() == 36,
          "PhotoTest did not return 36 amplitudes");

  const auto a0 = emi::MassModels::ConstantWidthBreitWigner(
      resonances[0], mass);
  Require(Near(Find(amplitudes, 'a', 0, 0, +1).value, 0.32270 * a0),
          "S wave does not match its resonance expression");

  const auto a2_1320 = emi::MassModels::ConstantWidthBreitWigner(
      resonances[2], mass);
  const auto a2_1700 = emi::MassModels::ConstantWidthBreitWigner(
      resonances[3], mass);
  const auto expectedD = -0.01000 * a2_1320 + 0.00400 * a2_1700;
  Require(Near(Find(amplitudes, 'a', 2, -1, +1).value, expectedD),
          "D wave does not contain the coherent two-resonance sum");

  const auto withBackground = PhotoTest({true, 1.0}).Evaluate(mass);
  for (const auto& amplitude : amplitudes) {
    const auto& key = amplitude.key;
    const auto difference =
        Find(withBackground, key.reflectivity, key.l, key.m, key.k).value -
        amplitude.value;
    if (key.l != 0) {
      Require(Near(difference, {0.0, 0.0}),
              "Background changed a non-S-wave amplitude");
    }
  }
  struct ExpectedBackground {
    char reflectivity;
    int k;
    std::complex<double> value;
  };
  const std::array<ExpectedBackground, 4> expectedBackgrounds{{
      {'a', +1, {+0.020, +0.020}},
      {'a', -1, {-0.012, +0.016}},
      {'b', +1, {+0.008, -0.006}},
      {'b', -1, {-0.005, -0.008}},
  }};
  for (const auto& background : expectedBackgrounds) {
    Require(Near(Find(withBackground, background.reflectivity, 0, 0,
                      background.k).value -
                     Find(amplitudes, background.reflectivity, 0, 0,
                          background.k).value,
                 background.value),
            "PhotoTest S-wave background is incorrect");
  }

  const auto noKMinus = PhotoTest({true, 0.0}).Evaluate(mass);
  for (const auto& amplitude : noKMinus) {
    if (amplitude.key.k == -1) {
      Require(amplitude.value == std::complex<double>{0.0, 0.0},
              "k-minus scale zero did not remove an amplitude");
    }
  }

  const auto first = withoutBackground.Evaluate(1.200);
  const auto atDifferentMass = withoutBackground.Evaluate(1.700);
  Require(!Near(Find(first, 'a', 0, 0, +1).value,
                Find(atDifferentMass, 'a', 0, 0, +1).value),
          "Changing mass did not reevaluate PhotoTest");
  const auto repeated = withoutBackground.Evaluate(1.200);
  Require(first.size() == repeated.size(),
          "Repeated PhotoTest evaluation changed set size");
  for (std::size_t i = 0; i < first.size(); ++i) {
    Require(Near(first[i].value, repeated[i].value),
            "PhotoTest retained stale mass-dependent state");
  }
}

emi::ModelConfig FullPhotoModel() {
  emi::ModelConfig model;
  model.SetWaves({
      {0, 0},
      {1, -1}, {1, 0}, {1, +1},
      {2, -2}, {2, -1}, {2, 0}, {2, +1}, {2, +2},
  });
  model.UseReflectivities(true, true);
  model.UseNucleonPolarization(emi::NucleonPolarization::Both);
  model.normalisationMoment = 2.0;
  return model;
}

void CheckGeneration(const std::filesystem::path& output) {
  emi::FixedMomentsConfig generation;
  generation.output = output;
  generation.photoproduction = true;
  generation.printAmplitudes = false;
  generation.model = FullPhotoModel();
  generation.mode = emi::FixedGenerationMode::PhotoTest;
  generation.massModel.massGeV = 1.306;
  generation.massModel.backgroundEnabled = true;
  generation.massModel.kMinusScale = 0.5;
  emi::GenerateFixedMoments(generation);

  std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "READ"));
  Require(file && !file->IsZombie(), "Could not read generated ROOT file");
  auto* tree = dynamic_cast<TTree*>(file->Get("genMoments"));
  Require(tree && tree->GetEntries() == 1,
          "Generated ROOT file has an invalid genMoments tree");
  tree->GetEntry(0);

  const auto rawAmplitudes = PhotoTest({true, 0.5}).Evaluate(1.306);
  for (const auto& amplitude : rawAmplitudes) {
    Require(tree->GetLeaf(AmplitudeName(amplitude, false).c_str()),
            "Generated file is missing a truth magnitude");
    Require(tree->GetLeaf(AmplitudeName(amplitude, true).c_str()),
            "Generated file is missing a truth phase");
  }
  Require(Near(LeafValue(*tree, "aphi_T_2_2_1"), 0.0),
          "Generated reference phase is not zero");
  Require(Near(LeafValue(*tree, "RH_0_0_0_0_0"), 2.0, 1e-10),
          "Generated zeroth moment is not normalized to two");

  auto* modelName = dynamic_cast<TObjString*>(file->Get("mass_model"));
  Require(modelName && std::string(modelName->GetString().Data()) == "PhotoTest",
          "Generated file is missing its mass-model name");
  auto* mass = dynamic_cast<TParameter<double>*>(
      file->Get("invariant_mass_GeV"));
  auto* scale = dynamic_cast<TParameter<double>*>(
      file->Get("amplitude_normalisation_scale"));
  auto* rawH000 = dynamic_cast<TParameter<double>*>(
      file->Get("raw_H_0_0_0"));
  Require(mass && Near(mass->GetVal(), 1.306),
          "Generated file has incorrect mass metadata");
  Require(scale && scale->GetVal() > 0.0,
          "Generated file has invalid normalization-scale metadata");
  Require(rawH000 && rawH000->GetVal() > 0.0,
          "Generated file has invalid raw-intensity metadata");

  file->Close();
  std::filesystem::remove(output);
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      throw std::invalid_argument("Expected a generated ROOT output path");
    }
    CheckMassModel();
    CheckGeneration(argv[1]);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "MassModelsTest: " << error.what() << '\n';
    return 1;
  }
}
