#pragma once

#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

namespace emi {

struct FitConfig {
  std::filesystem::path input;
  std::string tree = "expMoments";
  int bin = 0;
  std::filesystem::path output;
  double epsilon = 1.0;
  bool photoproduction = false;

  unsigned starts = 10000;
  unsigned workers = 0;
  std::uint32_t seed = 0;
  bool hesse = true;
  bool useNumericalGradients = false;
  bool verbose = true;

  FitConfig& SetStarts(unsigned value) { starts = value; return *this; }
  FitConfig& SetWorkers(unsigned value) { workers = value; return *this; }
  FitConfig& SetSeed(std::uint32_t value) { seed = value; return *this; }
  FitConfig& UseHesse(bool value = true) { hesse = value; return *this; }
  FitConfig& UseNumericalGradients(bool value = true) {
    useNumericalGradients = value;
    return *this;
  }
};

struct BootstrapConfig {
  FitConfig fit;
  unsigned toys = 1000;
  unsigned startsPerToy = 1000;
};

struct Wave {
  int l;
  int m;
};

struct ModelConfig {
  std::vector<Wave> waves = {{1, -1}, {1, 0}, {1, 1}};
  bool usePositiveReflectivity = true;
  bool useNegativeReflectivity = true;
  bool enforceLongitudinalParity = true;
  double normalisationMoment = 2.0;

  ModelConfig& SetWaves(std::initializer_list<Wave> values) {
    waves.assign(values);
    return *this;
  }
  ModelConfig& SetWaves(std::vector<Wave> values) {
    waves = std::move(values);
    return *this;
  }
  ModelConfig& AddWave(int l, int m) { waves.push_back({l, m}); return *this; }
  ModelConfig& UseReflectivities(bool positive, bool negative) {
    usePositiveReflectivity = positive;
    useNegativeReflectivity = negative;
    return *this;
  }
  ModelConfig& EnforceLongitudinalParity(bool value = true) {
    enforceLongitudinalParity = value;
    return *this;
  }
};

} // namespace emi
