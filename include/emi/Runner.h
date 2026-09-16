#pragma once

// High-level public API.  Each function owns a complete analysis operation and
// writes a self-describing ROOT file, while Config.h carries all user choices.

#include "emi/Config.h"

#include <filesystem>
#include <string>

namespace emi {

// Fit one selected input bin from many random starts. The output tree keeps all
// starts so callers can inspect alternative minima and convergence failures.
void RunFit(const FitConfig& config, const ModelConfig& model = {});

// Gaussian-fluctuate one input bin and retain the best start for each toy.
void RunBootstrap(const BootstrapConfig& config, const ModelConfig& model = {});

// Rebuild the published moment ROOT inputs embedded in the corresponding C++
// source. An empty output chooses the tutorial-local canonical input path.
void MakeLeptoproductionMoments(const std::string& dataset,
                                const std::filesystem::path& output = {},
                                const std::string& tree = "expMoments");
void MakePhotoproductionMoments(const std::string& dataset,
                               const std::filesystem::path& output = {},
                               const std::string& tree = "expMoments");

// Evaluate the generation mode and amplitudes supplied by UserSettings.h.
void GenerateFixedMoments(const FixedMomentsConfig& config);

// Generate the internal deterministic example used by closure-test scripts;
// this overload deliberately ignores the editable FixedMoments configuration.
void GenerateFixedMoments(const std::filesystem::path& output,
                          double epsilon = 0.8,
                          bool printAmplitudes = true,
                          const ModelConfig& model = {},
                          bool photoproduction = false);

void GenerateRandomMoments(const std::filesystem::path& output,
                           unsigned events,
                           double epsilon = 1.0,
                           std::uint32_t seed = 0,
                           const ModelConfig& model = {});

} // namespace emi
