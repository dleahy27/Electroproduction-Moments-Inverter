#pragma once

#include "emi/Config.h"

#include <filesystem>
#include <string>

namespace emi {

void RunFit(const FitConfig& config, const ModelConfig& model = {});
void RunBootstrap(const BootstrapConfig& config, const ModelConfig& model = {});

void MakeLeptoproductionMoments(const std::string& dataset,
                                const std::filesystem::path& output = {},
                                const std::string& tree = "expMoments");
void MakePhotoproductionMoments(const std::string& dataset,
                               const std::filesystem::path& output = {},
                               const std::string& tree = "expMoments");

void GenerateFixedMoments(const std::filesystem::path& output,
                          double epsilon = 0.8,
                          bool printValues = true,
                          const ModelConfig& model = {});
void GenerateRandomMoments(const std::filesystem::path& output,
                           unsigned events,
                           double epsilon = 1.0,
                           std::uint32_t seed = 0,
                           const ModelConfig& model = {});

} // namespace emi
