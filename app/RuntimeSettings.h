#pragma once

// Boundary between the CLI and ROOT's run-time C++ interpreter.  Keeping the
// interpreter details private prevents them leaking into the public library.
// The path may be omitted; RuntimeSettings.C then searches the current checkout
// before falling back to the source directory recorded by CMake.

#include "emi/Config.h"

#include <filesystem>

namespace emi::runtime {

UserSettingsConfig LoadUserSettings(const std::filesystem::path& requestedPath = {});

} // namespace emi::runtime
