#pragma once

// Boundary between the CLI and ROOT's run-time C++ interpreter.  Keeping the
// interpreter details private prevents them leaking into the public library.

#include "emi/Config.h"

#include <filesystem>

namespace emi::runtime {

UserSettingsConfig LoadUserSettings(
    const std::filesystem::path& requestedPath = {});

} // namespace emi::runtime
