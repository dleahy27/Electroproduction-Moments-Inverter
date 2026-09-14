#pragma once

#include "emi/Config.h"

#include <filesystem>

namespace emi::runtime {

UserSettingsConfig LoadUserSettings(
    const std::filesystem::path& requestedPath = {});

} // namespace emi::runtime
