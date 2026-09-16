// Load a user settings header at run time through ROOT's Cling interpreter.
// This preserves a compiled C++ configuration API while letting tutorial
// users change waves and bounds without rebuilding the executable.
#include "RuntimeSettings.h"

#include "TInterpreter.h"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>

#ifndef EMI_SOURCE_DIR
#define EMI_SOURCE_DIR "."
#endif

namespace emi::runtime {
namespace {

std::string QuoteForInclude(const std::filesystem::path& path) {
  // The path is inserted into generated C++ source, so backslashes and quotes
  // must be escaped according to C++ string-literal rules before Cling sees it.
  std::string value = std::filesystem::absolute(path).string();
  std::string escaped;
  escaped.reserve(value.size());
  for (char character : value) {
    if (character == '\\' || character == '"') escaped += '\\';
    escaped += character;
  }
  return escaped;
}

std::filesystem::path FindSettingsFile(
    const std::filesystem::path& requestedPath) {
  if (!requestedPath.empty()) return requestedPath;

  const std::filesystem::path local = "app/UserSettings.h";
  // Prefer the current checkout so copied study configurations work naturally;
  // fall back to the source directory recorded when this binary was built.
  if (std::filesystem::exists(local)) return local;
  return std::filesystem::path(EMI_SOURCE_DIR) / "app/UserSettings.h";
}

} // namespace

UserSettingsConfig LoadUserSettings(
    const std::filesystem::path& requestedPath) {
  const std::filesystem::path settingsPath = FindSettingsFile(requestedPath);
  if (!std::filesystem::is_regular_file(settingsPath)) {
    throw std::runtime_error(
        "Could not find user settings file: " + settingsPath.string());
  }

  gInterpreter->AddIncludePath(
      (std::filesystem::path(EMI_SOURCE_DIR) / "include").c_str());
  // Declare compiles the user header in ROOT's Cling interpreter. The header's
  // inline functions then have the same public Config.h types as the executable.
  const std::string declaration =
      "#include \"" + QuoteForInclude(settingsPath) + "\"";
  if (!gInterpreter->Declare(declaration.c_str())) {
    throw std::runtime_error(
        "Could not compile user settings: " + settingsPath.string());
  }

  TInterpreter::EErrorCode error = TInterpreter::kNoError;
  const Longptr_t address = gInterpreter->Calc(
      "new emi::UserSettingsConfig(emi::user::Settings())", &error);
  // Calc returns an integer large enough to hold a pointer. Move the aggregate
  // into ordinary compiled C++ storage, then delete the Cling-created object so
  // no interpreter-owned lifetime escapes this boundary.
  if (error != TInterpreter::kNoError || address == 0) {
    throw std::runtime_error(
        "User settings must define emi::user::Settings()");
  }

  auto* loaded = reinterpret_cast<UserSettingsConfig*>(address);
  UserSettingsConfig settings = std::move(*loaded);
  delete loaded;
  return settings;
}

} // namespace emi::runtime
