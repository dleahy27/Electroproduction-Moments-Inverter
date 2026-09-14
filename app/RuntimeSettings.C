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
  const std::string declaration =
      "#include \"" + QuoteForInclude(settingsPath) + "\"";
  if (!gInterpreter->Declare(declaration.c_str())) {
    throw std::runtime_error(
        "Could not compile user settings: " + settingsPath.string());
  }

  TInterpreter::EErrorCode error = TInterpreter::kNoError;
  const Longptr_t address = gInterpreter->Calc(
      "new emi::UserSettingsConfig(emi::user::Settings())", &error);
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
