#include "biodynamo.h"
#include <string>
#include "command_line_options.h"
#include "cpptoml/cpptoml.h"

namespace bdm {

/// Return only the executable name given the path
/// @param path path and filename of the executable
/// e.g. `executable`, `./executable`, './build/executable'
/// @return `executabl`
std::string ExtractExecutableName(const char* path) {
  std::string s(path);
  auto pos = s.find_last_of("/");
  if (pos == std::string::npos) {
    return s;
  } else {
    return s.substr(pos + 1, s.length() - 1);
  }
}

void InitializeBiodynamo(int argc, const char** argv) {
  // Removing this line causes an unexplainable segfault due to setting the
  // gErrorIngoreLevel global parameter of ROOT. We need to log at least one
  // thing before setting that parameter.
  Param::executable_name_ = ExtractExecutableName(argv[0]);
  auto options = bdm::DefaultSimulationOptionParser(argc, argv);
  constexpr auto kConfigFile = "bdm.toml";
  constexpr auto kConfigFileParentDir = "../bdm.toml";
  if (FileExists(kConfigFile)) {
    auto config = cpptoml::parse_file(kConfigFile);
    Param::AssignFromConfig(config);
  } else if (FileExists(kConfigFileParentDir)) {
    auto config = cpptoml::parse_file(kConfigFileParentDir);
    Param::AssignFromConfig(config);
  } else {
    Log::Warning("InitializeBiodynamo",
                 "Config file %s not found in `.` or `../` directory.",
                 kConfigFile);
  }
  if (options.backup_file_ != "") {
    Param::backup_file_ = options.backup_file_;
    Param::restore_file_ = options.restore_file_;
  }
}

void InitializeBiodynamo(const std::string& executable_name) {
  const char* argv[1] = {executable_name.c_str()};
  InitializeBiodynamo(1, argv);
}

}  // namespace bdm
