#include "emi/Runner.h"
#include "RuntimeSettings.h"

#include <exception>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace {

class Arguments {
public:
  Arguments(int argc, char** argv, int first) {
    for (int i = first; i < argc; ++i) {
      std::string key = argv[i];
      if (key.rfind("--", 0) != 0) {
        throw std::invalid_argument("Unexpected argument: " + key);
      }
      key.erase(0, 2);
      if (i + 1 < argc && std::string_view(argv[i + 1]).rfind("--", 0) != 0) {
        values_[key] = argv[++i];
      } else {
        flags_.insert(key);
      }
    }
  }

  bool Has(const std::string& key) const { return flags_.count(key) || values_.count(key); }

  std::string Get(const std::string& key, std::string fallback = {}) const {
    const auto found = values_.find(key);
    return found == values_.end() ? std::move(fallback) : found->second;
  }

  unsigned Unsigned(const std::string& key, unsigned fallback) const {
    return Has(key) ? static_cast<unsigned>(std::stoul(Get(key))) : fallback;
  }

  int Integer(const std::string& key, int fallback) const {
    return Has(key) ? std::stoi(Get(key)) : fallback;
  }

  double Real(const std::string& key, double fallback) const {
    return Has(key) ? std::stod(Get(key)) : fallback;
  }

  void RequireOnly(std::initializer_list<std::string_view> valueOptions,
                   std::initializer_list<std::string_view> flagOptions) const {
    auto contains = [](std::string_view key, const auto& options) {
      for (std::string_view option : options) {
        if (key == option) return true;
      }
      return false;
    };
    for (const auto& entry : values_) {
      const auto& key = entry.first;
      if (!contains(key, valueOptions)) {
        throw std::invalid_argument("Option --" + key +
                                    " is not valid for this command");
      }
    }
    for (const auto& key : flags_) {
      if (contains(key, valueOptions)) {
        throw std::invalid_argument("Option --" + key + " requires a value");
      }
      if (!contains(key, flagOptions)) {
        throw std::invalid_argument("Option --" + key +
                                    " is not valid for this command");
      }
    }
  }

private:
  std::unordered_map<std::string, std::string> values_;
  std::unordered_set<std::string> flags_;
};

void PrintHelp() {
  std::cout << R"(Electroproduction Moments Inverter

Usage:
  emi fit --input FILE --output FILE [options]
  emi bootstrap --input FILE --output FILE [options]
  emi make-lepto [--dataset e_rho] [--output FILE]
  emi make-photo [--dataset gluex] [--output FILE]
  emi generate-fixed [options]    (also accepts: emi --generate-fixed)
  emi generate-random [--output FILE] [--events N] [--seed N]

Fit options:
  --tree NAME             override the configured input tree
  --bin N                 override the configured input bin
  --epsilon VALUE         override virtual-photon polarisation
  --starts N              override the configured random-start count
  --workers N             override processes; 0 uses available CPUs
  --seed N                override seed; 0 selects a random seed
  --photo                 use the photoproduction model
  --electro               use the electroproduction model
  --polarization MODE     nucleon mode: initial, recoil, or both (default: none)
  --k-gauge REF           SO(2) real-wave reference, e.g. b:T:1:1
  --no-hesse              skip the Hessian uncertainty calculation
  --numerical-gradients    use Minuit's numerical derivatives for cross-checks
  --quiet                 suppress progress output
  --verbose               enable progress output

Bootstrap-only options:
  --toys N                Gaussian bootstrap samples (default: 1000)
  --starts-per-toy N      minimisations per sample (default: 1000)
  --hesse                 calculate a Hessian for every sample (off by default)

Synthetic-generation options:
  generate-fixed          use FixedMoments() from app/UserSettings.h
  --output FILE           override FixedMoments().output
  --epsilon VALUE         override FixedMoments().epsilon
  --seed N                override FixedMoments().seed
  --suppression VALUE     override suppression in a dominant-K mode
  --mass VALUE            PhotoTest invariant mass in GeV
  --k-minus-scale VALUE   scale the complete PhotoTest k=- sector
  --no-background         disable the PhotoTest S-wave background
  --quiet                 do not print generated amplitudes
  --verbose               print generated amplitudes

Configuration:
  --settings FILE         runtime C++ settings file (default: app/UserSettings.h)
)";
}

emi::FitConfig ReadFit(const Arguments& args, emi::FitConfig fit) {
  if (args.Has("input")) fit.input = args.Get("input");
  if (args.Has("output")) fit.output = args.Get("output");
  if (args.Has("tree")) fit.tree = args.Get("tree");
  if (args.Has("bin")) fit.bin = args.Integer("bin", fit.bin);
  if (args.Has("epsilon")) fit.epsilon = args.Real("epsilon", fit.epsilon);
  if (args.Has("starts")) fit.starts = args.Unsigned("starts", fit.starts);
  if (args.Has("workers")) fit.workers = args.Unsigned("workers", fit.workers);
  if (args.Has("seed")) fit.seed = args.Unsigned("seed", fit.seed);
  if (args.Has("photo") && args.Has("electro")) {
    throw std::invalid_argument("Choose either --photo or --electro, not both");
  }
  if (args.Has("photo")) fit.photoproduction = true;
  if (args.Has("electro")) fit.photoproduction = false;
  if (args.Has("no-hesse")) fit.hesse = false;
  if (args.Has("numerical-gradients")) fit.useNumericalGradients = true;
  if (args.Has("quiet")) fit.verbose = false;
  if (args.Has("verbose")) fit.verbose = true;
  return fit;
}

emi::NucleonPolarization ReadPolarization(const Arguments& args,
                                          emi::NucleonPolarization fallback) {
  if (!args.Has("polarization")) return fallback;
  const std::string value = args.Get("polarization");
  if (value == "initial") return emi::NucleonPolarization::Initial;
  if (value == "recoil") return emi::NucleonPolarization::Recoil;
  if (value == "both") return emi::NucleonPolarization::Both;
  if (value == "none") return emi::NucleonPolarization::None;
  throw std::invalid_argument(
      "Polarization must be one of: none, initial, recoil, both");
}

emi::ModelConfig ReadGauge(const Arguments& args, emi::ModelConfig model) {
  model.UseNucleonPolarization(
      ReadPolarization(args, model.nucleonPolarization));
  if (args.Has("k-gauge")) {
    const std::string value = args.Get("k-gauge");
    std::istringstream input(value);
    char reflectivity = '\0';
    char orientation = '\0';
    char separator1 = '\0';
    char separator2 = '\0';
    char separator3 = '\0';
    int l = -1;
    int m = 0;
    if (!(input >> reflectivity >> separator1 >> orientation >> separator2 >> l >>
          separator3 >> m) || separator1 != ':' || separator2 != ':' ||
        separator3 != ':' || input.peek() != std::char_traits<char>::eof()) {
      throw std::invalid_argument(
          "k-gauge must have the form reflectivity:orientation:l:m");
    }
    model.SetKMixingGauge(reflectivity, orientation, l, m);
  }
  return model;
}

emi::FixedMomentsConfig ReadFixedGeneration(
    const Arguments& args, emi::FixedMomentsConfig generation) {
  args.RequireOnly({"settings", "output", "epsilon", "seed", "suppression",
                    "mass", "k-minus-scale"},
                   {"quiet", "verbose", "no-background"});
  if (args.Has("quiet") && args.Has("verbose")) {
    throw std::invalid_argument("Choose either --quiet or --verbose, not both");
  }
  if (args.Has("output")) generation.output = args.Get("output");
  if (args.Has("epsilon")) {
    generation.epsilon = args.Real("epsilon", generation.epsilon);
  }
  if (args.Has("seed")) {
    generation.seed = args.Unsigned("seed", generation.seed);
  }
  if (args.Has("suppression")) {
    const bool dominantK =
        generation.mode == emi::FixedGenerationMode::KPositiveDominant ||
        generation.mode == emi::FixedGenerationMode::KReflectivitySplit;
    if (!dominantK) {
      throw std::invalid_argument(
          "--suppression applies only to KPositiveDominant or "
          "KReflectivitySplit in FixedMoments()");
    }
    generation.suppression =
        args.Real("suppression", generation.suppression);
  }
  const bool photoTest =
      generation.mode == emi::FixedGenerationMode::PhotoTest;
  if (!photoTest &&
      (args.Has("mass") || args.Has("k-minus-scale") ||
       args.Has("no-background"))) {
    throw std::invalid_argument(
        "--mass, --k-minus-scale, and --no-background apply only to "
        "PhotoTest in FixedMoments()");
  }
  if (args.Has("mass")) {
    generation.massModel.massGeV = args.Real("mass", 0.0);
  }
  if (args.Has("k-minus-scale")) {
    generation.massModel.kMinusScale = args.Real(
        "k-minus-scale", generation.massModel.kMinusScale);
  }
  if (args.Has("no-background")) {
    generation.massModel.backgroundEnabled = false;
  }
  if (args.Has("quiet")) generation.printAmplitudes = false;
  if (args.Has("verbose")) generation.printAmplitudes = true;
  return generation;
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2 || std::string_view(argv[1]) == "--help" ||
        std::string_view(argv[1]) == "help") {
      PrintHelp();
      return 0;
    }

    const std::string command = std::string_view(argv[1]) == "--generate-fixed"
                                    ? "generate-fixed"
                                    : argv[1];
    const Arguments args(argc, argv, 2);
    auto loadSettings = [&] {
      return emi::runtime::LoadUserSettings(args.Get("settings"));
    };

    if (command == "fit") {
      auto settings = loadSettings();
      emi::RunFit(ReadFit(args, settings.fit),
                  ReadGauge(args, settings.model));

    } else if (command == "bootstrap") {
      auto settings = loadSettings();
      emi::BootstrapConfig bootstrap = settings.bootstrap;
      bootstrap.fit = ReadFit(args, bootstrap.fit);
      bootstrap.fit.starts = args.Unsigned("starts-per-toy", bootstrap.startsPerToy);
      if (args.Has("hesse")) bootstrap.fit.hesse = true;
      bootstrap.toys = args.Unsigned("toys", bootstrap.toys);
      bootstrap.startsPerToy = bootstrap.fit.starts;
      emi::RunBootstrap(bootstrap, ReadGauge(args, settings.model));

    } else if (command == "make-lepto") {
      emi::MakeLeptoproductionMoments(args.Get("dataset", "e_rho"),
                                      args.Get("output"),
                                      args.Get("tree", "expMoments"));
    } else if (command == "make-photo") {
      emi::MakePhotoproductionMoments(args.Get("dataset", "gluex"),
                                      args.Get("output"),
                                      args.Get("tree", "expMoments"));
    } else if (command == "generate-fixed") {
      emi::GenerateFixedMoments(
          ReadFixedGeneration(args, loadSettings().fixedMoments));
    } else if (command == "generate-example") {
      const auto settings = loadSettings();
      emi::GenerateFixedMoments(args.Get("output"),
                                args.Real("epsilon", 0.8), !args.Has("quiet"),
                                ReadGauge(args, settings.model),
                                args.Has("photo"));
    } else if (command == "generate-random") {
      const auto settings = loadSettings();
      emi::GenerateRandomMoments(args.Get("output"),
                                 args.Unsigned("events", 1000000),
                                 args.Real("epsilon", 1.0),
                                 args.Unsigned("seed", 0),
                                 ReadGauge(args, settings.model));
    } else {
      throw std::invalid_argument("Unknown command: " + command);
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
