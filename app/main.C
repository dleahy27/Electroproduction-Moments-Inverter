#include "emi/Runner.h"
#include "UserSettings.h"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

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
  emi generate-fixed [--output FILE] [--epsilon VALUE]
  emi generate-random [--output FILE] [--events N] [--seed N]

Fit options:
  --tree NAME             input tree (default: expMoments)
  --bin N                 input bin (default: 0)
  --epsilon VALUE         virtual-photon polarisation (default: 1)
  --starts N              random starts (default: 10000)
  --workers N             processes; 0 uses available CPUs (default: 0)
  --seed N                reproducible random seed; 0 is random (default: 0)
  --photo                 use the photoproduction model
  --no-hesse              skip the Hessian uncertainty calculation
  --numerical-gradients    use Minuit's numerical derivatives for cross-checks
  --quiet                 suppress progress output

Bootstrap-only options:
  --toys N                Gaussian bootstrap samples (default: 1000)
  --starts-per-toy N      minimisations per sample (default: 1000)
  --hesse                 calculate a Hessian for every sample (off by default)
)";
}

emi::FitConfig ReadFit(const Arguments& args) {
  emi::FitConfig fit;
  emi::user::ConfigureFit(fit);
  fit.input = args.Get("input");
  fit.output = args.Get("output");
  fit.tree = args.Get("tree", "expMoments");
  fit.bin = args.Integer("bin", 0);
  fit.epsilon = args.Real("epsilon", 1.0);
  fit.starts = args.Unsigned("starts", fit.starts);
  fit.workers = args.Unsigned("workers", fit.workers);
  fit.seed = args.Unsigned("seed", fit.seed);
  fit.photoproduction = args.Has("photo");
  if (args.Has("no-hesse")) fit.hesse = false;
  if (args.Has("numerical-gradients")) fit.useNumericalGradients = true;
  fit.verbose = !args.Has("quiet");
  return fit;
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2 || std::string_view(argv[1]) == "--help" ||
        std::string_view(argv[1]) == "help") {
      PrintHelp();
      return 0;
    }

    const std::string command = argv[1];
    const Arguments args(argc, argv, 2);
    if (command == "fit") {
      emi::RunFit(ReadFit(args), emi::user::Model());
    } else if (command == "bootstrap") {
      emi::BootstrapConfig bootstrap;
      bootstrap.fit = ReadFit(args);
      emi::user::ConfigureBootstrap(bootstrap);
      bootstrap.fit.starts = args.Unsigned("starts-per-toy", bootstrap.startsPerToy);
      if (args.Has("hesse")) bootstrap.fit.hesse = true;
      bootstrap.toys = args.Unsigned("toys", bootstrap.toys);
      bootstrap.startsPerToy = bootstrap.fit.starts;
      emi::RunBootstrap(bootstrap, emi::user::Model());
    } else if (command == "make-lepto") {
      emi::MakeLeptoproductionMoments(args.Get("dataset", "e_rho"),
                                      args.Get("output"),
                                      args.Get("tree", "expMoments"));
    } else if (command == "make-photo") {
      emi::MakePhotoproductionMoments(args.Get("dataset", "gluex"),
                                      args.Get("output"),
                                      args.Get("tree", "expMoments"));
    } else if (command == "generate-fixed") {
      emi::GenerateFixedMoments(args.Get("output"),
                                args.Real("epsilon", 0.8), !args.Has("quiet"),
                                emi::user::GenerationModel());
    } else if (command == "generate-random") {
      emi::GenerateRandomMoments(args.Get("output"),
                                 args.Unsigned("events", 1000000),
                                 args.Real("epsilon", 1.0),
                                 args.Unsigned("seed", 0),
                                 emi::user::GenerationModel());
    } else {
      throw std::invalid_argument("Unknown command: " + command);
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
