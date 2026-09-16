// Validate and read a moment bin from a ROOT TTree.  Branch discovery accepts
// the project's indexed moment convention, then maps the observed covariance
// into the exact ordering used by the selected physics model.
#include "Detail.h"

#include "TBranch.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace emi::detail {

namespace {

bool TryReadArrayBranchElement(TTree* t,
                                      const char* branchName,
                                      int bin,
                                      double& value) {
  TBranch* branch = t->GetBranch(branchName);
  if (!branch) {
    value = 0.0;
    return false;
  }

  TLeaf* leaf = branch->GetLeaf(branchName);
  if (!leaf && branch->GetListOfLeaves() && branch->GetListOfLeaves()->GetEntries() > 0) {
    leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
  }
  if (!leaf) {
    value = 0.0;
    return false;
  }

  t->GetEntry(0);
  const int nData = leaf->GetNdata();
  // Published inputs store all kinematic bins in one array branch, whereas
  // generated closure files use scalar branches. Supporting both shapes lets
  // the same fit path consume experimental and synthetic data.
  if (nData <= 1) {
    if (bin != 0) {
      value = 0.0;
      return false;
    }
    value = leaf->GetValue(0);
    return true;
  }

  if (bin < 0 || bin >= nData) {
    value = 0.0;
    return false;
  }

  value = leaf->GetValue(bin);
  return true;
}

struct MomentBranch {
  int alpha = 0;
  int beta = 0;
  int delta = 0;
  int L = 0;
  int M = 0;
  bool mixed04 = false;
  bool tensor = false;
  std::string valueName;
};

bool ParseMomentBranch(const std::string& name, const InternalConfig& cfg,
                       MomentBranch& moment) {
  // Four spellings are accepted: raw/mixed photon responses, each with or
  // without explicit target and recoil Pauli indices. The selected
  // polarization determines which tensor components are observable.
  static const std::regex mixedPattern(R"(^(R?H04)_([0-9]+)_([0-9]+)$)");
  static const std::regex mixedTensorPattern(
      R"(^(R?H04)_([0-3])_([0-3])_([0-9]+)_([0-9]+)$)");
  static const std::regex rawPattern(R"(^(R?H)_([0-9]+)_([0-9]+)_([0-9]+)$)");
  static const std::regex rawTensorPattern(
      R"(^(R?H)_([0-9]+)_([0-3])_([0-3])_([0-9]+)_([0-9]+)$)");
  std::smatch match;
  const bool explicitPolarization = cfg.nucleonPolarization != NucleonPolarization::None;
  const bool initial = cfg.nucleonPolarization == NucleonPolarization::Initial ||
                       cfg.nucleonPolarization == NucleonPolarization::Both;
  const bool recoil = cfg.nucleonPolarization == NucleonPolarization::Recoil ||
                      cfg.nucleonPolarization == NucleonPolarization::Both;
  auto allowed = [&](int beta, int delta) {
    return (initial || beta == 0) && (recoil || delta == 0);
  };

  if (!cfg.photoproduction &&
      std::regex_match(name, match, mixedTensorPattern)) {
    const int beta = std::stoi(match[2]);
    const int delta = std::stoi(match[3]);
    if (explicitPolarization ? !allowed(beta, delta)
                             : (beta != 0 || delta != 0)) {
      return false;
    }
    moment = {0, beta, delta, std::stoi(match[4]), std::stoi(match[5]),
              true, true, name};
    return true;
  }
  if (!cfg.photoproduction && std::regex_match(name, match, mixedPattern)) {
    moment = {0, 0, 0, std::stoi(match[2]), std::stoi(match[3]),
              true, false, name};
    return true;
  }

  if (std::regex_match(name, match, rawTensorPattern)) {
    const int alpha = std::stoi(match[2]);
    const int beta = std::stoi(match[3]);
    const int delta = std::stoi(match[4]);
    if ((explicitPolarization ? !allowed(beta, delta)
                              : (beta != 0 || delta != 0)) ||
        (cfg.photoproduction ? alpha > 3 : (alpha == 0 || alpha == 4 || alpha > 8))) {
      return false;
    }
    moment = {alpha, beta, delta, std::stoi(match[5]), std::stoi(match[6]),
              false, true, name};
    return true;
  }
  if (!std::regex_match(name, match, rawPattern)) return false;

  const int alpha = std::stoi(match[2]);
  if (cfg.photoproduction ? alpha > 3 : (alpha == 0 || alpha == 4 || alpha > 8)) {
    return false;
  }
  moment = {alpha, 0, 0, std::stoi(match[3]), std::stoi(match[4]),
            false, false, name};
  return true;
}

} // namespace

std::vector<ObservedMoment> ReadObservedMoments(const InternalConfig& cfg) {
  std::unique_ptr<TFile> fin(TFile::Open(cfg.momentsFile.c_str(), "READ"));
  if (!fin || fin->IsZombie()) throw std::runtime_error("Failed to open file " + cfg.momentsFile);
  TTree* t = dynamic_cast<TTree*>(fin->Get(cfg.momentsTree.c_str()));
  if (!t) throw std::runtime_error("Could not find tree '" + cfg.momentsTree + "'");
  if (t->GetEntries() < 1) {
    throw std::runtime_error("Tree '" + cfg.momentsTree + "' is empty");
  }
  if (cfg.bin < 0) throw std::runtime_error("Bin index must be nonnegative");

  std::map<std::tuple<int, int, int, int, int, bool>, MomentBranch> branches;
  // The map key is physical rather than textual. It collapses equivalent H/RH
  // and unpolarized/tensor spellings before any data are read.
  TObjArray* branchList = t->GetListOfBranches();
  for (int index = 0; branchList && index < branchList->GetEntries(); ++index) {
    const std::string name = branchList->At(index)->GetName();
    MomentBranch candidate;
    if (!ParseMomentBranch(name, cfg, candidate)) continue;
    const auto key = std::make_tuple(candidate.alpha, candidate.beta, candidate.delta,
                                     candidate.L, candidate.M, candidate.mixed04);
    auto found = branches.find(key);
    auto priority = [](const MomentBranch& branch) {
      // Prefer experimental RH names over internal H names. For an
      // unpolarized fit, prefer alpha,L,M branches when both spellings exist,
      // but accept alpha,0,0,L,M from a polarized file as a fallback.
      return (branch.valueName.rfind("RH", 0) == 0 ? 2 : 0) +
             (branch.tensor ? 0 : 1);
    };
    if (found == branches.end() || priority(candidate) > priority(found->second)) {
      branches[key] = std::move(candidate);
    }
  }

  std::vector<ObservedMoment> obs;
  obs.reserve(branches.size());
  for (const auto& entry : branches) {
    // Missing or non-positive uncertainties cannot contribute a finite
    // standardized residual, so they are reported and excluded.
    const MomentBranch& branch = entry.second;
    const std::string errName = branch.valueName + "_err";
    double val = 0.0;
    double sig = 0.0;
    if (!TryReadArrayBranchElement(t, branch.valueName.c_str(), cfg.bin, val)) {
      if (cfg.verbose) {
        std::cout << "Skipping " << branch.valueName
                  << " because bin " << cfg.bin
                  << " is not present in that branch." << std::endl;
      }
      continue;
    }
    if (!TryReadArrayBranchElement(t, errName.c_str(), cfg.bin, sig) ||
        !(sig > 0.0) || !std::isfinite(sig)) {
      if (cfg.verbose) {
        std::cout << "Skipping " << branch.valueName
                  << " because its uncertainty is not positive\n";
      }
      continue;
    }

    ObservedMoment m;
    m.alpha = branch.alpha;
    m.beta = branch.beta;
    m.delta = branch.delta;
    m.L = branch.L;
    m.M = branch.M;
    m.value = val;
    m.sigma = sig;
    m.isMixed04 = branch.mixed04;
    m.name = branch.valueName;
    obs.push_back(std::move(m));
  }
  return obs;
}

} // namespace emi::detail
