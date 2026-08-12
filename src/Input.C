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
  if (nData <= 1) {
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
  int L = 0;
  int M = 0;
  bool mixed04 = false;
  std::string valueName;
};

bool ParseMomentBranch(const std::string& name, bool photoproduction,
                       MomentBranch& moment) {
  static const std::regex mixedPattern(R"(^(R?H04)_([0-9]+)_([0-9]+)$)");
  static const std::regex rawPattern(R"(^(R?H)_([0-9]+)_([0-9]+)_([0-9]+)$)");
  std::smatch match;
  if (!photoproduction && std::regex_match(name, match, mixedPattern)) {
    moment = {0, std::stoi(match[2]), std::stoi(match[3]), true, name};
    return true;
  }
  if (!std::regex_match(name, match, rawPattern)) return false;

  const int alpha = std::stoi(match[2]);
  if (photoproduction ? alpha > 3 : (alpha == 0 || alpha == 4 || alpha > 8)) {
    return false;
  }
  moment = {alpha, std::stoi(match[3]), std::stoi(match[4]), false, name};
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
  if (cfg.bin < 0 || cfg.bin >= kMaximumInputBins) throw std::runtime_error("Bin out of range");

  std::map<std::tuple<int, int, int, bool>, MomentBranch> branches;
  TObjArray* branchList = t->GetListOfBranches();
  for (int index = 0; branchList && index < branchList->GetEntries(); ++index) {
    const std::string name = branchList->At(index)->GetName();
    MomentBranch candidate;
    if (!ParseMomentBranch(name, cfg.photoproduction, candidate)) continue;
    const auto key = std::make_tuple(candidate.alpha, candidate.L,
                                     candidate.M, candidate.mixed04);
    auto found = branches.find(key);
    const bool isExperimentalName = name.rfind("RH", 0) == 0;
    if (found == branches.end() || isExperimentalName) branches[key] = std::move(candidate);
  }

  std::vector<ObservedMoment> obs;
  obs.reserve(branches.size());
  for (const auto& entry : branches) {
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
