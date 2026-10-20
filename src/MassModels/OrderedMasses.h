#pragma once

// Continuous bounded base/gap coordinates for shared resonance masses.
// Kept independent of ROOT so the mapping, Jacobian and covariance can be
// checked without starting a minimizer or interpreter.
#include "emi/Config.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace emi::MassModels {

struct FitCoordinate {
  std::string name;
  MassParameter control;
};

using JacobianRow = std::vector<std::pair<int, double>>;
struct PhysicalPoint {
  std::vector<double> values;
  std::vector<JacobianRow> jacobian;
};

class OrderedMasses {
public:
  struct Node {
    enum class Kind { Fixed, Base, Gap, Fraction, LogBase, LogGap };
    std::string group, resonance;
    int physical = -1, previous = -1, fit = -1;
    double initial = 0.0, lower = 0.0, upper = 0.0;
    Kind kind = Kind::Fixed;
  };

  std::vector<FitCoordinate> physicalCoordinates, fitCoordinates;
  std::vector<Node> nodes;
  double minimumGap = 1.e-6;

  OrderedMasses() = default;
  OrderedMasses(const MassModelConfig& model,
                const std::vector<FitCoordinate>& coordinates)
      : physicalCoordinates(coordinates), minimumGap(model.minimumMassGapGeV) {
    if (!std::isfinite(minimumGap) || minimumGap <= 0.0)
      throw std::invalid_argument("minimumMassGapGeV must be finite and positive");
    std::unordered_map<std::string, int> physicalIndex;
    for (size_t i = 0; i < coordinates.size(); ++i)
      physicalIndex.emplace(coordinates[i].name, static_cast<int>(i));
    std::vector<bool> positive(coordinates.size(), false);
    for (const auto& resonance : model.resonances) {
      for (const auto& suffix : {"_mass", "_width"}) {
        const auto found = physicalIndex.find(resonance.name + suffix);
        if (found != physicalIndex.end()) positive[found->second] = true;
      }
      for (const auto& channel : resonance.channels) {
        const auto found = physicalIndex.find(resonance.name + "_" + channel.name + "_g2");
        if (found != physicalIndex.end()) positive[found->second] = true;
      }
    }
    struct Record { const MassResonanceConfig* resonance; int index; };
    std::map<std::string, std::vector<Record>> groups;
    std::unordered_map<std::string, const MassResonanceConfig*> shared;
    for (const auto& resonance : model.resonances) {
      const auto inserted = shared.emplace(resonance.name, &resonance);
      if (!inserted.second) {
        const auto& previous = *inserted.first->second;
        if (previous.wave.l != resonance.wave.l ||
            previous.massOrderingGroup != resonance.massOrderingGroup)
          throw std::invalid_argument("Inconsistent shared resonance ordering sector: " + resonance.name);
        continue;
      }
      const std::string group = resonance.massOrderingGroup.empty()
          ? "l" + std::to_string(resonance.wave.l) : resonance.massOrderingGroup;
      // Reserve a branch-safe alphabet for diagnostic coordinate names.
      if (group.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_") != std::string::npos)
        throw std::invalid_argument("Mass ordering group must use letters, digits and underscores: " + group);
      groups[group].push_back({&resonance, physicalIndex.at(resonance.name + "_mass")});
    }
    std::vector<bool> transformed(coordinates.size(), false);
    for (auto& item : groups) {
      auto& records = item.second;
      if (records.size() < 2)
        continue;
      const int l = records.front().resonance->wave.l;
      for (const auto& record : records)
        if (record.resonance->wave.l != l)
          throw std::invalid_argument("Mass ordering group mixes different angular momenta: " + item.first);
      // Sort complete resonance identities once; all other coordinates remain
      // attached to those identities, never to a dynamically sorted mass.
      std::sort(records.begin(), records.end(), [](const auto& a, const auto& b) {
        return a.resonance->poleMass.initial < b.resonance->poleMass.initial;
      });
      for (size_t i = 1; i < records.size(); ++i) {
        const double gap = records[i].resonance->poleMass.initial -
                           records[i - 1].resonance->poleMass.initial;
        if (gap == 0.0)
          throw std::invalid_argument("Ambiguous equal starting masses in ordering group " + item.first);
        if (gap < minimumGap)
          throw std::invalid_argument("Starting mass gap is below minimumMassGapGeV in group " + item.first);
      }
      // Future upper bounds (including fixed masses) limit earlier masses.
      // This is the exact feasible envelope, not independent gap bounds that
      // would silently cut away part of the physical domain.
      std::vector<double> upper(records.size());
      for (size_t reverse = records.size(); reverse > 0; --reverse) {
        const size_t i = reverse - 1;
        const auto& control = records[i].resonance->poleMass;
        upper[i] = control.fixed ? control.initial :
            (control.bounded ? control.upper : std::numeric_limits<double>::infinity());
        if (i + 1 < records.size())
          upper[i] = std::min(upper[i], upper[i + 1] - minimumGap);
        const double lower = control.fixed ? control.initial : (control.bounded ? control.lower : 0.0);
        if (lower > upper[i] || control.initial > upper[i])
          throw std::invalid_argument("Infeasible ordered mass bounds in group " + item.first);
      }
      bool previousFixed = true;
      for (size_t i = 0; i < records.size(); ++i) {
        const auto& record = records[i];
        const auto& control = record.resonance->poleMass;
        Node node;
        node.group = item.first;
        node.resonance = record.resonance->name;
        node.physical = record.index;
        node.previous = i ? records[i - 1].index : -1;
        node.initial = control.initial;
        node.lower = control.fixed ? control.initial : (control.bounded ? control.lower : 0.0);
        node.upper = upper[i];
        transformed[node.physical] = true;
        if (control.fixed || node.lower == node.upper) {
          node.kind = Node::Kind::Fixed;
          node.initial = node.lower;
        } else if (!i) {
          node.kind = std::isfinite(node.upper) ? Node::Kind::Base : Node::Kind::LogBase;
        } else {
          node.kind = !std::isfinite(node.upper) ? Node::Kind::LogGap :
              (previousFixed ? Node::Kind::Gap : Node::Kind::Fraction);
        }
        previousFixed = node.kind == Node::Kind::Fixed;
        nodes.push_back(node);
      }
    }
    direct_.assign(coordinates.size(), -1);
    logarithmic_.assign(coordinates.size(), false);
    for (size_t i = 0; i < coordinates.size(); ++i) {
      if (!transformed[i]) {
        direct_[i] = static_cast<int>(fitCoordinates.size());
        auto coordinate = coordinates[i];
        if (!coordinate.control.bounded && !coordinate.control.fixed && positive[i]) {
          if (coordinate.control.initial <= 0.0)
            throw std::invalid_argument("Unbounded positive parameter needs a positive start: " + coordinate.name);
          logarithmic_[i] = true;
          coordinate.name = "positive__" + coordinate.name;
          coordinate.control.initial = std::log(coordinate.control.initial);
        }
        fitCoordinates.push_back(std::move(coordinate));
      }
    }
    for (auto& node : nodes) {
      if (node.kind == Node::Kind::Fixed)
        continue;
      MassParameter control;
      control.fixed = false;
      std::string name = "ordered__" + node.group + "__";
      if (node.kind == Node::Kind::LogBase) {
        name += "log_base_mass";
        control = {std::log(node.initial - node.lower), 0.0, 0.0, false, false};
      } else if (node.kind == Node::Kind::Base) {
        name += "base_mass";
        control = {node.initial, node.lower, node.upper, false};
      } else {
        const double previous = coordinates[node.previous].control.initial;
        const double lowerMass = std::max(node.lower, previous + minimumGap);
        if (previous + minimumGap == previous)
          throw std::invalid_argument("minimumMassGapGeV is below mass floating-point resolution");
        name += "gap_" + node.resonance;
        if (node.kind == Node::Kind::LogGap) {
          name += "_log";
          if (node.initial <= lowerMass)
            throw std::invalid_argument("Unbounded mass gap must start above minimumMassGapGeV");
          control = {std::log(node.initial - lowerMass), 0.0, 0.0, false, false};
        } else if (node.kind == Node::Kind::Gap) {
          control = {node.initial - previous, lowerMass - previous,
                     node.upper - previous, false};
        } else {
          name += "_fraction";
          const double span = node.upper - lowerMass;
          control = {span > 0 ? (node.initial - lowerMass) / span : 0.5, 0.0, 1.0, false};
        }
      }
      if (std::any_of(fitCoordinates.begin(), fitCoordinates.end(),
                      [&](const auto& c) { return c.name == name; }))
        throw std::invalid_argument("Duplicate ordered coordinate: " + name);
      node.fit = static_cast<int>(fitCoordinates.size());
      fitCoordinates.push_back({name, control});
    }
  }

  int DirectFitIndex(int physical) const { return direct_.at(physical); }

  PhysicalPoint Evaluate(const double* fit) const {
    PhysicalPoint result;
    result.values.resize(physicalCoordinates.size());
    result.jacobian.resize(physicalCoordinates.size());
    for (size_t i = 0; i < direct_.size(); ++i) {
      if (direct_[i] >= 0) {
        result.values[i] = logarithmic_[i] ? std::exp(fit[direct_[i]]) : fit[direct_[i]];
        result.jacobian[i].push_back({direct_[i], logarithmic_[i] ? result.values[i] : 1.0});
      }
    }
    for (const auto& node : nodes) {
      auto& value = result.values[node.physical];
      auto& row = result.jacobian[node.physical];
      switch (node.kind) {
      case Node::Kind::Fixed:
        value = node.initial;
        break;
      case Node::Kind::LogBase: {
        const double displacement = std::exp(fit[node.fit]);
        value = node.lower + displacement;
        row.push_back({node.fit, displacement});
        break;
      }
      case Node::Kind::LogGap: {
        const double previousLower = result.values[node.previous] + minimumGap;
        const bool activePrevious = previousLower >= node.lower;
        const double displacement = std::exp(fit[node.fit]);
        value = (activePrevious ? previousLower : node.lower) + displacement;
        if (activePrevious) row = result.jacobian[node.previous];
        row.push_back({node.fit, displacement});
        break;
      }
      case Node::Kind::Base:
        value = fit[node.fit];
        row.push_back({node.fit, 1.0});
        break;
      case Node::Kind::Gap:
        value = result.values[node.previous] + fit[node.fit];
        row = result.jacobian[node.previous];
        row.push_back({node.fit, 1.0});
        break;
      case Node::Kind::Fraction: {
        const double previousLower = result.values[node.previous] + minimumGap;
        const bool activePrevious = previousLower >= node.lower;
        const double lower = activePrevious ? previousLower : node.lower;
        const double fraction = fit[node.fit];
        value = lower + fraction * (node.upper - lower);
        if (activePrevious)
          for (const auto& entry : result.jacobian[node.previous])
            row.push_back({entry.first, (1.0 - fraction) * entry.second});
        row.push_back({node.fit, node.upper - lower});
        break;
      }
      }
    }
    return result;
  }

  // Reject numerical overflow/underflow of free log coordinates before a
  // line-shape helper is called. There are no finite optimization windows.
  bool Valid(const double* fit) const { return Valid(Evaluate(fit)); }

  bool Valid(const PhysicalPoint& point) const {
    for (size_t i = 0; i < point.values.size(); ++i)
      if (!std::isfinite(point.values[i]) || (logarithmic_[i] && point.values[i] <= 0.0))
        return false;
    for (const auto& node : nodes)
      if (point.values[node.physical] <= 0.0 ||
          (node.previous >= 0 && point.values[node.physical] <= point.values[node.previous]))
        return false;
    return true;
  }

  // Full row-major covariance; sparse Jacobian multiplication retains all
  // mass/width/channel/production cross-correlations.
  static std::vector<double> Propagate(const PhysicalPoint& point,
                                       const std::vector<double>& covariance,
                                       size_t fitSize) {
    const size_t n = point.values.size();
    std::vector<double> result(n * n, 0.0);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        for (const auto& a : point.jacobian[i])
          for (const auto& b : point.jacobian[j])
            result[i * n + j] += a.second * covariance[a.first * fitSize + b.first] * b.second;
    return result;
  }

private:
  std::vector<int> direct_;
  std::vector<bool> logarithmic_;
};

} // namespace emi::MassModels
