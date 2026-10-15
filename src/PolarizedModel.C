// Build polarized moment tensors from complex partial-wave amplitudes.  The
// explicit reflectivity, photon-helicity k, and target-orientation labels make
// the interference terms auditable against the angular-momentum formulae.
// Initial-, recoil-, and double-polarization observables expose different
// combinations of k=+1 and k=-1.  This changes which continuous basis rotations
// are observable and which gauge convention makes amplitude coordinates unique.
#include "Detail.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <tuple>

namespace emi::detail {

namespace {

struct AmplitudeComponent {
  // One physical helicity amplitude becomes a two-term reflectivity sum.  The
  // coefficient stores signs introduced by the basis transformation and parity.
  char reflectivity = 'a';
  char orientation = 'T';
  int k = 1;
  int l = 0;
  int m = 0;
  double coefficient = 1.0;
};

double PhotonResponse(int alpha, int row, int column, bool& carriesI) {
  // Non-zero entries of Appendix D in photon-helicity order (+1, 0, -1).
  // carriesI means that the returned real coefficient multiplies i.
  // Separating that factor lets the final real moment choose sine or cosine
  // without complex arithmetic in the minimizer hot path.
  carriesI = false;
  switch (alpha) {
  case 0:
    return (row == column && row != 1) ? 1.0 : 0.0;
  case 1:
    return ((row == 0 && column == 2) || (row == 2 && column == 0)) ? 1.0 : 0.0;
  case 2:
    carriesI = true;
    if (row == 0 && column == 2)
      return -1.0;
    if (row == 2 && column == 0)
      return 1.0;
    return 0.0;
  case 3:
    if (row == column && row == 0)
      return 1.0;
    if (row == column && row == 2)
      return -1.0;
    return 0.0;
  case 4:
    return (row == 1 && column == 1) ? 2.0 : 0.0;
  case 5:
    if ((row == 0 && column == 1) || (row == 1 && column == 0))
      return kInvSqrt2;
    if ((row == 1 && column == 2) || (row == 2 && column == 1))
      return -kInvSqrt2;
    return 0.0;
  case 6:
    carriesI = true;
    if ((row == 0 && column == 1) || (row == 2 && column == 1))
      return -kInvSqrt2;
    if ((row == 1 && column == 0) || (row == 1 && column == 2))
      return kInvSqrt2;
    return 0.0;
  case 7:
    return std::abs(row - column) == 1 ? kInvSqrt2 : 0.0;
  case 8:
    carriesI = true;
    if ((row == 0 && column == 1) || (row == 1 && column == 2))
      return -kInvSqrt2;
    if ((row == 1 && column == 0) || (row == 2 && column == 1))
      return kInvSqrt2;
    return 0.0;
  default:
    return 0.0;
  }
}

int PolarizationParity(int alpha, int beta, int delta) {
  // Combining the photon, target, and recoil transformations determines
  // whether this tensor component is even or odd under reflection in the
  // production plane. Nucleon spin is an axial vector: its in-plane x and z
  // components (Pauli indices 1 and 3) change sign, while the normal y
  // component (index 2) does not. Treating indices 2 and 3 as odd would remove
  // the physical interference between the k=+1 and k=-1 amplitudes.
  const int photon = (alpha == 2 || alpha == 3 || alpha == 6 || alpha == 7) ? -1 : 1;
  const int target = (beta == 1 || beta == 3) ? -1 : 1;
  const int recoil = (delta == 1 || delta == 3) ? -1 : 1;
  return photon * target * recoil;
}

struct SpinReplacement {
  int targetHelicity = 1;
  int recoilHelicity = 1;
  double coefficient = 1.0;
  bool carriesI = false;
};

struct SpinReplacementRule {
  int targetSign = 1;
  int kSign = 1;
  int coefficientSign = 1;
  bool coefficientHasTarget = false;
  bool coefficientHasK = false;
  bool carriesI = false;
};

// Rows are beta and columns are delta. Each rule is the corresponding entry
// C_{beta delta} from the reduced-k target--recoil replacement table.
constexpr std::array<SpinReplacementRule, 16> kSpinReplacementRules{{
    {+1, +1, +1, false, false, false}, // (0,0):       1
    {+1, -1, -1, false, false, false}, // (0,1):      -1
    {+1, -1, -1, false, true, true},   // (0,2):    -i k
    {+1, +1, +1, false, true, false},  // (0,3):     k
    {-1, -1, -1, true, true, false},   // (1,0):  -e k
    {-1, +1, +1, true, false, false},  // (1,1):   e
    {-1, +1, +1, true, false, true},   // (1,2): i e
    {-1, -1, -1, true, false, false},  // (1,3):  -e
    {-1, -1, +1, true, true, true},    // (2,0): i e k
    {-1, +1, -1, true, true, true},    // (2,1): -i e k
    {-1, +1, +1, true, false, false},  // (2,2):   e
    {-1, -1, +1, true, false, true},   // (2,3): i e
    {+1, +1, +1, false, false, false}, // (3,0):       1
    {+1, -1, -1, false, false, false}, // (3,1):      -1
    {+1, -1, -1, false, true, true},   // (3,2):    -i k
    {+1, +1, +1, false, true, false},  // (3,3):     k
}};

SpinReplacement ReplaceTargetRecoil(int beta, int delta, int targetHelicity, int k) {
  // Target--recoil replacement table in the reduced-k basis. Here epsilon is
  // targetHelicity and k is +1 for non-flip or -1 for flip. The table fixes
  // epsilon', k', and C_{beta delta}; keeping it explicit prevents an index
  // transpose or Pauli-matrix sign convention from being hidden in the loops.
  const auto& rule = kSpinReplacementRules[4 * beta + delta];
  const int targetPrime = rule.targetSign * targetHelicity;
  const int kPrime = rule.kSign * k;

  SpinReplacement replacement;
  replacement.targetHelicity = targetPrime;
  replacement.recoilHelicity = targetPrime * kPrime;
  replacement.coefficient = static_cast<double>(rule.coefficientSign);
  if (rule.coefficientHasTarget)
    replacement.coefficient *= targetHelicity;
  if (rule.coefficientHasK)
    replacement.coefficient *= k;
  replacement.carriesI = rule.carriesI;
  return replacement;
}

std::array<AmplitudeComponent, 2> PhysicalAmplitude(const InternalConfig& cfg,
                                                    int photonHelicity, int l, int m,
                                                    int targetHelicity,
                                                    int recoilHelicity) {
  // Section III reflectivity inversion followed by the nucleon parity relation.
  const char orientation = photonHelicity == 0 ? 'L' : 'T';
  const int reducedM = photonHelicity == -1 ? -m : m;
  const int k = targetHelicity == recoilHelicity ? 1 : -1;
  const double mParity = (std::abs(m) & 1) ? -1.0 : 1.0;

  std::array<AmplitudeComponent, 2> components;
  for (size_t i = 0; i < components.size(); ++i) {
    const int reflectivity = i == 0 ? 1 : -1;
    auto& component = components[i];
    component.reflectivity = reflectivity > 0 ? 'a' : 'b';
    component.orientation = orientation;
    component.k = k;
    component.l = l;
    component.m = reducedM;

    if (photonHelicity == -1) {
      component.coefficient *= -reflectivity * mParity;
    }
    if (targetHelicity < 0) {
      component.coefficient *= recoilHelicity < 0 ? reflectivity : -reflectivity;
    }
    if (orientation == 'L' && cfg.enforceLongitudinalParity && reducedM < 0) {
      component.coefficient *= reflectivity * mParity;
      component.m = -reducedM;
    }
  }
  return components;
}

} // namespace

std::vector<MomentModel>
BuildPolarizedMomentModels(const InternalConfig& cfg,
                           const std::unordered_map<long long, int>& paramIndex,
                           const std::unordered_set<std::string>& neededMoments,
                           std::vector<PhasePair>& phasePairs) {
  // alpha labels photon response, beta target polarization, and delta recoil
  // polarization.  Their allowed ranges are reduced to zero for an unmeasured
  // spin, avoiding identically absent tensor components.
  std::vector<MomentModel> models;
  const int alphaMax = cfg.photoproduction ? 3 : 8;
  const int betaMax = cfg.nucleonPolarization == NucleonPolarization::Recoil ? 0 : 3;
  const int deltaMax = cfg.nucleonPolarization == NucleonPolarization::Initial ? 0 : 3;
  int maximumL = 0;
  for (const auto& wave : cfg.waves)
    maximumL = std::max(maximumL, wave.l);
  models.reserve((alphaMax + 1) * (betaMax + 1) * (deltaMax + 1) * (2 * maximumL + 1) *
                 (maximumL + 1));

  std::unordered_map<long long, int> phasePairLookup;
  // Polarized tensors contain many repeated phase differences. Share each
  // pair so its sine and cosine are evaluated only once per objective call.
  auto parameterIndex = [&](const AmplitudeComponent& component, bool phase) {
    const auto found =
        paramIndex.find(MakeParameterKey(component.reflectivity, component.orientation,
                                         component.k, component.l, component.m, phase));
    return found == paramIndex.end() ? -1 : found->second;
  };
  auto phasePairIndex = [&](int first, int second) {
    const long long key =
        (static_cast<long long>(first) << 32) | static_cast<unsigned int>(second);
    const auto found = phasePairLookup.find(key);
    if (found != phasePairLookup.end())
      return found->second;
    const int index = static_cast<int>(phasePairs.size());
    phasePairs.push_back({first, second});
    phasePairLookup.emplace(key, index);
    return index;
  };

  for (int alpha = 0; alpha <= alphaMax; ++alpha) {
    for (int beta = 0; beta <= betaMax; ++beta) {
      for (int delta = 0; delta <= deltaMax; ++delta) {
        const int parity = PolarizationParity(alpha, beta, delta);
        const double responseSign =
            beta == 0 && delta == 0 && (alpha == 0 || alpha == 4) ? 1.0 : -1.0;
        for (int L = 0; L <= 2 * maximumL; ++L) {
          for (int M = 0; M <= L; ++M) {
            MomentModel model;
            model.alpha = alpha;
            model.beta = beta;
            model.delta = delta;
            model.L = L;
            model.M = M;
            model.name = MakeMomentName(cfg, alpha, beta, delta, L, M);
            if ((!neededMoments.empty() && !neededMoments.count(model.name)) ||
                (parity < 0 && M == 0))
              continue;

            auto addTerm = [&](const AmplitudeComponent& first,
                               const AmplitudeComponent& second, double coefficient,
                               bool coefficientCarriesI) {
              int mag1 = parameterIndex(first, false);
              int mag2 = parameterIndex(second, false);
              int phi1 = parameterIndex(first, true);
              int phi2 = parameterIndex(second, true);
              if (mag1 < 0 || mag2 < 0 || phi1 < 0 || phi2 < 0)
                return;

              auto emit = [&](double value, TrigKind trig) {
                if (std::abs(value) < 1e-14)
                  return;
                int termMag1 = mag1;
                int termMag2 = mag2;
                int termPhi1 = phi1;
                int termPhi2 = phi2;
                if (trig == TrigKind::kSin && termPhi1 == termPhi2)
                  return;
                if (termPhi1 > termPhi2) {
                  std::swap(termMag1, termMag2);
                  std::swap(termPhi1, termPhi2);
                  if (trig == TrigKind::kSin)
                    value = -value;
                }
                model.terms.push_back({value, termMag1, termMag2,
                                       phasePairIndex(termPhi1, termPhi2), trig,
                                       termPhi1 == termPhi2});
              };

              // The parity relation chooses the real or imaginary part of the
              // amplitude bilinear. All factors of i have already been reduced
              // into the real coefficient using i^2 = -1.
              if (parity > 0) {
                emit(coefficientCarriesI ? -coefficient : coefficient,
                     coefficientCarriesI ? TrigKind::kSin : TrigKind::kCos);
              } else {
                emit(coefficient,
                     coefficientCarriesI ? TrigKind::kCos : TrigKind::kSin);
              }
            };

            for (const auto& firstWave : cfg.waves) {
              for (const auto& secondWave : cfg.waves) {
                // Each observable is a weighted sum of bilinears A_i A_j*.
                // The angular factor couples the two partial waves to (L,M);
                // the inner loops then trace over photon and nucleon spins.
                const double angularFactor =
                    ClebschGordan(secondWave.l, L, firstWave.l, 0, 0, 0) *
                    ClebschGordan(secondWave.l, L, firstWave.l, secondWave.m, M,
                                  firstWave.m) *
                    std::sqrt((2.0 * secondWave.l + 1.0) / (2.0 * firstWave.l + 1.0));
                if (angularFactor == 0.0)
                  continue;

                constexpr std::array<int, 3> photonHelicity = {1, 0, -1};
                for (int photon1 = 0; photon1 < 3; ++photon1) {
                  for (int photon2 = 0; photon2 < 3; ++photon2) {
                    bool photonCarriesI = false;
                    const double photonWeight =
                        PhotonResponse(alpha, photon1, photon2, photonCarriesI);
                    if (photonWeight == 0.0)
                      continue;
                    for (int target = -1; target <= 1; target += 2) {
                      for (int recoil = -1; recoil <= 1; recoil += 2) {
                        const int k = target == recoil ? 1 : -1;
                        const auto replacement =
                            ReplaceTargetRecoil(beta, delta, target, k);
                        const auto first =
                            PhysicalAmplitude(cfg, photonHelicity[photon1], firstWave.l,
                                              firstWave.m, target, recoil);
                        const auto second = PhysicalAmplitude(
                            cfg, photonHelicity[photon2], secondWave.l, secondWave.m,
                            replacement.targetHelicity, replacement.recoilHelicity);
                        double coefficient = 0.5 * responseSign * angularFactor *
                                             photonWeight * replacement.coefficient;
                        const int imaginaryFactors =
                            static_cast<int>(photonCarriesI) +
                            static_cast<int>(replacement.carriesI);
                        if (imaginaryFactors == 2)
                          coefficient = -coefficient;
                        const bool coefficientCarriesI = (imaginaryFactors & 1) != 0;
                        for (const auto& firstComponent : first) {
                          for (const auto& secondComponent : second) {
                            addTerm(firstComponent, secondComponent,
                                    coefficient * firstComponent.coefficient *
                                        secondComponent.coefficient,
                                    coefficientCarriesI);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }

            std::sort(model.terms.begin(), model.terms.end(),
                      [](const Term& lhs, const Term& rhs) {
                        return std::tie(lhs.idxMag1, lhs.idxMag2, lhs.phasePairIdx,
                                        lhs.trig) < std::tie(rhs.idxMag1, rhs.idxMag2,
                                                             rhs.phasePairIdx,
                                                             rhs.trig);
                      });
            size_t write = 0;
            // Several spin paths can produce the same algebraic bilinear.
            // Sorting and coalescing them reduces work in every chi-square call.
            for (const auto& term : model.terms) {
              if (write > 0) {
                auto& previous = model.terms[write - 1];
                if (previous.idxMag1 == term.idxMag1 &&
                    previous.idxMag2 == term.idxMag2 &&
                    previous.phasePairIdx == term.phasePairIdx &&
                    previous.trig == term.trig) {
                  previous.coeff += term.coeff;
                  continue;
                }
              }
              model.terms[write++] = term;
            }
            model.terms.resize(write);
            model.terms.erase(std::remove_if(model.terms.begin(), model.terms.end(),
                                             [](const Term& term) {
                                               return std::abs(term.coeff) < 1e-13;
                                             }),
                              model.terms.end());
            if (!model.terms.empty())
              models.push_back(std::move(model));
          }
        }
      }
    }
  }
  return models;
}

} // namespace emi::detail
