#include "Detail.h"

#include "Math/SpecFuncMathMore.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace emi::detail {

namespace {

double ClebschGordan(int l1, int l2, int l3, int m1, int m2, int m3) {
  const double sign = ((l1 - l2 + m3) & 1) ? -1.0 : 1.0;
  return sign * std::sqrt(2.0 * l3 + 1.0) *
         ROOT::Math::wigner_3j(2 * l1, 2 * l2, 2 * l3,
                               2 * m1, 2 * m2, -2 * m3);
}

static std::string MString(int m) { return (m < 0) ? "m"+std::to_string(-m) : std::to_string(m); }
static std::string MagName(char refl, char orient, int l, int m) { return std::string(1, refl) + '_' + orient + '_' + std::to_string(l) + '_' + MString(m); }
static std::string PhiName(char refl, char orient, int l, int m) { return std::string(1, refl) + "phi_" + orient + '_' + std::to_string(l) + '_' + MString(m); }

} // namespace

ParameterLabel ParseParameterLabel(const std::string& name) {
  ParameterLabel out;
  if (name.size() < 5) return out;
  out.reflectivity = name[0];
  out.phase = (name.find("phi_") != std::string::npos);
  out.orientation = (name.find("_L_") != std::string::npos || name.find("phi_L_") != std::string::npos) ? 'L' : 'T';
  const size_t last = name.find_last_of('_');
  if (last == std::string::npos || last + 1 >= name.size()) return out;
  const size_t prev = name.find_last_of('_', last - 1);
  if (prev == std::string::npos || prev + 1 >= last) return out;
  out.l = std::stoi(name.substr(prev + 1, last - prev - 1));
  const std::string mstr = name.substr(last + 1);
  out.m = (!mstr.empty() && mstr[0] == 'm')
          ? -std::stoi(mstr.substr(1))
          : std::stoi(mstr);
  out.valid = true;
  return out;
}

long long MakeParameterKey(char refl, char orient, int l, int m, bool isPhase) {
  const long long reflBit = (refl == 'b');
  const long long orientBit = (orient == 'L');
  const long long phaseBit = isPhase;
  const long long mEnc = static_cast<long long>(m + 32);
  return ( ((((phaseBit << 1) | reflBit) << 1) | orientBit) << 12) | (static_cast<long long>(l) << 6) | mEnc;
}

static int LongitudinalParitySign(int refl, int absM) {
  const int mParity = (absM & 1) ? -1 : +1;
  return refl * mParity;
}

std::vector<Parameter> BuildParameters(const InternalConfig& cfg) {
  if (cfg.waves.empty()) throw std::runtime_error("At least one (l,m) wave must be selected");
  if (!cfg.usePositiveReflectivity && !cfg.useNegativeReflectivity) {
    throw std::runtime_error("At least one reflectivity must be selected");
  }

  std::vector<Parameter> pars;
  pars.reserve(cfg.waves.size() * 8);
  std::unordered_set<std::string> names;

  auto add = [&](char refl, char orient, int l, int m, bool isPhase) {
    if (orient == 'L' && cfg.enforceLongitudinalParity) m = std::abs(m);

    Parameter p;
    p.name = isPhase ? PhiName(refl, orient, l, m) : MagName(refl, orient, l, m);
    if (!names.insert(p.name).second) return;
    p.init = 0.0;
    p.step = isPhase ? 0.6 : 0.2;
    if (cfg.photoproduction)
    {
      p.low = 0.0; // set to 0 for photoproduction due to ambiguity i.e. one appears in both sides
    } else
    {
      p.low = isPhase ? -kPi : 0.0;
    }
    p.high = isPhase ? kPi : cfg.magnitudeMax;
    p.phase = isPhase;

    // For longitudinal m=0, negative reflectivity is odd under the parity relation
    // A_m = eps (-1)^m A_-m and therefore vanishes.
    if (cfg.enforceLongitudinalParity && orient == 'L' && refl == 'b' && m == 0) {
      p.init = 0.0;
      p.fixed = true;
      p.low = 0.0;
      p.high = 0.0;
      p.step = 0.0;
    }

    pars.push_back(std::move(p));
  };

  for (const auto& wave : cfg.waves) {
    if (wave.l < 0 || std::abs(wave.m) > wave.l) {
      throw std::runtime_error("Invalid selected wave (l,m) = (" +
                               std::to_string(wave.l) + "," +
                               std::to_string(wave.m) + ")");
    }
    for (char reflectivity : {'a', 'b'}) {
      if (reflectivity == 'a' && !cfg.usePositiveReflectivity) continue;
      if (reflectivity == 'b' && !cfg.useNegativeReflectivity) continue;
      add(reflectivity, 'T', wave.l, wave.m, false);
      add(reflectivity, 'L', wave.l, wave.m, false);
      add(reflectivity, 'T', wave.l, wave.m, true);
      add(reflectivity, 'L', wave.l, wave.m, true);
    }
  }

  auto fixReferencePhase = [&](char reflectivity) {
    Parameter* best = nullptr;
    int bestL = -1;
    int bestM = -999;
    for (auto& parameter : pars) {
      const auto label = ParseParameterLabel(parameter.name);
      if (!parameter.phase || parameter.fixed || !label.valid ||
          label.reflectivity != reflectivity || label.orientation != 'T') continue;
      if (label.l > bestL || (label.l == bestL && label.m > bestM)) {
        best = &parameter;
        bestL = label.l;
        bestM = label.m;
      }
    }
    if (!best) return;
    best->init = 0.0;
    best->fixed = true;
    best->low = 0.0;
    best->high = 0.0;
    best->step = 0.0;
  };

  if (cfg.usePositiveReflectivity) fixReferencePhase('a');
  if (cfg.useNegativeReflectivity) fixReferencePhase('b');

  if (cfg.photoproduction) {
    for (auto& p : pars) {
      const auto label = ParseParameterLabel(p.name);
      if (!label.valid || label.orientation != 'L') continue;
      p.init = 0.0;
      p.fixed = true;
      p.low = 0.0;
      p.high = 0.0;
      p.step = 0.0;
    }
  }

  return pars;
}

struct BruSelection {
  double coeff = 0.0;
  char reflectivity1 = 'a';
  char orientation1 = 'T';
  int l1 = 0;
  int m1 = 0;
  char reflectivity2 = 'a';
  char orientation2 = 'T';
  int l2 = 0;
  int m2 = 0;
  TrigKind trig = TrigKind::kCos;
};

static bool ResolveBruSelection(int reflsign, double factor,
                                int l, int m, int lpr, int mpr,
                                int alpha,
                                bool orientSwap,
                                BruSelection& out)
{
  out.reflectivity1 = (reflsign == -1) ? 'b' : 'a';
  out.reflectivity2 = (reflsign == -1) ? 'b' : 'a';

  out.orientation1 = 'T';
  out.orientation2 = 'T';

  if (alpha == 4) {
    out.orientation1 = 'L';
    out.orientation2 = 'L';
  }

  if (alpha >= 5) {
    out.orientation1 = 'L';
    out.orientation2 = 'T';
  }

  if (alpha >= 5 && orientSwap) {
    out.orientation1 = 'T';
    out.orientation2 = 'L';
  }

  // There are no independent negative-m longitudinal amplitudes.
    if(out.orientation1=='L'&&m<0){
      factor*=LongitudinalParitySign(reflsign,m);
      m=-m;
    }
    if(out.orientation2=='L'&&mpr<0){
      factor*=LongitudinalParitySign(reflsign,mpr);
      mpr=-mpr;
    }

  // rho^1 and rho^2 carry an overall reflectivity factor.

  if (alpha==1 || alpha==2)
  {
    out.coeff = reflsign * factor;
  }else
  {
    out.coeff = factor;
  }

  out.l1 = l;
  out.m1 = m;
  out.l2 = lpr;
  out.m2 = mpr;

  out.trig =
      (alpha == 3 || alpha == 7 || alpha == 8)
          ? TrigKind::kSin
          : TrigKind::kCos;

  return out.coeff != 0.0;
}



std::vector<MomentModel> BuildMomentModels(const InternalConfig& cfg,
                                                  const std::unordered_map<long long, int>& paramIndex,
                                                  const std::unordered_set<std::string>& neededMoments,
                                                  std::vector<PhasePair>& phasePairs) {
  std::vector<MomentModel> models;
  const int alphaMax = cfg.photoproduction ? 3 : 8;
  int maximumL = 0;
  for (const auto& wave : cfg.waves) maximumL = std::max(maximumL, wave.l);
  models.reserve((alphaMax + 1) * (2 * maximumL + 1) * (2 * maximumL + 2) / 2);

  std::unordered_map<long long, int> phasePairLookup;
  const auto& waves = cfg.waves;

  struct CGKey {
    int lpr, L, l, mpr, M, m;
    bool operator==(const CGKey& o) const { return lpr == o.lpr && L == o.L && l == o.l && mpr == o.mpr && M == o.M && m == o.m; }
  };
  struct CGKeyHash {
    size_t operator()(const CGKey& k) const {
      size_t h = 1469598103934665603ull;
      auto mix = [&](int v) { h ^= static_cast<size_t>(v + 32); h *= 1099511628211ull; };
      mix(k.lpr); mix(k.L); mix(k.l); mix(k.mpr); mix(k.M); mix(k.m);
      return h;
    }
  };
  std::unordered_map<CGKey, double, CGKeyHash> cgCache;
  cgCache.reserve(2048);

  auto getCG = [&](int lpr, int L, int l, int mpr, int M, int m) {
    const CGKey key{lpr, L, l, mpr, M, m};
    auto it = cgCache.find(key);
    if (it != cgCache.end()) return it->second;
    const double v = ClebschGordan(lpr, L, l, mpr, M, m);
    cgCache.emplace(key, v);
    return v;
  };

  auto getPhasePairIdx = [&](int idxPhi1, int idxPhi2) {
    const long long key = (static_cast<long long>(idxPhi1) << 32) | static_cast<unsigned int>(idxPhi2);
    auto it = phasePairLookup.find(key);
    if (it != phasePairLookup.end()) return it->second;
    const int idx = static_cast<int>(phasePairs.size());
    phasePairs.push_back({idxPhi1, idxPhi2});
    phasePairLookup.emplace(key, idx);
    return idx;
  };

  auto paramIdx = [&](char refl, char orient, int l, int m, bool isPhase) -> int {
    const auto it = paramIndex.find(MakeParameterKey(refl, orient, l, m, isPhase));
    return it == paramIndex.end() ? -1 : it->second;
  };

  auto emit = [&](MomentModel& mm, int reflsign, double factor, int l, int m,
                  int lpr, int mpr, int alpha, bool orientSwap) {
    if (reflsign > 0 && !cfg.usePositiveReflectivity) return;
    if (reflsign < 0 && !cfg.useNegativeReflectivity) return;
    BruSelection sel;
    if (!ResolveBruSelection(reflsign, factor, l, m, lpr, mpr, alpha,
                             orientSwap, sel)) return;

    Term t;
    t.coeff = sel.coeff;
    t.idxMag1 = paramIdx(sel.reflectivity1, sel.orientation1, sel.l1, sel.m1, false);
    t.idxMag2 = paramIdx(sel.reflectivity2, sel.orientation2, sel.l2, sel.m2, false);
    const int idxPhi1 = paramIdx(sel.reflectivity1, sel.orientation1, sel.l1, sel.m1, true);
    const int idxPhi2 = paramIdx(sel.reflectivity2, sel.orientation2, sel.l2, sel.m2, true);
    if (t.idxMag1 < 0 || t.idxMag2 < 0 || idxPhi1 < 0 || idxPhi2 < 0) return;
    t.phasePairIdx = getPhasePairIdx(idxPhi1, idxPhi2);
    t.trig = sel.trig;
    t.ignorePhase = (sel.orientation1 == sel.orientation2 && sel.l1 == sel.l2 && sel.m1 == sel.m2);
    mm.terms.push_back(t);
  };

  for (int alpha = 0; alpha <= alphaMax; ++alpha) {
    for (int L = 0; L <= 2 * maximumL; ++L) {
      for (int M = 0; M <= L; ++M)
      {
        MomentModel mm;
        mm.alpha = alpha;
        mm.L = L;
        mm.M = M;
        mm.name =
            std::string("H_")
            + std::to_string(alpha) + "_"
            + std::to_string(L) + "_"
            + std::to_string(M);

        if (!neededMoments.empty() && neededMoments.find(mm.name) == neededMoments.end()) continue;

        if ((alpha == 2 || alpha == 3 ||
             alpha == 6 || alpha == 7) &&
            M == 0) {
          continue;
        }

        for (const auto& w1 : waves)
        {
          const int il = w1.l;
          const int im = w1.m;

          for (const auto& w2 : waves)
          {
            const int ilpr = w2.l;
            const int impr = w2.m;

            const double CM =
                getCG(ilpr, L, il, impr, M, im);
            if (CM == 0.0) continue;

            const double C0 =
                getCG(ilpr, L, il, 0, 0, 0);
            if (C0 == 0.0) continue;

            double ccfactor =
                CM * C0
                * std::sqrt(
                      (2.0 * ilpr + 1.0)
                      / (2.0 * il + 1.0));

            if (ccfactor == 0.0) continue;

            const int mmprimesign =
                ((std::abs(im - impr) & 1) ? -1 : 1);

            const int mprimesign =
                ((std::abs(impr) & 1) ? -1 : 1);

            const int msign =
                ((std::abs(im) & 1) ? -1 : 1);

            if (alpha == 0) {
              emit(mm, +1, ccfactor,
                   il, im, ilpr, impr, 0, false);

              emit(mm, +1, mmprimesign * ccfactor,
                   il, -im, ilpr, -impr, 0, false);

              if (cfg.useNegativeReflectivity) {
                emit(mm, -1, ccfactor,
                     il, im, ilpr, impr, 0, false);

                emit(mm, -1, mmprimesign * ccfactor,
                     il, -im, ilpr, -impr, 0, false);
              }
            }

            // alpha = 1 and 2 negatives cancel (-ve from moment definition)
            else if (alpha == 1) {
              emit(mm, +1, msign * ccfactor,
                   il, -im, ilpr, impr, 1, false);

              emit(mm, +1, mprimesign * ccfactor,
                   il, im, ilpr, -impr, 1, false);

              if (cfg.useNegativeReflectivity) {
                emit(mm, -1, msign * ccfactor,
                     il, -im, ilpr, impr, 1, false);

                emit(mm, -1, mprimesign * ccfactor,
                     il, im, ilpr, -impr, 1, false);
              }
            }

            else if (alpha == 2) {
              emit(mm, +1, msign * ccfactor,
                   il, -im, ilpr, impr, 2, false);

              emit(mm, +1, -mprimesign * ccfactor,
                   il, im, ilpr, -impr, 2, false);

              if (cfg.useNegativeReflectivity) {
                emit(mm, -1, msign * ccfactor,
                     il, -im, ilpr, impr, 2, false);

                emit(mm, -1, -mprimesign * ccfactor,
                     il, im, ilpr, -impr, 2, false);
              }
            }

            // -ve factor from moment def same for rest other than 4
            else if (alpha == 3) {
              ccfactor *= -1;
              emit(mm, +1, ccfactor,
                   il, im, ilpr, impr, 3, false);

              emit(mm, +1, -mmprimesign * ccfactor,
                   il, -im, ilpr, -impr, 3, false);

              if (cfg.useNegativeReflectivity) {
                emit(mm, -1, ccfactor,
                     il, im, ilpr, impr, 3, false);

                emit(mm, -1, -mmprimesign * ccfactor,
                     il, -im, ilpr, -impr, 3, false);
              }
            }

            else if (alpha == 4) {
              const double f = 2.0 * ccfactor;

              emit(mm, +1, f,
                   il, im, ilpr, impr, 4, false);

              if (cfg.useNegativeReflectivity) {
                emit(mm, -1, f,
                     il, im, ilpr, impr, 4, false);
              }
            }

            else if (alpha == 5) {
              const double f =
                  -ccfactor / std::sqrt(2.0);

              int refl = +1;

              emit(mm, refl, f,
                   il, im, ilpr, impr, 5, false);

              emit(mm, refl, f,
                   il, im, ilpr, impr, 5, true);

              emit(mm, refl, mmprimesign * f,
                   il, -im, ilpr, -impr, 5, false);

              emit(mm, refl, mmprimesign * f,
                   il, -im, ilpr, -impr, 5, true);

              if (cfg.useNegativeReflectivity) {
                refl = -1;

                emit(mm, refl, f,
                   il, im, ilpr, impr, 5, false);

                emit(mm, refl, f,
                     il, im, ilpr, impr, 5, true);

                emit(mm, refl, mmprimesign * f,
                     il, -im, ilpr, -impr, 5, false);

                emit(mm, refl, mmprimesign * f,
                     il, -im, ilpr, -impr, 5, true);
              }
            }

            else if (alpha == 6) {
              const double f =
                  -ccfactor / std::sqrt(2.0);

              int refl = +1;

              emit(mm, refl, f,
                   il, im, ilpr, impr, 6, false);

              emit(mm, refl, -f,
                   il, im, ilpr, impr, 6, true);

              emit(mm, refl, -mmprimesign * f,
                   il, -im, ilpr, -impr, 6, false);

              emit(mm, refl, mmprimesign * f,
                   il, -im, ilpr, -impr, 6, true);

              if (cfg.useNegativeReflectivity) {
                refl = -1;

                emit(mm, refl, f,
                   il, im, ilpr, impr, 6, false);

                emit(mm, refl, -f,
                     il, im, ilpr, impr, 6, true);

                emit(mm, refl, -mmprimesign * f,
                     il, -im, ilpr, -impr, 6, false);

                emit(mm, refl, mmprimesign * f,
                     il, -im, ilpr, -impr, 6, true);
              }
            }

            else if (alpha == 7) {
              const double f =
                  -ccfactor / std::sqrt(2.0);

              int refl = +1;

              emit(mm, refl, f,
                   il, im, ilpr, impr, 7, false);

              emit(mm, refl, f,
                   il, im, ilpr, impr, 7, true);

              emit(mm, refl, -mmprimesign * f,
                   il, -im, ilpr, -impr, 7, false);

              emit(mm, refl, -mmprimesign * f,
                   il, -im, ilpr, -impr, 7, true);

              if (cfg.useNegativeReflectivity) {
                refl = -1;

                emit(mm, refl, f,
                   il, im, ilpr, impr, 7, false);

                emit(mm, refl, f,
                     il, im, ilpr, impr, 7, true);

                emit(mm, refl, -mmprimesign * f,
                     il, -im, ilpr, -impr, 7, false);

                emit(mm, refl, -mmprimesign * f,
                     il, -im, ilpr, -impr, 7, true);
              }
            }

            // no negative as it cancels with i^2 factor
            else if (alpha == 8)
            {
              const double f =
                  ccfactor / std::sqrt(2.0);

              int refl = +1;

              emit(mm, refl, f,
                   il, im, ilpr, impr, 8, false);

              emit(mm, refl, -f,
                   il, im, ilpr, impr, 8, true);

              emit(mm, refl, mmprimesign * f,
                   il, -im, ilpr, -impr, 8, false);

              emit(mm, refl, -mmprimesign * f,
                   il, -im, ilpr, -impr, 8, true);

              if (cfg.useNegativeReflectivity) {
                refl = -1;

                emit(mm, refl, f,
                   il, im, ilpr, impr, 8, false);

                emit(mm, refl, -f,
                     il, im, ilpr, impr, 8, true);

                emit(mm, refl, mmprimesign * f,
                     il, -im, ilpr, -impr, 8, false);

                emit(mm, refl, -mmprimesign * f,
                     il, -im, ilpr, -impr, 8, true);
              }
            }
          }
        }

        if (!mm.terms.empty()) models.push_back(std::move(mm));
      }
    }
  }
  return models;
}
} // namespace emi::detail
