#include "emi/Runner.h"

#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TObjString.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <filesystem>

namespace emi {
namespace {


constexpr double kSqrt6Over5  = 0.48989794855663561963945681494118;
constexpr double kSqrt12Over5 = 0.69282032302755092063339055356909;

struct ValErr2 {
  double v{0.0};
  double stat{0.0};
  double syst{0.0};
};

inline ValErr2 VE(double v, double stat, double syst) {
  return ValErr2{v, stat, syst};
}

inline double Comb(const ValErr2& x) {
  return std::sqrt(x.stat * x.stat + x.syst * x.syst);
}

struct SDMEsTable {
  ValErr2 r00_04, re_r10_04, r1m1_04;
  ValErr2 r1m1_1, re_r10_1, im_r10_2, r00_1, im_r10_3, r00_8;
  ValErr2 r11_5, r1m1_5, im_r1m1_6, im_r1m1_7, r11_8, r1m1_8;
  ValErr2 r11_1, im_r1m1_3, im_r1m1_2, re_r10_5, im_r10_6, im_r10_7, re_r10_8, r00_5;
};

struct OutputArrays {
  static constexpr int kMaxBins = 64;
  int Nbins{0};
  double Q2[kMaxBins]{};
  double RH04_0_0[kMaxBins]{}, RH04_0_0_err[kMaxBins]{};
  double RH04_1_0[kMaxBins]{}, RH04_1_0_err[kMaxBins]{}, RH04_1_1[kMaxBins]{}, RH04_1_1_err[kMaxBins]{}, RH04_2_0[kMaxBins]{}, RH04_2_0_err[kMaxBins]{}, RH04_2_1[kMaxBins]{}, RH04_2_1_err[kMaxBins]{}, RH04_2_2[kMaxBins]{}, RH04_2_2_err[kMaxBins]{};
  double RH_1_0_0[kMaxBins]{}, RH_1_0_0_err[kMaxBins]{}, RH_1_1_0[kMaxBins]{}, RH_1_1_0_err[kMaxBins]{}, RH_1_1_1[kMaxBins]{}, RH_1_1_1_err[kMaxBins]{}, RH_1_2_0[kMaxBins]{}, RH_1_2_0_err[kMaxBins]{}, RH_1_2_1[kMaxBins]{}, RH_1_2_1_err[kMaxBins]{}, RH_1_2_2[kMaxBins]{}, RH_1_2_2_err[kMaxBins]{};
  double RH_2_1_0[kMaxBins]{}, RH_2_1_0_err[kMaxBins]{}, RH_2_1_1[kMaxBins]{}, RH_2_1_1_err[kMaxBins]{}, RH_2_2_1[kMaxBins]{}, RH_2_2_1_err[kMaxBins]{}, RH_2_2_2[kMaxBins]{}, RH_2_2_2_err[kMaxBins]{};
  double RH_3_1_0[kMaxBins]{}, RH_3_1_0_err[kMaxBins]{}, RH_3_1_1[kMaxBins]{}, RH_3_1_1_err[kMaxBins]{}, RH_3_2_1[kMaxBins]{}, RH_3_2_1_err[kMaxBins]{}, RH_3_2_2[kMaxBins]{}, RH_3_2_2_err[kMaxBins]{};
  double RH_5_0_0[kMaxBins]{}, RH_5_0_0_err[kMaxBins]{}, RH_5_1_0[kMaxBins]{}, RH_5_1_0_err[kMaxBins]{}, RH_5_1_1[kMaxBins]{}, RH_5_1_1_err[kMaxBins]{}, RH_5_2_0[kMaxBins]{}, RH_5_2_0_err[kMaxBins]{}, RH_5_2_1[kMaxBins]{}, RH_5_2_1_err[kMaxBins]{}, RH_5_2_2[kMaxBins]{}, RH_5_2_2_err[kMaxBins]{};
  double RH_6_1_0[kMaxBins]{}, RH_6_1_0_err[kMaxBins]{}, RH_6_1_1[kMaxBins]{}, RH_6_1_1_err[kMaxBins]{}, RH_6_2_1[kMaxBins]{}, RH_6_2_1_err[kMaxBins]{}, RH_6_2_2[kMaxBins]{}, RH_6_2_2_err[kMaxBins]{};
  double RH_7_1_0[kMaxBins]{}, RH_7_1_0_err[kMaxBins]{}, RH_7_1_1[kMaxBins]{}, RH_7_1_1_err[kMaxBins]{}, RH_7_2_1[kMaxBins]{}, RH_7_2_1_err[kMaxBins]{}, RH_7_2_2[kMaxBins]{}, RH_7_2_2_err[kMaxBins]{};
  double RH_8_0_0[kMaxBins]{}, RH_8_0_0_err[kMaxBins]{}, RH_8_1_0[kMaxBins]{}, RH_8_1_0_err[kMaxBins]{}, RH_8_1_1[kMaxBins]{}, RH_8_1_1_err[kMaxBins]{}, RH_8_2_0[kMaxBins]{}, RH_8_2_0_err[kMaxBins]{}, RH_8_2_1[kMaxBins]{}, RH_8_2_1_err[kMaxBins]{}, RH_8_2_2[kMaxBins]{}, RH_8_2_2_err[kMaxBins]{};
};

struct DatasetSpec {
  std::string key;
  std::string title;
  std::vector<double> q2;
  std::vector<SDMEsTable> bins;
};

static void ResizeOutputs(OutputArrays& out, std::size_t n) {
  if (n > static_cast<std::size_t>(OutputArrays::kMaxBins)) {
    throw std::runtime_error("Too many bins for fixed-size output arrays");
  }
  out.Nbins = static_cast<int>(n);
}

static void FillBin(OutputArrays& out, int i, double q2, const SDMEsTable& p) {
  out.Q2[i] = q2;

  const double s_r00_04    = Comb(p.r00_04);
  const double s_re_r10_04 = Comb(p.re_r10_04);
  const double s_r1m1_04   = Comb(p.r1m1_04);

  const double s_r1m1_1    = Comb(p.r1m1_1);
  const double s_re_r10_1  = Comb(p.re_r10_1);
  const double s_r00_1     = Comb(p.r00_1);
  const double s_r11_1     = Comb(p.r11_1);

  const double s_im_r10_2   = Comb(p.im_r10_2);
  const double s_im_r1m1_2  = Comb(p.im_r1m1_2);
  const double s_im_r10_3   = Comb(p.im_r10_3);
  const double s_im_r1m1_3  = Comb(p.im_r1m1_3);

  const double s_r00_5      = Comb(p.r00_5);
  const double s_r11_5      = Comb(p.r11_5);
  const double s_re_r10_5   = Comb(p.re_r10_5);
  const double s_r1m1_5     = Comb(p.r1m1_5);

  const double s_im_r10_6   = Comb(p.im_r10_6);
  const double s_im_r1m1_6  = Comb(p.im_r1m1_6);
  const double s_im_r10_7   = Comb(p.im_r10_7);
  const double s_im_r1m1_7  = Comb(p.im_r1m1_7);

  const double s_r00_8      = Comb(p.r00_8);
  const double s_r11_8      = Comb(p.r11_8);
  const double s_re_r10_8   = Comb(p.re_r10_8);
  const double s_r1m1_8     = Comb(p.r1m1_8);

  out.RH04_0_0[i] = 2;                                             out.RH04_0_0_err[i] = 0.001;
  out.RH04_2_0[i] = 2 * 0.2 * (3.0 * p.r00_04.v - 1.0);            out.RH04_2_0_err[i] = 2 * 0.6 * s_r00_04;
  out.RH04_2_1[i] = 2 * kSqrt12Over5 * p.re_r10_04.v;              out.RH04_2_1_err[i] = 2 * kSqrt12Over5 * s_re_r10_04;
  out.RH04_2_2[i] =  -2 * kSqrt6Over5 * p.r1m1_04.v;                out.RH04_2_2_err[i] = 2 * kSqrt6Over5 * s_r1m1_04;

  out.RH_1_0_0[i] =  -2 * (2.0 * p.r11_1.v + p.r00_1.v);            out.RH_1_0_0_err[i] = 2 * std::sqrt(4.0 * s_r11_1 * s_r11_1 + s_r00_1 * s_r00_1);
  out.RH_1_2_0[i] = 2 * 0.4 * (p.r11_1.v - p.r00_1.v);             out.RH_1_2_0_err[i] = 2 * 0.4 * std::sqrt(s_r11_1 * s_r11_1 + s_r00_1 * s_r00_1);
  out.RH_1_2_1[i] = -2 *  kSqrt12Over5 * p.re_r10_1.v;              out.RH_1_2_1_err[i] = 2 * kSqrt12Over5 * s_re_r10_1;
  out.RH_1_2_2[i] = 2 * kSqrt6Over5 * p.r1m1_1.v;                 out.RH_1_2_2_err[i] = 2 * kSqrt6Over5 * s_r1m1_1;

  out.RH_2_2_1[i] =  -2 * kSqrt12Over5 * p.im_r10_2.v;              out.RH_2_2_1_err[i] = 2 * kSqrt12Over5 * s_im_r10_2;
  out.RH_2_2_2[i] = 2 * kSqrt6Over5 * p.im_r1m1_2.v;              out.RH_2_2_2_err[i] = 2 * kSqrt6Over5 * s_im_r1m1_2;
  out.RH_3_2_1[i] =  -2 * kSqrt12Over5 * p.im_r10_3.v;              out.RH_3_2_1_err[i] = 2 * kSqrt12Over5 * s_im_r10_3;
  out.RH_3_2_2[i] = 2 * kSqrt6Over5 * p.im_r1m1_3.v;              out.RH_3_2_2_err[i] = 2 * kSqrt6Over5 * s_im_r1m1_3;

  out.RH_5_0_0[i] =  -2 * (2.0 * p.r11_5.v + p.r00_5.v);                out.RH_5_0_0_err[i] = 2 * std::sqrt(4.0 * s_r11_5 * s_r11_5 + s_r00_5 * s_r00_5);
  out.RH_5_2_0[i] = 2 * 0.4 * (p.r11_5.v - p.r00_5.v);                 out.RH_5_2_0_err[i] = 2 * 0.4 * std::sqrt(s_r11_5 * s_r11_5 + s_r00_5 * s_r00_5);
  out.RH_5_2_1[i] =  -2 * kSqrt12Over5 * p.re_r10_5.v;                  out.RH_5_2_1_err[i] = 2 * kSqrt12Over5 * s_re_r10_5;
  out.RH_5_2_2[i] = 2 * kSqrt6Over5 * p.r1m1_5.v;                     out.RH_5_2_2_err[i] = 2 * kSqrt6Over5 * s_r1m1_5;

  out.RH_6_2_1[i] =  -2 * kSqrt12Over5 * p.im_r10_6.v;              out.RH_6_2_1_err[i] = 2 * kSqrt12Over5 * s_im_r10_6;
  out.RH_6_2_2[i] = 2 * kSqrt6Over5 * p.im_r1m1_6.v;              out.RH_6_2_2_err[i] = 2 * kSqrt6Over5 * s_im_r1m1_6;
  out.RH_7_2_1[i] =  -2 * kSqrt12Over5 * p.im_r10_7.v;              out.RH_7_2_1_err[i] = 2 * kSqrt12Over5 * s_im_r10_7;
  out.RH_7_2_2[i] = 2 * kSqrt6Over5 * p.im_r1m1_7.v;              out.RH_7_2_2_err[i] = 2 * kSqrt6Over5 * s_im_r1m1_7;

  out.RH_8_0_0[i] =  -2 * (2.0 * p.r11_8.v + p.r00_8.v);                out.RH_8_0_0_err[i] = 2 * std::sqrt(4.0 * s_r11_8 * s_r11_8 + s_r00_8 * s_r00_8);
  out.RH_8_2_0[i] = 2 * 0.4 * (p.r11_8.v - p.r00_8.v);                 out.RH_8_2_0_err[i] = 2 * 0.4 * std::sqrt(s_r11_8 * s_r11_8 + s_r00_8 * s_r00_8);
  out.RH_8_2_1[i] =  -2 * kSqrt12Over5 * p.re_r10_8.v;                  out.RH_8_2_1_err[i] = 2 * kSqrt12Over5 * s_re_r10_8;
  out.RH_8_2_2[i] = 2 * kSqrt6Over5 * p.r1m1_8.v;                     out.RH_8_2_2_err[i] = 2 * kSqrt6Over5 * s_r1m1_8;
}

static std::string NormalizeKey(std::string key) {
  std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) { return std::tolower(c); });
  if (key == "rho" || key == "e-rho" || key == "electron_rho" || key == "hermes" || key == "hermes_rho") return "e_rho";
  if (key == "phi" || key == "e-phi" || key == "electron_phi" || key == "hermes_phi" || key == "e_phi_h" || key == "phi_h" || key == "phi_hydrogen") return "e_phi";
  if (key == "e-omega" || key == "electron_omega" || key == "hermes_omega") return "e_omega";
  if (key == "mu-rho" || key == "muon_rho" || key == "m_rho") return "mu_rho";
  if (key == "omega" || key == "mu-omega" || key == "muon_omega" || key == "m_omega") return "mu_omega";
  return key;
}

static std::string DefaultOutFile(const std::string& datasetKey) {
  return  datasetKey + "_moments.root";
}

static DatasetSpec GetDatasetSpec(const std::string& requestedKey) {
  const std::string key = NormalizeKey(requestedKey);

  if (key == "e_rho") {
    DatasetSpec ds;
    ds.key = "e_rho";
    ds.title = "HERMES electroproduction rho (electron beam)";
    ds.q2 = {0.82, 1.19, 1.66, 3.06};
    ds.bins = {
      SDMEsTable{
        VE(0.349,0.026,0.061), VE(0.028,0.028,0.020), VE(-0.024,0.013,0.021),
        VE(0.283,0.023,0.049), VE(-0.037,0.044,0.032), VE(0.023,0.019,0.007), VE(-0.054,0.039,0.013), VE(0.002,0.041,0.008), VE(0.022,0.079,0.026),
        VE(-0.015,0.010,0.007), VE(0.009,0.011,0.019), VE(-0.011,0.010,0.013), VE(-0.003,0.078,0.021), VE(0.019,0.053,0.007), VE(0.013,0.062,0.008),
        VE(-0.039,0.017,0.018), VE(0.021,0.051,0.010), VE(-0.294,0.019,0.038), VE(0.151,0.028,0.026), VE(-0.149,0.015,0.010), VE(0.079,0.068,0.011), VE(0.040,0.043,0.011), VE(0.121,0.038,0.039)
      },
      SDMEsTable{
        VE(0.368,0.018,0.011), VE(0.029,0.007,0.003), VE(-0.014,0.010,0.010),
        VE(0.262,0.018,0.024), VE(-0.043,0.012,0.006), VE(0.022,0.012,0.018), VE(0.011,0.032,0.018), VE(-0.041,0.026,0.005), VE(0.040,0.084,0.014),
        VE(-0.011,0.006,0.006), VE(0.008,0.007,0.006), VE(0.002,0.007,0.007), VE(0.023,0.056,0.013), VE(0.056,0.045,0.004), VE(0.072,0.053,0.011),
        VE(-0.034,0.013,0.013), VE(0.000,0.033,0.004), VE(-0.255,0.016,0.022), VE(0.171,0.007,0.000), VE(-0.167,0.007,0.003), VE(0.092,0.038,0.010), VE(0.020,0.031,0.008), VE(0.094,0.017,0.017)
      },
      SDMEsTable{
        VE(0.397,0.017,0.018), VE(0.035,0.007,0.011), VE(-0.019,0.010,0.003),
        VE(0.274,0.019,0.024), VE(-0.036,0.012,0.012), VE(0.005,0.012,0.024), VE(0.007,0.031,0.009), VE(-0.074,0.025,0.005), VE(0.054,0.086,0.011),
        VE(-0.008,0.006,0.011), VE(-0.013,0.007,0.003), VE(0.002,0.007,0.004), VE(-0.005,0.055,0.010), VE(0.051,0.044,0.006), VE(-0.018,0.054,0.004),
        VE(-0.023,0.013,0.008), VE(-0.031,0.032,0.007), VE(-0.239,0.017,0.011), VE(0.161,0.006,0.004), VE(-0.167,0.006,0.005), VE(0.039,0.036,0.004), VE(0.074,0.034,0.002), VE(0.057,0.015,0.019)
      },
      SDMEsTable{
        VE(0.454,0.014,0.011), VE(0.026,0.007,0.003), VE(0.001,0.009,0.007),
        VE(0.204,0.017,0.012), VE(-0.009,0.013,0.010), VE(0.022,0.013,0.008), VE(0.037,0.034,0.002), VE(0.048,0.024,0.006), VE(0.010,0.085,0.016),
        VE(-0.021,0.006,0.016), VE(0.020,0.007,0.008), VE(-0.010,0.007,0.007), VE(-0.109,0.047,0.004), VE(-0.002,0.035,0.005), VE(0.004,0.045,0.014),
        VE(-0.018,0.012,0.010), VE(-0.026,0.028,0.005), VE(-0.197,0.017,0.012), VE(0.141,0.006,0.008), VE(-0.156,0.006,0.010), VE(0.187,0.034,0.018), VE(0.098,0.032,0.005), VE(0.151,0.015,0.007)
      }
    };
    return ds;
  }


  if (key == "e_phi") {
    DatasetSpec ds;
    ds.key = "e_phi";
    ds.title = "HERMES electroproduction phi (hydrogen target only)";
    ds.q2 = {1.20, 1.70, 4.50};
    ds.bins = {
      SDMEsTable{
        VE(0.26587,0.01747,0.00898), VE(-0.01791,0.01363,0.02109), VE(-0.00560,0.01949,0.00554),
        VE(0.31918,0.02780,0.00435), VE(0.01965,0.01914,0.03653), VE(-0.00778,0.01975,0.02628), VE(0.04024,0.02424,0.00360), VE(0.00772,0.06678,0.00674), VE(0.14508,0.13646,0.04381),
        VE(-0.01891,0.01186,0.02587), VE(0.01249,0.01561,0.00930), VE(-0.02733,0.01480,0.00972), VE(-0.12475,0.15702,0.01144), VE(-0.13716,0.10893,0.02097), VE(-0.24102,0.13534,0.01771),
        VE(-0.02139,0.02424,0.01095), VE(-0.08304,0.10327,0.01463), VE(-0.27705,0.02854,0.01162), VE(0.17033,0.00955,0.01217), VE(-0.14678,0.00914,0.00393), VE(0.02572,0.08588,0.02175), VE(0.08476,0.07824,0.01492), VE(-0.01768,0.01362,0.02586)
      },
      SDMEsTable{
        VE(0.34252,0.01965,0.00919), VE(-0.00440,0.01221,0.01444), VE(0.02598,0.01865,0.00826),
        VE(0.28006,0.02684,0.01108), VE(-0.00920,0.01618,0.04668), VE(-0.01632,0.01745,0.00879), VE(-0.00281,0.02999,0.00474), VE(0.09648,0.05523,0.01392), VE(-0.20985,0.14228,0.00837),
        VE(-0.00042,0.01119,0.02339), VE(0.00177,0.01421,0.01114), VE(-0.02196,0.01398,0.00752), VE(-0.07973,0.11871,0.00798), VE(0.00648,0.09232,0.01549), VE(-0.15893,0.11723,0.00675),
        VE(0.02931,0.02305,0.01632), VE(-0.08855,0.07884,0.00406), VE(-0.26642,0.02632,0.02976), VE(0.16208,0.00860,0.00875), VE(-0.16801,0.00867,0.00655), VE(0.16662,0.06566,0.01911), VE(0.10157,0.07364,0.01133), VE(0.02859,0.01744,0.02096)
      },
      SDMEsTable{
        VE(0.40301,0.02022,0.01855), VE(-0.01809,0.01275,0.01803), VE(-0.02619,0.01746,0.01435),
        VE(0.26413,0.02369,0.02038), VE(0.00418,0.01848,0.02072), VE(-0.00443,0.01820,0.02672), VE(0.01879,0.03502,0.02007), VE(-0.04420,0.05054,0.01140), VE(0.29582,0.14240,0.02704),
        VE(0.01220,0.01083,0.02672), VE(-0.01883,0.01326,0.01451), VE(-0.01407,0.01256,0.00680), VE(-0.08306,0.10176,0.01025), VE(-0.15341,0.08244,0.00616), VE(0.08515,0.10236,0.01747),
        VE(-0.00178,0.02187,0.01718), VE(0.01682,0.06746,0.01041), VE(-0.23127,0.02351,0.00530), VE(0.14914,0.00968,0.00675), VE(-0.18073,0.00884,0.00701), VE(-0.05928,0.06387,0.01189), VE(0.16047,0.06218,0.01116), VE(0.02543,0.01886,0.01851)
      }
    };
    return ds;
  }

  if (key == "e_omega") {
    DatasetSpec ds;
    ds.key = "e_omega";
    ds.title = "HERMES leptoproduction omega (electron beam)";
    ds.q2 = {1.28, 2.00, 4.00};
    ds.bins = {
      SDMEsTable{
        VE(0.164,0.034,0.022), VE(0.005,0.021,0.004), VE(-0.004,0.032,0.000),
        VE(-0.032,0.050,0.032), VE(-0.005,0.032,0.013), VE(0.012,0.030,0.012), VE(0.009,0.049,0.011),
        VE(0.044,0.096,0.008), VE(-0.147,0.210,0.039),
        VE(-0.074,0.020,0.021), VE(-0.047,0.024,0.007),
        VE(0.070,0.025,0.013), VE(-0.326,0.223,0.058), VE(0.276,0.171,0.049), VE(-0.507,0.212,0.093),
        VE(0.063,0.040,0.015), VE(0.074,0.153,0.013),
        VE(0.172,0.048,0.027), VE(0.038,0.016,0.018),
        VE(-0.062,0.015,0.012), VE(0.163,0.139,0.030), VE(0.088,0.143,0.021), VE(0.031,0.029,0.001)
      },
      SDMEsTable{
        VE(0.166,0.030,0.044), VE(-0.060,0.020,0.011), VE(-0.023,0.031,0.003),
        VE(-0.175,0.049,0.037), VE(-0.090,0.031,0.012), VE(0.042,0.030,0.003), VE(0.039,0.049,0.013),
        VE(0.047,0.076,0.009), VE(0.035,0.196,0.026),
        VE(-0.050,0.020,0.012), VE(-0.078,0.025,0.021),
        VE(-0.015,0.024,0.017), VE(-0.161,0.198,0.030), VE(-0.120,0.155,0.021), VE(-0.026,0.188,0.005),
        VE(-0.037,0.041,0.012), VE(-0.110,0.131,0.021),
        VE(0.133,0.050,0.043), VE(0.022,0.015,0.010),
        VE(-0.069,0.012,0.014), VE(-0.006,0.125,0.009), VE(0.078,0.137,0.028), VE(0.029,0.025,0.012)
      },
      SDMEsTable{
        VE(0.179,0.031,0.036), VE(0.016,0.019,0.022), VE(0.008,0.031,0.014),
        VE(-0.314,0.053,0.090), VE(0.073,0.034,0.016), VE(0.036,0.034,0.016), VE(-0.032,0.053,0.015),
        VE(0.073,0.076,0.018), VE(-0.197,0.171,0.045),
        VE(-0.070,0.021,0.029), VE(0.008,0.025,0.009), VE(0.043,0.026,0.026), VE(0.046,0.204,0.023),
        VE(-0.312,0.144,0.080), VE(0.185,0.178,0.063),
        VE(0.003,0.044,0.012), VE(0.088,0.124,0.024),
        VE(0.163,0.057,0.029), VE(0.053,0.015,0.022), VE(-0.046,0.014,0.013), VE(0.170,0.128,0.042),
        VE(0.280,0.119,0.067), VE(0.068,0.027,0.016)
      }
    };
    return ds;
  }

  if (key == "mu_omega") {
    DatasetSpec ds;
    ds.key = "mu_omega";
    ds.title = "COMPASS Leptoproduction omega (muon beam)";
    ds.q2 = {1.16, 1.64, 3.58};
    ds.bins = {
      SDMEsTable{
        VE(0.313,0.027,0.027), VE(0.020,0.017,0.016), VE(0.034,0.025,0.034),
        VE(-0.029,0.037,0.015), VE(-0.080,0.023,0.033), VE(0.064,0.025,0.012), VE(-0.002,0.041,0.028), VE(-0.026,0.080,0.187), VE(-0.064,0.213,0.404),
        VE(0.009,0.014,0.035), VE(-0.002,0.019,0.015), VE(0.021,0.018,0.016), VE(0.159,0.161,0.216), VE(0.092,0.123,0.120), VE(-0.023,0.163,0.451),
        VE(-0.024,0.028,0.028), VE(0.167,0.106,0.229), VE(-0.027,0.037,0.025), VE(0.120,0.012,0.008), VE(-0.102,0.012,0.009), VE(0.079,0.127,0.140), VE(0.164,0.119,0.183), VE(0.059,0.021,0.022)
      },
      SDMEsTable{
        VE(0.322,0.029,0.046), VE(0.014,0.017,0.016), VE(0.009,0.026,0.043),
        VE(-0.011,0.037,0.036), VE(-0.074,0.025,0.055), VE(0.066,0.025,0.018), VE(-0.090,0.043,0.073), VE(0.175,0.087,0.109), VE(0.178,0.209,0.166),
        VE(-0.028,0.015,0.041), VE(-0.028,0.020,0.052), VE(0.038,0.019,0.024), VE(-0.006,0.186,0.079), VE(-0.260,0.132,0.126), VE(-0.186,0.169,0.191),
        VE(-0.013,0.029,0.046), VE(0.082,0.119,0.162), VE(0.036,0.037,0.047), VE(0.086,0.013,0.013), VE(-0.075,0.012,0.012), VE(0.061,0.136,0.070), VE(-0.131,0.112,0.116), VE(0.130,0.023,0.039)
      },
      SDMEsTable{
        VE(0.415,0.038,0.052), VE(0.016,0.023,0.013), VE(0.024,0.032,0.016),
        VE(-0.101,0.048,0.090), VE(-0.080,0.035,0.045), VE(0.069,0.032,0.046), VE(-0.144,0.063,0.087), VE(-0.038,0.109,0.062), VE(0.224,0.245,0.101),
        VE(-0.037,0.019,0.026), VE(-0.026,0.024,0.020), VE(0.008,0.024,0.017), VE(0.330,0.246,0.150), VE(-0.188,0.163,0.159), VE(0.231,0.203,0.300),
        VE(-0.042,0.038,0.024), VE(0.024,0.161,0.131), VE(0.103,0.049,0.065), VE(0.101,0.017,0.050), VE(-0.088,0.016,0.033), VE(-0.150,0.167,0.182), VE(0.340,0.153,0.144), VE(0.219,0.031,0.077)
      }
    };
    return ds;
  }

  if (key == "mu_rho") {
    DatasetSpec ds;
    ds.key = "mu_rho";
    ds.title = "COMPASS Leptoproduction rho (muon beam)";
    ds.q2 = {1.14, 1.60, 2.80, 6.02};
    ds.bins = {
      SDMEsTable{
        VE(0.4080,0.0056,0.0243), VE(0.0452,0.0034,0.0058), VE(-0.0213,0.0044,0.0055),
        VE(0.2781,0.0058,0.0088), VE(-0.0521,0.0044,0.0049), VE(0.0505,0.0043,0.0031), VE(-0.0217,0.0092,0.0109), VE(0.0144,0.0104,0.0027), VE(0.0095,0.0302,0.0269),
        VE(0.0014,0.0026,0.0040), VE(0.0017,0.0032,0.0027), VE(0.0006,0.0031,0.0025), VE(-0.0079,0.0215,0.0444), VE(0.0227,0.0163,0.0310), VE(-0.0154,0.0206,0.0209),
        VE(-0.0252,0.0051,0.0083), VE(0.0038,0.0134,0.0110), VE(-0.2763,0.0060,0.0083), VE(0.1774,0.0023,0.0042), VE(-0.1695,0.0021,0.0033), VE(0.0230,0.0148,0.0158), VE(0.0253,0.0147,0.0111), VE(0.1150,0.0050,0.0080)
      },
      SDMEsTable{
        VE(0.4749,0.0055,0.0201), VE(0.0431,0.0034,0.0060), VE(-0.0074,0.0042,0.0030),
        VE(0.2337,0.0057,0.0074), VE(-0.0439,0.0045,0.0055), VE(0.0508,0.0045,0.0041), VE(-0.0441,0.0097,0.0103), VE(-0.0068,0.0105,0.0084), VE(-0.0041,0.0304,0.0294),
        VE(0.0009,0.0025,0.0017), VE(0.0079,0.0031,0.0054), VE(-0.0074,0.0031,0.0027), VE(0.0063,0.0206,0.0086), VE(0.0168,0.0156,0.0140), VE(-0.0105,0.0195,0.0177),
        VE(-0.0099,0.0049,0.0069), VE(0.0279,0.0131,0.0192), VE(-0.2300,0.0059,0.0045), VE(0.1726,0.0023,0.0025), VE(-0.1591,0.0023,0.0038), VE(0.0482,0.0152,0.0082), VE(0.0400,0.0146,0.0059), VE(0.1419,0.0052,0.0122)
      },
      SDMEsTable{
        VE(0.5490,0.0085,0.0193), VE(0.0508,0.0052,0.0069), VE(-0.0081,0.0064,0.0073),
        VE(0.2437,0.0089,0.0061), VE(-0.0713,0.0070,0.0014), VE(0.0612,0.0069,0.0060), VE(-0.0532,0.0156,0.0420), VE(0.0209,0.0161,0.0114), VE(-0.0498,0.0477,0.0385),
        VE(0.0014,0.0037,0.0055), VE(0.0087,0.0047,0.0042), VE(0.0003,0.0046,0.0029), VE(-0.0400,0.0314,0.0156), VE(0.0397,0.0243,0.0562), VE(-0.0575,0.0309,0.0569),
        VE(-0.0157,0.0074,0.0101), VE(0.0051,0.0205,0.0094), VE(-0.2586,0.0089,0.0165), VE(0.1938,0.0036,0.0083), VE(-0.1829,0.0035,0.0067), VE(0.0851,0.0230,0.0343), VE(0.0466,0.0231,0.0365), VE(0.1950,0.0081,0.0213)
      },
      SDMEsTable{
        VE(0.5319,0.0183,0.0555), VE(0.0358,0.0110,0.0089), VE(0.0059,0.0136,0.0057),
        VE(0.1647,0.0193,0.0220), VE(-0.0613,0.0150,0.0234), VE(0.0628,0.0151,0.0284), VE(-0.0419,0.0326,0.0254), VE(0.0212,0.0380,0.0694), VE(0.2174,0.1065,0.1142),
        VE(0.0316,0.0080,0.0230), VE(-0.0096,0.0100,0.0087), VE(-0.0067,0.0102,0.0033), VE(0.0716,0.0755,0.0571), VE(-0.0800,0.0546,0.0492), VE(-0.1683,0.0698,0.0418),
        VE(-0.0122,0.0159,0.0146), VE(0.0702,0.0495,0.0476), VE(-0.1450,0.0199,0.0228), VE(0.1562,0.0078,0.0145), VE(-0.1513,0.0077,0.0077), VE(0.0296,0.0553,0.0073), VE(0.0471,0.0548,0.0625), VE(0.2021,0.0167,0.0406)
      }
    };
    return ds;
  }

  throw std::runtime_error(
    "Unknown dataset key '" + requestedKey +
    "'. Supported keys: e_rho, e_phi, e_omega, mu_rho or mu_omega."
  );
}

static void BuildOutputs(const DatasetSpec& ds, OutputArrays& out) {
  if (ds.q2.size() != ds.bins.size()) {
    throw std::runtime_error("Dataset '" + ds.key + "' has mismatched q2/bin sizes");
  }

  ResizeOutputs(out, ds.q2.size());
  for (std::size_t i = 0; i < ds.q2.size(); ++i) {
    FillBin(out, static_cast<int>(i), ds.q2[i], ds.bins[i]);
  }
}

} // namespace

void MakeLeptoproductionMoments(const std::string& dataset,
                                const std::filesystem::path& output,
                                const std::string& treeName) {

  const DatasetSpec ds = GetDatasetSpec(dataset.empty() ? "e_rho" : dataset);
  const std::filesystem::path outPath = output.empty()
      ? std::filesystem::path("InputFiles/Experiment") / DefaultOutFile(ds.key)
      : output;
  if (!outPath.parent_path().empty()) {
    std::filesystem::create_directories(outPath.parent_path());
  }

  OutputArrays out;
  BuildOutputs(ds, out);

  TFile fout(outPath.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    throw std::runtime_error("Failed to open output file: " + outPath.string());
  }

  TTree tree(treeName.c_str(), (std::string("Experimental moments and uncertainties vs Q2 for ") + ds.title).c_str());

  tree.Branch("Nbins", &out.Nbins, "Nbins/I");
  auto br = [&](const char* name, double* arr) {
    tree.Branch(name, arr, (std::string(name) + "[Nbins]/D").c_str());
  };

  br("Q2", out.Q2);
  br("RH04_0_0", out.RH04_0_0); br("RH04_0_0_err", out.RH04_0_0_err);
  br("RH04_2_0", out.RH04_2_0); br("RH04_2_0_err", out.RH04_2_0_err);
  br("RH04_2_1", out.RH04_2_1); br("RH04_2_1_err", out.RH04_2_1_err);
  br("RH04_2_2", out.RH04_2_2); br("RH04_2_2_err", out.RH04_2_2_err);

  br("RH_1_0_0", out.RH_1_0_0); br("RH_1_0_0_err", out.RH_1_0_0_err);
  br("RH_1_2_0", out.RH_1_2_0); br("RH_1_2_0_err", out.RH_1_2_0_err);
  br("RH_1_2_1", out.RH_1_2_1); br("RH_1_2_1_err", out.RH_1_2_1_err);
  br("RH_1_2_2", out.RH_1_2_2); br("RH_1_2_2_err", out.RH_1_2_2_err);
  br("RH_2_2_1", out.RH_2_2_1); br("RH_2_2_1_err", out.RH_2_2_1_err);
  br("RH_2_2_2", out.RH_2_2_2); br("RH_2_2_2_err", out.RH_2_2_2_err);
  br("RH_3_2_1", out.RH_3_2_1); br("RH_3_2_1_err", out.RH_3_2_1_err);
  br("RH_3_2_2", out.RH_3_2_2); br("RH_3_2_2_err", out.RH_3_2_2_err);

  br("RH_5_0_0", out.RH_5_0_0); br("RH_5_0_0_err", out.RH_5_0_0_err);
  br("RH_5_2_0", out.RH_5_2_0); br("RH_5_2_0_err", out.RH_5_2_0_err);
  br("RH_5_2_1", out.RH_5_2_1); br("RH_5_2_1_err", out.RH_5_2_1_err);
  br("RH_5_2_2", out.RH_5_2_2); br("RH_5_2_2_err", out.RH_5_2_2_err);

  br("RH_6_2_1", out.RH_6_2_1); br("RH_6_2_1_err", out.RH_6_2_1_err);
  br("RH_6_2_2", out.RH_6_2_2); br("RH_6_2_2_err", out.RH_6_2_2_err);
  br("RH_7_2_1", out.RH_7_2_1); br("RH_7_2_1_err", out.RH_7_2_1_err);
  br("RH_7_2_2", out.RH_7_2_2); br("RH_7_2_2_err", out.RH_7_2_2_err);

  br("RH_8_0_0", out.RH_8_0_0); br("RH_8_0_0_err", out.RH_8_0_0_err);
  br("RH_8_2_0", out.RH_8_2_0); br("RH_8_2_0_err", out.RH_8_2_0_err);
  br("RH_8_2_1", out.RH_8_2_1); br("RH_8_2_1_err", out.RH_8_2_1_err);
  br("RH_8_2_2", out.RH_8_2_2); br("RH_8_2_2_err", out.RH_8_2_2_err);

  TObjString datasetKey(ds.key.c_str());
  datasetKey.Write("dataset_key");
  TObjString datasetTitle(ds.title.c_str());
  datasetTitle.Write("dataset_title");
  TParameter<int>("NbinsMeta", static_cast<int>(ds.q2.size())).Write();

  tree.Fill();
  tree.Write();
  fout.Close();

  std::cout << "Saved dataset '" << ds.key << "' (" << ds.title << ") to " << outPath << std::endl;
  std::cout << "Q2 bins: ";
  for (std::size_t i = 0; i < ds.q2.size(); ++i) {
    std::cout << ds.q2[i] << (i + 1 == ds.q2.size() ? "" : ", ");
  }
  std::cout << std::endl;
}

} // namespace emi
