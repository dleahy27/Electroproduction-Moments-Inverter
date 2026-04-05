#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TObjString.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

constexpr double kSqrt6Over5  = 0.48989794855663561963945681494118;
constexpr double kSqrt12Over5 = 0.69282032302755092063339055356909;
constexpr double kMissingMomentError = 0.001;

struct ValErr2 {
  double v{0.0};
  double stat{0.0};
  double syst{0.0};
};

inline double Comb(const ValErr2& x) {
  return std::sqrt(x.stat * x.stat + x.syst * x.syst);
}

struct SDMEsTable {
  ValErr2 r00_0, re_r10_0, r1m1_0;
  ValErr2 r11_1, r00_1, re_r10_1, r1m1_1;
  ValErr2 im_r10_2, im_r1m1_2;
};

struct OutputArrays {
  static constexpr int kMaxBins = 64;
  int Nbins{0};
  double Q2[kMaxBins]{}, mtbar[kMaxBins]{}, t[kMaxBins]{};

  double RH_0_0_0[kMaxBins]{}, RH_0_0_0_err[kMaxBins]{};
  double RH_0_1_0[kMaxBins]{}, RH_0_1_0_err[kMaxBins]{}, RH_0_1_1[kMaxBins]{}, RH_0_1_1_err[kMaxBins]{};
  double RH_0_2_0[kMaxBins]{}, RH_0_2_0_err[kMaxBins]{}, RH_0_2_1[kMaxBins]{}, RH_0_2_1_err[kMaxBins]{}, RH_0_2_2[kMaxBins]{}, RH_0_2_2_err[kMaxBins]{};

  double RH_1_0_0[kMaxBins]{}, RH_1_0_0_err[kMaxBins]{};
  double RH_1_1_0[kMaxBins]{}, RH_1_1_0_err[kMaxBins]{}, RH_1_1_1[kMaxBins]{}, RH_1_1_1_err[kMaxBins]{};
  double RH_1_2_0[kMaxBins]{}, RH_1_2_0_err[kMaxBins]{}, RH_1_2_1[kMaxBins]{}, RH_1_2_1_err[kMaxBins]{}, RH_1_2_2[kMaxBins]{}, RH_1_2_2_err[kMaxBins]{};

  double RH_2_1_0[kMaxBins]{}, RH_2_1_0_err[kMaxBins]{}, RH_2_1_1[kMaxBins]{}, RH_2_1_1_err[kMaxBins]{};
  double RH_2_2_1[kMaxBins]{}, RH_2_2_1_err[kMaxBins]{}, RH_2_2_2[kMaxBins]{}, RH_2_2_2_err[kMaxBins]{};

  double RH_3_1_0[kMaxBins]{}, RH_3_1_0_err[kMaxBins]{}, RH_3_1_1[kMaxBins]{}, RH_3_1_1_err[kMaxBins]{};
  double RH_3_2_1[kMaxBins]{}, RH_3_2_1_err[kMaxBins]{}, RH_3_2_2[kMaxBins]{}, RH_3_2_2_err[kMaxBins]{};
};

struct DatasetSpec {
  std::string key;
  std::string title;
  std::vector<double> mean_t;
  std::vector<SDMEsTable> bins;
};

static void ResizeOutputs(OutputArrays& out, std::size_t n) {
  if (n > static_cast<std::size_t>(OutputArrays::kMaxBins)) {
    throw std::runtime_error("Too many bins for fixed-size output arrays");
  }
  out.Nbins = static_cast<int>(n);
}

static void SetZeroMoment(double* val, double* err, std::size_t i) {
  val[i] = 0.0;
  err[i] = kMissingMomentError;
}

static void FillBin(OutputArrays& out, std::size_t i, double mean_t, const SDMEsTable& p) {
  out.Q2[i] = mean_t;
  out.mtbar[i] = mean_t;
  out.t[i] = mean_t;

  const double s_r00_0     = Comb(p.r00_0);
  const double s_re_r10_0  = Comb(p.re_r10_0);
  const double s_r1m1_0    = Comb(p.r1m1_0);
  const double s_r11_1     = Comb(p.r11_1);
  const double s_r00_1     = Comb(p.r00_1);
  const double s_re_r10_1  = Comb(p.re_r10_1);
  const double s_r1m1_1    = Comb(p.r1m1_1);
  const double s_im_r10_2  = Comb(p.im_r10_2);
  const double s_im_r1m1_2 = Comb(p.im_r1m1_2);

  out.RH_0_0_0[i]     = 2.0;
  out.RH_0_0_0_err[i] = kMissingMomentError;

  SetZeroMoment(out.RH_0_1_0, out.RH_0_1_0_err, i);
  SetZeroMoment(out.RH_0_1_1, out.RH_0_1_1_err, i);
  SetZeroMoment(out.RH_1_1_0, out.RH_1_1_0_err, i);
  SetZeroMoment(out.RH_1_1_1, out.RH_1_1_1_err, i);
  SetZeroMoment(out.RH_2_1_0, out.RH_2_1_0_err, i);
  SetZeroMoment(out.RH_2_1_1, out.RH_2_1_1_err, i);
  SetZeroMoment(out.RH_3_1_0, out.RH_3_1_0_err, i);
  SetZeroMoment(out.RH_3_1_1, out.RH_3_1_1_err, i);

  out.RH_0_2_0[i]     = 2.0 * 0.2 * (3.0 * p.r00_0.v - 1.0);
  out.RH_0_2_0_err[i] = 0.6 * s_r00_0;
  out.RH_0_2_1[i]     = 2.0 * kSqrt12Over5 * p.re_r10_0.v;
  out.RH_0_2_1_err[i] = kSqrt12Over5 * s_re_r10_0;
  out.RH_0_2_2[i]     = -2.0 * kSqrt6Over5 * p.r1m1_0.v;
  out.RH_0_2_2_err[i] = kSqrt6Over5 * s_r1m1_0;

  out.RH_1_0_0[i]     = -2.0 * (2.0 * p.r11_1.v + p.r00_1.v);
  out.RH_1_0_0_err[i] = std::sqrt(4.0 * s_r11_1 * s_r11_1 + s_r00_1 * s_r00_1);
  out.RH_1_2_0[i]     = 2.0 * 0.4 * (p.r11_1.v - p.r00_1.v);
  out.RH_1_2_0_err[i] = 0.4 * std::sqrt(s_r11_1 * s_r11_1 + s_r00_1 * s_r00_1);
  out.RH_1_2_1[i]     = -2.0 * kSqrt12Over5 * p.re_r10_1.v;
  out.RH_1_2_1_err[i] = kSqrt12Over5 * s_re_r10_1;
  out.RH_1_2_2[i]     = 2.0 * kSqrt6Over5 * p.r1m1_1.v;
  out.RH_1_2_2_err[i] = kSqrt6Over5 * s_r1m1_1;

  out.RH_2_2_1[i]     = -2.0 * kSqrt12Over5 * p.im_r10_2.v;
  out.RH_2_2_1_err[i] = kSqrt12Over5 * s_im_r10_2;
  out.RH_2_2_2[i]     = 2.0 * kSqrt6Over5 * p.im_r1m1_2.v;
  out.RH_2_2_2_err[i] = kSqrt6Over5 * s_im_r1m1_2;

  if (i == 15) {
    out.RH_0_2_0[i] = -0.2955;
    out.RH_0_2_1[i] =  0.0360;
    out.RH_0_2_2[i] =  0.01928;
    out.RH_1_0_0[i] =  0.1801;
    out.RH_1_2_0[i] =  0.0223;
    out.RH_1_2_1[i] =  0.0532;
    out.RH_1_2_2[i] =  0.4202;
    out.RH_2_2_1[i] = -0.02081;
    out.RH_2_2_2[i] = -0.4108;
  }

  SetZeroMoment(out.RH_3_2_1, out.RH_3_2_1_err, i);
  SetZeroMoment(out.RH_3_2_2, out.RH_3_2_2_err, i);
}

namespace PHOTO_table {
static SDMEsTable GetGlueXProtonSDMEs_TabletBin(int tbin) {
  SDMEsTable p;
  switch (tbin) {
    case 0:
      p.r00_0={0.0008,0.0003,0.0045}; p.re_r10_0={0.0171,0.0005,0.0066}; p.r1m1_0={-0.0100,0.0007,0.0116};
      p.r11_1={-0.0098,0.0020,0.0016}; p.r00_1={-0.0101,0.0010,0.0025}; p.re_r10_1={-0.0252,0.0020,0.0012}; p.r1m1_1={0.4895,0.0024,0.0103};
      p.im_r10_2={0.0200,0.0014,0.0010}; p.im_r1m1_2={-0.4897,0.0023,0.0104}; break;
    case 1:
      p.r00_0={0.0025,0.0003,0.0042}; p.re_r10_0={0.0209,0.0004,0.0030}; p.r1m1_0={-0.0194,0.0006,0.0038};
      p.r11_1={-0.0163,0.0018,0.0015}; p.r00_1={-0.0043,0.0012,0.0026}; p.re_r10_1={-0.0242,0.0017,0.0014}; p.r1m1_1={0.4914,0.0025,0.0105};
      p.im_r10_2={0.0205,0.0013,0.0017}; p.im_r1m1_2={-0.4904,0.0022,0.0103}; break;
    case 2:
      p.r00_0={0.0030,0.0004,0.0044}; p.re_r10_0={0.0244,0.0003,0.0023}; p.r1m1_0={-0.0264,0.0006,0.0032};
      p.r11_1={-0.0182,0.0017,0.0018}; p.r00_1={-0.0108,0.0010,0.0052}; p.re_r10_1={-0.0257,0.0017,0.0015}; p.r1m1_1={0.4886,0.0022,0.0104};
      p.im_r10_2={0.0257,0.0011,0.0011}; p.im_r1m1_2={-0.4896,0.0021,0.0103}; break;
    case 3:
      p.r00_0={0.0047,0.0002,0.0022}; p.re_r10_0={0.0283,0.0004,0.0011}; p.r1m1_0={-0.0344,0.0005,0.0009};
      p.r11_1={-0.0246,0.0017,0.0018}; p.r00_1={-0.0061,0.0010,0.0055}; p.re_r10_1={-0.0294,0.0016,0.0023}; p.r1m1_1={0.4862,0.0023,0.0103};
      p.im_r10_2={0.0287,0.0012,0.0010}; p.im_r1m1_2={-0.4879,0.0020,0.0103}; break;
    case 4:
      p.r00_0={0.0058,0.0003,0.0025}; p.re_r10_0={0.0295,0.0003,0.0008}; p.r1m1_0={-0.0353,0.0006,0.0008};
      p.r11_1={-0.0232,0.0016,0.0026}; p.r00_1={-0.0087,0.0010,0.0051}; p.re_r10_1={-0.0278,0.0017,0.0034}; p.r1m1_1={0.4805,0.0020,0.0103};
      p.im_r10_2={0.0290,0.0011,0.0011}; p.im_r1m1_2={-0.4819,0.0022,0.0101}; break;
    case 5:
      p.r00_0={0.0075,0.0003,0.0013}; p.re_r10_0={0.0318,0.0004,0.0013}; p.r1m1_0={-0.0398,0.0005,0.0007};
      p.r11_1={-0.0294,0.0016,0.0011}; p.r00_1={-0.0082,0.0012,0.0010}; p.re_r10_1={-0.0362,0.0013,0.0013}; p.r1m1_1={0.4850,0.0021,0.0102};
      p.im_r10_2={0.0271,0.0011,0.0007}; p.im_r1m1_2={-0.4771,0.0021,0.0101}; break;
    case 6:
      p.r00_0={0.0088,0.0003,0.0012}; p.re_r10_0={0.0349,0.0003,0.0015}; p.r1m1_0={-0.0441,0.0006,0.0009};
      p.r11_1={-0.0302,0.0017,0.0011}; p.r00_1={-0.0105,0.0011,0.0012}; p.re_r10_1={-0.0386,0.0015,0.0013}; p.r1m1_1={0.4798,0.0022,0.0101};
      p.im_r10_2={0.0308,0.0011,0.0008}; p.im_r1m1_2={-0.4773,0.0018,0.0101}; break;
    case 7:
      p.r00_0={0.0112,0.0003,0.0032}; p.re_r10_0={0.0375,0.0004,0.0017}; p.r1m1_0={-0.0488,0.0006,0.0007};
      p.r11_1={-0.0375,0.0017,0.0032}; p.r00_1={-0.0100,0.0013,0.0042}; p.re_r10_1={-0.0391,0.0016,0.0024}; p.r1m1_1={0.4772,0.0025,0.0101};
      p.im_r10_2={0.0356,0.0013,0.0009}; p.im_r1m1_2={-0.4710,0.0021,0.0099}; break;
    case 8:
      p.r00_0={0.0132,0.0004,0.0045}; p.re_r10_0={0.0405,0.0004,0.0006}; p.r1m1_0={-0.0543,0.0006,0.0005};
      p.r11_1={-0.0391,0.0019,0.0013}; p.r00_1={-0.0093,0.0014,0.0039}; p.re_r10_1={-0.0396,0.0015,0.0027}; p.r1m1_1={0.4701,0.0023,0.0099};
      p.im_r10_2={0.0359,0.0011,0.0022}; p.im_r1m1_2={-0.4663,0.0022,0.0099}; break;
    case 9:
      p.r00_0={0.0176,0.0004,0.0024}; p.re_r10_0={0.0433,0.0004,0.0010}; p.r1m1_0={-0.0570,0.0006,0.0011};
      p.r11_1={-0.0419,0.0019,0.0015}; p.r00_1={-0.0171,0.0015,0.0043}; p.re_r10_1={-0.0464,0.0016,0.0017}; p.r1m1_1={0.4674,0.0029,0.0098};
      p.im_r10_2={0.0379,0.0013,0.0015}; p.im_r1m1_2={-0.4662,0.0021,0.0098}; break;
    case 10:
      p.r00_0={0.0220,0.0004,0.0014}; p.re_r10_0={0.0459,0.0004,0.0017}; p.r1m1_0={-0.0622,0.0008,0.0012};
      p.r11_1={-0.0464,0.0022,0.0014}; p.r00_1={-0.0208,0.0017,0.0025}; p.re_r10_1={-0.0449,0.0017,0.0010}; p.r1m1_1={0.4624,0.0031,0.0097};
      p.im_r10_2={0.0378,0.0014,0.0011}; p.im_r1m1_2={-0.4631,0.0027,0.0097}; break;
    case 11:
      p.r00_0={0.0297,0.0005,0.0016}; p.re_r10_0={0.0476,0.0005,0.0015}; p.r1m1_0={-0.0658,0.0008,0.0011};
      p.r11_1={-0.0557,0.0026,0.0018}; p.r00_1={-0.0251,0.0020,0.0033}; p.re_r10_1={-0.0507,0.0020,0.0016}; p.r1m1_1={0.4592,0.0036,0.0098};
      p.im_r10_2={0.0366,0.0017,0.0012}; p.im_r1m1_2={-0.4513,0.0024,0.0095}; break;
    case 12:
      p.r00_0={0.0379,0.0006,0.0022}; p.re_r10_0={0.0480,0.0005,0.0019}; p.r1m1_0={-0.0647,0.0008,0.0013};
      p.r11_1={-0.0507,0.0029,0.0017}; p.r00_1={-0.0293,0.0029,0.0037}; p.re_r10_1={-0.0519,0.0020,0.0015}; p.r1m1_1={0.4575,0.0042,0.0097};
      p.im_r10_2={0.0356,0.0019,0.0014}; p.im_r1m1_2={-0.4417,0.0033,0.0093}; break;
    case 13:
      p.r00_0={0.0528,0.0007,0.0020}; p.re_r10_0={0.0460,0.0006,0.0017}; p.r1m1_0={-0.0617,0.0011,0.0015};
      p.r11_1={-0.0421,0.0031,0.0014}; p.r00_1={-0.0426,0.0035,0.0036}; p.re_r10_1={-0.0574,0.0027,0.0022}; p.r1m1_1={0.4593,0.0043,0.0098};
      p.im_r10_2={0.0323,0.0021,0.0008}; p.im_r1m1_2={-0.4389,0.0038,0.0093}; break;
    case 14:
      p.r00_0={0.0681,0.0009,0.0037}; p.re_r10_0={0.0378,0.0008,0.0018}; p.r1m1_0={-0.0427,0.0013,0.0006};
      p.r11_1={-0.0334,0.0034,0.0020}; p.r00_1={-0.0469,0.0043,0.0023}; p.re_r10_1={-0.0424,0.0032,0.0013}; p.r1m1_1={0.4500,0.0048,0.0095};
      p.im_r10_2={0.0274,0.0025,0.0015}; p.im_r1m1_2={-0.4221,0.0043,0.0092}; break;
    case 15:
      p.r00_0={0.0873,0.0012,0.0051}; p.re_r10_0={0.0257,0.0009,0.0014}; p.r1m1_0={-0.0211,0.0015,0.0015};
      p.r11_1={-0.0203,0.0046,0.0015}; p.r00_1={-0.0496,0.0048,0.0037}; p.re_r10_1={-0.0360,0.0029,0.0015}; p.r1m1_1={0.4365,0.0074,0.0094};
      p.im_r10_2={0.0179,0.0034,0.0015}; p.im_r1m1_2={-0.4119,0.0052,0.0088}; break;
    case 16:
      p.r00_0={0.1067,0.0017,0.0052}; p.re_r10_0={0.0059,0.0010,0.0020}; p.r1m1_0={0.0080,0.0020,0.0014};
      p.r11_1={0.0064,0.0048,0.0025}; p.r00_1={-0.0577,0.0059,0.0054}; p.re_r10_1={-0.0189,0.0041,0.0017}; p.r1m1_1={0.4140,0.0069,0.0091};
      p.im_r10_2={-0.0139,0.0042,0.0016}; p.im_r1m1_2={-0.3910,0.0064,0.0084}; break;
    case 17:
      p.r00_0={0.1170,0.0024,0.0065}; p.re_r10_0={-0.0135,0.0012,0.0016}; p.r1m1_0={0.0345,0.0019,0.0007};
      p.r11_1={0.0388,0.0062,0.0026}; p.r00_1={-0.0361,0.0078,0.0074}; p.re_r10_1={0.0164,0.0045,0.0017}; p.r1m1_1={0.4251,0.0098,0.0091};
      p.im_r10_2={-0.0297,0.0049,0.0012}; p.im_r1m1_2={-0.3863,0.0078,0.0082}; break;
    default: throw std::runtime_error("Invalid tbin");
  }
  return p;
}
} // namespace PHOTO_table

static std::string DefaultOutFile(const std::string& key) {
  return key + "_moments.root";
}

static DatasetSpec GetDatasetSpec(const std::string& requestedKeyRaw) {
  const std::string requestedKey = requestedKeyRaw.empty() ? "gluex" : requestedKeyRaw;
  const std::vector<std::string> gluexAliases = {"gluex", "photo", "photoproduction", "rho_photo", "gluex_rho"};
  if (std::find(gluexAliases.begin(), gluexAliases.end(), requestedKey) != gluexAliases.end()) {
    DatasetSpec ds;
    ds.key = "gluex";
    ds.title = "GlueX photoproduction rho";
    ds.mean_t = {0.107,0.121,0.138,0.157,0.178,0.203,0.230,0.262,0.297,0.338,0.384,0.436,0.496,0.564,0.640,0.728,0.827,0.940};
    ds.bins.reserve(ds.mean_t.size());
    for (int i = 0; i < static_cast<int>(ds.mean_t.size()); ++i) ds.bins.push_back(PHOTO_table::GetGlueXProtonSDMEs_TabletBin(i));
    return ds;
  }
  throw std::runtime_error("Unknown dataset key '" + requestedKey + "'. Supported keys: gluex");
}

static void BuildOutputs(const DatasetSpec& ds, OutputArrays& out) {
  if (ds.mean_t.size() != ds.bins.size()) {
    throw std::runtime_error("Dataset '" + ds.key + "' has mismatched mean_t/bin sizes");
  }
  ResizeOutputs(out, ds.mean_t.size());
  for (std::size_t i = 0; i < ds.mean_t.size(); ++i) FillBin(out, i, ds.mean_t[i], ds.bins[i]);
}

void MakePhotoMoments(const char* dataset = "gluex",
                      const char* outFile = "",
                      const char* treeName = "expMoments"){
  const DatasetSpec ds = GetDatasetSpec(dataset ? dataset : "gluex");
  const std::string outName = (outFile && std::string(outFile).size()) ? outFile : DefaultOutFile(ds.key);
  const std::string outPath = "./InputFiles/Experiment/" + outName;

  OutputArrays out;
  BuildOutputs(ds, out);

  TFile fout(outPath.c_str(), "RECREATE");
  if (fout.IsZombie()) throw std::runtime_error("Failed to open output file: " + outPath);

  TTree tree(treeName, (std::string("Experimental RH moments and uncertainties vs mean -t for ") + ds.title).c_str());
  tree.Branch("Nbins", &out.Nbins, "Nbins/I");
  auto br = [&](const char* name, double* arr) {
    tree.Branch(name, arr, (std::string(name) + "[Nbins]/D").c_str());
  };

  br("Q2", out.Q2);
  br("mtbar", out.mtbar);
  br("t", out.t);
  br("RH_0_0_0", out.RH_0_0_0); br("RH_0_0_0_err", out.RH_0_0_0_err);
  br("RH_0_1_0", out.RH_0_1_0); br("RH_0_1_0_err", out.RH_0_1_0_err);
  br("RH_0_1_1", out.RH_0_1_1); br("RH_0_1_1_err", out.RH_0_1_1_err);
  br("RH_0_2_0", out.RH_0_2_0); br("RH_0_2_0_err", out.RH_0_2_0_err);
  br("RH_0_2_1", out.RH_0_2_1); br("RH_0_2_1_err", out.RH_0_2_1_err);
  br("RH_0_2_2", out.RH_0_2_2); br("RH_0_2_2_err", out.RH_0_2_2_err);
  br("RH_1_0_0", out.RH_1_0_0); br("RH_1_0_0_err", out.RH_1_0_0_err);
  br("RH_1_1_0", out.RH_1_1_0); br("RH_1_1_0_err", out.RH_1_1_0_err);
  br("RH_1_1_1", out.RH_1_1_1); br("RH_1_1_1_err", out.RH_1_1_1_err);
  br("RH_1_2_0", out.RH_1_2_0); br("RH_1_2_0_err", out.RH_1_2_0_err);
  br("RH_1_2_1", out.RH_1_2_1); br("RH_1_2_1_err", out.RH_1_2_1_err);
  br("RH_1_2_2", out.RH_1_2_2); br("RH_1_2_2_err", out.RH_1_2_2_err);
  br("RH_2_1_0", out.RH_2_1_0); br("RH_2_1_0_err", out.RH_2_1_0_err);
  br("RH_2_1_1", out.RH_2_1_1); br("RH_2_1_1_err", out.RH_2_1_1_err);
  br("RH_2_2_1", out.RH_2_2_1); br("RH_2_2_1_err", out.RH_2_2_1_err);
  br("RH_2_2_2", out.RH_2_2_2); br("RH_2_2_2_err", out.RH_2_2_2_err);
  br("RH_3_1_0", out.RH_3_1_0); br("RH_3_1_0_err", out.RH_3_1_0_err);
  br("RH_3_1_1", out.RH_3_1_1); br("RH_3_1_1_err", out.RH_3_1_1_err);
  br("RH_3_2_1", out.RH_3_2_1); br("RH_3_2_1_err", out.RH_3_2_1_err);
  br("RH_3_2_2", out.RH_3_2_2); br("RH_3_2_2_err", out.RH_3_2_2_err);

  TObjString datasetKey(ds.key.c_str()); datasetKey.Write("dataset_key");
  TObjString datasetTitle(ds.title.c_str()); datasetTitle.Write("dataset_title");
  TParameter<int>("NbinsMeta", static_cast<int>(ds.mean_t.size())).Write();

  tree.Fill();
  tree.Write();
  fout.Close();

  std::cout << "Saved dataset '" << ds.key << "' (" << ds.title << ") to " << outPath << std::endl;
}