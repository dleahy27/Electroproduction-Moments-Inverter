#!/usr/bin/env python3

"""Compare numerical, analytical, and BruFit minimizer performance."""

import matplotlib.pyplot as plt
import numpy as np
import ROOT

from analysis_utils import PROJECT_DIR, save_figure, style_axis, wrap_phase


# Settings: edit the timing and external BruFit paths when needed.
NUMERIC_TIMES = PROJECT_DIR / "numeric_times.csv"
ANALYTIC_TIMES = PROJECT_DIR / "analytic_times.csv"
ANALYTIC_FIT = PROJECT_DIR / "OutputFiles/gradient_analytic_0.root"
NUMERIC_FIT = PROJECT_DIR / "OutputFiles/gradient_numeric_0.root"
EMI_FIT = PROJECT_DIR / "OutputFiles/e_rho_fit_0.root"
BRUFIT_DIR = (
    PROJECT_DIR.parents[1]
    / "brufit/tutorials/ElectroAmps/TwoSpin0AmpsFromMoments"
)
BRUFIT_TIMES = BRUFIT_DIR / "brufit_times.csv"
BRUFIT_FIT = BRUFIT_DIR / "resultsGivenMoments_eRhoFirstBin.root"
OUTPUT_DIR = PROJECT_DIR / "AnalysisScripts/outputs/performance"
SHOW_PLOTS = False


def slope_through_origin(x, y):
    """Least-squares slope for the timing model t = slope * starts."""
    return np.dot(x, y) / np.dot(x, x)


def dataframe_columns(frame):
    return {str(name) for name in frame.GetColumnNames()}


def best_row(path, tree_name, objective):
    """Use RDataFrame to select the lowest finite accepted minimum."""
    frame = ROOT.RDataFrame(tree_name, str(path))
    names = dataframe_columns(frame)
    condition = f"std::isfinite({objective})"
    if "fit_ok" in names:
        condition += " && fit_ok"
    if "status" in names:
        condition += " && status == 0"
    frame = frame.Filter(condition)
    minimum = frame.Min(objective).GetValue()
    columns = sorted(
        name for name in names
        if name.startswith(("a_", "b_", "aphi_", "bphi_"))
        and not name.startswith("err__")
    )
    values = frame.Filter(f"{objective} == {minimum:.17g}").Range(1).AsNumpy(
        columns=[objective] + columns
    )
    return minimum, {name: data[0] for name, data in values.items()}


# Timing comparison. This section runs when all three timing tables exist.
if NUMERIC_TIMES.exists() and ANALYTIC_TIMES.exists() and BRUFIT_TIMES.exists():
    numeric_starts, numeric_times = np.loadtxt(
        NUMERIC_TIMES, delimiter=",", unpack=True, skiprows=1
    )
    analytic_starts, analytic_times = np.loadtxt(
        ANALYTIC_TIMES, delimiter=",", unpack=True, skiprows=1
    )
    brufit_starts, brufit_times = np.loadtxt(
        BRUFIT_TIMES, delimiter=",", unpack=True, skiprows=1
    )
    timing_sets = [
        (numeric_starts, numeric_times, "EMI numerical", "C0"),
        (analytic_starts, analytic_times, "EMI analytical", "C1"),
        (brufit_starts, brufit_times, "BruFit", "C2"),
    ]
    model_starts = np.logspace(1, 6, 1001)
    figure, axis = plt.subplots(figsize=(8, 5))
    for starts, times, label, colour in timing_sets:
        slope = slope_through_origin(starts, times)
        axis.plot(starts, times, "o", color=colour, label=label)
        axis.plot(model_starts, slope * model_starts, color=colour)
        print(f"{label}: t = {slope:.7g} N_starts")
    axis.set(xlabel=r"$N_{\mathrm{starts}}$", ylabel="Time (s)")
    axis.legend(frameon=False)
    style_axis(axis)
    save_figure(
        figure, OUTPUT_DIR / "numeric_v_analytic.png", show=SHOW_PLOTS
    )

# Compare the recorded status of analytical and numerical minimizations.
status_columns = ["fit_ok", "hesse_ok", "status", "cov_status", "chi2"]
status_summaries = []
for label, path in (("Analytical", ANALYTIC_FIT), ("Numerical", NUMERIC_FIT)):
    data = ROOT.RDataFrame("fitResults", str(path)).AsNumpy(status_columns)
    data = {name: np.asarray(values) for name, values in data.items()}
    finite = np.isfinite(data["chi2"])
    good_minimum = finite & data["fit_ok"] & (data["status"] == 0)
    good_hessian = good_minimum & data["hesse_ok"] & (data["cov_status"] >= 2)
    counts = [len(data["chi2"]), np.count_nonzero(data["fit_ok"]),
              np.count_nonzero(good_minimum), np.count_nonzero(good_hessian)]
    status_summaries.append(counts)
    print(label, dict(zip(["total", "fit_ok", "minimum", "Hessian"], counts)))

status_summaries = np.asarray(status_summaries)
x = np.arange(4)
figure, axis = plt.subplots(figsize=(8, 5))
axis.bar(x - 0.18, status_summaries[0], 0.36, label="Analytical")
axis.bar(x + 0.18, status_summaries[1], 0.36, label="Numerical")
axis.set_xticks(x, ["Total", "fit_ok", "Good minimum", "Good Hessian"])
axis.set_ylabel("Number of starts")
axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
axis.legend(frameon=False)
style_axis(axis)
save_figure(
    figure, OUTPUT_DIR / "analytic_v_numeric_status.png", show=SHOW_PLOTS
)

# Print the best analytical and numerical amplitudes branch by branch.
analytic_minimum, analytic = best_row(ANALYTIC_FIT, "fitResults", "chi2")
numeric_minimum, numeric = best_row(NUMERIC_FIT, "fitResults", "chi2")
print(f"Analytical minimum chi2: {analytic_minimum:.12g}")
print(f"Numerical minimum chi2:  {numeric_minimum:.12g}")
for name in sorted(analytic.keys() & numeric.keys()):
    if name == "chi2":
        continue
    difference = numeric[name] - analytic[name]
    if name.startswith(("aphi_", "bphi_")):
        difference = wrap_phase(difference)
    print(f"{name:<20} {analytic[name]:14.7g} {numeric[name]:14.7g} {difference:14.7g}")

# Print the common best-fit EMI and BruFit parameters.
emi_minimum, emi = best_row(EMI_FIT, "fitResults", "chi2")
brufit_minimum, brufit = best_row(BRUFIT_FIT, "PartialWaves", "val")
print(f"EMI minimum chi2: {emi_minimum:.12g}")
print(f"BruFit minimum val: {brufit_minimum:.12g}")
for name in sorted(emi.keys() & brufit.keys()):
    if name in {"chi2", "val"}:
        continue
    difference = brufit[name] - emi[name]
    if name.startswith(("aphi_", "bphi_")):
        difference = wrap_phase(difference)
    print(f"{name:<20} {emi[name]:14.7g} {brufit[name]:14.7g} {difference:14.7g}")
