#!/usr/bin/env python3
"""Validate the effective single-k PhotoTest mass-dependent fits."""

import matplotlib.pyplot as plt
import numpy as np
from tutorials.analysis_utils import save_figure, style_axis, range_errors
from tutorials.photo_test_m.analysis import load_mass_scan, pole_names, coupling_name, coupling_sector_name
from tutorials.paths import PROJECT_DIR

CASE = PROJECT_DIR / "tutorials/photo_test_m/single_k_test/output"
PLOTS = CASE / "plots"


def intensity_plot(cases):
    """Overlay the unpolarized intensity and its raw prediction residual."""
    figure, (top, residual) = plt.subplots(
        2, 1, sharex=True, figsize=(8, 6), gridspec_kw={"height_ratios": [3, 1]}
    )
    for case in cases:
        tag = case["tag"]
        label = f"s={case['asymmetry']:+.2f}"
        intensity = case["intensity"]
        truth = case["moments"]
        prediction = case["predictions"]
        if not np.array_equal(prediction["mass_GeV"], truth["mass_GeV"]):
            raise RuntimeError(f"Prediction mass grid does not match truth for {tag}")
        model = prediction["model__" + intensity]
        top.plot(truth["mass_GeV"], truth[intensity], "x", ms=4, label=f"{label} truth")
        colour = top.lines[-1].get_color()
        top.plot(truth["mass_GeV"], model, "-", color=colour, lw=1.4, label=f"{label} model")
        residual.plot(truth["mass_GeV"], model - truth[intensity], "-", color=colour, lw=1.2)
    top.set_ylabel(r"$H^{0}(00)$")
    residual.set(xlabel="mass [GeV]", ylabel="model − input")
    top.legend(fontsize=6, ncol=2)
    style_axis(top)
    style_axis(residual)
    save_figure(figure, PLOTS / "intensity_and_residuals.pdf")


def pole_plot(case):
    """Compare the accepted global minimum with the injected resonance poles."""
    tag = case["tag"]
    poles = [name for name in case["coordinates"] if name.endswith(("_mass", "_width"))]
    values = case["fit"]
    best = 0  # Selected rows are sorted by chi2, with the best fit first.
    x = np.arange(len(poles))
    figure, axis = plt.subplots(figsize=(max(8, 0.85 * len(poles)), 4.5))
    axis.scatter(
        x, [case["truth"][name] for name in poles], marker="x", s=52, label="truth"
    )
    axis.errorbar(
        x,
        [values[name][best] for name in poles],
        yerr=range_errors(
            [values[name][best] for name in poles],
            np.column_stack([values[name] for name in poles]),
        ),
        fmt="_", capsize=3,
        ms=17,
        mew=2,
        label="best fit; best 5% range",
    )
    axis.set(xticks=x, xticklabels=pole_names(poles), ylabel="value [GeV]", ylim=(0,2.55), title=f"s={case['asymmetry']:+.2f}")
    axis.tick_params(axis="x", rotation=45)
    axis.legend()
    style_axis(axis)
    save_figure(figure, PLOTS / f"{tag}_masses_widths_5.pdf")


def coupling_plot(case):
    """Show effective k=0 coordinates beside both distinct injected sectors."""
    tag = case["tag"]
    coupling_names = [
        name for name in case["coordinates"] if name.endswith(("_mag", "_phase"))
    ]
    values = case["fit"]
    best = 0
    figure, axes = plt.subplots(
        2, 1, figsize=(max(10, 0.32 * len(coupling_names)), 7), sharex=True
    )
    for axis, suffix, ylabel in zip(
        axes, ("_mag", "_phase"), ("coupling magnitude", "coupling phase [rad]")
    ):
        selected = [name for name in coupling_names if name.endswith(suffix)]
        x = np.arange(len(selected))
        plus = [case["truth"][coupling_sector_name(name, 1)] for name in selected]
        minus = [case["truth"][coupling_sector_name(name, -1)] for name in selected]
        fitted = [values[name][best] for name in selected]
        axis.scatter(x, plus, marker="^", s=26, label=r"input $k=+1$")
        axis.scatter(x, minus, marker="v", s=26, label=r"input $k=-1$")
        axis.errorbar(
            x, fitted,
            yerr=range_errors(fitted, np.column_stack([values[name] for name in selected])),
            fmt="o", ms=4, capsize=3,
            label=r"effective $k$ fit; best 5% range",
        )

        axis.set_ylabel(ylabel)

        style_axis(axis)
    axes[0].legend(fontsize=7, ncol=3)
    axes[0].set_ylim(-0.05,2)
    axes[0].set_title(f"s={case['asymmetry']:+.2f}")

    axes[1].set_ylim(-3.2,3.2)
    axes[-1].set(
        xticks=x, xticklabels=coupling_name(selected),xlabel="coupling", title=f"s={case['asymmetry']:+.2f}"
    )
    axes[-1].tick_params(axis="x", rotation=90)
    save_figure(figure, PLOTS / f"{tag}_couplings_5.pdf")
def main():
    cases = load_mass_scan(CASE)
    intensity_plot(cases)
    for case in cases:
        pole_plot(case)
        coupling_plot(case)


if __name__ == "__main__":
    main()
