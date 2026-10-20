#!/usr/bin/env python3
"""Validate the effective single-k PhotoTest mass-dependent fits."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tutorials.analysis_utils import save_figure, style_axis, range_errors
from tutorials.photo_test_m.analysis import load_mass_scan, pole_names, coupling_name, coupling_sector_name, moment_indices
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
        curve = case["curves"]
        top.plot(truth["mass_GeV"], truth[intensity], "x", ms=4, label=f"{label} truth")
        colour = top.lines[-1].get_color()
        top.plot(curve["mass_GeV"], curve["model__" + intensity], "-",
                 color=colour, lw=1.4, label=f"{label} model")
        residual.plot(truth["mass_GeV"], model - truth[intensity], "-", color=colour, lw=1.2)
    top.set_ylabel(r"$H^{0}(00)$")
    residual.set(xlabel="mass [GeV]", ylabel="model − input")
    top.legend(fontsize=6, ncol=2)
    style_axis(top)
    style_axis(residual)
    save_figure(figure, PLOTS / "intensity_and_residuals.pdf")


def moment_plots(case):
    """One file per (L,M), showing exact points and densely evaluated BWs."""
    groups = {}
    for name in case["moment_names"]:
        alpha, l, m = moment_indices(name)
        groups.setdefault((l, m), {})[alpha] = name
    truth, curves = case["moments"], case["curves"]
    order = np.argsort(curves["mass_GeV"])
    for (l, m), moments in sorted(groups.items()):
        alphas = (0, 1) if m == 0 else (0, 1, 2, 3)
        missing = set(alphas) - moments.keys()
        if missing:
            raise RuntimeError(f"Missing alpha={sorted(missing)} for L={l}, M={m}")
        figure, axes = plt.subplots(1 if m == 0 else 2, 2, figsize=(10, 3.8 if m == 0 else 7),
                                    constrained_layout=True, squeeze=False)
        for axis, alpha in zip(axes.flat, alphas):
            name = moments[alpha]
            axis.plot(truth["mass_GeV"], truth[name], "o", ms=4, color="C0",
                      linestyle="none", label="generated points", zorder=3)
            axis.plot(curves["mass_GeV"][order], curves["model__" + name][order],
                      color="C1", lw=1.5, label="best-fit mass model")
            axis.set(xlabel="mass [GeV]", ylabel=rf"$H^{{{alpha}}}({l}{m})$")
            style_axis(axis)
        axes.flat[0].legend(fontsize=8)
        figure.suptitle(f"s={case['asymmetry']:+.2f}; L={l}, M={m}")
        save_figure(figure, PLOTS / "moments" / f"{case['tag']}_moments_L{l}_M{m}.pdf")


def pole_plot(case):
    """Compare the best accepted start with the injected resonance parameters."""
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
    axes[0].set_ylim(-0.05, 2)
    axes[0].set_title(f"s={case['asymmetry']:+.2f}")

    axes[1].set_ylim(-3.2, 3.2)
    axes[-1].set(
        xticks=x, xticklabels=coupling_name(selected),xlabel="coupling", title=f"s={case['asymmetry']:+.2f}"
    )
    axes[-1].tick_params(axis="x", rotation=90)
    save_figure(figure, PLOTS / f"{tag}_couplings_5.pdf")


def main():
    global PLOTS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE,
                        help="directory containing scan.root")
    options = parser.parse_args()
    PLOTS = options.output / "plots"
    cases = load_mass_scan(options.output)
    intensity_plot(cases)
    for case in cases:
        pole_plot(case)
        coupling_plot(case)
        moment_plots(case)


if __name__ == "__main__":
    main()
