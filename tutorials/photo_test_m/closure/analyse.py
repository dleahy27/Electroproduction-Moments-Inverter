#!/usr/bin/env python3
"""Plot every moment and validate mass, width, magnitude, and closure residuals.

The single-wave phase is unobservable; its gauge value is deliberately excluded
from the recovery check. Select the best accepted start, as the C++ writer does.
"""

import argparse
import math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from tutorials.analysis_utils import save_figure, style_axis, range_errors
from tutorials.photo_test_m.analysis import load_closure
from tutorials.paths import PROJECT_DIR


def plot_coordinates(data, plots, check):
    """Plot injected and recovered coordinates; check observable parameters."""
    coordinates, results, best = data["coordinates"], data["results"], data["best"]
    for suffix, ylabel, filename in (
        ("_mass", "mass [GeV]", "masses"),
        ("_width", "width [GeV]", "widths"),
        ("_mag", "coupling magnitude", "magnitudes"),
        ("_phase", "coupling phase [rad]", "phases"),
    ):
        selected = [name for name in coordinates if name.endswith(suffix)]
        injected = np.array([data["truth"][name] for name in selected])
        fitted = np.array([results[name][best] for name in selected])
        x = np.arange(len(selected))
        figure, axis = plt.subplots(figsize=(max(7, 0.35 * len(selected)), 4))
        axis.scatter(x, injected, marker="x", label="injected")
        axis.errorbar(
            x, fitted,
            yerr=range_errors(fitted, np.column_stack(
                [results[name][data["top"]] for name in selected]
            )),
            fmt="_", ms=11, capsize=3, label="best fit; best 0.1% range",
        )
        axis.set(xticks=x, xticklabels=selected, ylabel=ylabel)
        axis.tick_params(axis="x", rotation=90)
        axis.legend()
        style_axis(axis)
        save_figure(figure, plots / f"{filename}.pdf")
        for name, target, value in zip(selected, injected, fitted):
            print(f"{name}: {value:.9g} (injected {target:.9g})")
        # In the full two-k truth-start closure both models use the same phase
        # gauge. The one-wave closure deliberately uses a different gauge.
        if check and suffix == "_phase" and any("_k1_" in name for name in selected):
            difference = np.angle(np.exp(1j * (fitted - injected)))
            if np.max(np.abs(difference)) > 1.0e-4:
                raise AssertionError("Full PhotoTest coupling phases do not close")
        if (
            check
            and suffix != "_phase"
            and not np.allclose(fitted, injected, atol=1.0e-5, rtol=1.0e-4)
        ):
            raise AssertionError(f"Recovered {suffix} coordinates do not close")


def plot_moments(spectra, data, model, plots):
    """Plot every common moment in pages of nine and return the largest pull."""
    max_pull = 0.0
    for first in range(0, len(spectra), 9):
        page = spectra[first : first + 9]
        figure, axes = plt.subplots(
            math.ceil(len(page) / 3),
            3,
            figsize=(12, 3.2 * math.ceil(len(page) / 3)),
            squeeze=False,
        )
        for axis, name in zip(axes.flat, page):
            prediction = model["model__" + name]
            axis.errorbar(
                data["mass_GeV"],
                data[name],
                yerr=data[name + "_err"],
                fmt="o",
                ms=3,
                label="input",
            )
            axis.plot(model["mass_GeV"], prediction, label="fit")
            axis.set(xlabel="mass [GeV]", ylabel=name)
            style_axis(axis)
            max_pull = max(
                max_pull, np.max(np.abs(prediction - data[name]) / data[name + "_err"])
            )
        for axis in axes.flat[len(page) :]:
            axis.set_visible(False)
        axes.flat[0].legend()
        save_figure(figure, plots / f"input_vs_model_{first // 9 + 1:02d}.pdf")
    return max_pull


def plot_intensity(spectra, data, model, plots):
    """Show the intensity and its fitted-minus-input residual."""
    intensity = "RH_0_0_0" if "RH_0_0_0" in spectra else "RH_0_0_0_0_0"
    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    axes[0].plot(data["mass_GeV"], data[intensity], "x", label="input")
    axes[0].plot(model["mass_GeV"], model["model__" + intensity], label="fit")
    axes[0].set_ylabel(intensity)
    axes[0].legend()
    axes[1].plot(data["mass_GeV"], model["model__" + intensity] - data[intensity])
    axes[1].set(xlabel="mass [GeV]", ylabel="model − input")
    save_figure(figure, plots / "intensity_and_residuals.pdf")


def plot_starts(results, plots):
    """Keep accepted and failed starts visible in the minimizer landscape."""
    figure, axis = plt.subplots(figsize=(7, 4))
    valid = results["fit_ok"].astype(bool)
    for accepted, label in ((True, "accepted"), (False, "failed")):
        mask = valid == accepted
        axis.scatter(results["start"][mask], results["chi2"][mask], label=label)
    axis.set(xlabel="global start index", ylabel="chi2")
    axis.legend()
    save_figure(figure, plots / "start_landscape.pdf")


def analyse(case, check=False):
    """Load once, draw each plot family, and optionally enforce closure."""
    plots = case / "output/plots"
    data = load_closure(case)
    results, best = data["results"], data["best"]
    print(f"Best start: {results['start'][best]}; chi2 = {results['chi2'][best]:.8g}")
    plot_coordinates(data, plots, check)
    moments, model, spectra = data["moments"], data["predictions"], data["spectra"]
    if not np.array_equal(moments["mass_GeV"], model["mass_GeV"]):
        raise AssertionError("Prediction mass grid differs from input")
    max_pull = plot_moments(spectra, moments, model, plots)
    plot_intensity(spectra, moments, model, plots)
    plot_starts(results, plots)
    print(f"Largest absolute moment pull: {max_pull:.8g}")
    if check and (not np.isfinite(max_pull) or max_pull > 1.0e-3):
        raise AssertionError("Moment predictions do not close")
    print(f"Plots: {plots}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="fail if parameter or moment closure fails"
    )
    parser.add_argument(
        "--case", type=Path, default=PROJECT_DIR / "tutorials/photo_test_m/closure"
    )
    options = parser.parse_args()
    analyse(options.case, options.check)


if __name__ == "__main__":
    main()
