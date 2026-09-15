#!/usr/bin/env python3

"""Plot bootstrap estimates of R = sigma_L / sigma_T."""

import matplotlib.pyplot as plt
import numpy as np
import ROOT

from analysis_utils import (
    PROJECT_DIR, MAGNITUDE_PATTERN, r_expression, save_figure, style_axis,
)


# Settings: each channel gives its file stem, Q2 bins, and published values.
OUTPUT_FILES = PROJECT_DIR / "OutputFiles"
OUTPUT_DIR = PROJECT_DIR / "AnalysisScripts/outputs/R_analysis_bootstrap"
TREE_NAME = "PartialWaves"
SHOW_PLOT = False

CHANNELS = [
    {
        "stem": "e_rho", "title": r"$e^-\rho^0$",
        "q2": np.array([0.82, 1.19, 1.66, 3.06]),
        "published": np.array([[0.649, 0.188], [0.671, 0.065],
                               [0.794, 0.083], [1.068, 0.087]]),
        "source": "Published HERMES",
    },
    {
        "stem": "mu_rho", "title": r"$\mu^-\rho^0$",
        "q2": np.array([1.14, 1.60, 2.80, 6.02]),
        "published": np.array([[0.724, 0.078], [0.930, 0.081],
                               [1.227, 0.122], [1.200, 0.729]]),
        "source": "Published COMPASS",
    },
    {
        "stem": "e_omega", "title": r"$e^-\omega$",
        "q2": np.array([1.28, 2.00, 4.00]),
        "published": np.array([[0.226, 0.066], [0.237, 0.091],
                               [0.283, 0.091]]),
        "source": "Published HERMES",
    },
    {
        "stem": "mu_omega", "title": r"$\mu^-\omega$",
        "q2": np.array([1.16, 1.64, 3.58]),
        "published": np.array([[0.475, 0.085], [0.495, 0.124],
                               [0.739, 0.197]]),
        "source": "Published COMPASS",
    },
]


# RDataFrame evaluates R directly and returns only its NumPy values.
results = {}
for channel in CHANNELS:
    rows = []
    for bin_index, q2 in enumerate(channel["q2"]):
        path = OUTPUT_FILES / f"{channel['stem']}_bootstrap_{bin_index}.root"
        if not path.exists():
            continue
        frame = ROOT.RDataFrame(TREE_NAME, str(path))
        names = {str(name) for name in frame.GetColumnNames()}
        if "R" not in names:
            magnitudes = [name for name in names if MAGNITUDE_PATTERN.fullmatch(name)]
            frame = frame.Define("R", r_expression(magnitudes))
        values = np.asarray(frame.AsNumpy(columns=["R"])["R"])
        values = values[np.isfinite(values)]
        rows.append((bin_index, q2, np.mean(values), np.std(values, ddof=1)))
        print(
            f"{channel['stem']} bin {bin_index}: "
            f"R = {np.mean(values):.6g} +/- {np.std(values, ddof=1):.6g}"
        )
    results[channel["stem"]] = rows

# Retain only channels for which at least one bootstrap file was found.
channels = [channel for channel in CHANNELS if results[channel["stem"]]]
figure = plt.figure(figsize=(6 * len(channels), 9))
grid = figure.add_gridspec(
    2, len(channels), height_ratios=[3.2, 1.25], hspace=0.08, wspace=0.12
)
top_axes = [figure.add_subplot(grid[0, index]) for index in range(len(channels))]
difference_axes = [
    figure.add_subplot(grid[1, index]) for index in range(len(channels))
]

for index, channel in enumerate(channels):
    rows = results[channel["stem"]]
    bins = np.array([row[0] for row in rows])
    q2 = np.array([row[1] for row in rows])
    measured = np.array([row[2] for row in rows])
    measured_error = np.array([row[3] for row in rows])
    published = channel["published"][bins]
    axis = top_axes[index]
    difference_axis = difference_axes[index]

    axis.errorbar(
        q2, measured, yerr=measured_error, fmt="o", markersize=8,
        capsize=4, label="Minimizer",
    )
    axis.errorbar(
        q2, published[:, 0], yerr=published[:, 1], fmt="s", markersize=7,
        capsize=4, label=channel["source"],
    )

    difference = published[:, 0] - measured
    difference_error = np.hypot(published[:, 1], measured_error)
    difference_axis.errorbar(
        q2, difference, yerr=difference_error, fmt="o", capsize=4,
        color="black",
    )
    difference_axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1)

    axis.set_title(channel["title"])
    axis.set_ylim(0, 2)
    axis.tick_params(labelbottom=False)
    axis.legend(frameon=False)
    difference_axis.set_xlabel(r"$Q^2\ [\mathrm{GeV}^2]$")
    difference_axis.set_ylim(-1, 1)
    for current_axis in (axis, difference_axis):
        style_axis(current_axis)
    if index:
        axis.tick_params(labelleft=False)
        difference_axis.tick_params(labelleft=False)

top_axes[0].set_ylabel(r"$R$")
difference_axes[0].set_ylabel(r"$\Delta R$")
save_figure(
    figure, OUTPUT_DIR / "R_scaling.pdf",
    show=SHOW_PLOT, transparent=True,
)
