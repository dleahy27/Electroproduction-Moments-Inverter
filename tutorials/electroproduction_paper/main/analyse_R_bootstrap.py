#!/usr/bin/env python3

"""Plot bootstrap estimates of R = sigma_L / sigma_T.

Each point is the mean and sample width over accepted bootstrap toys. Older
files without an R branch derive the same ratio lazily from amplitude squares.
"""

import matplotlib.pyplot as plt
import numpy as np
import ROOT

from tutorials.analysis_utils import (
    PROJECT_DIR, MAGNITUDE_PATTERN, r_expression, save_figure, style_axis,
)
from tutorials.electroproduction_paper.datasets import CHANNELS, MAIN_CHANNEL_KEYS


# Settings: channel values are shared with the Hessian analysis.
OUTPUT_FILES = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/root"
OUTPUT_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/plots/R_bootstrap"
TREE_NAME = "PartialWaves"
SHOW_PLOT = False

CHANNELS_TO_PLOT = [CHANNELS[key] for key in MAIN_CHANNEL_KEYS]


# RDataFrame evaluates R directly and returns only its NumPy values.
results = {}
for channel in CHANNELS_TO_PLOT:
    rows = []
    for bin_index, q2 in enumerate(channel.q2):
        path = OUTPUT_FILES / f"{channel.key}_bootstrap_{bin_index}.root"
        if not path.exists():
            continue
        frame = ROOT.RDataFrame(TREE_NAME, str(path))
        names = {str(name) for name in frame.GetColumnNames()}
        # Reject failed toys before deriving R; failed amplitudes can be finite
        # and would otherwise contaminate the bootstrap width.
        if "valid" in names:
            frame = frame.Filter("valid != 0")
        if "status" in names:
            frame = frame.Filter("status == 0")
        if "R" not in names:
            magnitudes = [name for name in names if MAGNITUDE_PATTERN.fullmatch(name)]
            frame = frame.Define("R", r_expression(magnitudes))
        values = np.asarray(frame.AsNumpy(columns=["R"])["R"])
        values = values[np.isfinite(values)]
        rows.append((bin_index, q2, np.mean(values), np.std(values, ddof=1)))
        print(
            f"{channel.key} bin {bin_index}: "
            f"R = {np.mean(values):.6g} +/- {np.std(values, ddof=1):.6g}"
        )
    results[channel.key] = rows

# Retain only channels for which at least one bootstrap file was found.
channels = [channel for channel in CHANNELS_TO_PLOT if results[channel.key]]
figure = plt.figure(figsize=(6 * len(channels), 9))
grid = figure.add_gridspec(
    2, len(channels), height_ratios=[3.2, 1.25], hspace=0.08, wspace=0.12
)
top_axes = [figure.add_subplot(grid[0, index]) for index in range(len(channels))]
difference_axes = [
    figure.add_subplot(grid[1, index]) for index in range(len(channels))
]

for index, channel in enumerate(channels):
    rows = results[channel.key]
    bins = np.array([row[0] for row in rows])
    q2 = np.array([row[1] for row in rows])
    measured = np.array([row[2] for row in rows])
    measured_error = np.array([row[3] for row in rows])
    published = np.asarray(channel.published_r)[bins]
    axis = top_axes[index]
    difference_axis = difference_axes[index]

    axis.errorbar(
        q2, measured, yerr=measured_error, fmt="o", markersize=8,
        capsize=4, label="Minimizer",
    )
    axis.errorbar(
        q2, published[:, 0], yerr=published[:, 1], fmt="s", markersize=7,
        capsize=4, label=channel.published_r_source,
    )

    difference = published[:, 0] - measured
    difference_error = np.hypot(published[:, 1], measured_error)
    difference_axis.errorbar(
        q2, difference, yerr=difference_error, fmt="o", capsize=4,
        color="black",
    )
    difference_axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1)

    axis.set_title(channel.title)
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
