#!/usr/bin/env python3

"""Plot best-fit Hessian estimates of R = sigma_L / sigma_T.

The selected row is the lowest accepted chi-square start; its R uncertainty is
the covariance-propagated value written by the C++ fitter.
"""

import matplotlib.pyplot as plt
import numpy as np

from tutorials.analysis_utils import PROJECT_DIR, best_index, branches, read_tree, save_figure, style_axis
from tutorials.electroproduction_paper.datasets import CHANNELS, MAIN_CHANNEL_KEYS


# Settings: channel values are shared with the corresponding bootstrap analysis.
OUTPUT_FILES = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/root"
OUTPUT_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/plots/R_hessian"
TREE_NAME = "fitResults"
SHOW_PLOT = False

CHANNELS_TO_PLOT = [CHANNELS[key] for key in MAIN_CHANNEL_KEYS]


# Read R and its propagated Hessian uncertainty from the global minimum.
results = {}
for channel in CHANNELS_TO_PLOT:
    stem = channel.key
    rows = []
    for bin_index, q2 in enumerate(channel.q2):
        path = OUTPUT_FILES / f"{stem}_fit_{bin_index}.root"
        if not path.exists():
            continue
        names = set(branches(path, TREE_NAME))
        columns = ["chi2", "R", "err__R"]
        columns += [name for name in ("fit_ok", "status") if name in names]
        data = read_tree(path, TREE_NAME, columns)
        entry = best_index(data)
        rows.append((bin_index, q2, data["R"][entry], data["err__R"][entry]))
        print(
            f"{stem} bin {bin_index}: R = {data['R'][entry]:.6g} "
            f"+/- {data['err__R'][entry]:.6g}"
        )
    results[stem] = rows

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

    axis.errorbar(q2, measured, yerr=measured_error, fmt="o", capsize=4,
                  label="Minimizer")
    axis.errorbar(q2, published[:, 0], yerr=published[:, 1], fmt="s",
                  capsize=4, label=channel.published_r_source)
    difference_axis.errorbar(
        q2, published[:, 0] - measured,
        yerr=np.hypot(published[:, 1], measured_error),
        fmt="o", capsize=4, color="black",
    )
    difference_axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1)
    axis.set(title=channel.title, ylim=(0, 2))
    axis.tick_params(labelbottom=False)
    axis.legend(frameon=False)
    difference_axis.set(xlabel=r"$Q^2\ [\mathrm{GeV}^2]$", ylim=(-1, 1))
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
