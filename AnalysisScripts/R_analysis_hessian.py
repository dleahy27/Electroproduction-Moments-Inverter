#!/usr/bin/env python3

"""Plot Hessian estimates of R = sigma_L / sigma_T."""

import matplotlib.pyplot as plt
import numpy as np

from analysis_utils import PROJECT_DIR, best_index, branches, read_tree, save_figure, style_axis


# Settings: values match the corresponding bootstrap R analysis.
OUTPUT_FILES = PROJECT_DIR / "OutputFiles"
OUTPUT_DIR = PROJECT_DIR / "AnalysisScripts/outputs/R_analysis_hessian"
TREE_NAME = "fitResults"
SHOW_PLOT = False

CHANNELS = [
    ("e_rho", r"$e^-\rho^0$", np.array([0.82, 1.19, 1.66, 3.06]),
     np.array([[0.649, 0.188], [0.671, 0.065], [0.794, 0.083], [1.068, 0.087]]),
     "Published HERMES"),
    ("mu_rho", r"$\mu^-\rho^0$", np.array([1.14, 1.60, 2.80, 6.02]),
     np.array([[0.724, 0.078], [0.930, 0.081], [1.227, 0.122], [1.200, 0.729]]),
     "Published COMPASS"),
    ("e_omega", r"$e^-\omega$", np.array([1.28, 2.00, 4.00]),
     np.array([[0.226, 0.066], [0.237, 0.091], [0.283, 0.091]]),
     "Published HERMES"),
    ("mu_omega", r"$\mu^-\omega$", np.array([1.16, 1.64, 3.58]),
     np.array([[0.475, 0.085], [0.495, 0.124], [0.739, 0.197]]),
     "Published COMPASS"),
]


# Read R and its propagated Hessian uncertainty from the global minimum.
results = {}
for stem, _, q2_values, _, _ in CHANNELS:
    rows = []
    for bin_index, q2 in enumerate(q2_values):
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

channels = [channel for channel in CHANNELS if results[channel[0]]]
figure = plt.figure(figsize=(6 * len(channels), 9))
grid = figure.add_gridspec(
    2, len(channels), height_ratios=[3.2, 1.25], hspace=0.08, wspace=0.12
)
top_axes = [figure.add_subplot(grid[0, index]) for index in range(len(channels))]
difference_axes = [
    figure.add_subplot(grid[1, index]) for index in range(len(channels))
]

for index, (stem, title, _, published_all, source) in enumerate(channels):
    rows = results[stem]
    bins = np.array([row[0] for row in rows])
    q2 = np.array([row[1] for row in rows])
    measured = np.array([row[2] for row in rows])
    measured_error = np.array([row[3] for row in rows])
    published = published_all[bins]
    axis = top_axes[index]
    difference_axis = difference_axes[index]

    axis.errorbar(q2, measured, yerr=measured_error, fmt="o", capsize=4,
                  label="Minimizer")
    axis.errorbar(q2, published[:, 0], yerr=published[:, 1], fmt="s",
                  capsize=4, label=source)
    difference_axis.errorbar(
        q2, published[:, 0] - measured,
        yerr=np.hypot(published[:, 1], measured_error),
        fmt="o", capsize=4, color="black",
    )
    difference_axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1)
    axis.set(title=title, ylim=(0, 2))
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
