#!/usr/bin/env python3

"""Compare generated and fitted double-polarization amplitudes.

Both nucleon spin indices are observed, removing the continuous spin-basis
ambiguity left in a single-polarization measurement.
"""

import re

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from tutorials.analysis_utils import (
    PROJECT_DIR, best_index, branches, complex_amplitude, phase_branch,
    read_tree, save_figure, style_axis, wave_label,
)


# Settings: normally only the input and fit filenames need changing.
CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/fixed_double"
TRUTH_FILE = CASE_DIR / "output/generated/fixed_both_moments.root"
FIT_FILE = CASE_DIR / "output/root/fixed_both_fit.root"
TRUTH_TREE = "genMoments"
FIT_TREE = "fitResults"
SHOW_PLOT = False

OUTPUT_DIR = (
    CASE_DIR / "output/plots"
)


# Find amplitude magnitude branches that are present in both files.
pattern = re.compile(r"^[ab]_[TL]_\d+_m?\d+_(?:1|m1)$")
truth_branches = set(branches(TRUTH_FILE, TRUTH_TREE))
fit_branches = set(branches(FIT_FILE, FIT_TREE))
magnitudes = sorted(
    name for name in truth_branches & fit_branches if pattern.fullmatch(name)
)
phases = [phase_branch(name) for name in magnitudes]

# The generated tree has one entry. The fit tree contains many random starts.
truth = read_tree(TRUTH_FILE, TRUTH_TREE, magnitudes + phases)
status = [name for name in ("chi2", "fit_ok", "status") if name in fit_branches]
fit = read_tree(FIT_FILE, FIT_TREE, status + magnitudes + phases)
entry = best_index(fit)

truth_amplitudes = {
    name: complex_amplitude(truth[name][0], truth[phase_branch(name)][0])
    for name in magnitudes
}
fit_amplitudes = {
    name: complex_amplitude(fit[name][entry], fit[phase_branch(name)][entry])
    for name in magnitudes
}

# Print the numerical closure before drawing it.
print(f"Best fit entry: {entry}, chi2 = {fit['chi2'][entry]:.7g}")
print(
    f"{'amplitude':<22} {'truth Re':>11} {'truth Im':>11} "
    f"{'fit Re':>11} {'fit Im':>11} {'|difference|':>13}"
)
for name in magnitudes:
    generated = truth_amplitudes[name]
    fitted = fit_amplitudes[name]
    print(
        f"{name:<22} {generated.real:11.6f} {generated.imag:11.6f} "
        f"{fitted.real:11.6f} {fitted.imag:11.6f} "
        f"{abs(fitted - generated):13.6g}"
    )

# Show each reflectivity and photon-polarization sector on its own Argand axis.
groups = [
    (("a", "T"), "Natural transverse"),
    (("b", "T"), "Unnatural transverse"),
    (("a", "L"), "Natural longitudinal"),
    (("b", "L"), "Unnatural longitudinal"),
]
groups = [
    (key, title) for key, title in groups
    if any(tuple(name.split("_")[:2]) == key for name in magnitudes)
]

figure, axes = plt.subplots(
    len(groups), 1, figsize=(9, 6 * len(groups)),
    constrained_layout=True, squeeze=False,
)
axes = axes[:, 0]
extent = 0.1

for axis, (group, title) in zip(axes, groups):
    selected = [
        name for name in magnitudes if tuple(name.split("_")[:2]) == group
    ]
    for colour_index, name in enumerate(selected):
        generated = truth_amplitudes[name]
        fitted = fit_amplitudes[name]
        colour = f"C{colour_index % 10}"
        axis.plot(
            fitted.real, fitted.imag, "o", color=colour,
            markersize=11, label=wave_label("_".join(name.split("_")[:-1])),
        )
        axis.plot(
            generated.real, generated.imag, "*", color=colour,
            markeredgecolor="black", markersize=17,
        )
        extent = max(extent, abs(generated.real), abs(generated.imag),
                     abs(fitted.real), abs(fitted.imag))

    axis.set_title(title)
    axis.set_ylabel("Im")
    axis.set_aspect("equal", adjustable="box")
    style_axis(axis, grid=True)
    handles, _ = axis.get_legend_handles_labels()
    handles.extend([
        Line2D([], [], marker="o", color="black", linestyle="none",
               label="Minimizer"),
        Line2D([], [], marker="*", color="white", markeredgecolor="black",
               linestyle="none", label="Set amplitude"),
    ])
    axis.legend(handles=handles, frameon=False, ncol=3)

limit = 1.15 * extent
for axis in axes:
    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
axes[-1].set_xlabel("Re")

figure.suptitle("Polarized fixed-amplitude closure")
save_figure(
    figure, OUTPUT_DIR / f"{FIT_FILE.stem}_argand.pdf",
    show=SHOW_PLOT, dpi=600, transparent=True,
)
