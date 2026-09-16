#!/usr/bin/env python3

"""Compare a fixed generated point with fitted amplitudes and moments.

The generated tree resolves k=+1 and k=-1, while the unpolarized fit contains
one effective amplitude. Moment closure is therefore the direct observable
comparison; the Argand plot also exposes how the two truth sectors combine.
"""

import re

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from tutorials.analysis_utils import (
    PROJECT_DIR, best_index, branches, complex_amplitude, phase_branch,
    read_tree, save_figure, style_axis, wave_label,
)


# Settings: edit these paths to select a different closure test.
CASE_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main"
TRUTH_FILE = CASE_DIR / "output/generated/fixed_photomoments.root"
FIT_FILE = CASE_DIR / "output/root/fixed_photomoments_unpolarized_fit.root"
TRUTH_TREE = "genMoments"
FIT_TREE = "fitResults"
SHOW_PLOTS = False

OUTPUT_DIR = CASE_DIR / "output/plots/fixed_amplitudes"


# The fit contains one magnitude per wave. The generated file contains the
# corresponding k=+1 and k=-1 amplitudes.
fit_branches = set(branches(FIT_FILE, FIT_TREE))
truth_branches = set(branches(TRUTH_FILE, TRUTH_TREE))
wave_pattern = re.compile(r"^[ab]_[TL]_\d+_(?:m?\d+)$")
waves = sorted(name for name in fit_branches if wave_pattern.fullmatch(name))
fit_phases = [phase_branch(wave) for wave in waves]
truth_columns = [
    branch
    for wave in waves
    for sector in ("1", "m1")
    for branch in (f"{wave}_{sector}", f"{phase_branch(wave)}_{sector}")
]

# Read the truth entry and select the converged fit with the lowest chi squared.
truth = read_tree(TRUTH_FILE, TRUTH_TREE, truth_columns)
status = [name for name in ("chi2", "fit_ok", "status") if name in fit_branches]
fit = read_tree(FIT_FILE, FIT_TREE, status + waves + fit_phases)
entry = best_index(fit)

fit_amplitudes = {
    wave: complex_amplitude(fit[wave][entry], fit[phase_branch(wave)][entry])
    for wave in waves
}
truth_amplitudes = {
    wave: [
        complex_amplitude(
            truth[f"{wave}_{sector}"][0],
            truth[f"{phase_branch(wave)}_{sector}"][0],
        )
        for sector in ("1", "m1")
    ]
    for wave in waves
}

print(f"Best fit entry: {entry}, chi2 = {fit['chi2'][entry]:.7g}")

# Plot the fitted amplitude and the two generated k sectors in the Argand plane.
groups = [
    (("a", "T"), "Natural transverse"),
    (("b", "T"), "Unnatural transverse"),
    (("a", "L"), "Natural longitudinal"),
    (("b", "L"), "Unnatural longitudinal"),
]
groups = [
    (key, title) for key, title in groups
    if any(tuple(wave.split("_")[:2]) == key for wave in waves)
]
figure, axes = plt.subplots(
    len(groups), 1, figsize=(8, 6 * len(groups)),
    sharex=True, constrained_layout=True, squeeze=False,
)
axes = axes[:, 0]

for axis, (group, title) in zip(axes, groups):
    selected = [wave for wave in waves if tuple(wave.split("_")[:2]) == group]
    for colour_index, wave in enumerate(selected):
        colour = f"C{colour_index % 10}"
        fitted = fit_amplitudes[wave]
        generated = truth_amplitudes[wave]
        axis.plot(
            fitted.real, fitted.imag, "o", color=colour,
            markersize=12, label=wave_label(wave),
        )
        axis.plot(
            [value.real for value in generated],
            [value.imag for value in generated],
            "*", color=colour, markeredgecolor="black", markersize=16,
        )

    axis.set(xlim=(-1, 1), ylim=(-1, 1), ylabel="Im", title=title)
    axis.set_aspect("equal", adjustable="box")
    style_axis(axis, grid=True)
    handles, _ = axis.get_legend_handles_labels()
    handles.append(
        Line2D([], [], marker="*", color="white", markeredgecolor="black",
               linestyle="none", label=r"Set amplitude ($k=\pm1$)")
    )
    axis.legend(handles=handles, frameon=False, loc="upper left")
axes[-1].set_xlabel("Re")

save_figure(
    figure, OUTPUT_DIR / "fixed_amplitudes.pdf",
    show=SHOW_PLOTS, dpi=600, transparent=True,
)


# Compare fitted moments with the unpolarized projection of the generated
# tensor moments. RH_alpha_0_0_L_M maps to H_alpha_L_M; alpha=0 is H04.
truth_moment_pattern = re.compile(r"^RH_(\d+)_0_0_(\d+)_(m?\d+)$")
moment_pairs = []
for truth_name in truth_branches:
    match = truth_moment_pattern.fullmatch(truth_name)
    if not match:
        continue
    alpha, ell, projection = match.groups()
    fit_name = f"H04_{ell}_{projection}" if alpha == "0" else f"H_{alpha}_{ell}_{projection}"
    if fit_name in fit_branches:
        moment_pairs.append((truth_name, fit_name))
moment_pairs.sort()

truth_moments = read_tree(
    TRUTH_FILE, TRUTH_TREE, [truth_name for truth_name, _ in moment_pairs]
)
fit_moments = read_tree(
    FIT_FILE, FIT_TREE, [fit_name for _, fit_name in moment_pairs]
)
generated = np.array([truth_moments[name][0] for name, _ in moment_pairs])
fitted = np.array([fit_moments[name][entry] for _, name in moment_pairs])
labels = [fit_name for _, fit_name in moment_pairs]
x = np.arange(len(labels))

figure, axes = plt.subplots(
    2, 1, figsize=(max(14, 0.35 * len(labels)), 9), sharex=True,
    constrained_layout=True, gridspec_kw={"height_ratios": [3, 1]},
)
axes[0].plot(x, fitted, "o", label="Minimizer")
axes[0].plot(x, generated, "x", markersize=8, label="Truth")
axes[0].axhline(0.0, color="0.5", linewidth=1)
axes[0].set_ylabel("Moment value")
axes[0].legend(frameon=False)

axes[1].plot(x, fitted - generated, "o", color="black")
axes[1].axhline(0.0, color="black", linestyle="--", linewidth=1)
axes[1].set(xlabel="Moment", ylabel=r"$\Delta$")
axes[1].set_xticks(x, labels, rotation=90)
for axis in axes:
    style_axis(axis)

save_figure(
    figure, OUTPUT_DIR / "fixed_moments.png",
    show=SHOW_PLOTS, dpi=300,
)
