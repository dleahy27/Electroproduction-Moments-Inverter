#!/usr/bin/env python3

"""Plot best-fit experimental amplitudes with Hessian uncertainties.

The covariance errors on magnitude and phase are propagated to approximate
Cartesian error bars for the Argand representation.
"""

import matplotlib.pyplot as plt
import numpy as np

from tutorials.analysis_utils import (
    PROJECT_DIR, MAGNITUDE_PATTERN, best_index, branches, phase_branch,
    read_tree, save_figure, style_axis, wave_label,
)
from tutorials.electroproduction_paper.datasets import CHANNELS, MAIN_CHANNEL_KEYS


# Settings: physics metadata lives in the family-level datasets module.
OUTPUT_FILES = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/root"
OUTPUT_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/plots/experimental_hessian"
TREE_NAME = "fitResults"
SHOW_PLOTS = False
CHANNELS_TO_PLOT = [CHANNELS[key] for key in MAIN_CHANNEL_KEYS]
GROUPS = [
    (("a", "T"), "Natural transverse"),
    (("b", "T"), "Unnatural transverse"),
    (("a", "L"), "Natural longitudinal"),
    (("b", "L"), "Unnatural longitudinal"),
]


def cartesian_errors(magnitude, phase, magnitude_error, phase_error):
    """Propagate independent magnitude and phase errors into Re and Im."""
    real_error = np.hypot(
        np.cos(phase) * magnitude_error,
        magnitude * np.sin(phase) * phase_error,
    )
    imaginary_error = np.hypot(
        np.sin(phase) * magnitude_error,
        magnitude * np.cos(phase) * phase_error,
    )
    return real_error, imaginary_error


for channel in CHANNELS_TO_PLOT:
    stem = channel.key
    channel_title = channel.title
    q2_values = np.asarray(channel.q2)
    invariant_mass = channel.invariant_mass
    minima = []
    channel_dir = OUTPUT_DIR / stem

    for bin_index, q2 in enumerate(q2_values):
        path = OUTPUT_FILES / f"{stem}_fit_{bin_index}.root"
        if not path.exists():
            continue
        names = set(branches(path, TREE_NAME))
        waves = sorted(
            name for name in names
            if MAGNITUDE_PATTERN.fullmatch(name)
            and phase_branch(name) in names
            and f"err__{name}" in names
            and f"err__{phase_branch(name)}" in names
        )
        status = ["chi2"] + [
            name for name in ("fit_ok", "status") if name in names
        ]
        columns = status + [
            branch for wave in waves
            for branch in (
                wave, phase_branch(wave),
                f"err__{wave}", f"err__{phase_branch(wave)}",
            )
        ]
        data = read_tree(path, TREE_NAME, columns)
        entry = best_index(data)
        minima.append((q2, waves, data, entry))

        # Plot the best complex amplitudes for this Q2 bin.
        figure, axes = plt.subplots(
            4, 1, figsize=(8, 24), sharex=True, constrained_layout=True
        )
        axes[0].text(
            0.04, 0.94,
            rf"$\langle Q^2\rangle={q2:g}\,\mathrm{{GeV}}^2$\n"
            rf"$\langle W\rangle={invariant_mass:.2f}\,\mathrm{{GeV}}$",
            transform=axes[0].transAxes, va="top",
        )
        for axis, (group, group_title) in zip(axes, GROUPS):
            selected = [wave for wave in waves if tuple(wave.split("_")[:2]) == group]
            for colour_index, wave in enumerate(selected):
                magnitude = data[wave][entry]
                phase = data[phase_branch(wave)][entry]
                real = magnitude * np.cos(phase)
                imaginary = magnitude * np.sin(phase)
                errors = cartesian_errors(
                    magnitude, phase, data[f"err__{wave}"][entry],
                    data[f"err__{phase_branch(wave)}"][entry],
                )
                axis.errorbar(
                    real, imaginary, xerr=errors[0], yerr=errors[1],
                    fmt="o", capsize=4, color=f"C{colour_index % 10}",
                    label=wave_label(wave),
                )
            axis.set(xlim=(-1, 1), ylim=(-1, 1), ylabel="Im",
                     title=f"{channel_title} — {group_title}")
            axis.set_aspect("equal", adjustable="box")
            axis.legend(frameon=False)
            style_axis(axis, grid=True)
        axes[-1].set_xlabel("Re")
        save_figure(
            figure, channel_dir / f"{stem}_fit_{bin_index}_argand.pdf",
            show=SHOW_PLOTS, transparent=True,
        )

    if not minima:
        continue

    # Plot every common amplitude magnitude against Q2.
    common_waves = sorted(set.intersection(*[set(result[1]) for result in minima]))
    figure, axes = plt.subplots(
        4, 4, figsize=(22, 18), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for axis, wave in zip(axes.flat, common_waves):
        q2 = np.array([result[0] for result in minima])
        values = np.array([result[2][wave][result[3]] for result in minima])
        errors = np.array([
            result[2][f"err__{wave}"][result[3]] for result in minima
        ])
        axis.errorbar(q2, values, yerr=errors, fmt="o", capsize=4)
        axis.set_ylim(-0.05, 1.0)
        if len(q2) > 1:
            axis.set_xscale("log")
        axis.set_title(wave_label(wave))
        style_axis(axis)
    for axis in axes.flat[len(common_waves):]:
        axis.set_visible(False)
    figure.suptitle(f"{channel_title} amplitude scaling (Hessian)")
    figure.supylabel("Magnitude")
    figure.supxlabel(r"$Q^2\ [\mathrm{GeV}^2]$")
    save_figure(
        figure, channel_dir / f"{stem}_scaling.png", show=SHOW_PLOTS
    )
