#!/usr/bin/env python3

"""Plot GlueX photoproduction amplitudes as functions of momentum transfer."""

import matplotlib.pyplot as plt
import numpy as np
import ROOT

from analysis_utils import (
    PROJECT_DIR, circular_mean, circular_std, phase_branch,
    save_figure, style_axis, wave_label,
)


# Settings: the t values correspond to ROOT-file bin indices 0 through 17.
OUTPUT_FILES = PROJECT_DIR / "OutputFiles"
OUTPUT_DIR = PROJECT_DIR / "AnalysisScripts/outputs/gluex_analysis"
SHOW_PLOTS = False
T_VALUES = np.array([
    0.107, 0.121, 0.138, 0.157, 0.178, 0.203,
    0.230, 0.262, 0.297, 0.338, 0.384, 0.436,
    0.496, 0.564, 0.640, 0.728, 0.827, 0.940,
])
WAVES = [
    "a_T_1_1", "a_T_1_0", "a_T_1_m1", "a_T_0_0",
    "b_T_1_1", "b_T_1_0", "b_T_1_m1", "b_T_0_0",
]


def read_columns(path, tree_name, select_minimum=False):
    """Read the required amplitudes, optionally near the best objective."""
    frame = ROOT.RDataFrame(tree_name, str(path))
    names = {str(name) for name in frame.GetColumnNames()}
    objective = "chi2" if "chi2" in names else "log_val"
    if select_minimum:
        minimum = frame.Min(objective).GetValue()
        frame = frame.Filter(f"{objective} <= {minimum + 0.01:.17g}")
    columns = [branch for wave in WAVES for branch in (wave, phase_branch(wave))]
    return {name: np.asarray(values) for name, values in frame.AsNumpy(columns).items()}


bootstrap_means = []
bootstrap_errors = []
available_t = []

for bin_index, momentum_transfer in enumerate(T_VALUES):
    fit_path = OUTPUT_FILES / f"gluex_fit_{bin_index}.root"
    bootstrap_path = OUTPUT_FILES / f"gluex_bootstrap_{bin_index}.root"

    if fit_path.exists():
        data = read_columns(fit_path, "fitResults", select_minimum=True)
        figure, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
        for row, reflectivity in enumerate(("a", "b")):
            for wave in [name for name in WAVES if name.startswith(reflectivity)]:
                magnitude = data[wave]
                phase = data[phase_branch(wave)]
                amplitude = magnitude * np.exp(1j * phase)
                axes[row, 0].plot(phase, magnitude, "o", label=wave_label(wave))
                axes[row, 1].plot(
                    amplitude.real, amplitude.imag, "o", label=wave_label(wave)
                )
            axes[row, 0].set(xlim=(-np.pi, np.pi), ylim=(-0.05, 1),
                             ylabel="Magnitude")
            axes[row, 1].set(xlim=(-1, 1), ylim=(-1, 1), ylabel="Im")
            axes[row, 1].set_aspect("equal", adjustable="box")
            axes[row, 1].legend(frameon=False)
            for axis in axes[row]:
                style_axis(axis, grid=True)
        axes[-1, 0].set_xlabel("Phase (rad)")
        axes[-1, 1].set_xlabel("Re")
        figure.suptitle(rf"GlueX amplitudes at $-\bar{{t}}={momentum_transfer:g}$")
        save_figure(
            figure, OUTPUT_DIR / f"gluex_fit_{bin_index}.png", show=SHOW_PLOTS
        )

    if not bootstrap_path.exists():
        continue
    data = read_columns(bootstrap_path, "PartialWaves")
    means = np.array([np.mean(data[wave]) for wave in WAVES])
    errors = np.array([np.std(data[wave]) for wave in WAVES])
    bootstrap_means.append(means)
    bootstrap_errors.append(errors)
    available_t.append(momentum_transfer)

    # Bootstrap means and widths in magnitude-phase and Argand coordinates.
    figure, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    for row, reflectivity in enumerate(("a", "b")):
        for colour_index, wave in enumerate(
            name for name in WAVES if name.startswith(reflectivity)
        ):
            magnitude = data[wave]
            phase = data[phase_branch(wave)]
            amplitude = magnitude * np.exp(1j * phase)
            colour = f"C{colour_index % 10}"
            axes[row, 0].errorbar(
                circular_mean(phase), np.mean(magnitude),
                xerr=circular_std(phase), yerr=np.std(magnitude),
                fmt="o", capsize=3, color=colour, label=wave_label(wave),
            )
            axes[row, 1].errorbar(
                np.mean(amplitude.real), np.mean(amplitude.imag),
                xerr=np.std(amplitude.real), yerr=np.std(amplitude.imag),
                fmt="o", capsize=3, color=colour, label=wave_label(wave),
            )
        axes[row, 0].set(xlim=(-np.pi, np.pi), ylim=(-0.05, 1),
                         ylabel="Magnitude")
        axes[row, 1].set(xlim=(-1, 1), ylim=(-1, 1), ylabel="Im")
        axes[row, 1].set_aspect("equal", adjustable="box")
        axes[row, 1].legend(frameon=False)
        for axis in axes[row]:
            style_axis(axis, grid=True)
    axes[-1, 0].set_xlabel("Phase (rad)")
    axes[-1, 1].set_xlabel("Re")
    figure.suptitle(rf"GlueX bootstrap at $-\bar{{t}}={momentum_transfer:g}$")
    save_figure(
        figure, OUTPUT_DIR / f"gluex_bootstrap_{bin_index}.png", show=SHOW_PLOTS
    )

# The original scaling figure contains the transverse amplitudes used above.
if bootstrap_means:
    means = np.asarray(bootstrap_means)
    errors = np.asarray(bootstrap_errors)
    figure, axes = plt.subplots(
        2, 4, figsize=(16, 8), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for index, (axis, wave) in enumerate(zip(axes.flat, WAVES)):
        axis.errorbar(available_t, means[:, index], yerr=errors[:, index],
                      fmt="o", capsize=3)
        axis.set_xscale("log")
        axis.set_ylim(-0.05, 1)
        axis.set_title(wave_label(wave))
        style_axis(axis)
    figure.supylabel("Magnitude")
    figure.supxlabel(r"$-\bar{t}\ [\mathrm{GeV}^2]$")
    save_figure(
        figure, OUTPUT_DIR / "gluex_scaling.pdf",
        show=SHOW_PLOTS, dpi=800, transparent=True,
    )
