#!/usr/bin/env python3

"""Plot bootstrap amplitudes for the experimental channels."""

import matplotlib.pyplot as plt
import numpy as np
import ROOT

from analysis_utils import (
    PROJECT_DIR, circular_mean, circular_std, phase_branch,
    save_figure, style_axis, wave_label,
)


# Settings: add or remove channels here without changing the analysis below.
OUTPUT_FILES = PROJECT_DIR / "OutputFiles"
OUTPUT_DIR = PROJECT_DIR / "AnalysisScripts/outputs/experimental_analysis_bootstrap"
TREE_NAME = "PartialWaves"
SHOW_PLOTS = False
CHANNELS = [
    ("e_rho", r"$e^-\rho^0$", np.array([0.82, 1.19, 1.66, 3.06]), 4.8),
    ("mu_rho", r"$\mu^-\rho^0$", np.array([1.14, 1.60, 2.80, 6.02]), 9.9),
    ("e_omega", r"$e^-\omega$", np.array([1.28, 2.00, 4.00]), 4.8),
    ("mu_omega", r"$\mu^-\omega$", np.array([1.16, 1.64, 3.58]), 7.6),
    ("e_phi", r"$e^-\phi$", np.array([1.2, 1.7, 4.5]), np.sqrt(21.89)),
]
GROUPS = [
    (("a", "T"), "Natural transverse"),
    (("b", "T"), "Unnatural transverse"),
    (("a", "L"), "Natural longitudinal"),
    (("b", "L"), "Unnatural longitudinal"),
]


def read_amplitudes(path):
    """Read valid bootstrap amplitudes and add parity-related L,-1 waves."""
    frame = ROOT.RDataFrame(TREE_NAME, str(path))
    names = {str(name) for name in frame.GetColumnNames()}
    if "valid" in names:
        frame = frame.Filter("valid != 0")

    definitions = [
        ("a_L_1_m1", "-a_L_1_1"),
        ("b_L_1_m1", "b_L_1_1"),
        ("aphi_L_1_m1", "aphi_L_1_1"),
        ("bphi_L_1_m1", "bphi_L_1_1"),
    ]
    for name, expression in definitions:
        if name not in names and expression.lstrip("-") in names:
            frame = frame.Define(name, expression)
            names.add(name)

    magnitudes = sorted(
        name for name in names
        if name.startswith(("a_", "b_"))
        and len(name.split("_")) == 4
        and phase_branch(name) in names
    )
    columns = [branch for wave in magnitudes for branch in (wave, phase_branch(wave))]
    arrays = frame.AsNumpy(columns=columns)
    return magnitudes, {name: np.asarray(values) for name, values in arrays.items()}


# Store mean magnitudes for the Q2-scaling panels after drawing each bin.
for stem, channel_title, q2_values, invariant_mass in CHANNELS:
    bin_results = []
    channel_dir = OUTPUT_DIR / stem

    for bin_index, q2 in enumerate(q2_values):
        path = OUTPUT_FILES / f"{stem}_bootstrap_{bin_index}.root"
        if not path.exists():
            continue
        waves, data = read_amplitudes(path)
        bin_results.append((q2, waves, data))

        # Argand means and Cartesian standard deviations summarize each bootstrap.
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
                amplitude = data[wave] * np.exp(1j * data[phase_branch(wave)])
                axis.errorbar(
                    np.mean(amplitude.real), np.mean(amplitude.imag),
                    xerr=np.std(amplitude.real), yerr=np.std(amplitude.imag),
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
            figure, channel_dir / f"{stem}_bootstrap_{bin_index}_argand.pdf",
            show=SHOW_PLOTS, transparent=True,
        )

        # Print the same magnitude and phase summaries numerically.
        with (channel_dir / f"{stem}_bootstrap_{bin_index}.txt").open("w") as stream:
            print(f"Q2 = {q2:g} GeV^2", file=stream)
            for wave in waves:
                print(
                    f"{wave:14s} magnitude = {np.mean(data[wave]):.7g} "
                    f"+/- {np.std(data[wave]):.7g}, phase = "
                    f"{circular_mean(data[phase_branch(wave)]):.7g} +/- "
                    f"{circular_std(data[phase_branch(wave)]):.7g}",
                    file=stream,
                )

    if not bin_results:
        continue

    # One panel per amplitude shows its magnitude as a function of Q2.
    common_waves = sorted(set.intersection(*[set(result[1]) for result in bin_results]))
    figure, axes = plt.subplots(
        4, 4, figsize=(22, 18), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for axis, wave in zip(axes.flat, common_waves):
        q2 = np.array([result[0] for result in bin_results])
        means = np.array([np.mean(result[2][wave]) for result in bin_results])
        errors = np.array([np.std(result[2][wave]) for result in bin_results])
        axis.errorbar(q2, means, yerr=errors, fmt="o", capsize=4)
        axis.set_ylim(-0.05, 1.0)
        if len(q2) > 1:
            axis.set_xscale("log")
        axis.set_title(wave_label(wave))
        style_axis(axis)
    for axis in axes.flat[len(common_waves):]:
        axis.set_visible(False)
    figure.suptitle(f"Amplitude scaling for {channel_title}")
    figure.supylabel(r"Amplitude magnitude $|A|$")
    figure.supxlabel(r"$Q^2\ [\mathrm{GeV}^2]$")
    save_figure(
        figure, channel_dir / f"{stem}_scaling.png", show=SHOW_PLOTS
    )
