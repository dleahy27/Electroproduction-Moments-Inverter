"""Shared amplitude checks and plots for the fixed-polarization tutorials.

Both fixed studies store the same physical quantities; only the measured
nucleon-polarization information differs. Keeping their plotting machinery in
one place ensures branch discovery, best-fit selection, and visual conventions
cannot drift between the single- and double-polarization examples.
"""

import re

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from tutorials.analysis_utils import (
    best_index,
    branches,
    complex_amplitude,
    phase_branch,
    read_tree,
    save_figure,
    style_axis,
    wave_label,
)

# ROOT identifiers use ``m1`` for -1. Capturing the base wave and spin sector
# separately lets colour describe the wave while marker direction describes k.
AMPLITUDE_PATTERN = re.compile(r"^(?P<wave>[ab]_[TL]_\d+_m?\d+)_(?P<sector>1|m1)$")
SECTORS = ("1", "m1")
SECTOR_MARKERS = {"1": "^", "m1": "v"}
# Generated magnitudes below this numerical tolerance represent amplitudes
# that are absent from the physical model, including parity-forbidden waves.
ZERO_MAGNITUDE_TOLERANCE = 1.0e-12
# Exact synthetic moments should permit chi-square values very close to zero.
# Minima inside this interval are equivalent for selecting the branch nearest
# the known truth in this closure-only diagnostic.
CHI2_EQUIVALENCE_TOLERANCE = 1.0e-6
GROUPS = [
    (("a", "T"), "Natural transverse"),
    (("b", "T"), "Unnatural transverse"),
    (("a", "L"), "Natural longitudinal"),
    (("b", "L"), "Unnatural longitudinal"),
]


def _amplitude_schema(truth_file, truth_tree, fit_file, fit_tree):
    """Return matched magnitude branches after validating the ROOT schemas."""
    truth_names = set(branches(truth_file, truth_tree))
    fit_names = set(branches(fit_file, fit_tree))
    truth_magnitudes = {
        name for name in truth_names if AMPLITUDE_PATTERN.fullmatch(name)
    }
    fit_magnitudes = {name for name in fit_names if AMPLITUDE_PATTERN.fullmatch(name)}

    # An intersection would silently discard missing results and could produce
    # a plausible but incomplete figure, so require exact schema agreement.
    if truth_magnitudes != fit_magnitudes:
        raise RuntimeError(
            "Truth and fit amplitude branches do not match: "
            f"missing from fit={sorted(truth_magnitudes - fit_magnitudes)}, "
            f"missing from truth={sorted(fit_magnitudes - truth_magnitudes)}"
        )
    if not truth_magnitudes:
        raise RuntimeError("No explicit k=+1/k=-1 amplitude branches were found")

    magnitudes = sorted(
        truth_magnitudes,
        key=lambda name: (
            AMPLITUDE_PATTERN.fullmatch(name).group("wave"),
            AMPLITUDE_PATTERN.fullmatch(name).group("sector") == "m1",
        ),
    )
    missing_phases = [
        phase_branch(name)
        for name in magnitudes
        if phase_branch(name) not in truth_names or phase_branch(name) not in fit_names
    ]
    if missing_phases:
        raise RuntimeError(f"Missing paired phase branches: {missing_phases}")

    base_waves = sorted(
        {AMPLITUDE_PATTERN.fullmatch(name).group("wave") for name in magnitudes}
    )
    for wave in base_waves:
        present = {
            AMPLITUDE_PATTERN.fullmatch(name).group("sector")
            for name in magnitudes
            if AMPLITUDE_PATTERN.fullmatch(name).group("wave") == wave
        }
        if present != set(SECTORS):
            raise RuntimeError(
                f"{wave} does not contain both k sectors: {sorted(present)}"
            )
    return magnitudes, base_waves, truth_names, fit_names


def _add_legend(axis, waves, colours):
    """Keep wave labels outside the complex-amplitude data area."""
    handles = [
        Line2D([], [], color=colours[wave], linewidth=3, label=wave_label(wave))
        for wave in waves
    ]
    axis.legend(handles=handles, frameon=False, loc="upper left",
                bbox_to_anchor=(1.03, 1.0), borderaxespad=0)


def _finish_figure(figure, axes, title, reference, rotation):
    """Explain marker and gauge conventions without covering any amplitudes."""
    handles = [
        Line2D([], [], marker=marker, color="black", linestyle="none",
               markersize=9, label=label)
        for marker, label in [
            ("^", r"Generated $k=+1$"),
            ("v", r"Generated $k=-1$"),
            ("o", "Fit (truth-closest minimum)"),
        ]
    ]
    figure.legend(handles=handles, loc="lower center", ncol=3,
                  bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=14)
    figure.suptitle(title, fontsize=22)
    box = dict(boxstyle="round,pad=0.5", facecolor="white",
               edgecolor="0.7", linewidth=0.8)
    # Reserve a footer for the common marker key; panel legends and gauge boxes
    # are included in constrained layout and in the tight PDF bounding box.
    figure.get_layout_engine().set(rect=(0, 0.035, 1, 0.965))


def _reference_wave(base_waves, reflectivity):
    """Mirror EMI's highest-l, highest-m transverse non-flip reference."""
    candidates = [wave for wave in base_waves
                  if wave.startswith(f"{reflectivity}_T_")]
    return max(candidates, key=lambda wave: (
        int(wave.split("_")[2]), int(wave.split("_")[3].replace("m", "-"))
    ))


def analyse_fixed_closure(
    truth_file,
    fit_file,
    output_dir,
    *,
    truth_tree="genMoments",
    fit_tree="fitResults",
    show_plot=False,
    polarization="both",
    k_gauge=None,
):
    """Compare generated and fitted amplitudes in Cartesian and polar form."""
    if polarization not in ("both", "initial", "recoil"):
        raise ValueError(f"Unknown nucleon polarization: {polarization}")
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix",
        "font.size": 18, "axes.titlesize": 18, "axes.labelsize": 24,
        "xtick.labelsize": 17, "ytick.labelsize": 17,
        "legend.fontsize": 16, "axes.linewidth": 1.8,
        "lines.linewidth": 2.2,
    })
    magnitudes, base_waves, truth_names, fit_names = _amplitude_schema(
        truth_file, truth_tree, fit_file, fit_tree
    )
    reference = _reference_wave(base_waves, "a")
    rotation = None
    if polarization != "both":
        rotation = (k_gauge.replace(":", "_").replace("_-", "_m")
                    if k_gauge else _reference_wave(base_waves, "b"))
        if rotation not in base_waves:
            raise ValueError(f"Spin gauge wave absent from input: {rotation}")
    study_title = {
        "both": "Initial and recoil polarization",
        "initial": "Initial nucleon polarization",
        "recoil": "Recoil nucleon polarization",
    }[polarization]
    phases = [phase_branch(name) for name in magnitudes]
    # Generated observable branches carry a leading R, while reconstructed fit
    # moments omit it. These pairs test physical closure independently of which
    # ambiguous amplitude solution the minimizer finds.
    moment_pairs = sorted(
        (truth_name, truth_name[1:])
        for truth_name in truth_names
        if truth_name.startswith("RH") and truth_name[1:] in fit_names
    )

    # The truth tree has one point; every fit row is an independent random
    # start. Include available quality flags when choosing the best minimum.
    truth = read_tree(
        truth_file,
        truth_tree,
        magnitudes + phases + [truth_name for truth_name, _ in moment_pairs],
    )
    status = [name for name in ("chi2", "fit_ok", "status") if name in fit_names]
    fit = read_tree(
        fit_file,
        fit_tree,
        status + magnitudes + phases + [fit_name for _, fit_name in moment_pairs],
    )
    best_entry = best_index(fit)
    truth_amplitudes = {
        name: complex_amplitude(truth[name][0], truth[phase_branch(name)][0])
        for name in magnitudes
    }
    # Truth may select among observationally equivalent minima in a synthetic
    # closure test. Experimental-data analyses cannot and do not use this step.
    truth_vector = np.asarray([truth_amplitudes[name] for name in magnitudes])
    fit_matrix = np.column_stack(
        [complex_amplitude(fit[name], fit[phase_branch(name)]) for name in magnitudes]
    )
    amplitude_rms = np.sqrt(
        np.mean(np.abs(fit_matrix - truth_vector[None, :]) ** 2, axis=1)
    )
    accepted = np.isfinite(fit["chi2"])
    if "fit_ok" in fit:
        accepted &= fit["fit_ok"].astype(bool)
    elif "status" in fit:
        accepted &= fit["status"] == 0
    equivalent = accepted & (
        fit["chi2"] <= fit["chi2"][best_entry] + CHI2_EQUIVALENCE_TOLERANCE
    )
    equivalent_entries = np.flatnonzero(equivalent)
    entry = equivalent_entries[np.argmin(amplitude_rms[equivalent_entries])]
    fit_amplitudes = {
        name: complex_amplitude(fit[name][entry], fit[phase_branch(name)][entry])
        for name in magnitudes
    }

    print(f"\nClosure: {truth_file.stem}")
    print(
        f"Lowest accepted chi2: entry {best_entry}, "
        f"chi2 = {fit['chi2'][best_entry]:.7g}, "
        f"amplitude RMS = {amplitude_rms[best_entry]:.7g}"
    )
    print(
        f"Truth-closest equivalent minimum: entry {entry}, "
        f"chi2 = {fit['chi2'][entry]:.7g}, "
        f"amplitude RMS = {amplitude_rms[entry]:.7g} "
        f"({len(equivalent_entries)} minima within "
        f"Delta chi2 <= {CHI2_EQUIVALENCE_TOLERANCE:g})"
    )
    if moment_pairs:
        best_moment_residual = max(
            abs(fit[fit_name][best_entry] - truth[truth_name][0])
            for truth_name, fit_name in moment_pairs
        )
        selected_moment_residual = max(
            abs(fit[fit_name][entry] - truth[truth_name][0])
            for truth_name, fit_name in moment_pairs
        )
        print(
            "Maximum moment residual: "
            f"lowest-chi2={best_moment_residual:.7g}, "
            f"truth-closest={selected_moment_residual:.7g}"
        )
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

    # The truth model decides whether an amplitude physically exists. A fit
    # value is deliberately hidden with its zero truth partner: displaying
    # numerical noise for a parity-forbidden branch would suggest a wave that
    # is not part of the model. Keep such branches in the table above so the
    # stored zero remains available as a diagnostic.
    plotted_magnitudes = {
        name
        for name in magnitudes
        if abs(truth_amplitudes[name]) > ZERO_MAGNITUDE_TOLERANCE
    }
    plotted_waves = [
        wave
        for wave in base_waves
        if any(f"{wave}_{sector}" in plotted_magnitudes for sector in SECTORS)
    ]
    if not plotted_waves:
        raise RuntimeError(
            "All generated amplitudes are zero; there is nothing to plot"
        )

    groups = [
        (key, title)
        for key, title in GROUPS
        if any(tuple(wave.split("_")[:2]) == key for wave in plotted_waves)
    ]
    # Reuse the same wave colours in transverse/longitudinal and natural/
    # unnatural panels, as in the other fixed-amplitude Argands.
    projections = sorted({tuple(wave.split("_")[2:]) for wave in plotted_waves},
                         key=lambda pair: (int(pair[0]), int(pair[1].replace("m", "-"))))
    colours = {wave: f"C{projections.index(tuple(wave.split('_')[2:])) % 10}"
               for wave in plotted_waves}

    # Cartesian Argand representation. Draw the fit circle larger and first so
    # its rim remains visible behind the triangle when closure is exact.
    figure, axes = plt.subplots(
        len(groups),
        1,
        figsize=(11, 6 * len(groups)),
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes[:, 0]
    extent = 0.1
    for axis, (group, title) in zip(axes, groups):
        selected = [
            wave for wave in plotted_waves if tuple(wave.split("_")[:2]) == group
        ]
        for wave in selected:
            for sector in SECTORS:
                name = f"{wave}_{sector}"
                if name not in plotted_magnitudes:
                    continue
                generated = truth_amplitudes[name]
                fitted = fit_amplitudes[name]
                axis.plot(
                    fitted.real,
                    fitted.imag,
                    "o",
                    color=colours[wave],
                    markeredgewidth=1.8,
                    markersize=13,
                    zorder=3,
                )
                axis.plot(
                    generated.real,
                    generated.imag,
                    SECTOR_MARKERS[sector],
                    color=colours[wave],
                    markeredgecolor="black",
                    markersize=9,
                    zorder=4,
                )
                extent = max(
                    extent,
                    abs(generated.real),
                    abs(generated.imag),
                    abs(fitted.real),
                    abs(fitted.imag),
                )
        axis.set_title(title.title(), pad=12)
        axis.set_xlabel("Re")
        axis.set_ylabel("Im")
        axis.axhline(0, color="0.5", linewidth=0.8, zorder=0)
        axis.axvline(0, color="0.5", linewidth=0.8, zorder=0)
        axis.set_aspect("equal", adjustable="box")
        style_axis(axis, grid=True)
        _add_legend(axis, selected, colours)

    limit = 1.15 * extent
    for axis in axes:
        axis.set_xlim(-limit, limit)
        axis.set_ylim(-limit, limit)
    axes[-1].set_xlabel("Re")
    _finish_figure(figure, axes, study_title, reference, rotation)
    save_figure(
        figure,
        output_dir / f"{fit_file.stem}_argand.pdf",
        show=show_plot,
        dpi=600,
        transparent=True,
    )

    # Polar coordinates separate radial and angular discrepancies. Phase is
    # undefined at exactly zero magnitude, so those points stay in Argand only.
    figure, axes = plt.subplots(
        len(groups),
        1,
        figsize=(11, 6 * len(groups)),
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes[:, 0]
    for axis, (group, title) in zip(axes, groups):
        selected = [
            wave for wave in plotted_waves if tuple(wave.split("_")[:2]) == group
        ]
        for wave in selected:
            for sector in SECTORS:
                name = f"{wave}_{sector}"
                if name not in plotted_magnitudes:
                    continue
                generated = truth_amplitudes[name]
                fitted = fit_amplitudes[name]
                if abs(fitted) > 1.0e-12:
                    axis.plot(
                        np.angle(fitted),
                        abs(fitted),
                        "o",
                        color=colours[wave],
                        markeredgecolor="black",
                        markersize=12,
                        zorder=3,
                    )
                if abs(generated) > 1.0e-12:
                    axis.plot(
                        np.angle(generated),
                        abs(generated),
                        SECTOR_MARKERS[sector],
                        color=colours[wave],
                        markeredgecolor="black",
                        markersize=9,
                        zorder=4,
                    )
        axis.set_title(title.title(), pad=12)
        axis.set_xlim(-np.pi, np.pi)
        axis.set_ylabel("Magnitude")
        axis.set_xticks(
            [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
            [r"$-\pi$", r"$-\pi/2$", "$0$", r"$\pi/2$", r"$\pi$"],
        )
        style_axis(axis, grid=True)
        _add_legend(axis, selected, colours)

    axes[-1].set_xlabel("Phase (rad)")
    _finish_figure(figure, axes, study_title, reference, rotation)
    save_figure(
        figure,
        output_dir / f"{fit_file.stem}_magnitude_phase.pdf",
        show=show_plot,
        dpi=600,
        transparent=True,
    )
