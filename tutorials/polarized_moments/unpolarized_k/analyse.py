#!/usr/bin/env python3

"""Study how an unpolarized fit combines two generated k sectors.

The observable reference is their incoherent quadrature magnitude. Repeated
generated points reveal bias and spread when that two-sector truth is fitted
with one unpolarized amplitude per wave.
"""

import re

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from tutorials.analysis_utils import (
    PROJECT_DIR, best_index, branches, complex_amplitude, phase_branch,
    read_tree, save_figure, style_axis, wave_label,
)


# Settings: BATCH=None selects the most recently modified common batch.
BATCH = None
REPEAT_TO_PLOT = 26
SUPPRESSION_INDEX_TO_PLOT = 0
SHOW_PLOTS = False

CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/unpolarized_k"
INPUT_BASE = CASE_DIR / "output/generated"
FIT_BASE = CASE_DIR / "output/root"
OUTPUT_DIR = CASE_DIR / "output/plots"


def amplitude_ratio(magnitudes):
    """Calculate R from the fitted or quadrature amplitude magnitudes."""
    transverse = sum(
        value**2 for wave, value in magnitudes.items()
        if wave.split("_")[1] == "T"
    )
    longitudinal = sum(
        (1 if wave.split("_")[3] == "0" else 2) * value**2
        for wave, value in magnitudes.items()
        if wave.split("_")[1] == "L"
    )
    return longitudinal / transverse


# Locate matching generated and fitted batch directories.
input_batches = {path.name: path for path in INPUT_BASE.iterdir() if path.is_dir()}
fit_batches = {path.name: path for path in FIT_BASE.iterdir() if path.is_dir()}
common_batches = sorted(input_batches.keys() & fit_batches.keys())
if BATCH is None:
    BATCH = max(
        common_batches,
        key=lambda name: max(
            (path.stat().st_mtime for path in input_batches[name].glob("*.root")),
            default=0,
        ),
    )
truth_dir = input_batches[BATCH]
fit_dir = fit_batches[BATCH]
output_dir = OUTPUT_DIR / BATCH

file_pattern = re.compile(r"^fixed_N(\d+)(?:_M(\d+))?_moments\.root$")
truth_amplitude_pattern = re.compile(r"^([ab]_[TL]_\d+_m?\d+)_(1|m1)$")

# Each list stores one value per repetition, grouped by suppression index and wave.
signed_amplitude_errors = {}
signed_R_errors = {}
mean_magnitudes = {}
fit_records = []

for truth_file in sorted(truth_dir.glob("fixed_N*_moments.root")):
    match = file_pattern.fullmatch(truth_file.name)
    if not match:
        continue
    repeat = int(match.group(1))
    suppression_index = int(match.group(2) or 0)
    fit_file = fit_dir / truth_file.name.replace(
        "_moments.root", "_unpolarized_fit.root"
    )
    if not fit_file.exists():
        continue

    truth_names = set(branches(truth_file, "genMoments"))
    fit_names = set(branches(fit_file, "fitResults"))
    waves = sorted({
        found.group(1)
        for name in truth_names
        if (found := truth_amplitude_pattern.fullmatch(name))
        and found.group(1) in fit_names
    })

    truth_columns = [
        branch
        for wave in waves
        for sector in ("1", "m1")
        for branch in (f"{wave}_{sector}", f"{phase_branch(wave)}_{sector}")
    ]
    fit_columns = ["chi2"]
    fit_columns += [name for name in ("fit_ok", "status") if name in fit_names]
    fit_columns += [branch for wave in waves for branch in (wave, phase_branch(wave))]

    truth = read_tree(truth_file, "genMoments", truth_columns)
    fit = read_tree(fit_file, "fitResults", fit_columns)
    entry = best_index(fit, require_fit_ok=False)

    generated_magnitudes = {
        wave: np.hypot(truth[f"{wave}_1"][0], truth[f"{wave}_m1"][0])
        for wave in waves
    }
    fitted_magnitudes = {wave: abs(fit[wave][entry]) for wave in waves}

    for wave in waves:
        key = (suppression_index, wave)
        mean_magnitudes.setdefault(key, []).append(
            (generated_magnitudes[wave], fitted_magnitudes[wave])
        )
        signed_amplitude_errors.setdefault(key, []).append(
            fitted_magnitudes[wave] - generated_magnitudes[wave]
        )

    generated_R = amplitude_ratio(generated_magnitudes)
    fitted_R = amplitude_ratio(fitted_magnitudes)
    signed_R_errors.setdefault(suppression_index, []).append(
        fitted_R - generated_R
    )
    fit_records.append({
        "repeat": repeat,
        "suppression": suppression_index,
        "truth_file": truth_file,
        "fit_file": fit_file,
        "entry": entry,
        "chi2": fit["chi2"][entry],
        "waves": waves,
        "truth": truth,
        "fit": fit,
    })


def mean_and_error(values):
    """Return the sample mean and standard error over repeated fits."""
    values = np.asarray(values)
    error = np.std(values, ddof=1) / np.sqrt(len(values)) if len(values) > 1 else 0.0
    return np.mean(values), error


# Print the quadrature and fitted magnitudes averaged over repetitions.
print(f"Batch: {BATCH}")
for key in sorted(mean_magnitudes):
    values = np.asarray(mean_magnitudes[key])
    print(
        f"M={key[0]:2d} {key[1]:14s} "
        f"generated={values[:, 0].mean():.6g} fit={values[:, 1].mean():.6g}"
    )

# Plot signed amplitude and R errors across suppression settings.
figure, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
for wave in sorted({wave for _, wave in signed_amplitude_errors}):
    settings = sorted(
        setting for setting, name in signed_amplitude_errors if name == wave
    )
    values = [
        mean_and_error(signed_amplitude_errors[(setting, wave)])
        for setting in settings
    ]
    axes[0].errorbar(
        settings, [value[0] for value in values],
        yerr=[value[1] for value in values], marker="o", capsize=3,
        label=wave_label(wave),
    )

if signed_R_errors:
    settings = sorted(signed_R_errors)
    values = [mean_and_error(signed_R_errors[setting]) for setting in settings]
    axes[1].errorbar(
        settings, [value[0] for value in values],
        yerr=[value[1] for value in values], marker="o", capsize=3,
        color="black",
    )
else:
    axes[1].text(0.5, 0.5, "No longitudinal amplitudes",
                 ha="center", va="center", transform=axes[1].transAxes)
for axis, title in zip(
    axes,
    [
        r"Amplitude signed error: $|A_{\rm fit}|-|A_{\rm quad}|$",
        r"$R$ signed error: $R_{\rm fit}-R_{\rm generated}$",
    ],
):
    axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1)
    axis.set(
        xlabel="Suppression/dominance index M",
        ylabel="Fit - generated",
        title=title,
    )
    style_axis(axis, grid=True)
axes[0].legend(frameon=False, ncol=2)
save_figure(
    figure, output_dir / "signed_errors.pdf", show=SHOW_PLOTS
)

# Select one result for direct complex-amplitude plots.
selected = next(
    record for record in fit_records
    if record["repeat"] == REPEAT_TO_PLOT
    and record["suppression"] == SUPPRESSION_INDEX_TO_PLOT
)
wave_rows = []
for wave in selected["waves"]:
    fit_value = complex_amplitude(
        selected["fit"][wave][selected["entry"]],
        selected["fit"][phase_branch(wave)][selected["entry"]],
    )
    truth_values = [
        complex_amplitude(
            selected["truth"][f"{wave}_{sector}"][0],
            selected["truth"][f"{phase_branch(wave)}_{sector}"][0],
        )
        for sector in ("1", "m1")
    ]
    wave_rows.append((wave, fit_value, truth_values))

groups = [
    (("a", "T"), "Natural transverse"),
    (("b", "T"), "Unnatural transverse"),
    (("a", "L"), "Natural longitudinal"),
    (("b", "L"), "Unnatural longitudinal"),
]

figure, axes = plt.subplots(
    4, 1, figsize=(8, 24), sharex=True, constrained_layout=True
)
for axis, (group, title) in zip(axes, groups):
    for colour_index, (wave, fitted, generated) in enumerate(wave_rows):
        if tuple(wave.split("_")[:2]) != group:
            continue
        colour = f"C{colour_index % 10}"
        axis.plot(fitted.real, fitted.imag, "o", color=colour,
                  markersize=12, label=wave_label(wave))
        axis.plot(generated[0].real, generated[0].imag, "^", color=colour,
                  markersize=11)
        axis.plot(generated[1].real, generated[1].imag, "v", color=colour,
                  markersize=11)
    axis.set(xlim=(-1, 1), ylim=(-1, 1), ylabel="Im", title=title)
    axis.set_aspect("equal", adjustable="box")
    style_axis(axis, grid=True)
    handles, _ = axis.get_legend_handles_labels()
    handles.extend([
        Line2D([], [], marker="^", color="black", linestyle="none",
               label=r"Generated $k=+1$"),
        Line2D([], [], marker="v", color="black", linestyle="none",
               label=r"Generated $k=-1$"),
    ])
    axis.legend(handles=handles, loc="upper left", frameon=False)
axes[-1].set_xlabel("Re")
save_figure(figure, output_dir / "selected_argand.pdf", show=SHOW_PLOTS)

# Show the same selected result as magnitude against phase.
figure, axes = plt.subplots(
    4, 1, figsize=(8, 24), sharex=True, sharey=True,
    constrained_layout=True,
)
for axis, (group, title) in zip(axes, groups):
    for colour_index, (wave, fitted, generated) in enumerate(wave_rows):
        if tuple(wave.split("_")[:2]) != group:
            continue
        colour = f"C{colour_index % 10}"
        axis.plot(np.angle(fitted), abs(fitted), "o", color=colour,
                  markersize=12, label=wave_label(wave))
        axis.plot(np.angle(generated[0]), abs(generated[0]), "^",
                  color=colour, markersize=11)
        if abs(generated[1]) > 1.0e-12:
            axis.plot(np.angle(generated[1]), abs(generated[1]), "v",
                      color=colour, markersize=11)
    axis.set(xlim=(-np.pi, np.pi), ylabel="Magnitude", title=title)
    axis.set_xticks(
        [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
        [r"$-\pi$", r"$-\pi/2$", "$0$", r"$\pi/2$", r"$\pi$"],
    )
    style_axis(axis, grid=True)
    if axis.get_legend_handles_labels()[0]:
        axis.legend(loc="upper left", frameon=False)
axes[-1].set_xlabel("Phase (rad)")
save_figure(figure, output_dir / "selected_magnitude_phase.pdf", show=SHOW_PLOTS)
