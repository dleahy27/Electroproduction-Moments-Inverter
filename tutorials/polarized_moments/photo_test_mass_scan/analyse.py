#!/usr/bin/env python3

"""Plot the generated and fitted PhotoTest mass scan.

The CSV manifest is the authoritative join between mass, k scale, truth file,
and fit file. Arrays are indexed by scale, mass, and wave so every plot uses
the same ordering and missing scan points remain explicit NaNs.
"""

import csv
import re

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import ROOT

from tutorials.analysis_utils import (
    PROJECT_DIR,
    best_index,
    complex_amplitude,
    phase_branch,
    save_figure,
    style_axis,
    wave_label,
)


# Settings: edit these values before running the analysis.
CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/photo_test_mass_scan"
MANIFEST_FILE = CASE_DIR / "output/root/scan_points.csv"
OUTPUT_DIR = CASE_DIR / "output/plots"
SHOW_PLOTS = False
FIGURE_DPI = 220
POINT_OFFSET_FRACTION = 0.06


ROOT.gROOT.SetBatch(not SHOW_PLOTS)


def wave_sort_key(name):
    """Sort the transverse waves by reflectivity, l, and m."""
    reflectivity, _, ell, projection = name.split("_")
    return reflectivity, int(ell), int(projection.replace("m", "-"))


def moment_sort_key(pair):
    """Sort moments by alpha, L, and M."""
    fields = pair[1].replace("H04_", "0_").replace("H_", "").split("_")
    return tuple(int(field.replace("m", "-")) for field in fields)


def offset_masses(masses, number_of_scales):
    """Separate scale points while retaining mass as the x variable."""
    if len(masses) < 2 or number_of_scales < 2:
        return np.tile(masses, (number_of_scales, 1))
    offsets = np.linspace(
        -POINT_OFFSET_FRACTION / 2,
        POINT_OFFSET_FRACTION / 2,
        number_of_scales,
    )
    return masses[None, :] + offsets[:, None]


def read_scan_point(point, waves, moment_pairs):
    """Read one generated point and its lowest-chi-square accepted fit."""
    truth_path = PROJECT_DIR / point["truth_file"]
    fit_path = PROJECT_DIR / point["fit_file"]
    if not truth_path.is_file() or not fit_path.is_file():
        return None

    truth_columns = [
        branch
        for wave in waves
        for sector in ("1", "m1")
        for branch in (f"{wave}_{sector}", f"{phase_branch(wave)}_{sector}")
    ]
    truth_columns += [truth_name for truth_name, _ in moment_pairs]

    fit_frame = ROOT.RDataFrame("fitResults", str(fit_path))
    fit_names = {str(name) for name in fit_frame.GetColumnNames()}
    status_columns = [
        name for name in ("chi2", "fit_ok", "status") if name in fit_names
    ]
    fit_columns = status_columns + [
        branch for wave in waves for branch in (wave, phase_branch(wave))
    ]
    fit_columns += [fit_name for _, fit_name in moment_pairs]

    truth = ROOT.RDataFrame("genMoments", str(truth_path)).AsNumpy(
        columns=truth_columns
    )
    fit = fit_frame.AsNumpy(columns=fit_columns)
    entry = best_index(fit)

    plus = np.empty(len(waves), dtype=complex)
    minus = np.empty(len(waves), dtype=complex)
    fitted = np.empty(len(waves), dtype=complex)

    # Each reflectivity has its own fixed phase in the unpolarized fit. Rotate
    # the generated amplitudes into those same two conventions.
    rotations = {}
    for reflectivity in ("a", "b"):
        reference = f"{reflectivity}phi_T_2_2_1"
        rotations[reflectivity] = np.exp(-1j * truth[reference][0])

    for index, wave in enumerate(waves):
        phase = phase_branch(wave)
        fitted[index] = complex_amplitude(
            fit[wave][entry], fit[phase][entry]
        )
        plus[index] = complex_amplitude(
            truth[f"{wave}_1"][0], truth[f"{phase}_1"][0]
        ) * rotations[wave[0]]
        minus[index] = complex_amplitude(
            truth[f"{wave}_m1"][0], truth[f"{phase}_m1"][0]
        ) * rotations[wave[0]]

    generated_moments = np.array([
        truth[truth_name][0] for truth_name, _ in moment_pairs
    ])
    fitted_moments = np.array([
        fit[fit_name][entry] for _, fit_name in moment_pairs
    ])
    return plus, minus, fitted, generated_moments, fitted_moments


def scale_handles(scales, colours):
    """Make a compact colour legend shared by the scan figures."""
    return [
        Line2D(
            [], [], marker="o", linestyle="none", color=colour,
            label=fr"$s={scale:.3g}$",
        )
        for scale, colour in zip(scales, colours)
    ]


def finish_mass_axis(axis, title, ylabel):
    """Apply the common mass-scan labels and styling."""
    axis.set(xlabel=r"Invariant mass $M$ (GeV)", ylabel=ylabel, title=title)
    #axis.set_xscale("log")
    style_axis(axis, grid=True)


# Read the scan definition written by the runner.
with MANIFEST_FILE.open(newline="") as stream:
    points = list(csv.DictReader(stream))

masses = np.array(sorted({float(point["mass_GeV"]) for point in points}))
scales = np.array(sorted({float(point["k_minus_scale"]) for point in points}))
mass_index = {mass: index for index, mass in enumerate(masses)}
scale_index = {scale: index for index, scale in enumerate(scales)}

# Discover the common unpolarized amplitudes and moments from one scan point.
first_point = next(
    point for point in points
    if (PROJECT_DIR / point["truth_file"]).is_file()
    and (PROJECT_DIR / point["fit_file"]).is_file()
)
first_truth = ROOT.RDataFrame(
    "genMoments", str(PROJECT_DIR / first_point["truth_file"])
)
first_fit = ROOT.RDataFrame(
    "fitResults", str(PROJECT_DIR / first_point["fit_file"])
)
truth_names = {str(name) for name in first_truth.GetColumnNames()}
fit_names = {str(name) for name in first_fit.GetColumnNames()}

wave_pattern = re.compile(r"^[ab]_T_\d+_(?:m?\d+)$")
waves = sorted(
    (name for name in fit_names if wave_pattern.fullmatch(name)),
    key=wave_sort_key,
)

# RH_alpha_0_0_L_M is the unpolarized projection H_alpha_L_M. The alpha=0
# moments use the conventional H04 name in the unpolarized fit output.
truth_moment_pattern = re.compile(r"^RH_(\d+)_0_0_(\d+)_(m?\d+)$")
moment_pairs = []
for truth_name in truth_names:
    match = truth_moment_pattern.fullmatch(truth_name)
    if not match:
        continue
    alpha, ell, projection = match.groups()
    fit_name = (
        f"H04_{ell}_{projection}"
        if alpha == "0"
        else f"H_{alpha}_{ell}_{projection}"
    )
    if fit_name in fit_names:
        moment_pairs.append((truth_name, fit_name))
moment_pairs.sort(key=moment_sort_key)

# Store all results as arrays indexed by scale, mass, and observable.
amplitude_shape = (len(scales), len(masses), len(waves))
moment_shape = (len(scales), len(masses), len(moment_pairs))
truth_plus = np.full(amplitude_shape, np.nan + 1j * np.nan, dtype=complex)
truth_minus = np.full(amplitude_shape, np.nan + 1j * np.nan, dtype=complex)
fitted_amplitudes = np.full(
    amplitude_shape, np.nan + 1j * np.nan, dtype=complex
)
truth_moments = np.full(moment_shape, np.nan)
fitted_moments = np.full(moment_shape, np.nan)

read_points = 0
for point in points:
    result = read_scan_point(point, waves, moment_pairs)
    if result is None:
        continue
    scale_row = scale_index[float(point["k_minus_scale"])]
    mass_column = mass_index[float(point["mass_GeV"])]
    truth_plus[scale_row, mass_column] = result[0]
    truth_minus[scale_row, mass_column] = result[1]
    fitted_amplitudes[scale_row, mass_column] = result[2]
    truth_moments[scale_row, mass_column] = result[3]
    fitted_moments[scale_row, mass_column] = result[4]
    read_points += 1


# Use colour for scale and a small horizontal displacement for coincident
# nominal masses. The displacement is multiplicative because the mass grid is
# logarithmic.
colours = plt.colormaps["viridis"](np.linspace(0.08, 0.92, len(scales)))
plot_masses = offset_masses(masses, len(scales))
#log_spacing = np.median(np.diff(np.log(masses))) if len(masses) > 1 else 0.0
#marker_gap = 0.0
colour_legend = scale_handles(scales, colours)
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 11,
})


# Plot each input/output moment comparison in its own file.
moment_dir = OUTPUT_DIR / "moments"
for index, (_, moment) in enumerate(moment_pairs):
    figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=True)
    for row, colour in enumerate(colours):
        axis.scatter(
            plot_masses[row], truth_moments[row, :, index],
            marker="*", s=46, color=colour, zorder=4, edgecolors="black",
        )
        axis.scatter(
            plot_masses[row], fitted_moments[row, :, index],
            marker="o", s=60, color=colour, zorder=3,
        )
    axis.axhline(0.0, color="0.45", linewidth=0.8)
    finish_mass_axis(axis, moment, "Moment value")
    marker_legend = [
        Line2D([], [], marker="*", linestyle="none", color="black",
               markersize=10, label="Input moment"),
        Line2D([], [], marker="o", linestyle="none", color="black",
               label="Fitted moment"),
    ]
    axis.legend(handles=colour_legend + marker_legend, frameon=False, ncol=2)
    save_figure(
        figure, moment_dir / f"{moment}.pdf",
        show=SHOW_PLOTS, dpi=FIGURE_DPI,
    )


# Plot each amplitude magnitude and phase in separate mass-scan figures.
amplitude_dir = OUTPUT_DIR / "amplitudes"
error_dir = OUTPUT_DIR / "signed_errors"
for index, wave in enumerate(waves):
    title = wave_label(wave)

    figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=True)
    for row, colour in enumerate(colours):
        axis.scatter(
            plot_masses[row],
            np.abs(truth_plus[row, :, index]),
            marker="^", s=46, color=colour, zorder=4,edgecolors="black",
        )
        axis.scatter(
            plot_masses[row], np.abs(truth_minus[row, :, index]),
            marker="v", s=46, color=colour, zorder=4,edgecolors="black",
        )
        axis.scatter(
            plot_masses[row],
            np.abs(fitted_amplitudes[row, :, index]),
            marker="o", s=60, color=colour, zorder=3,
        )
    finish_mass_axis(axis, title, "Amplitude magnitude")
    amplitude_markers = [
        Line2D([], [], marker="^", linestyle="none", color="black",
               label=r"Generated $k=+1$"),
        Line2D([], [], marker="v", linestyle="none", color="black",
               label=r"Generated $k=-1$"),
        Line2D([], [], marker="o", linestyle="none", color="black",
               label="Unpolarized fit"),
    ]
    axis.legend(
        handles=colour_legend + amplitude_markers, frameon=False, ncol=2
    )
    save_figure(
        figure, amplitude_dir / f"{wave}_magnitude.pdf",
        show=SHOW_PLOTS, dpi=FIGURE_DPI,
    )

    figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=True)
    for row, colour in enumerate(colours):
        plus_phase = np.angle(truth_plus[row, :, index])
        minus_phase = np.angle(truth_minus[row, :, index])
        minus_phase[np.abs(truth_minus[row, :, index]) < 1.0e-12] = np.nan
        axis.scatter(
            plot_masses[row], plus_phase,
            marker="^", s=46, color=colour, zorder=4, edgecolors="black",
        )
        axis.scatter(
            plot_masses[row], minus_phase,
            marker="v", s=46, color=colour, zorder=4, edgecolors="black",
        )
        axis.scatter(
            plot_masses[row],
            np.angle(fitted_amplitudes[row, :, index]),
            marker="o", s=60, color=colour, zorder=3,
        )
    finish_mass_axis(axis, title, "Phase (rad)")
    axis.set_ylim(-np.pi, np.pi)
    axis.set_yticks(
        [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
        [r"$-\pi$", r"$-\pi/2$", "$0$", r"$\pi/2$", r"$\pi$"],
    )
    axis.legend(
        handles=colour_legend + amplitude_markers, frameon=False, ncol=2
    )
    save_figure(
        figure, amplitude_dir / f"{phase_branch(wave)}.pdf",
        show=SHOW_PLOTS, dpi=FIGURE_DPI,
    )

# Plot the signed amplitude error as a function of the k=-1 scaling.
# Make one figure per mass bin, with natural (a) and unnatural (b)
# linear and squared-amplitude differences above their respective totals.
for mass_index, mass in enumerate(masses):

    quadrature = np.hypot(
        np.abs(truth_plus[:, mass_index, :]),
        np.abs(truth_minus[:, mass_index, :]),
    )
    signed_error = (
        np.abs(fitted_amplitudes[:, mass_index, :]) - quadrature
    )
    signed_squared_error = (
        np.abs(fitted_amplitudes[:, mass_index, :]) ** 2 - quadrature ** 2
    )

    figure, axes = plt.subplots(
        2, 4,
        figsize=(28, 10),
        sharex=True,
        constrained_layout=True,
    )

    (
        natural_axis,
        unnatural_axis,
        natural_squared_axis,
        unnatural_squared_axis,
    ) = axes[0]
    (
        natural_total_axis,
        unnatural_total_axis,
        natural_squared_total_axis,
        unnatural_squared_total_axis,
    ) = axes[1]
    unnatural_axis.sharey(natural_axis)
    unnatural_squared_axis.sharey(natural_squared_axis)

    natural_waves = [
        (index, wave)
        for index, wave in enumerate(waves)
        if wave.startswith("a_")
    ]
    unnatural_waves = [
        (index, wave)
        for index, wave in enumerate(waves)
        if wave.startswith("b_")
    ]

    natural_colours = plt.colormaps["tab10"](
        np.linspace(0.0, 0.9, max(len(natural_waves), 1))
    )
    unnatural_colours = plt.colormaps["tab10"](
        np.linspace(0.0, 0.9, max(len(unnatural_waves), 1))
    )

    for colour, (wave_index, wave) in zip(
        natural_colours, natural_waves
    ):
        natural_axis.plot(
            scales,
            signed_error[:, wave_index],
            marker="o",
            color=colour,
            label=wave_label(wave),
        )
        natural_squared_axis.plot(
            scales,
            signed_squared_error[:, wave_index],
            marker="o",
            color=colour,
            label=wave_label(wave),
        )

    for colour, (wave_index, wave) in zip(
        unnatural_colours, unnatural_waves
    ):
        unnatural_axis.plot(
            scales,
            signed_error[:, wave_index],
            marker="o",
            color=colour,
            label=wave_label(wave),
        )
        unnatural_squared_axis.plot(
            scales,
            signed_squared_error[:, wave_index],
            marker="o",
            color=colour,
            label=wave_label(wave),
        )

    for total_axis, exchange_waves, errors in (
        (natural_total_axis, natural_waves, signed_error),
        (unnatural_total_axis, unnatural_waves, signed_error),
        (natural_squared_total_axis, natural_waves, signed_squared_error),
        (unnatural_squared_total_axis, unnatural_waves, signed_squared_error),
    ):
        wave_indices = [index for index, _ in exchange_waves]
        total_axis.scatter(
            scales,
            np.sum(errors[:, wave_indices], axis=1),
            marker="o",
            color="tab:blue",
        )

    for axis in axes.flat:
        axis.axhline(
            0.0,
            color="black",
            linestyle="--",
            linewidth=0.9,
        )
        axis.set_xlabel(r"Scaling factor (relative to $k=+1$)")
        # A symmetric-log transform keeps the physically useful zero-scale
        # point while still separating values that span several decades.
        axis.set_xscale("symlog", linthresh=0.05)
        style_axis(axis, grid=True)

    natural_axis.set_title("Natural exchange — amplitude difference")
    unnatural_axis.set_title("Unnatural exchange — amplitude difference")
    natural_squared_axis.set_title(
        "Natural exchange — squared-amplitude difference"
    )
    unnatural_squared_axis.set_title(
        "Unnatural exchange — squared-amplitude difference"
    )

    natural_total_axis.set_title("Natural exchange — summed difference")
    unnatural_total_axis.set_title("Unnatural exchange — summed difference")
    natural_total_axis.set_ylabel("Summed signed error")

    natural_squared_total_axis.set_title(
        "Natural exchange — summed squared-amplitude difference"
    )
    unnatural_squared_total_axis.set_title(
        "Unnatural exchange — summed squared-amplitude difference"
    )
    natural_squared_total_axis.set_ylabel(
        r"Summed $|A_{\rm fit}|^2-|A_{\rm quad}|^2$"
    )

    natural_axis.set_ylabel(
        r"$|A_{\rm fit}|-\sqrt{|A_{+}|^2+|A_{-}|^2}$"
    )
    natural_squared_axis.set_ylabel(
        r"$|A_{\rm fit}|^2-|A_{\rm quad}|^2$"
    )

    natural_axis.legend(frameon=False, ncol=2)
    unnatural_axis.legend(frameon=False, ncol=2)
    natural_squared_axis.legend(frameon=False, ncol=2)
    unnatural_squared_axis.legend(frameon=False, ncol=2)

    figure.suptitle(
        rf"Invariant mass $M={mass:.3f}$ GeV"
    )

    save_figure(
        figure,
        error_dir / f"signed_error_mass_{mass:.3f}.pdf",
        show=SHOW_PLOTS,
        dpi=FIGURE_DPI,
    )


# Each moment produces one figure; each wave produces magnitude and phase
# figures; each mass point produces one signed-error summary figure.
figure_count = len(moment_pairs) + 2 * len(waves) + len(masses)
print(f"Read {read_points} of {len(points)} scan points")
print(f"Wrote {figure_count} figures beneath {OUTPUT_DIR}")
