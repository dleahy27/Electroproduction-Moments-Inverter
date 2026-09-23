#!/usr/bin/env python3

"""Plot the generated and fitted PhotoTest mass scan.

The CSV manifest is the authoritative join between mass, k asymmetry, truth
file, and fit file. Arrays are indexed by asymmetry, mass, and wave so every plot uses
the same ordering and missing scan points remain explicit NaNs.

Truth and fit are compared in both amplitude and observable moment space.  The
latter remains meaningful when amplitude coordinates depend on a gauge choice.
"""

import csv
import gc
import re

import matplotlib

matplotlib.use("Agg")

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
FIGURE_DPI = 100
POINT_OFFSET_FRACTION = 0.06
CORRELATION_PREFIX = "corr__"


ROOT.gROOT.SetBatch(not SHOW_PLOTS)


def wave_sort_key(name):
    """Sort the transverse waves by reflectivity, l, and m."""
    reflectivity, _, ell, projection = name.split("_")
    return reflectivity, int(ell), int(projection.replace("m", "-"))


def moment_sort_key(pair):
    """Sort moments by alpha, L, and M."""
    fields = pair[1].replace("H04_", "0_").replace("H_", "").split("_")
    return tuple(int(field.replace("m", "-")) for field in fields)


def offset_masses(masses, number_of_asymmetries):
    """Separate coincident scale markers while retaining mass as x.

    The displacement is purely graphical and is expressed in GeV.  It does not
    modify the mass used to generate or fit any scan point.
    """
    if len(masses) < 2 or number_of_asymmetries < 2:
        return np.tile(masses, (number_of_asymmetries, 1))
    offsets = np.linspace(
        -POINT_OFFSET_FRACTION / 2,
        POINT_OFFSET_FRACTION / 2,
        number_of_asymmetries,
    )
    return masses[None, :] + offsets[:, None]


def read_scan_point(point, waves, moment_pairs, correlation_parameters):
    """Read one generated point and its lowest-chi-square accepted fit.

    Returns complex k=+1 truth, k=-1 truth, fitted amplitudes, generated moments,
    and fitted moments in the caller-provided ordering.  ``None`` marks an
    incomplete point so a resumed scan can still be analysed without shifting
    the indices of every later mass.
    """
    truth_path = PROJECT_DIR / point["truth_file"]
    fit_path = PROJECT_DIR / point["fit_file"]
    if not truth_path.is_file() or not fit_path.is_file():
        return None

    # Assemble the exact branch set before starting either ROOT event loop.  A
    # magnitude and phase are needed to reconstruct each complex amplitude.
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
    fit_columns = [
        branch
        for wave in waves
        for branch in (wave, phase_branch(wave))
    ]
    fit_columns += [
        fit_name for _, fit_name in moment_pairs
    ]
    fit_columns += sorted(
        name for name in fit_names if name.startswith(CORRELATION_PREFIX)
    )

    truth = ROOT.RDataFrame("genMoments", str(truth_path)).AsNumpy(
        columns=truth_columns
    )
    # Selecting the minimum needs only a few small columns.  Loading hundreds
    # of amplitude, error, moment, and correlation branches for every random
    # start would otherwise dominate memory use for a large scan.
    fit_quality = fit_frame.AsNumpy(columns=status_columns)
    selected_entry = int(best_index(fit_quality))
    fit = fit_frame.Range(selected_entry, selected_entry + 1).AsNumpy(
        columns=fit_columns
    )
    entry = 0

    # Preallocate fixed-shape arrays so all scan points can later be inserted
    # into the common [asymmetry, mass, wave] data cube without Python object arrays.
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

    # Moments require no phase alignment: they are bilinear observables and are
    # invariant under the common rotations applied above for amplitude plots.
    generated_moments = np.array([
        truth[truth_name][0] for truth_name, _ in moment_pairs
    ])
    fitted_moments = np.array([
        fit[fit_name][entry] for _, fit_name in moment_pairs
    ])
    correlation = read_correlation_matrix(
        fit,
        entry,
        correlation_parameters,
    )

    result = (
        plus,
        minus,
        fitted,
        generated_moments,
        fitted_moments,
        correlation,
    )

    del truth
    del fit
    del fit_quality
    del fit_frame

    gc.collect()

    return result


def asymmetry_handles(asymmetries, colours):
    """Make a compact colour legend shared by the scan figures."""
    return [
        Line2D(
            [], [], marker="o", linestyle="none", color=colour,
            label=fr"$s={asymmetry:.3g}$",
        )
        for asymmetry, colour in zip(asymmetries, colours)
    ]


def finish_mass_axis(axis, title, ylabel):
    """Apply the common mass-scan labels and styling."""
    axis.set(xlabel=r"Invariant mass $M$ (GeV)", ylabel=ylabel, title=title)
    #axis.set_xscale("log")
    style_axis(axis, grid=True)


def save_scan_figure(figure, path):
    """Save one scan figure and promptly collect its cyclic artist objects."""
    save_figure(
        figure, path, show=SHOW_PLOTS, dpi=FIGURE_DPI
    )
    plt.close(figure)
    del figure
    gc.collect()


def read_correlation_matrix(fit, entry, parameters):
    """Read pairwise correlation branches from one selected fit-result row."""
    matrix = np.full((len(parameters), len(parameters)), np.nan)
    np.fill_diagonal(matrix, 1.0)
    indices = {name: index for index, name in enumerate(parameters)}
    for branch, values in fit.items():
        if not branch.startswith(CORRELATION_PREFIX):
            continue
        pair = branch.removeprefix(CORRELATION_PREFIX).split("__", maxsplit=1)
        if len(pair) != 2:
            continue
        name_x, name_y = pair
        if name_x not in indices or name_y not in indices:
            continue
        value = float(values[entry])
        if not np.isfinite(value) or abs(value) > 1.0 + 1.0e-10:
            continue
        x_index = indices[name_x]
        y_index = indices[name_y]
        matrix[x_index, y_index] = np.clip(value, -1.0, 1.0)
        matrix[y_index, x_index] = matrix[x_index, y_index]
    return matrix


def parameter_label(name):
    """Return a compact label for an amplitude magnitude or phase parameter."""
    if "phi_" in name:
        magnitude_name = name[0] + name[name.index("_"):]
        return rf"$\phi$({wave_label(magnitude_name)})"
    return wave_label(name)


# Read the scan definition written by the runner.  The manifest, rather than a
# filename parser, is the authoritative mapping from physical coordinates to
# files and remains reliable if the filename formatting changes.
with MANIFEST_FILE.open(newline="") as stream:
    points = list(csv.DictReader(stream))

masses = np.array(sorted({float(point["mass_GeV"]) for point in points}))
asymmetries = np.array(sorted({float(point["k_asymmetry"]) for point in points}))
# Dictionaries turn manifest values into dense array indices.  Values originate
# from the same decimal strings, so exact lookup is safe within this process.
mass_index = {mass: index for index, mass in enumerate(masses)}
asymmetry_index = {
    asymmetry: index for index, asymmetry in enumerate(asymmetries)
}

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
    (
        name
        for name in fit_names
        if wave_pattern.fullmatch(name)
        and phase_branch(name) in fit_names
    ),
    key=wave_sort_key,
)
if not waves:
    raise RuntimeError("The fit tree contains no amplitude magnitude/phase pairs")
correlation_parameters = [
    parameter
    for wave in waves
    for parameter in (wave, phase_branch(wave))
]

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

# Store all results as arrays indexed by asymmetry, mass, and observable.
amplitude_shape = (len(asymmetries), len(masses), len(waves))
moment_shape = (len(asymmetries), len(masses), len(moment_pairs))
truth_plus = np.full(amplitude_shape, np.nan + 1j * np.nan, dtype=complex)
truth_minus = np.full(amplitude_shape, np.nan + 1j * np.nan, dtype=complex)
fitted_amplitudes = np.full(
    amplitude_shape, np.nan + 1j * np.nan, dtype=complex
)
truth_moments = np.full(moment_shape, np.nan)
fitted_moments = np.full(moment_shape, np.nan)
correlation_shape = (
    len(asymmetries), len(masses),
    len(correlation_parameters), len(correlation_parameters),
)
correlations = np.full(correlation_shape, np.nan)

read_points = 0
for point in points:
    result = read_scan_point(
        point, waves, moment_pairs, correlation_parameters
    )
    if result is None:
        continue
    asymmetry_row = asymmetry_index[float(point["k_asymmetry"])]
    mass_column = mass_index[float(point["mass_GeV"])]
    truth_plus[asymmetry_row, mass_column] = result[0]
    truth_minus[asymmetry_row, mass_column] = result[1]
    fitted_amplitudes[asymmetry_row, mass_column] = result[2]
    truth_moments[asymmetry_row, mass_column] = result[3]
    fitted_moments[asymmetry_row, mass_column] = result[4]
    correlations[asymmetry_row, mass_column] = result[5]
    read_points += 1


# Use colour for asymmetry and a small horizontal displacement for coincident
# nominal masses. The displacement is multiplicative because the mass grid is
# logarithmic.
colours = plt.colormaps["viridis"](np.linspace(0.08, 0.92, len(asymmetries)))
plot_masses = offset_masses(masses, len(asymmetries))
#log_spacing = np.median(np.diff(np.log(masses))) if len(masses) > 1 else 0.0
#marker_gap = 0.0
colour_legend = asymmetry_handles(asymmetries, colours)
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 11,
})


# Plot each input/output moment comparison in its own file.
moment_dir = OUTPUT_DIR / "moments"
for index, (_, moment) in enumerate(moment_pairs):
    figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=False)
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
    save_scan_figure(figure, moment_dir / f"{moment}.pdf")


# Plot each amplitude magnitude and phase in separate mass-scan figures.
amplitude_dir = OUTPUT_DIR / "amplitudes"
error_dir = OUTPUT_DIR / "signed_errors"
for index, wave in enumerate(waves):
    title = wave_label(wave)

    figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=False)
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
    save_scan_figure(figure, amplitude_dir / f"{wave}_magnitude.pdf")

    figure, axis = plt.subplots(figsize=(16, 8), constrained_layout=False)
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
    save_scan_figure(figure, amplitude_dir / f"{phase_branch(wave)}.pdf")


# Plot the Hessian correlation matrix stored with each selected fit result.
correlation_dir = OUTPUT_DIR / "correlations"
correlation_figure_count = 0
labels = [parameter_label(name) for name in correlation_parameters]
for asymmetry_row, asymmetry in enumerate(asymmetries):
    for mass_column, mass in enumerate(masses):
        matrix = correlations[asymmetry_row, mass_column]
        off_diagonal = ~np.eye(len(matrix), dtype=bool)
        if not np.any(np.isfinite(matrix[off_diagonal])):
            continue
        figure, axis = plt.subplots(figsize=(18, 16), constrained_layout=False)
        image = axis.imshow(
            matrix, vmin=-1.0, vmax=1.0, cmap="coolwarm", interpolation="none"
        )
        axis.set_xticks(np.arange(len(labels)), labels, rotation=90)
        axis.set_yticks(np.arange(len(labels)), labels)
        axis.set_title(
            rf"Amplitude correlations: $M={mass:.3f}$ GeV, $s={asymmetry:.3g}$"
        )
        figure.colorbar(image, ax=axis, label="Correlation coefficient")
        save_scan_figure(
            figure,
            correlation_dir
            / f"asymmetry_{asymmetry_row:02d}_mass_{mass_column:02d}_correlation.pdf",
        )
        correlation_figure_count += 1


# Plot the signed amplitude error as a function of the k-sector asymmetry.
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
        constrained_layout=False,
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
            asymmetries,
            signed_error[:, wave_index],
            marker="o",
            color=colour,
            label=wave_label(wave),
        )
        natural_squared_axis.plot(
            asymmetries,
            signed_squared_error[:, wave_index],
            marker="o",
            color=colour,
            label=wave_label(wave),
        )

    for colour, (wave_index, wave) in zip(
        unnatural_colours, unnatural_waves
    ):
        unnatural_axis.plot(
            asymmetries,
            signed_error[:, wave_index],
            marker="o",
            color=colour,
            label=wave_label(wave),
        )
        unnatural_squared_axis.plot(
            asymmetries,
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
            asymmetries,
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
        axis.set_xlabel(r"K-sector asymmetry $s$")
        axis.set_xlim(-1.05, 1.05)
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

    save_scan_figure(
        figure,
        error_dir / f"signed_error_mass_{mass:.3f}.pdf",
    )



# Each moment produces one figure; each wave produces magnitude and phase
# figures; each mass point produces one signed-error summary figure.
figure_count = (
    len(moment_pairs) + 2 * len(waves) + len(masses)
    + correlation_figure_count
)
print(f"Read {read_points} of {len(points)} scan points")
print(f"Wrote {figure_count} figures beneath {OUTPUT_DIR}")
