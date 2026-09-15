#!/usr/bin/env python3

"""Plot amplitude errors and Argand diagrams for the PhotoTest mass scan."""

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import ROOT


# Settings: edit these values before running the analysis.
MANIFEST_FILE = "OutputFiles/photo_test_mass_scan/scan_points.csv"
SHOW_PLOTS = False
ARGAND_LIMIT = 1.0
ERROR_FIGURE_DPI = 180
ARGAND_FIGURE_DPI = 600


ROOT.gROOT.SetBatch(not SHOW_PLOTS)
PROJECT_DIR = Path(__file__).resolve().parents[1]
manifest_path = PROJECT_DIR / MANIFEST_FILE


def wave_sort_key(name):
    """Sort amplitudes by reflectivity, orbital angular momentum, and m."""
    reflectivity, _, ell, projection = name.split("_")
    return reflectivity, int(ell), int(projection.replace("m", "-"))


def complex_amplitude(magnitude, phase):
    """Construct a complex amplitude from a magnitude and phase."""
    magnitude = float(magnitude)
    phase = float(phase)
    if magnitude < 0.0:
        magnitude = -magnitude
        phase += np.pi
    return magnitude * np.exp(1j * phase)


def read_scan_point(point, waves):
    """Read the truth amplitudes and best finite fit for one scan point."""
    truth_path = PROJECT_DIR / point["truth_file"]
    fit_path = PROJECT_DIR / point["fit_file"]
    if not truth_path.is_file() or not fit_path.is_file():
        return None

    truth_frame = ROOT.RDataFrame("genMoments", str(truth_path))
    fit_frame = ROOT.RDataFrame("fitResults", str(fit_path))
    truth_branches = {str(name) for name in truth_frame.GetColumnNames()}
    fit_branches = {str(name) for name in fit_frame.GetColumnNames()}

    plus = np.full(len(waves), np.nan + 1j * np.nan, dtype=complex)
    minus = np.full(len(waves), np.nan + 1j * np.nan, dtype=complex)
    fitted = np.full(len(waves), np.nan + 1j * np.nan, dtype=complex)

    # List only the branches required for amplitudes shared by the two trees.
    fit_columns = ["chi2"]
    truth_columns = ["aphi_T_2_2_1", "bphi_T_2_2_1"]
    for index, wave in enumerate(waves):
        phase_wave = wave[0] + "phi" + wave[1:]
        needed_fit = {wave, phase_wave}
        needed_truth = {
            wave + "_1", phase_wave + "_1",
            wave + "_m1", phase_wave + "_m1",
        }
        if not needed_fit <= fit_branches or not needed_truth <= truth_branches:
            continue
        fit_columns.extend([wave, phase_wave])
        truth_columns.extend(sorted(needed_truth))

    # RDataFrame reads the selected columns; NumPy chooses the lowest finite chi2.
    truth = truth_frame.AsNumpy(columns=list(dict.fromkeys(truth_columns)))
    fit = fit_frame.AsNumpy(columns=list(dict.fromkeys(fit_columns)))
    finite = np.flatnonzero(np.isfinite(fit["chi2"]))
    if len(finite) == 0:
        return None
    best_entry = finite[np.argmin(fit["chi2"][finite])]
    chi2 = fit["chi2"][best_entry]
    rotations = {
        reflectivity: np.exp(
            -1j * truth[reflectivity + "phi_T_2_2_1"][0]
        )
        for reflectivity in ("a", "b")
    }

    for index, wave in enumerate(waves):
        phase_wave = wave[0] + "phi" + wave[1:]
        if wave not in fit:
            continue
        fitted[index] = complex_amplitude(
            fit[wave][best_entry], fit[phase_wave][best_entry]
        )
        rotation = rotations[wave[0]]
        plus[index] = complex_amplitude(
            truth[wave + "_1"][0], truth[phase_wave + "_1"][0]
        ) * rotation
        minus[index] = complex_amplitude(
            truth[wave + "_m1"][0], truth[phase_wave + "_m1"][0]
        ) * rotation
    return plus, minus, fitted, chi2


def wave_label(name):
    """Convert a_T_1_0 to a compact P-wave plot label."""
    reflectivity, polarization, ell, projection = name.split("_")
    orbital = ("S", "P", "D")[int(ell)]
    sign = "+" if reflectivity == "a" else "-"
    projection = projection.replace("m", "-")
    return orbital + "$^" + sign + r"_{\mathrm{" + polarization + "}" + projection + "}$"


def style_argand_axis(axis):
    """Apply the large-format style used by the fixed-amplitude figures."""
    axis.tick_params(which="major", labelsize=22, width=2.0, length=8, pad=6)
    axis.tick_params(which="minor", width=1.6, length=5)
    for spine in axis.spines.values():
        spine.set_linewidth(2.0)


# Read the scan layout from the runner's CSV manifest.
with manifest_path.open(newline="") as stream:
    points = list(csv.DictReader(stream))

masses = np.array(sorted({float(point["mass_GeV"]) for point in points}))
scales = np.array(sorted({float(point["k_minus_scale"]) for point in points}))
mass_index = {mass: index for index, mass in enumerate(masses)}
scale_index = {scale: index for index, scale in enumerate(scales)}

# Discover the transverse waves from the first available fit file.
first_fit_path = next(
    PROJECT_DIR / point["fit_file"]
    for point in points
    if (PROJECT_DIR / point["fit_file"]).is_file()
)
first_fit = ROOT.RDataFrame("fitResults", str(first_fit_path))
amplitude_pattern = re.compile(r"^[ab]_T_[0-9]+_(?:m[0-9]+|[0-9]+)$")
waves = sorted(
    (
        str(name) for name in first_fit.GetColumnNames()
        if amplitude_pattern.match(str(name))
    ),
    key=wave_sort_key,
)

# Store every scan point in arrays indexed by scale, mass, and wave.
shape = (len(scales), len(masses), len(waves))
truth_plus = np.full(shape, np.nan + 1j * np.nan, dtype=complex)
truth_minus = np.full(shape, np.nan + 1j * np.nan, dtype=complex)
fitted = np.full(shape, np.nan + 1j * np.nan, dtype=complex)
best_chi2 = np.full((len(scales), len(masses)), np.nan)

for point in points:
    result = read_scan_point(point, waves)
    if result is None:
        continue
    row = scale_index[float(point["k_minus_scale"])]
    column = mass_index[float(point["mass_GeV"])]
    (
        truth_plus[row, column],
        truth_minus[row, column],
        fitted[row, column],
        best_chi2[row, column],
    ) = result


# Shared plotting settings and output directory.
figure_dir = manifest_path.parent / "figures"
figure_dir.mkdir(parents=True, exist_ok=True)
scale_colours = plt.cm.plasma(np.linspace(0.08, 0.92, len(scales)))
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 11,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "savefig.dpi": ERROR_FIGURE_DPI,
})

# Plot absolute and relative magnitude errors for each fitted wave.
for index, wave in enumerate(waves):
    reflectivity, _, ell, projection = wave.split("_")
    orbital = ("S", "P", "D")[int(ell)]
    sign = "+" if reflectivity == "a" else "-"
    title = rf"${orbital}_{{{projection.replace('m', '-')}}}^{{{sign}}}$ transverse"

    generated_magnitude = np.sqrt(
        np.abs(truth_plus[:, :, index]) ** 2
        + np.abs(truth_minus[:, :, index]) ** 2
    )
    absolute_error = np.abs(np.abs(fitted[:, :, index]) - generated_magnitude)
    relative_error = np.divide(
        absolute_error,
        generated_magnitude,
        out=np.full_like(absolute_error, np.nan),
        where=generated_magnitude > 1.0e-14,
    )

    figure, axes = plt.subplots(
        2, 1, figsize=(7.2, 7.2), sharex=True, constrained_layout=True
    )
    for row, (scale, colour) in enumerate(zip(scales, scale_colours)):
        axes[0].plot(
            masses, absolute_error[row], marker="o", markersize=3.5,
            color=colour, label=fr"$s={scale:.3f}$",
        )
        axes[1].plot(
            masses, relative_error[row], marker="o", markersize=3.5,
            color=colour,
        )
    axes[0].set_ylabel(r"$\left||A_{\rm fit}|-\sqrt{|T_+|^2+|T_-|^2}\right|$")
    axes[1].set_ylabel("Relative absolute error")
    axes[1].set_xlabel(r"Invariant mass $M$ (GeV)")
    axes[1].set_yscale("log")
    axes[0].set_title(title)
    axes[0].legend(frameon=False, ncol=3)
    figure.savefig(figure_dir / f"photo_test_{wave}_errors.pdf")
    if SHOW_PLOTS:
        plt.show()
    plt.close(figure)


# Plot all waves at each mass, with one column per k-minus scale.
poster_style = {
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 36,
    "axes.titlesize": 42,
    "axes.labelsize": 40,
    "lines.linewidth": 4.0,
    "axes.linewidth": 2.8,
}
wave_colours = plt.colormaps["tab20"](np.linspace(0.0, 1.0, len(waves)))
reflectivity_groups = [("a", "Natural Transverse"), ("b", "Unnatural Transverse")]

for column, mass in enumerate(masses):
    with plt.rc_context(poster_style):
        figure, axes = plt.subplots(
            2, len(scales), figsize=(6 * len(scales), 12),
            sharex=True, sharey=True, constrained_layout=True, squeeze=False,
        )
        for row, (reflectivity, group_title) in enumerate(reflectivity_groups):
            selected = [
                index for index, wave in enumerate(waves)
                if wave.startswith(reflectivity + "_T_")
            ]
            for scale_row, scale in enumerate(scales):
                axis = axes[row, scale_row]
                for index in selected:
                    colour = wave_colours[index]
                    axis.plot(
                        fitted[scale_row, column, index].real,
                        fitted[scale_row, column, index].imag,
                        color=colour, marker="o", linestyle="none",
                        markersize=12, label=wave_label(waves[index]), zorder=3,
                    )
                    for truth in (truth_plus, truth_minus):
                        value = truth[scale_row, column, index]
                        axis.plot(
                            value.real, value.imag, color=colour, marker="*",
                            linestyle="none", markersize=10,
                            markeredgecolor="black", markeredgewidth=0.6,
                            zorder=4,
                        )

                axis.set(xlim=(-ARGAND_LIMIT, ARGAND_LIMIT), ylim=(-ARGAND_LIMIT, ARGAND_LIMIT))
                axis.set_aspect("equal", adjustable="box")
                if scale_row == 0:
                    axis.set_ylabel(group_title + "\nIm", fontsize=20)
                if row == 1:
                    axis.set_xlabel("Re")
                if row == 0:
                    axis.set_title(fr"$s={scale:.3f}$")
                style_argand_axis(axis)

        # Use one legend entry per wave plus entries explaining the markers.
        handles = [
            Line2D(
                [], [], color=wave_colours[index], marker="o",
                linestyle="none", label=wave_label(wave),
            )
            for index, wave in enumerate(waves)
        ]
        handles.extend([
            Line2D([], [], color="black", marker="o", linestyle="none", label="Minimizer"),
            Line2D(
                [], [], color="black", marker="*", markerfacecolor="white",
                linestyle="none", label=r"Set Amp ($k=\pm1$)",
            ),
        ])
        figure.legend(
            handles=handles, loc="outside lower center", frameon=False,
            fontsize=14, ncol=10, reverse=True,
        )
        figure.suptitle(fr"Invariant Mass $M={mass:.3f}$ GeV")
        output_name = f"photo_test_mass_bin_{column:02d}_{mass:.3f}_GeV_argand.pdf"
        figure.savefig(
            figure_dir / output_name,
            dpi=ARGAND_FIGURE_DPI,
            bbox_inches="tight",
            transparent=True,
        )
        if SHOW_PLOTS:
            plt.show()
        plt.close(figure)

print(f"Read {np.isfinite(best_chi2).sum()} of {best_chi2.size} scan points")
print(f"Wrote {len(waves) + len(masses)} figures to {figure_dir}")
