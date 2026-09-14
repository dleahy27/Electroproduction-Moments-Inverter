#!/usr/bin/env python3

"""Plot amplitude errors and Argand diagrams for the PhotoTest mass scan."""

import argparse
import csv
import re
from pathlib import Path

try:
    import ROOT
except (ImportError, RuntimeError) as error:
    raise SystemExit(
        "PyROOT is required. Load the same ROOT environment used to build EMI "
        "and run this script with the Python version for which ROOT was built.\n"
        f"Original import error: {error}"
    ) from error

ROOT.PyConfig.IgnoreCommandLineOptions = True

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import numpy as np


parser = argparse.ArgumentParser(
    description="Analyse the PhotoTest two-k truth versus single-k fits."
)
parser.add_argument("--manifest", type=Path, default=None)
parser.add_argument("--show", action="store_true")
args = parser.parse_args()
ROOT.gROOT.SetBatch(not args.show)

project_dir = Path(__file__).resolve().parents[1]
manifest = args.manifest.resolve() if args.manifest else (
    project_dir / "OutputFiles" / "photo_test_mass_scan" / "scan_points.csv"
)
if not manifest.is_file():
    raise SystemExit(
        f"Missing {manifest}; run scripts/run-photo-test-mass-scan.py first"
    )

with manifest.open(newline="") as stream:
    points = list(csv.DictReader(stream))
if not points:
    raise SystemExit(f"The manifest is empty: {manifest}")

masses = np.array(sorted({float(point["mass_GeV"]) for point in points}))
scales = np.array(sorted({float(point["k_minus_scale"]) for point in points}))
mass_lookup = {value: index for index, value in enumerate(masses)}
scale_lookup = {value: index for index, value in enumerate(scales)}

available_fit_files = [
    project_dir / point["fit_file"]
    for point in points
    if (project_dir / point["fit_file"]).is_file()
]
if not available_fit_files:
    raise SystemExit("None of the fit files listed in the manifest exists")
first_fit_file = available_fit_files[0]
first_fit_root = ROOT.TFile.Open(str(first_fit_file))
if not first_fit_root or first_fit_root.IsZombie():
    raise SystemExit(f"Could not open {first_fit_file}")
first_fit_tree = first_fit_root.Get("fitResults")
if not first_fit_tree:
    raise SystemExit(f"Missing fitResults in {first_fit_file}")
fit_branch_names = {
    branch.GetName() for branch in first_fit_tree.GetListOfBranches()
}
amplitude_pattern = re.compile(r"^[ab]_T_[0-9]+_(?:m[0-9]+|[0-9]+)$")
waves = sorted(
    (name for name in fit_branch_names if amplitude_pattern.match(name)),
    key=lambda name: (
        name.split("_")[0],
        int(name.split("_")[2]),
        int(name.split("_")[3].replace("m", "-")),
    ),
)
first_fit_root.Close()
if not waves:
    raise SystemExit("No fitted transverse amplitudes were found")

shape = (len(scales), len(masses), len(waves))
truth_plus = np.full(shape, np.nan + 1j * np.nan, dtype=complex)
truth_minus = np.full(shape, np.nan + 1j * np.nan, dtype=complex)
fitted = np.full(shape, np.nan + 1j * np.nan, dtype=complex)
best_chi2 = np.full((len(scales), len(masses)), np.nan)

missing = []
for point in points:
    mass = float(point["mass_GeV"])
    scale = float(point["k_minus_scale"])
    mass_index = mass_lookup[mass]
    scale_index = scale_lookup[scale]
    truth_file = project_dir / point["truth_file"]
    fit_file = project_dir / point["fit_file"]
    if not truth_file.is_file() or not fit_file.is_file():
        missing.append((truth_file, fit_file))
        continue

    truth_root = ROOT.TFile.Open(str(truth_file))
    fit_root = ROOT.TFile.Open(str(fit_file))
    if (not truth_root or truth_root.IsZombie()
            or not fit_root or fit_root.IsZombie()):
        missing.append((truth_file, fit_file))
        if truth_root:
            truth_root.Close()
        if fit_root:
            fit_root.Close()
        continue

    truth_tree = truth_root.Get("genMoments")
    fit_tree = fit_root.Get("fitResults")
    if not truth_tree or not fit_tree or truth_tree.GetEntries() < 1:
        missing.append((truth_file, fit_file))
        truth_root.Close()
        fit_root.Close()
        continue

    truth_tree.GetEntry(0)
    fit_names = {branch.GetName() for branch in fit_tree.GetListOfBranches()}
    truth_names = {
        branch.GetName() for branch in truth_tree.GetListOfBranches()
    }
    finite_entries = []
    for entry in range(fit_tree.GetEntries()):
        fit_tree.GetEntry(entry)
        chi2 = float(getattr(fit_tree, "chi2"))
        if np.isfinite(chi2):
            finite_entries.append((chi2, entry))
    if not finite_entries:
        missing.append((truth_file, fit_file))
        truth_root.Close()
        fit_root.Close()
        continue

    chi2, entry = min(finite_entries)
    fit_tree.GetEntry(entry)
    best_chi2[scale_index, mass_index] = chi2

    # The single-k fit fixes the highest a and b phases independently. Rotate
    # each generated reflectivity into the same plotting convention. This does
    # not change the quadrature truth magnitude used for the error plots.
    truth_rotations = {}
    for reflectivity in ("a", "b"):
        reference_phase = float(
            getattr(truth_tree, reflectivity + "phi_T_2_2_1")
        )
        truth_rotations[reflectivity] = np.exp(-1j * reference_phase)

    for wave_index, wave in enumerate(waves):
        phase_wave = wave[0] + "phi" + wave[1:]
        required_fit = {wave, phase_wave}
        required_truth = {
            wave + "_1", phase_wave + "_1",
            wave + "_m1", phase_wave + "_m1",
        }
        if not required_fit <= fit_names or not required_truth <= truth_names:
            continue

        fit_magnitude = float(getattr(fit_tree, wave))
        fit_phase = float(getattr(fit_tree, phase_wave))
        if fit_magnitude < 0.0:
            fit_magnitude = -fit_magnitude
            fit_phase += np.pi
        fitted[scale_index, mass_index, wave_index] = (
            fit_magnitude * np.exp(1j * fit_phase)
        )

        rotation = truth_rotations[wave[0]]
        plus_magnitude = float(getattr(truth_tree, wave + "_1"))
        plus_phase = float(getattr(truth_tree, phase_wave + "_1"))
        minus_magnitude = float(getattr(truth_tree, wave + "_m1"))
        minus_phase = float(getattr(truth_tree, phase_wave + "_m1"))
        truth_plus[scale_index, mass_index, wave_index] = (
            plus_magnitude * np.exp(1j * plus_phase) * rotation
        )
        truth_minus[scale_index, mass_index, wave_index] = (
            minus_magnitude * np.exp(1j * minus_phase) * rotation
        )

    truth_root.Close()
    fit_root.Close()

if missing:
    print(f"Warning: {len(missing)} scan points could not be read")
if not np.isfinite(fitted.real).any():
    raise SystemExit("No finite fitted amplitudes were read")

figure_dir = manifest.parent / "figures"
figure_dir.mkdir(parents=True, exist_ok=True)
scale_colours = plt.cm.plasma(np.linspace(0.08, 0.92, len(scales)))
mass_colours = plt.cm.viridis(LogNorm(masses.min(), masses.max())(masses))

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 11,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "savefig.dpi": 180,
})

for wave_index, wave in enumerate(waves):
    reflectivity, _, ell_text, projection_text = wave.split("_")
    orbital = ("S", "P", "D")[int(ell_text)]
    projection = int(projection_text.replace("m", "-"))
    reflectivity_sign = "+" if reflectivity == "a" else "-"
    wave_title = (
        rf"${orbital}_{{{projection}}}^{{{reflectivity_sign}}}$ transverse"
    )

    generated_magnitude = np.sqrt(
        np.abs(truth_plus[:, :, wave_index]) ** 2
        + np.abs(truth_minus[:, :, wave_index]) ** 2
    )
    fitted_magnitude = np.abs(fitted[:, :, wave_index])
    absolute_error = np.abs(fitted_magnitude - generated_magnitude)
    relative_error = np.divide(
        absolute_error,
        generated_magnitude,
        out=np.full_like(absolute_error, np.nan),
        where=generated_magnitude > 1.0e-14,
    )

    figure, axes = plt.subplots(
        2, 1, figsize=(7.2, 7.2), sharex=True, constrained_layout=True
    )
    for scale_index, (scale, colour) in enumerate(zip(scales, scale_colours)):
        axes[0].plot(
            masses, absolute_error[scale_index], marker="o", markersize=3.5,
            color=colour, label=fr"$f={scale:.3f}$",
        )
        axes[1].plot(
            masses, relative_error[scale_index], marker="o", markersize=3.5,
            color=colour,
        )
    axes[0].set_ylabel(
        r"$\left||A_{\rm fit}|-\sqrt{|T_+|^2+|T_-|^2}\right|$"
    )
    axes[1].set_ylabel("Relative absolute error")
    axes[1].set_xlabel(r"Invariant mass $w$ (GeV)")
    axes[1].set_xscale("log")
    axes[0].set_title(wave_title)
    axes[0].legend(frameon=False, ncol=3)
    figure.savefig(figure_dir / f"errors_{wave}.pdf")
    if args.show:
        plt.show()
    plt.close(figure)

    figure, axes = plt.subplots(
        1, len(scales), figsize=(3.4 * len(scales), 3.7),
        sharex=True, sharey=True, constrained_layout=True,
    )
    if len(scales) == 1:
        axes = np.array([axes])

    all_values = np.concatenate([
        truth_plus[:, :, wave_index].ravel(),
        truth_minus[:, :, wave_index].ravel(),
        fitted[:, :, wave_index].ravel(),
    ])
    finite_values = all_values[
        np.isfinite(all_values.real) & np.isfinite(all_values.imag)
    ]
    if finite_values.size == 0:
        print(f"Warning: no finite points for {wave}; skipping its figures")
        continue
    limit = 1.08 * max(
        np.max(np.abs(finite_values.real)),
        np.max(np.abs(finite_values.imag)),
        1.0e-3,
    )

    for scale_index, (axis, scale) in enumerate(zip(axes, scales)):
        plus = truth_plus[scale_index, :, wave_index]
        minus = truth_minus[scale_index, :, wave_index]
        fit = fitted[scale_index, :, wave_index]
        axis.plot(plus.real, plus.imag, color="0.72", linewidth=0.8, zorder=1)
        axis.plot(minus.real, minus.imag, color="0.82", linewidth=0.8, zorder=1)
        axis.plot(fit.real, fit.imag, color="0.55", linewidth=0.8, zorder=1)
        axis.scatter(
            plus.real, plus.imag, c=mass_colours, marker="*", s=62,
            edgecolors="black", linewidths=0.35, zorder=3,
        )
        axis.scatter(
            minus.real, minus.imag, c=mass_colours, marker="*", s=34,
            edgecolors="white", linewidths=0.45, zorder=3,
        )
        axis.scatter(
            fit.real, fit.imag, c=mass_colours, marker="o", s=22,
            edgecolors="black", linewidths=0.35, zorder=4,
        )
        axis.axhline(0.0, color="0.45", linewidth=0.7)
        axis.axvline(0.0, color="0.45", linewidth=0.7)
        axis.set_xlim(-limit, limit)
        axis.set_ylim(-limit, limit)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(fr"$f={scale:.3f}$")
        axis.set_xlabel(r"Re $A$")
    axes[0].set_ylabel(r"Im $A$")
    figure.suptitle(wave_title)

    legend = [
        Line2D([], [], marker="*", linestyle="", markersize=10,
               markerfacecolor="0.55", markeredgecolor="black", label=r"truth $k=+$"),
        Line2D([], [], marker="*", linestyle="", markersize=7,
               markerfacecolor="0.55", markeredgecolor="white", label=r"truth $k=-$"),
        Line2D([], [], marker="o", linestyle="", markersize=6,
               markerfacecolor="0.55", markeredgecolor="black", label="single-k fit"),
    ]
    axes[-1].legend(handles=legend, frameon=False, loc="best", fontsize=9)
    scalar_map = plt.cm.ScalarMappable(
        norm=LogNorm(masses.min(), masses.max()), cmap="viridis"
    )
    figure.colorbar(
        scalar_map, ax=axes, location="bottom", shrink=0.5, pad=0.08,
        label=r"Invariant mass $w$ (GeV)",
    )
    figure.savefig(figure_dir / f"argand_{wave}.pdf")
    if args.show:
        plt.show()
    plt.close(figure)

print(f"Read {np.isfinite(best_chi2).sum()} of {best_chi2.size} scan points")
print(f"Wrote {2 * len(waves)} figures to {figure_dir}")
