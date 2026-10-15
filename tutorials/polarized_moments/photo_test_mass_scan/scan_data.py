"""Load a PhotoTest mass scan once; plotting modules receive named arrays.

The merged ROOT trees carry scan_case identifiers from the CSV manifest.
RDataFrame selects fit rows before large amplitude/covariance columns are read.
All array axes are [asymmetry, mass, observable]; absent points remain NaNs.
"""

import csv
import re
from dataclasses import dataclass, replace
from pathlib import Path
import numpy as np
from tutorials.analysis_utils import (
    branches,
    read_tree,
    read_best_rows,
    complex_amplitude,
    phase_branch,
    wave_label,
)

CORRELATION_PREFIX = "corr__"


def wave_sort_key(name):
    """Sort the transverse waves by reflectivity, l, and m."""
    reflectivity, _, ell, projection = name.split("_")
    return reflectivity, int(ell), int(projection.replace("m", "-"))


def moment_sort_key(pair):
    """Sort moments by alpha, L, and M."""
    fields = pair[1].replace("H04_", "0_").replace("H_", "").split("_")
    return tuple(int(field.replace("m", "-")) for field in fields)


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
        magnitude_name = name[0] + name[name.index("_") :]
        return rf"$\phi$({wave_label(magnitude_name)})"
    return wave_label(name)


@dataclass
class MassScan:
    """Physics arrays and their axis labels, shared by all scan plotters."""

    masses: np.ndarray
    asymmetries: np.ndarray
    waves: list
    moment_pairs: list
    truth_plus: np.ndarray
    truth_minus: np.ndarray
    fitted_amplitudes: np.ndarray
    truth_moments: np.ndarray
    fitted_moments: np.ndarray
    minimum_chi2: np.ndarray
    correlation_parameters: list
    correlations: np.ndarray
    unnormalised_intensities: np.ndarray
    moment_scales: np.ndarray
    read_points: int
    total_points: int
    amplitude_samples: np.ndarray
    moment_samples: np.ndarray


def load_scan(manifest, require_unnormalised=False):
    """Read truth, best-fit rows, and optional correlations from one ROOT file."""
    with Path(manifest).open(newline="") as stream:
        points = list(csv.DictReader(stream))
    if not points:
        raise RuntimeError("The scan manifest is empty")
    path = Path(manifest).with_name("scan.root")
    if not path.is_file():
        raise FileNotFoundError(
            "Merge existing results with photo_test_mass_scan.run --merge-only"
        )
    truth_names = set(branches(path, "genMoments"))
    fit_names = set(branches(path, "fitResults"))
    wave_pattern = re.compile(r"^[ab]_T_\d+_(?:m?\d+)$")
    waves = sorted(
        (
            name
            for name in fit_names
            if wave_pattern.fullmatch(name) and phase_branch(name) in fit_names
        ),
        key=wave_sort_key,
    )
    if not waves:
        raise RuntimeError("The fit tree contains no amplitude magnitude/phase pairs")
    moment_pairs = []
    for name in truth_names:
        match = re.fullmatch(r"RH_(\d+)_0_0_(\d+)_(m?\d+)", name)
        if match:
            alpha, ell, projection = match.groups()
            fit_name = (
                f"H04_{ell}_{projection}"
                if alpha == "0"
                else f"H_{alpha}_{ell}_{projection}"
            )
            if fit_name in fit_names:
                moment_pairs.append((name, fit_name))
    moment_pairs.sort(key=moment_sort_key)
    parameters = [branch for wave in waves for branch in (wave, phase_branch(wave))]
    raw_columns = [
        "raw_photo_intensity",
        "raw_photo_intensity_l0",
        "raw_photo_intensity_l1",
        "raw_photo_intensity_l2",
        "raw_photo_H_0_0_0",
        "normalised_photo_H_0_0_0",
    ]
    has_raw = set(raw_columns).issubset(truth_names)
    if require_unnormalised and not has_raw:
        raise RuntimeError(
            "Regenerate truth with the current EMI executable for unnormalised plots"
        )
    truth_columns = [
        f"{branch}_{sector}" for branch in parameters for sector in ("1", "m1")
    ]
    truth_columns += [name for name, _ in moment_pairs]
    if has_raw:
        truth_columns += raw_columns
    truth = read_tree(path, "genMoments", ["scan_case", *truth_columns])
    fit_columns = parameters + [name for _, name in moment_pairs]
    fit_columns += sorted(
        name for name in fit_names if name.startswith(CORRELATION_PREFIX)
    )
    fit_columns += [name for name in ("val", "chi2") if name in fit_names]
    fits, _ = read_best_rows(path, "fitResults", fit_columns, fraction=0.001)
    truth_rows = {int(case): index for index, case in enumerate(truth["scan_case"])}
    fit_rows = {}
    for index, case in enumerate(fits["scan_case"]):
        fit_rows.setdefault(int(case), []).append(index)

    masses = np.array(sorted({float(point["mass_GeV"]) for point in points}))
    asymmetries = np.array(sorted({float(point["k_asymmetry"]) for point in points}))
    mass_rows = {value: index for index, value in enumerate(masses)}
    asymmetry_rows = {value: index for index, value in enumerate(asymmetries)}
    shape = (len(asymmetries), len(masses))
    amplitude_shape = (*shape, len(waves))
    plus = np.full(amplitude_shape, np.nan + 1j * np.nan)
    minus = np.full(amplitude_shape, np.nan + 1j * np.nan)
    fitted = np.full(amplitude_shape, np.nan + 1j * np.nan)
    moments = np.full((*shape, len(moment_pairs)), np.nan)
    reconstructed = np.full_like(moments, np.nan)
    minimum_chi2 = np.full(shape, np.nan)
    correlations = np.full((*shape, len(parameters), len(parameters)), np.nan)
    intensities = np.full((*shape, 4), np.nan)
    scales = np.full(shape, np.nan)
    amplitude_samples = np.empty(shape, dtype=object)
    moment_samples = np.empty(shape, dtype=object)
    for cell in np.ndindex(shape):
        amplitude_samples[cell] = np.full((1, len(waves)), np.nan + 1j * np.nan)
        moment_samples[cell] = np.full((1, len(moment_pairs)), np.nan)
    completed = 0
    for case, point in enumerate(points):
        if case not in truth_rows or case not in fit_rows:
            continue
        t, selected = truth_rows[case], fit_rows[case]
        f = selected[0]
        row = asymmetry_rows[float(point["k_asymmetry"])]
        column = mass_rows[float(point["mass_GeV"])]
        # Align each reflectivity with the unpolarized fit's phase reference.
        rotations = {
            refl: np.exp(-1j * truth[f"{refl}phi_T_2_2_1"][t]) for refl in ("a", "b")
        }
        for index, wave in enumerate(waves):
            phase = phase_branch(wave)
            plus[row, column, index] = (
                complex_amplitude(truth[f"{wave}_1"][t], truth[f"{phase}_1"][t])
                * rotations[wave[0]]
            )
            minus[row, column, index] = (
                complex_amplitude(truth[f"{wave}_m1"][t], truth[f"{phase}_m1"][t])
                * rotations[wave[0]]
            )
            fitted[row, column, index] = complex_amplitude(
                fits[wave][f], fits[phase][f]
            )
        amplitude_samples[row, column] = np.column_stack([
            complex_amplitude(fits[wave][selected], fits[phase_branch(wave)][selected])
            for wave in waves
        ])
        moment_samples[row, column] = np.column_stack([
            fits[name][selected] for _, name in moment_pairs
        ])
        moments[row, column] = [truth[name][t] for name, _ in moment_pairs]
        reconstructed[row, column] = [fits[name][f] for _, name in moment_pairs]
        correlations[row, column] = read_correlation_matrix(fits, f, parameters)
        # Preserve the original plotted quantity: val is the stored objective.
        minimum_chi2[row, column] = fits["val" if "val" in fits else "chi2"][f]
        if has_raw:
            raw = float(truth["raw_photo_H_0_0_0"][t])
            normalized = float(truth["normalised_photo_H_0_0_0"][t])
            if (
                not np.isfinite(raw)
                or not np.isfinite(normalized)
                or raw <= 0
                or normalized <= 0
            ):
                raise RuntimeError(f"Invalid intensity scale for scan case {case}")
            scales[row, column] = raw / normalized
            intensities[row, column] = [truth[name][t] for name in raw_columns[:4]]
        completed += 1
    return MassScan(
        masses,
        asymmetries,
        waves,
        moment_pairs,
        plus,
        minus,
        fitted,
        moments,
        reconstructed,
        minimum_chi2,
        parameters,
        correlations,
        intensities,
        scales,
        completed,
        len(points),
        amplitude_samples,
        moment_samples,
    )


def unnormalised(scan):
    """Restore the common amplitude scale; bilinear moments use its square."""
    scale = np.sqrt(scan.moment_scales)[:, :, None]
    amplitude_samples = np.empty_like(scan.amplitude_samples)
    moment_samples = np.empty_like(scan.moment_samples)
    for cell in np.ndindex(scan.moment_scales.shape):
        amplitude_samples[cell] = scan.amplitude_samples[cell] * np.sqrt(scan.moment_scales[cell])
        moment_samples[cell] = scan.moment_samples[cell] * scan.moment_scales[cell]
    return replace(
        scan,
        amplitude_samples=amplitude_samples,
        moment_samples=moment_samples,
        truth_plus=scan.truth_plus * scale,
        truth_minus=scan.truth_minus * scale,
        fitted_amplitudes=scan.fitted_amplitudes * scale,
        truth_moments=scan.truth_moments * scan.moment_scales[:, :, None],
        fitted_moments=scan.fitted_moments * scan.moment_scales[:, :, None],
    )


def sample_errors(centre, samples, transform=lambda values: values):
    """Return [lower/upper, asymmetry, mass, observable] range bars."""
    from tutorials.analysis_utils import range_errors
    errors = np.full((2, *centre.shape), np.nan)
    for cell in np.ndindex(samples.shape):
        errors[(slice(None), *cell)] = range_errors(centre[cell], transform(samples[cell]))
    return errors


def intensity_bounds(scan):
    """Min/max of wave-summed intensities, evaluated within each selected fit."""
    lower = np.full(scan.amplitude_samples.shape, np.nan)
    upper = np.full_like(lower, np.nan)
    for cell in np.ndindex(lower.shape):
        intensity = np.sum(np.abs(scan.amplitude_samples[cell]) ** 2, axis=1)
        lower[cell], upper[cell] = np.min(intensity), np.max(intensity)
    return lower, upper


def comparison_errors(values, reference, lower, upper, reference_lower, reference_upper):
    """Range bars for differences/ratios of independently fitted scan cases."""
    difference = values - reference
    difference_errors = np.stack((
        difference - (lower - reference_upper),
        (upper - reference_lower) - difference,
    ))
    ratio = np.divide(values, reference, out=np.full_like(values, np.nan),
                      where=reference > 0.0)
    ratio_lower = np.divide(lower, reference_upper, out=np.full_like(values, np.nan),
                            where=reference_upper > 0.0)
    ratio_upper = np.divide(upper, reference_lower, out=np.full_like(values, np.nan),
                            where=reference_lower > 0.0)
    return difference_errors, np.stack((ratio - ratio_lower, ratio_upper - ratio))
