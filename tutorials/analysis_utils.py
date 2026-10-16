"""ROOT-to-NumPy and plotting helpers shared by the tutorial analyses.

The functions here contain mechanics that should behave identically in every
study.  Wave choices, kinematic values, and other physics assumptions remain in
the calling module so they are visible to a reader reproducing that analysis.

ROOT's ``RDataFrame`` is used only at the file boundary.  Converting selected
columns to NumPy arrays early keeps the subsequent statistics and plotting
steps explicit while avoiding repeated scans of the same tree.
"""

import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import ROOT

from .paths import PROJECT_DIR as PROJECT_DIR

MAGNITUDE_PATTERN = re.compile(r"^[ab]_[TL]_\d+_(?:m?\d+)$")
PHASE_PATTERN = re.compile(r"^[ab]phi_[TL]_\d+_(?:m?\d+)$")


def branches(path, tree_name):
    """Return the branch names stored in one ROOT tree.

    ``RDataFrame`` opens the file lazily.  Asking for column names reads the
    tree schema without materialising every entry, which is cheaper than
    exporting data merely to discover which optional diagnostics exist.
    """
    frame = ROOT.RDataFrame(tree_name, str(path))
    return [str(name) for name in frame.GetColumnNames()]


def read_tree(path, tree_name, columns, definitions=None):
    """Read selected scalar ROOT branches into a dictionary of NumPy arrays.

    A definition is a ``(name, C++ expression)`` pair evaluated by ROOT before
    conversion. ``AsNumpy`` then performs one event loop for all requested
    columns; calling it separately for each branch would reread the tree.
    """
    frame = ROOT.RDataFrame(tree_name, str(path))
    for name, expression in definitions or []:
        frame = frame.Define(name, expression)
    return {
        name: np.asarray(values)
        for name, values in frame.AsNumpy(columns=list(columns)).items()
    }


def best_indices(data, objective="chi2", require_fit_ok=True, fraction=0.001):
    """Return the lowest-objective fraction of accepted rows, best first.

    Round the selected count up and retain at least one fit (0.001 = 0.1%).

    When ``fit_ok`` is stored, it records the return value of Minuit's
    minimization call and is the authoritative convergence flag.  A non-zero
    Minuit status can still accompany a converged fit (for example when its
    covariance needed adjustment), so status is used only as a fallback for
    older trees that do not contain ``fit_ok``.
    """
    if not 0.0 < fraction <= 1.0:
        raise ValueError("The selected fit fraction must be in (0, 1]")
    good = np.isfinite(data[objective])
    if require_fit_ok and "fit_ok" in data:
        good &= data["fit_ok"].astype(bool)
    elif require_fit_ok and "status" in data:
        good &= data["status"] == 0
    indices = np.flatnonzero(good)
    if not len(indices):
        raise ValueError(f"No accepted finite rows were found for {objective!r}")
    order = np.argsort(data[objective][indices], kind="stable")
    count = max(1, int(np.ceil(fraction * len(indices))))
    return indices[order[:count]]


def best_index(data, objective="chi2", require_fit_ok=True):
    """Return the accepted row with the smallest finite objective value."""
    return best_indices(data, objective, require_fit_ok)[0]


def range_errors(centre, samples):
    """Asymmetric bars spanning sample minima/maxima along the first axis.

    These are minimization ranges, not statistical uncertainties. Phase callers
    supply their displayed phase convention before calculating the range.
    """
    centre = np.asarray(centre)
    samples = np.asarray(samples)
    return np.stack((centre - np.min(samples, axis=0),
                     np.max(samples, axis=0) - centre))


def phase_branch(magnitude_branch):
    """Return the phase branch paired with an amplitude magnitude branch.

    EMI inserts ``phi`` after the one-character reflectivity label; keeping that
    convention in one helper avoids ad-hoc string replacements in every plot.
    """
    return magnitude_branch[0] + "phi" + magnitude_branch[1:]


def input_moment_branch(fit_branch):
    """Map an output moment name to its experimental-input branch name.

    EMI omits the leading ``R`` when it writes a reconstructed response, so
    both ``H04_*`` and ``H_*`` map back to the published ``RH...`` branches.
    Keeping the rule here makes the Hessian and bootstrap comparisons agree.
    """
    return "R" + fit_branch


def wave_label(branch):
    """Convert a storage name such as ``a_T_1_0`` to a plot label.

    Orbital angular momentum l=0,1,2,... maps to S,P,D,...; ``a`` and ``b`` are
    positive and negative reflectivity, while ``T``/``L`` denotes photon
    orientation and the final field is the magnetic projection m.
    """
    reflectivity, polarization, ell, projection = branch.split("_")
    orbital = ("S", "P", "D", "F", "G")[int(ell)]
    sign = "+" if reflectivity == "a" else "-"
    projection = projection.replace("m", "-")
    return (
        orbital + "$^" + sign + r"_{\mathrm{" + polarization + "}" + projection + "}$"
    )


def circular_mean(angles):
    """Calculate the mean of periodic angles in radians.

    Averaging unit vectors avoids the false mean near zero obtained by directly
    averaging samples clustered on either side of the ``-pi/pi`` boundary.
    """
    return np.angle(np.mean(np.exp(1j * angles)))


def circular_std(angles):
    """Calculate circular standard deviation from the mean resultant length.

    The resultant is the length of the mean unit phasor.  It approaches one for
    concentrated phases and zero for phases spread around the full circle.
    """
    # Round-off can put a theoretically unit-length resultant just above one.
    resultant = np.clip(np.abs(np.mean(np.exp(1j * angles))), 0.0, 1.0)
    return np.sqrt(-2.0 * np.log(resultant))


def wrap_phase(angles):
    """Map angles to the principal interval using a unit-vector argument.

    The sine/cosine form handles arbitrarily large inputs and avoids branching
    separately for values above and below the chosen 2*pi interval.
    """
    return np.arctan2(np.sin(angles), np.cos(angles))


def complex_amplitude(magnitude, phase):
    """Construct a complex amplitude, allowing signed fitted magnitudes.

    A negative radial coordinate represents the same point as a positive
    magnitude with phase shifted by pi.  Canonicalising it here prevents an
    artificial reflection in Argand and phase plots.
    """
    magnitude = np.asarray(magnitude)
    phase = np.asarray(phase) + np.where(magnitude < 0.0, np.pi, 0.0)
    return np.abs(magnitude) * np.exp(1j * phase)


def r_expression(magnitudes):
    """Build ROOT C++ for the longitudinal/transverse ratio ``R``.

    Non-zero longitudinal projections occur as symmetry-related ``+m`` and
    ``-m`` terms in this branch convention, hence their factor of two.
    """
    transverse = [f"({name}*{name})" for name in magnitudes if "_T_" in name]
    longitudinal = []
    for name in magnitudes:
        if "_L_" not in name:
            continue
        projection = name.split("_")[3].replace("m", "-")
        weight = 1 if projection == "0" else 2
        longitudinal.append(f"{weight}*({name}*{name})")
    return "(" + "+".join(longitudinal) + ")/(" + "+".join(transverse) + ")"


def style_axis(axis, grid=False):
    """Use plain inward ticks suitable for scientific figures."""
    axis.minorticks_on()
    axis.tick_params(which="both", direction="in", top=True, right=True)
    axis.tick_params(which="major", length=7, width=1.4)
    axis.tick_params(which="minor", length=4, width=1.0)
    if grid:
        axis.grid(alpha=0.25)


def save_figure(figure, path, show=False, **kwargs):
    """Create parent directories, save, and release a Matplotlib figure.

    Closing figures matters for scans that produce hundreds of PDFs: otherwise
    Matplotlib retains every artist and memory usage grows throughout the job.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, bbox_inches="tight", **kwargs)
    if show:
        plt.show()
    plt.close(figure)


def is_mass_parameter(name):
    """Recognize concise mass parameters and names from older ROOT files."""
    if name.startswith("err__"):
        return False
    return name.startswith(("bw_", "bg_")) or bool(re.fullmatch(
        r".+_(?:mass|width)|.+_[ab](?:_L)?_\d+_m?\d+(?:_km?\d+)?_(?:mag|phase)", name
    ))


def truth_parameters(path):
    """Read injected mass coordinates together, opening the ROOT file once."""
    source = ROOT.TFile.Open(str(path), "READ")
    if not source or source.IsZombie():
        raise OSError(f"Cannot open truth file {path}")
    try:
        return {
            key.GetName(): float(source.Get(key.GetName()).GetVal())
            for key in source.GetListOfKeys()
            if is_mass_parameter(key.GetName())
        }
    finally:
        source.Close()


def read_best_rows(path, tree_name, columns, group="scan_case", require_fit_ok=True,
                   fraction=None):
    """Read quality columns once, then let RDataFrame filter the selected rows.

    Only the small quality table is used to find one minimum per scan case.
    Amplitudes, moments, and covariance columns are read for those rows alone.
    Return selected values and all quality rows (useful for start landscapes).
    With fraction set, read that fraction per case, sorted by case then chi2.
    """
    frame = ROOT.RDataFrame(tree_name, str(path)).Define("scan_entry", "rdfentry_")
    names = {str(name) for name in frame.GetColumnNames()}
    quality_columns = [group, "scan_entry", "chi2"]
    quality_columns += [
        name for name in ("fit_ok", "status", "start", "val") if name in names
    ]
    quality = frame.AsNumpy(columns=quality_columns)
    entries = []
    for case in np.unique(quality[group]):
        mask = quality[group] == case
        rows = {name: values[mask] for name, values in quality.items()}
        selected_indices = (
            [best_index(rows, require_fit_ok=require_fit_ok)]
            if fraction is None
            else best_indices(rows, require_fit_ok=require_fit_ok, fraction=fraction)
        )
        entries.extend(int(rows["scan_entry"][index]) for index in selected_indices)
    entry_list = ",".join(f"{entry}ULL" for entry in sorted(entries))
    selection = (
        "static const std::vector<ULong64_t> entries{" + entry_list + "}; "
        "return std::binary_search(entries.begin(), entries.end(), rdfentry_);"
    )
    selected = frame.Filter(selection).AsNumpy(
        columns=list(dict.fromkeys([group, *(["chi2"] if fraction is not None else []), *columns]))
    )
    if fraction is not None:
        order = np.lexsort((selected["chi2"], selected[group]))
        selected = {name: values[order] for name, values in selected.items()}
    return selected, quality
