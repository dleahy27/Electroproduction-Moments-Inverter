"""ROOT-to-NumPy and plotting helpers shared by the tutorial analyses.

The functions here contain mechanics that should behave identically in every
study.  Wave choices, kinematic values, and other physics assumptions remain in
the calling module so they are visible to a reader reproducing that analysis.
"""

import re
import matplotlib.pyplot as plt
import numpy as np
import ROOT

from .paths import PROJECT_DIR


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


def best_index(data, objective="chi2", require_fit_ok=True):
    """Return the accepted row with the smallest finite objective value.

    A Minuit invocation can finish with finite parameters while still carrying
    a failure status.  Applying all available quality flags prevents such a row
    from winning solely because its reported objective happens to be small.
    """
    good = np.isfinite(data[objective])
    if require_fit_ok and "fit_ok" in data:
        good &= data["fit_ok"].astype(bool)
    if "status" in data:
        good &= data["status"] == 0
    indices = np.flatnonzero(good)
    if not len(indices):
        raise ValueError(f"No accepted finite rows were found for {objective!r}")
    return indices[np.argmin(data[objective][indices])]


def phase_branch(magnitude_branch):
    """Return the phase branch paired with an amplitude magnitude branch."""
    return magnitude_branch[0] + "phi" + magnitude_branch[1:]


def input_moment_branch(fit_branch):
    """Map an output moment name to its experimental-input branch name.

    EMI omits the leading ``R`` when it writes a reconstructed response, so
    both ``H04_*`` and ``H_*`` map back to the published ``RH...`` branches.
    Keeping the rule here makes the Hessian and bootstrap comparisons agree.
    """
    return "R" + fit_branch


def wave_label(branch):
    """Convert a_T_1_0 to the compact scientific label P+_T0."""
    reflectivity, polarization, ell, projection = branch.split("_")
    orbital = ("S", "P", "D", "F", "G")[int(ell)]
    sign = "+" if reflectivity == "a" else "-"
    projection = projection.replace("m", "-")
    return (
        orbital + "$^" + sign
        + r"_{\mathrm{" + polarization + "}" + projection + "}$"
    )


def circular_mean(angles):
    """Calculate the mean of periodic angles in radians.

    Averaging unit vectors avoids the false mean near zero obtained by directly
    averaging samples clustered on either side of the ``-pi/pi`` boundary.
    """
    return np.angle(np.mean(np.exp(1j * angles)))


def circular_std(angles):
    """Calculate circular standard deviation from the mean resultant length."""
    # Round-off can put a theoretically unit-length resultant just above one.
    resultant = np.clip(np.abs(np.mean(np.exp(1j * angles))), 0.0, 1.0)
    return np.sqrt(-2.0 * np.log(resultant))


def wrap_phase(angles):
    """Map angles to ``[-pi, pi)`` using the argument of a unit vector."""
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
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, bbox_inches="tight", **kwargs)
    if show:
        plt.show()
    plt.close(figure)
