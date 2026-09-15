"""Small shared helpers for the standalone analysis scripts."""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ROOT


PROJECT_DIR = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = Path(__file__).resolve().parent

MAGNITUDE_PATTERN = re.compile(r"^[ab]_[TL]_\d+_(?:m?\d+)$")
PHASE_PATTERN = re.compile(r"^[ab]phi_[TL]_\d+_(?:m?\d+)$")


def branches(path, tree_name):
    """Return the branch names stored in one ROOT tree."""
    frame = ROOT.RDataFrame(tree_name, str(path))
    return [str(name) for name in frame.GetColumnNames()]


def read_tree(path, tree_name, columns, definitions=None):
    """Read selected scalar ROOT branches into a dictionary of NumPy arrays."""
    frame = ROOT.RDataFrame(tree_name, str(path))
    for name, expression in definitions or []:
        frame = frame.Define(name, expression)
    return {
        name: np.asarray(values)
        for name, values in frame.AsNumpy(columns=list(columns)).items()
    }


def best_index(data, objective="chi2", require_fit_ok=True):
    """Return the row with the smallest finite objective value."""
    good = np.isfinite(data[objective])
    if require_fit_ok and "fit_ok" in data:
        good &= data["fit_ok"].astype(bool)
    if "status" in data:
        good &= data["status"] == 0
    indices = np.flatnonzero(good)
    return indices[np.argmin(data[objective][indices])]


def phase_branch(magnitude_branch):
    """Return the phase branch paired with an amplitude magnitude branch."""
    return magnitude_branch[0] + "phi" + magnitude_branch[1:]


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
    """Calculate the mean of periodic angles in radians."""
    return np.angle(np.mean(np.exp(1j * angles)))


def circular_std(angles):
    """Calculate the circular standard deviation in radians."""
    resultant = np.abs(np.mean(np.exp(1j * angles)))
    return np.sqrt(-2.0 * np.log(resultant))


def wrap_phase(angles):
    """Map angles to the interval [-pi, pi)."""
    return np.arctan2(np.sin(angles), np.cos(angles))


def complex_amplitude(magnitude, phase):
    """Construct a complex amplitude, allowing signed fitted magnitudes."""
    magnitude = np.asarray(magnitude)
    phase = np.asarray(phase) + np.where(magnitude < 0.0, np.pi, 0.0)
    return np.abs(magnitude) * np.exp(1j * phase)


def r_expression(magnitudes):
    """Build the RDataFrame expression for sigma_L / sigma_T."""
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
    """Save a figure and close it unless interactive display is requested."""
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, bbox_inches="tight", **kwargs)
    if show:
        plt.show()
    plt.close(figure)
