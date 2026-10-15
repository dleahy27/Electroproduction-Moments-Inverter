#!/usr/bin/env python3

"""Compare one fixed electroproduction amplitude point with its fitted result.

The generated and fitted trees both use the unpolarized, single-k branch
convention. Truth amplitudes and response moments are read from the generated
ROOT file, so selecting another fixed test requires changing only TEST_NAME.
"""

import re

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from tutorials.analysis_utils import (
    PROJECT_DIR,
    best_index,
    branches,
    complex_amplitude,
    phase_branch,
    read_tree,
    save_figure,
    style_axis,
    wave_label,
)

# Select the fixed-amplitude test produced by run_fixed_amplitudes.py.
TEST_NAME = "large_b"
SHOW_PLOTS = False

CASE_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main"
TRUTH_FILE = CASE_DIR / f"output/generated/fixed_{TEST_NAME}_moments.root"
FIT_FILE = CASE_DIR / f"output/root/fixed_{TEST_NAME}_fit.root"
OUTPUT_DIR = CASE_DIR / "output/plots/fixed_amplitudes"
TRUTH_TREE = "genMoments"
FIT_TREE = "fitResults"


AMPLITUDE_PATTERN = re.compile(r"^[ab]_[TL]_\d+_(?:m?\d+)$")
MOMENT_PATTERN = re.compile(r"^RH_(\d+)_(\d+)_(m?\d+)$")
GROUPS = [
    (("a", "T"), "Natural Transverse"),
    (("b", "T"), "Unnatural Transverse"),
    (("a", "L"), "Natural Longitudinal"),
    (("b", "L"), "Unnatural Longitudinal"),
]


def parse_amplitude(name):
    """Return reflectivity, orientation, l, and m from an amplitude branch."""
    reflectivity, orientation, ell, projection = name.split("_")
    projection = -int(projection[1:]) if projection.startswith("m") else int(projection)
    return reflectivity, orientation, int(ell), projection


def keep_amplitude(name):
    """Apply the physical longitudinal conventions used by the legacy plot."""
    reflectivity, orientation, _, projection = parse_amplitude(name)
    if orientation == "L" and projection < 0:
        return False
    if reflectivity == "b" and orientation == "L" and projection == 0:
        return False
    return True


def moment_label(name):
    """Convert H_alpha_L_M into the compact label used in the notebook."""
    _, alpha, ell, projection = name.split("_")
    projection = projection.replace("m", "-")
    return rf"$H^{alpha}({ell}{projection})$"


def main():
    """Run the configured study when invoked as a Python module."""
    for path in (TRUTH_FILE, FIT_FILE):
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing fixed-amplitude input for TEST_NAME={TEST_NAME!r}: {path}"
            )

    truth_branches = set(branches(TRUTH_FILE, TRUTH_TREE))
    fit_branches = set(branches(FIT_FILE, FIT_TREE))

    # Single-k amplitudes have no trailing _1 or _m1 sector suffix. Require the
    # generated and fitted schemas to agree rather than silently dropping waves.
    truth_magnitudes = {
        name for name in truth_branches if AMPLITUDE_PATTERN.fullmatch(name)
    }
    fit_magnitudes = {
        name for name in fit_branches if AMPLITUDE_PATTERN.fullmatch(name)
    }
    if truth_magnitudes != fit_magnitudes:
        raise RuntimeError(
            "Truth and fit amplitude branches do not match: "
            f"missing from fit={sorted(truth_magnitudes - fit_magnitudes)}, "
            f"missing from truth={sorted(fit_magnitudes - truth_magnitudes)}"
        )
    if not truth_magnitudes:
        raise RuntimeError(
            "No single-k amplitude branches were found; expected names such as "
            "a_T_1_0 without a k-sector suffix"
        )

    all_magnitudes = sorted(truth_magnitudes)
    missing_phases = [
        phase_branch(name)
        for name in all_magnitudes
        if phase_branch(name) not in truth_branches
        or phase_branch(name) not in fit_branches
    ]
    if missing_phases:
        raise RuntimeError(f"Missing paired phase branches: {missing_phases}")

    magnitudes = [name for name in all_magnitudes if keep_amplitude(name)]
    if not magnitudes:
        raise RuntimeError("No physical amplitudes remain after parity filtering")
    phases = [phase_branch(name) for name in magnitudes]

    # Raw generated responses carry a leading R; fitted responses omit it. RH04 is
    # an experimentally mixed observable and is intentionally excluded so that the
    # nine panels show the separate alpha=0,...,8 responses from the notebook.
    truth_moment_names = {
        name for name in truth_branches if MOMENT_PATTERN.fullmatch(name)
    }
    fit_moment_names = {
        name for name in fit_branches if MOMENT_PATTERN.fullmatch("R" + name)
    }
    expected_fit_moments = {name[1:] for name in truth_moment_names}
    if expected_fit_moments != fit_moment_names:
        raise RuntimeError(
            "Truth and fit raw-moment branches do not match: "
            f"missing from fit={sorted(expected_fit_moments - fit_moment_names)}, "
            f"missing from truth={sorted(fit_moment_names - expected_fit_moments)}"
        )
    if not truth_moment_names:
        raise RuntimeError("No raw RH_alpha_L_M truth moments were found")

    moment_pairs = sorted(
        ((truth_name, truth_name[1:]) for truth_name in truth_moment_names),
        key=lambda pair: tuple(
            int(value.replace("m", "-")) for value in pair[1].split("_")[1:]
        ),
    )
    alphas = {
        int(MOMENT_PATTERN.fullmatch(truth_name).group(1))
        for truth_name, _ in moment_pairs
    }
    if alphas != set(range(9)):
        raise RuntimeError(
            "Expected raw moments for alpha=0,...,8; " f"found alpha={sorted(alphas)}"
        )

    truth_columns = magnitudes + phases + [truth_name for truth_name, _ in moment_pairs]
    status = [name for name in ("chi2", "fit_ok", "status") if name in fit_branches]
    if "chi2" not in status:
        raise RuntimeError("The fit tree has no chi2 branch for best-entry selection")
    fit_columns = (
        status + magnitudes + phases + [fit_name for _, fit_name in moment_pairs]
    )

    truth = read_tree(TRUTH_FILE, TRUTH_TREE, truth_columns)
    if len(truth[magnitudes[0]]) != 1:
        raise RuntimeError(
            f"Expected one generated truth entry, found {len(truth[magnitudes[0]])}"
        )
    fit = read_tree(FIT_FILE, FIT_TREE, fit_columns)
    entry = best_index(fit)

    truth_amplitudes = {
        name: complex_amplitude(truth[name][0], truth[phase_branch(name)][0])
        for name in magnitudes
    }
    fit_amplitudes = {
        name: complex_amplitude(fit[name][entry], fit[phase_branch(name)][entry])
        for name in magnitudes
    }

    print(
        f"Fixed test {TEST_NAME!r}: best fit entry {entry}, "
        f"chi2 = {fit['chi2'][entry]:.7g}"
    )
    print(f"Plotting {len(magnitudes)} amplitudes and {len(moment_pairs)} moments")

    # Match the legacy notebook's four stacked Argand panels. Colour identifies a
    # wave; marker shape distinguishes the minimizer from the generated point.
    plt.rcParams.update(
        {
            "font.family": "serif",
            "text.usetex": True,
            "mathtext.fontset": "stix",
            "font.size": 22,
            "axes.titlesize": 18,
            "axes.labelsize": 26,
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
            "legend.fontsize": 18,
            "lines.linewidth": 2.2,
            "axes.linewidth": 1.8,
        }
    )

    groups = [
        (key, title)
        for key, title in GROUPS
        if any(tuple(name.split("_")[:2]) == key for name in magnitudes)
    ]
    figure, axes = plt.subplots(
        len(groups),
        1,
        figsize=(8, 6 * len(groups)),
        sharex=True,
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes[:, 0]

    for group_index, (axis, (group, title)) in enumerate(zip(axes, groups)):
        selected = [name for name in magnitudes if tuple(name.split("_")[:2]) == group]
        for colour_index, name in enumerate(selected):
            colour = f"C{colour_index % 10}"
            fitted = fit_amplitudes[name]
            generated = truth_amplitudes[name]
            axis.plot(
                fitted.real,
                fitted.imag,
                "o",
                color=colour,
                markersize=12,
                label=wave_label(name),
            )
            axis.plot(
                generated.real,
                generated.imag,
                "*",
                color=colour,
                markeredgecolor="black",
                markersize=16,
            )

        axis.set(xlim=(-1, 1), ylim=(-1, 1), ylabel="Im", title=title)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xticks([-1, -0.5, 0, 0.5, 1])
        axis.set_yticks([-1, -0.5, 0, 0.5, 1])
        style_axis(axis, grid=True)
        handles, _ = axis.get_legend_handles_labels()
        if group_index == 0:
            handles.append(
                Line2D(
                    [],
                    [],
                    marker="*",
                    color="white",
                    markeredgecolor="black",
                    linestyle="none",
                    markersize=12,
                    label="Set Amp",
                )
            )
        axis.legend(handles=handles, frameon=False, loc="upper left", reverse=True)
    axes[-1].set_xlabel("Re")

    save_figure(
        figure,
        OUTPUT_DIR / f"fixed_{TEST_NAME}.pdf",
        show=SHOW_PLOTS,
        dpi=300,
        transparent=True,
    )

    # Plot every raw response moment in a separate alpha panel, as in the notebook.
    figure, axes = plt.subplots(3, 3, figsize=(18, 18), constrained_layout=True)
    axes = axes.ravel()

    for alpha, axis in enumerate(axes):
        alpha_pairs = [
            pair
            for pair in moment_pairs
            if int(MOMENT_PATTERN.fullmatch(pair[0]).group(1)) == alpha
        ]
        generated = np.array([truth[truth_name][0] for truth_name, _ in alpha_pairs])
        fitted = np.array([fit[fit_name][entry] for _, fit_name in alpha_pairs])
        labels = [moment_label(fit_name) for _, fit_name in alpha_pairs]
        x = np.arange(len(alpha_pairs))

        axis.scatter(x, fitted, s=90, label="Minimizer")
        axis.scatter(
            x,
            generated,
            s=90,
            marker="x",
            linewidths=2,
            color="tab:orange",
            zorder=4,
            label="Truth",
        )
        axis.axhline(0.0, color="0.45", linestyle=":", linewidth=1)
        axis.set_title(rf"$\alpha = {alpha}$")
        axis.set_xticks(x, labels, rotation=90, fontsize=20)
        axis.minorticks_off()
        axis.margins(x=0.05)
        # axis.grid(axis="y", alpha=0.25)
        if alpha == 0:
            axis.legend(frameon=False)

    figure.supxlabel("Moment", fontsize=28)
    figure.supylabel("Moment value", fontsize=28)

    save_figure(
        figure,
        OUTPUT_DIR / f"fixed_moments_{TEST_NAME}.pdf",
        show=SHOW_PLOTS,
        dpi=400,
    )


if __name__ == "__main__":
    main()
