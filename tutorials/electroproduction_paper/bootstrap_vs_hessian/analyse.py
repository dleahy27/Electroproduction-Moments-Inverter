#!/usr/bin/env python3

"""Compare bootstrap widths with Hessian errors for matching data bins.

The Hessian describes local curvature around one minimum, while the bootstrap
refits fluctuated datasets and can expose skewness or alternative solutions.
Phases therefore use circular statistics and every comparison retains the fit
quality filters stored in the ROOT trees.
"""

import re

import matplotlib.pyplot as plt
import numpy as np
import ROOT

from tutorials.analysis_utils import (
    PROJECT_DIR, MAGNITUDE_PATTERN, PHASE_PATTERN, best_index, branches,
    circular_mean, circular_std, r_expression, read_tree, save_figure,
    style_axis, wave_label, wrap_phase,
)
from tutorials.electroproduction_paper.datasets import CHANNELS


# Settings: PAIRS=None discovers matching *_bootstrap_i.root and *_fit_i.root.
PAIRS = None
BOOTSTRAP_TREE = "PartialWaves"
HESSIAN_TREE = "fitResults"
SHOW_PLOTS = False
MAIN_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main"
OUTPUT_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/bootstrap_vs_hessian/output/plots"
Q2_VALUES = {name: np.asarray(channel.q2) for name, channel in CHANNELS.items()}


def comparison_figure(names, bootstrap, bootstrap_error, hessian,
                      hessian_error, ylabel, title, output, phase=False,
                      x_values=None):
    """Draw values, value differences, and uncertainty differences."""
    x = np.arange(len(names)) if x_values is None else np.asarray(x_values)
    difference = bootstrap - hessian
    if phase:
        # A raw subtraction can report nearly 2*pi for neighbouring angles on
        # opposite sides of the conventional -pi/pi branch cut.
        difference = wrap_phase(difference)

    figure, axes = plt.subplots(
        3, 1, figsize=(max(10, 0.55 * len(names)), 9), sharex=True,
        constrained_layout=True, gridspec_kw={"height_ratios": [3, 1.3, 1.3]},
    )
    offset = 0.10 if x_values is None else 0.01 * max(np.ptp(x), 1.0)
    axes[0].errorbar(
        x - offset, bootstrap, yerr=bootstrap_error,
        fmt="o", capsize=3, label="Bootstrap mean and standard deviation",
    )
    axes[0].errorbar(
        x + offset, hessian, yerr=hessian_error,
        fmt="s", capsize=3, label="Best fit and Hessian error",
    )
    axes[0].set(ylabel=ylabel, title=title)
    axes[0].legend(frameon=False)
    axes[1].plot(x, difference, "o", color="black")
    axes[1].set_ylabel(r"$x_{\rm boot}-x_{\rm Hess}$")
    axes[2].plot(x, bootstrap_error - hessian_error, "s", color="black")
    axes[2].set_ylabel(r"$\sigma_{\rm boot}-\sigma_{\rm Hess}$")
    for axis in axes[1:]:
        axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1)
    for axis in axes:
        style_axis(axis, grid=True)
    if x_values is None:
        axes[2].set_xticks(x, names, rotation=50, ha="right")
    else:
        axes[2].set_xlabel(r"$Q^2\ [\mathrm{GeV}^2]$")
    save_figure(figure, output, show=SHOW_PLOTS)


# Discover pairs by replacing the bootstrap part of each filename.
if PAIRS is None:
    PAIRS = []
    for bootstrap_file in sorted((MAIN_DIR / "output/root").glob("*_bootstrap_*.root")):
        fit_file = bootstrap_file.with_name(
            bootstrap_file.name.replace("_bootstrap_", "_fit_")
        )
        if fit_file.exists():
            PAIRS.append((bootstrap_file, fit_file))

file_pattern = re.compile(r"(.+)_bootstrap_(\d+)\.root$")
R_RESULTS = {}

for bootstrap_file, fit_file in PAIRS:
    match = file_pattern.fullmatch(bootstrap_file.name)
    dataset, bin_text = match.groups()
    bin_index = int(bin_text)
    output_dir = OUTPUT_DIR / dataset / f"bin_{bin_index}"

    bootstrap_names = set(branches(bootstrap_file, BOOTSTRAP_TREE))
    fit_names = set(branches(fit_file, HESSIAN_TREE))
    magnitudes = sorted(
        name for name in bootstrap_names & fit_names
        if MAGNITUDE_PATTERN.fullmatch(name) and f"err__{name}" in fit_names
    )
    phases = sorted(
        name for name in bootstrap_names & fit_names
        if PHASE_PATTERN.fullmatch(name) and f"err__{name}" in fit_names
    )
    moments = sorted(
        name for name in bootstrap_names & fit_names
        if name.startswith(("H04_", "H_")) and f"err__{name}" in fit_names
    )

    # Older result files may predate the stored R branch. RDataFrame can
    # evaluate the same sum-of-squared-amplitudes expression without rewriting
    # those files.
    r_column = "R"
    bootstrap_frame = ROOT.RDataFrame(BOOTSTRAP_TREE, str(bootstrap_file))
    if "R" not in bootstrap_names:
        r_column = "derived_R"
        bootstrap_frame = bootstrap_frame.Define(r_column, r_expression(magnitudes))
    bootstrap_columns = magnitudes + phases + moments + [r_column]
    bootstrap_columns += [
        name for name in ("valid", "status") if name in bootstrap_names
    ]
    if "valid" in bootstrap_names:
        bootstrap_frame = bootstrap_frame.Filter("valid != 0")
    if "status" in bootstrap_names:
        bootstrap_frame = bootstrap_frame.Filter("status == 0")
    bootstrap = {
        name: np.asarray(values)
        for name, values in bootstrap_frame.AsNumpy(bootstrap_columns).items()
    }
    bootstrap["R"] = bootstrap.pop(r_column)

    fit_columns = ["chi2", "R", "err__R"]
    fit_columns += [name for name in ("fit_ok", "status") if name in fit_names]
    fit_columns += magnitudes + phases + moments
    fit_columns += [f"err__{name}" for name in magnitudes + phases + moments]
    fit = read_tree(fit_file, HESSIAN_TREE, fit_columns)
    entry = best_index(fit)
    print(
        f"{dataset} bin {bin_index}: {len(bootstrap['R'])} bootstrap samples; "
        f"Hessian entry {entry}, chi2={fit['chi2'][entry]:.7g}"
    )

    for names, ylabel, filename, phase in (
        (magnitudes, r"Amplitude magnitude $|A|$", "amplitudes.png", False),
        (phases, "Amplitude phase (rad)", "phases.png", True),
        (moments, "Moment value", "moments.png", False),
    ):
        # Phase samples live on a circle: values close to -pi and +pi are
        # neighbours and must not be averaged as ordinary real numbers.
        statistic = circular_mean if phase else np.mean
        uncertainty = circular_std if phase else (
            lambda values: np.std(values, ddof=1)
        )
        bootstrap_value = np.array([statistic(bootstrap[name]) for name in names])
        bootstrap_error = np.array([
            uncertainty(bootstrap[name]) for name in names
        ])
        hessian_value = np.array([fit[name][entry] for name in names])
        hessian_error = np.array([fit[f"err__{name}"][entry] for name in names])
        labels = [
            wave_label(name) if name in magnitudes else name for name in names
        ]
        comparison_figure(
            labels, bootstrap_value, bootstrap_error,
            hessian_value, hessian_error, ylabel,
            f"{dataset}, bin {bin_index}: bootstrap and Hessian",
            output_dir / filename, phase=phase,
        )

    R_RESULTS.setdefault(dataset, []).append((
        bin_index,
        np.mean(bootstrap["R"]), np.std(bootstrap["R"], ddof=1),
        fit["R"][entry], fit["err__R"][entry],
    ))

# R is clearest as one Q2-dependent figure per reaction.
for dataset, rows in R_RESULTS.items():
    rows.sort()
    bins = np.array([row[0] for row in rows])
    q2 = Q2_VALUES[dataset][bins] if dataset in Q2_VALUES else bins
    comparison_figure(
        [str(value) for value in q2],
        np.array([row[1] for row in rows]),
        np.array([row[2] for row in rows]),
        np.array([row[3] for row in rows]),
        np.array([row[4] for row in rows]),
        r"$R=\sigma_L/\sigma_T$",
        f"{dataset}: bootstrap and Hessian comparison of R",
        OUTPUT_DIR / dataset / "R.png",
        x_values=q2,
    )
