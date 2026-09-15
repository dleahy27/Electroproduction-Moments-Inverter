#!/usr/bin/env python3

"""Compare Hessian-fitted moments with the experimental input moments."""

import re

import matplotlib.pyplot as plt
import numpy as np
import ROOT

from analysis_utils import PROJECT_DIR, best_index, branches, read_tree, save_figure, style_axis


# Settings: set ROOT_FILES to a list of paths to override automatic discovery.
ROOT_FILES = None
TREE_NAME = "fitResults"
SHOW_PLOTS = False
OUTPUT_DIR = PROJECT_DIR / "AnalysisScripts/outputs/hessian_only_analysis"


def input_branch(fit_branch):
    """Map an output moment name to its experimental-input branch name."""
    return "R" + fit_branch


if ROOT_FILES is None:
    ROOT_FILES = sorted((PROJECT_DIR / "OutputFiles").glob("*_fit_*.root"))

file_pattern = re.compile(r"(.+)_fit_(\d+)\.root$")
for fit_file in ROOT_FILES:
    match = file_pattern.fullmatch(fit_file.name)
    if not match:
        continue
    dataset, bin_text = match.groups()
    bin_index = int(bin_text)
    input_file = PROJECT_DIR / "InputFiles/Experiment" / f"{dataset}_moments.root"
    if not input_file.exists():
        continue

    fit_names = set(branches(fit_file, TREE_NAME))
    input_names = set(branches(input_file, "expMoments"))
    moments = sorted(
        name for name in fit_names
        if name.startswith(("H04_", "H_"))
        and f"err__{name}" in fit_names
        and input_branch(name) in input_names
    )

    # Select the global minimum and retain its propagated moment errors.
    columns = ["chi2"]
    columns += [name for name in ("fit_ok", "status") if name in fit_names]
    columns += [branch for name in moments for branch in (name, f"err__{name}")]
    fit = read_tree(fit_file, TREE_NAME, columns)
    entry = best_index(fit)

    input_data = ROOT.RDataFrame("expMoments", str(input_file)).AsNumpy(
        columns=[input_branch(name) for name in moments]
        + [input_branch(name) + "_err" for name in moments]
    )
    experimental = np.array([
        input_data[input_branch(name)][0][bin_index] for name in moments
    ])
    experimental_error = np.array([
        input_data[input_branch(name) + "_err"][0][bin_index] for name in moments
    ])
    fitted = np.array([fit[name][entry] for name in moments])
    fitted_error = np.array([fit[f"err__{name}"][entry] for name in moments])

    # Compare the experimental and best-fit moments.
    x = np.arange(len(moments))
    figure, axis = plt.subplots(
        figsize=(max(10, 0.58 * len(moments)), 6), constrained_layout=True
    )
    axis.errorbar(
        x - 0.10, experimental, yerr=experimental_error,
        fmt="o", capsize=3, label="Experimental moment",
    )
    axis.errorbar(
        x + 0.10, fitted, yerr=fitted_error,
        fmt="s", capsize=3, label="Best fit (Hessian)",
    )
    axis.axhline(0.0, color="0.45", linestyle=":", linewidth=1)
    axis.set_xticks(x, moments, rotation=45, ha="right")
    axis.set(xlabel="Moment", ylabel="Moment value",
             title=f"{dataset}, bin {bin_index}: experimental and Hessian moments")
    axis.legend(frameon=False)
    style_axis(axis, grid=True)
    output = OUTPUT_DIR / dataset / f"{fit_file.stem}_moments.png"
    save_figure(figure, output, show=SHOW_PLOTS)
    print(
        f"{fit_file.name}: entry {entry}, chi2 = {fit['chi2'][entry]:.7g}; "
        f"wrote {output}"
    )
