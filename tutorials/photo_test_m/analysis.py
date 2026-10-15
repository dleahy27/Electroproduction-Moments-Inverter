"""Load mass-fit data once, keeping ROOT I/O separate from Matplotlib plots."""

import numpy as np
from tutorials.analysis_utils import (
    branches,
    best_index,
    best_indices,
    read_tree,
    read_best_rows,
    truth_parameters,
)


def pole_names(names):
    # names of form:
    # bw_meson_massvalue_mass/width
    plot_names = []
    for i,name in enumerate(names):
        strings = name.split("_")
        mass = strings[2]
        if strings[1] == "a0":
            meson = r"a_{0}("+mass+")"
        elif strings[1] == "pi1":
            meson = r"\pi_{1}("+mass+")"
        elif strings[1] == "a2":
            meson = r"a_{2}("+mass+")"
        else:
            raise ValueError("Meson not recognised for what was in the photo model.")

        if strings[-1]=="width":
            plot_names.append(r"$\Gamma_{"+meson+"}$")
        elif strings[-1]=="mass":
            plot_names.append(r"$m_{"+meson+"}$")

        else:
            raise ValueError("The string does not end with mass or width as expected.")
    return plot_names

def coupling_name(names):
    # names of form:
    # bw_meson_massvalue_a/b_l_m_k_phase/mag
    plot_names = []

    for i,name in enumerate(names):
        strings = name.split("_")
        mass = strings[2]
        refl_val = strings[3]
        m = strings[5]

        if refl_val == "a":
            refl = r"+"
        elif refl_val == "b":
            refl = r"-"
        else:
            raise ValueError("Reflectivity not recognised.")

        if strings[1] == "a0":
            meson = r"a_{0"
        elif strings[1] == "pi1":
            meson = r"\pi_{1"
        elif strings[1] == "a2":
            meson = r"a_{2"
        else:
            raise ValueError("Meson not recognised for what was in the photo model.")

        plot_names.append(r"$"+meson+m+"}^{"+refl+"}"+"("+mass+")$")

    return plot_names
def mass_coordinates(names):
    """Select physical BW coordinates, excluding status and uncertainty columns."""
    return [name for name in names if name.startswith("bw_")]


def load_closure(case):
    """Read the closure's fit table, truth parameters, and both mass spectra."""
    truth_file = case / "output/moments.root"
    fit_file = case / "output/fit.root"
    coordinates = mass_coordinates(branches(fit_file, "massFitResults"))
    results = read_tree(
        fit_file, "massFitResults", ["chi2", "ndf", "fit_ok", "start", *coordinates]
    )
    predictions = set(branches(fit_file, "massFitPredictions"))
    spectra = [
        name
        for name in branches(truth_file, "massMoments")
        if name.startswith("RH")
        and not name.endswith("_err")
        and "model__" + name in predictions
    ]
    return {
        "coordinates": coordinates,
        "results": results,
        "best": best_index(results),
        "top": best_indices(results),
        "truth": truth_parameters(truth_file),
        "spectra": spectra,
        "moments": read_tree(
            truth_file,
            "massMoments",
            ["mass_GeV", *spectra, *[name + "_err" for name in spectra]],
        ),
        "predictions": read_tree(
            fit_file,
            "massFitPredictions",
            ["mass_GeV", *["model__" + name for name in spectra]],
        ),
    }


def load_mass_scan(output):
    """Read merged scan trees once and return one small record per case.

    Start quality is retained for the landscape plots; physical fit coordinates
    are loaded for the best 0.1% of accepted fits per case by RDataFrame.
    """
    path = output / "scan.root"
    if not path.is_file():
        raise FileNotFoundError(
            "Run python -m tutorials.photo_test_m.single_k_test.run first"
        )
    coordinates = mass_coordinates(branches(path, "massFitResults"))
    fits, quality = read_best_rows(path, "massFitResults", coordinates, fraction=0.05)
    truth_names = branches(path, "massModelTruth")
    truth = read_tree(path, "massModelTruth", truth_names)
    moment_names = branches(path, "massMoments")
    intensity = "RH_0_0_0" if "RH_0_0_0" in moment_names else "RH_0_0_0_0_0"
    moments = read_tree(path, "massMoments", ["scan_case", "mass_GeV", intensity])
    predictions = read_tree(
        path, "massFitPredictions", ["scan_case", "mass_GeV", "model__" + intensity]
    )
    cases = []
    for index in np.unique(truth["scan_case"]):

        def rows(table):
            mask = table["scan_case"] == index
            return {name: values[mask] for name, values in table.items()}

        injected = rows(truth)
        cases.append(
            {
                "tag": f"s_{int(index):02d}",
                "asymmetry": float(injected["scan_k_asymmetry"][0]),
                "coordinates": coordinates,
                "fit": rows(fits),
                "quality": rows(quality),
                "intensity": intensity,
                "truth": {name: values[0] for name, values in injected.items()},
                "moments": rows(moments),
                "predictions": rows(predictions),
            }
        )
    return cases
