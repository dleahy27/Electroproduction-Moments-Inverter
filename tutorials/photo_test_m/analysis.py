"""Load mass-fit data once, keeping ROOT I/O separate from Matplotlib plots."""

import re
import numpy as np
from tutorials.analysis_utils import (
    branches,
    best_index,
    best_indices,
    read_tree,
    read_best_rows,
    truth_parameters,
    is_mass_parameter,
)


def resonance_label(name):
    """Matplotlib label for a named resonance, including custom channels."""
    match = re.fullmatch(r"(pi|[af])(\d+)_(\d+)", name)
    if match:
        family, spin, mass = match.groups()
        family = r"\pi" if family == "pi" else family
        return family + "_{" + spin + "}(" + mass + ")"
    return r"\mathrm{" + name.replace("_", r"\_") + "}"


def pole_names(names):
    labels = []
    for name in names:
        resonance, parameter = mass_coordinate_name(name).rsplit("_", 1)
        if parameter not in ("mass", "width"):
            raise ValueError(f"Not a pole mass or width: {name}")
        symbol = "m" if parameter == "mass" else r"\Gamma"
        labels.append("$" + symbol + "_{" + resonance_label(resonance) + "}$")
    return labels


def coupling_name(names):
    labels = []
    for name in names:
        match = re.fullmatch(
            r"(.+)_([ab])(?:_L)?_\d+_(m?\d+)(?:_km?\d+)?_(?:mag|phase)",
            mass_coordinate_name(name),
        )
        if not match:
            raise ValueError(f"Not a mass coupling: {name}")
        resonance, reflectivity, projection = match.groups()
        sign = "+" if reflectivity == "a" else "-"
        projection = projection.replace("m", "-")
        labels.append("$" + resonance_label(resonance) + "^{" + sign + "}_{" + projection + "}$")
    return labels


def coupling_sector_name(name, sector):
    """Add an explicit k sector to a single-k coupling name."""
    base, parameter = mass_coordinate_name(name).rsplit("_", 1)
    index = "m" + str(-sector) if sector < 0 else str(sector)
    return base + ("_k" + index if sector else "") + "_" + parameter


def mass_coordinate_name(name):
    """Normalize old bw/k0/minus spellings to concise mass coordinates."""
    name = name.removeprefix("bw_")
    name = re.sub(r"(^|_)(k?)-(\d+)(?=_|$)", r"\1\2m\3", name)
    return name.replace("_k0_", "_")


def mass_coordinate_table(table):
    """Normalize coordinate keys after reading the original ROOT branches."""
    normalized = {}
    for name, values in table.items():
        canonical = mass_coordinate_name(name)
        if canonical in normalized:
            raise ValueError(f"Duplicate mass coordinate alias: {canonical}")
        normalized[canonical] = values
    return normalized


def mass_coordinates(names):
    """Select resonance parameters, excluding background and error columns."""
    return [name for name in names if is_mass_parameter(name) and not name.startswith("bg_")]


def load_closure(case):
    """Read the closure's fit table, truth parameters, and both mass spectra."""
    truth_file = case / "output/moments.root"
    fit_file = case / "output/fit.root"
    coordinates = mass_coordinates(branches(fit_file, "massFitResults"))
    results = read_tree(
        fit_file, "massFitResults", ["chi2", "ndf", "fit_ok", "start", *coordinates]
    )
    results = mass_coordinate_table(results)
    coordinates = [mass_coordinate_name(name) for name in coordinates]
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
        "truth": mass_coordinate_table(truth_parameters(truth_file)),
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
    fits = mass_coordinate_table(fits)
    truth = mass_coordinate_table(truth)
    coordinates = [mass_coordinate_name(name) for name in coordinates]
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
