"""Integration checks for mass-fit input scaling and selection.

Run with a Python providing PyROOT after building EMI:
    python tests/mass_fit_input.py [path/to/emi]
All generated settings and ROOT files live in a temporary directory.
"""

from array import array
import math
from pathlib import Path
import subprocess
import sys
import tempfile

import ROOT

ROOT.gROOT.SetBatch(True)
PROJECT = Path(__file__).resolve().parents[1]


def main():
    executable = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else PROJECT / "build/emi"
    template = (PROJECT / "tutorials/photo_test_m/closure/settings.h").read_text()
    with tempfile.TemporaryDirectory(prefix="emi-mass-input-") as temporary:
        directory = Path(temporary)
        raw = directory / "raw.root"
        template = template.replace(
            "tutorials/photo_test_m/closure/output/moments.root", str(raw)
        ).replace("config.fit.starts = 8;", "config.fit.starts = 1;").replace(
            "config.fit.workers = 2;", "config.fit.workers = 1;"
        )
        count = 0

        def invoke(extra="", source=raw, error=None, command="mass-fit"):
            nonlocal count
            count += 1
            settings = directory / f"settings_{count}.h"
            settings.write_text(template.replace(
                "  config.model = TruthModel();", extra + "\n  config.model = TruthModel();"
            ))
            output = directory / f"result_{count}.root"
            arguments = [str(executable), command, "--settings", str(settings),
                         "--output", str(output)]
            if command == "mass-fit":
                arguments += ["--input", str(source)]
            result = subprocess.run(arguments, cwd=PROJECT, text=True, capture_output=True)
            if error is not None:
                assert result.returncode != 0 and error in result.stdout + result.stderr, result.stdout + result.stderr
                return
            assert result.returncode == 0, result.stdout + result.stderr
            return output

        generated = invoke(command="generate-mass-model")
        # Generation's --output override keeps the fixture outside the checkout.
        generated.rename(raw)
        with ROOT.TFile.Open(str(raw)) as file:
            tree = file.Get("massMoments")
            names = [branch.GetName() for branch in tree.GetListOfBranches()]
            rows = [{name: float(getattr(row, name)) for name in names} for row in tree]
        zeroth = "RH_0_0_0"
        assert zeroth in names

        def normalised_file(name, denominator=2.0, bad_yield=None, poison_first=False):
            path = directory / name
            with ROOT.TFile.Open(str(path), "RECREATE") as file:
                tree = ROOT.TTree("massMoments", "Normalised input")
                buffers = {key: array("d", [0.0]) for key in names + ["yield"]}
                for key, buffer in buffers.items():
                    tree.Branch(key, buffer, key + "/D")
                for index, row in enumerate(rows):
                    scale = row[zeroth] / denominator
                    for key in names:
                        buffers[key][0] = row[key] if key == "mass_GeV" else row[key] / scale
                        if poison_first and index == 0 and key.endswith("_err"):
                            buffers[key][0] = 0.0
                    buffers["yield"][0] = row[zeroth] if bad_yield is None else bad_yield
                    tree.Fill()
                tree.Write()
            return path

        def read_result(path):
            with ROOT.TFile.Open(str(path)) as file:
                tree = file.Get("massFitResults")
                tree.GetEntry(0)
                summary = {name: float(getattr(tree, name)) for name in (
                    "initial_chi2", "chi2", "ndf", "D2_mass", "D2_width", "D2_a_2_1_mag"
                )}
                predictions = file.Get("massFitPredictions")
                keys = [branch.GetName() for branch in predictions.GetListOfBranches()]
                values = [[float(getattr(row, key)) for key in keys] for row in predictions]
            return summary, values

        def equivalent(first, second):
            a, av = read_result(first)
            b, bv = read_result(second)
            for name in a:
                assert math.isclose(a[name], b[name], rel_tol=1e-8, abs_tol=1e-8), (name, a, b)
            assert len(av) == len(bv)
            for x, y in zip(av, bv):
                assert all(math.isclose(i, j, rel_tol=1e-8, abs_tol=1e-8) for i, j in zip(x, y))

        normalised = normalised_file("normalised.root")
        baseline = invoke()
        equivalent(baseline, invoke('  config.yieldBranch = "yield";', normalised))
        unit = normalised_file("unit.root", denominator=1.0)
        equivalent(baseline, invoke('  config.yieldBranch = "yield"; config.normalisation = 1.0;', unit))
        summary, predictions = read_result(baseline)
        assert len(predictions) == len(rows)
        assert abs(summary["D2_mass"] - 1.320) < 1e-5
        assert abs(summary["D2_width"] - 0.120) < 1e-5

        lower, upper = rows[8]["mass_GeV"], rows[30]["mass_GeV"]
        bounds = f"  config.minMassGeV = {lower!r}; config.maxMassGeV = {upper!r};"
        indices = "  config.bins = {" + ",".join(map(str, range(8, 31))) + "};"
        selected = invoke(bounds)
        equivalent(selected, invoke(indices))
        equivalent(selected, invoke(bounds + indices))
        equivalent(selected, invoke(bounds + ' config.yieldBranch = "yield";', normalised))
        selected_summary, selected_predictions = read_result(selected)
        assert len(selected_predictions) == 23
        assert selected_predictions[0][0] == lower and selected_predictions[-1][0] == upper
        assert selected_summary["ndf"] + 3 == (summary["ndf"] + 3) / len(rows) * 23
        equivalent(invoke("  config.bins = {8, 12, 20};"),
                   invoke("  config.bins = {20, 8, 12, 8};"))
        equivalent(invoke(f"  config.minMassGeV = {lower!r};"),
                   invoke("  config.bins = {" + ",".join(map(str, range(8, 41))) + "};"))
        equivalent(invoke(f"  config.maxMassGeV = {upper!r};"),
                   invoke("  config.bins = {" + ",".join(map(str, range(31))) + "};"))
        poisoned = normalised_file("poisoned.root", poison_first=True)
        equivalent(selected, invoke(bounds + ' config.yieldBranch = "yield";', poisoned))

        for selection, message in [
            ("config.bins = {-1};", "outside the input tree"),
            ("config.bins = {41};", "outside the input tree"),
            ("config.minMassGeV = 1.5; config.maxMassGeV = 1.1;", "ordered bounds"),
            ("config.minMassGeV = 4.0;", "no input entries"),
            ("config.bins = {0}; config.minMassGeV = 1.5;", "no input entries"),
            ('config.yieldBranch = "missing";', "Missing scalar branch"),
            ('config.yieldBranch = "";', "positive finite normalisation"),
            ('config.yieldBranch = "yield"; config.normalisation = 0;', "positive finite normalisation"),
        ]:
            invoke(selection, normalised, error=message)
        for value in (0.0, -1.0, float("nan"), float("inf")):
            invalid = normalised_file("invalid.root", bad_yield=value)
            invoke('config.yieldBranch = "yield";', invalid, error="Yield must be positive and finite")
        print("Mass-fit integration checks passed: scaling, selections, discovery, predictions, ndf, and invalid inputs.")


if __name__ == "__main__":
    main()
