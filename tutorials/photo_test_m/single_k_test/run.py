#!/usr/bin/env python3
"""Generate two-k PhotoTest truth and fit an effective k=0 mass model.

Edit ASYMMETRIES below and settings.h for fit and mass-grid defaults.
Run from the repository root with python -m tutorials.photo_test_m.single_k_test.run.
"""

import argparse
import subprocess
from pathlib import Path
from tutorials.paths import PROJECT_DIR

ASYMMETRIES = (0.0, 1.0)
CASE = PROJECT_DIR / "tutorials/photo_test_m/single_k_test"
EMI = PROJECT_DIR / "build/emi"


def main():
    from tutorials.root_io import merge_scan

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts", type=int, help="override starts per asymmetry")
    parser.add_argument("--workers", type=int, help="override worker processes")
    parser.add_argument("--max-calls", type=int, help="override each start's call budget")
    parser.add_argument("--settings", type=Path, default=CASE / "settings.h")
    parser.add_argument("--output", type=Path, default=CASE / "output")
    parser.add_argument("--asymmetries", type=float, nargs="+", default=ASYMMETRIES)
    options = parser.parse_args()
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    common = ["--settings", str(options.settings.resolve())]
    fit_options = []
    for flag in ("starts", "workers", "max-calls"):
        value = getattr(options, flag.replace("-", "_"))
        if value is not None:
            fit_options.extend(["--" + flag, str(value)])
    points = []
    for index, asymmetry in enumerate(options.asymmetries):
        print(f"[{index + 1}/{len(options.asymmetries)}] s={asymmetry:+.2f}", flush=True)
        truth = output / f"s_{index:02d}_truth.root"
        fit = output / f"s_{index:02d}_truth_fit.root"
        subprocess.run(
            [str(EMI), "generate-mass-model", *common, "--photo-test",
             "--k-asymmetry", str(asymmetry), "--polarization", "both",
             "--output", str(truth)],
            cwd=PROJECT_DIR, check=True,
        )
        # Use the model and phase references from settings.h. --photo-test would
        # rebuild the default model and discard those fit controls.
        subprocess.run(
            [str(EMI), "mass-fit", *common, *fit_options, "--polarization", "none",
             "--input", str(truth), "--output", str(fit)],
            cwd=PROJECT_DIR, check=True,
        )
        points.append({"truth_file": str(truth), "fit_file": str(fit),
                       "k_asymmetry": asymmetry})
    # Store scan coordinates alongside the results; analysis needs only scan.root.
    merge_scan(
        points, output / "scan.root",
        {"massMoments": "truth_file", "massFitResults": "fit_file",
         "massFitPredictions": "fit_file", "massFitCurves": "fit_file"},
        mass_truth_key="truth_file",
    )


if __name__ == "__main__":
    main()
