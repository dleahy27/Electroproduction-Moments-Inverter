#!/usr/bin/env python3
"""Generate two-k PhotoTest truth and fit an effective k=0 mass model.

Edit ASYMMETRIES below and settings.h for fit and mass-grid defaults.
Run from the repository root with python -m tutorials.photo_test_m.single_k_test.run.
"""

import subprocess
from tutorials.paths import PROJECT_DIR

ASYMMETRIES = (0.0, 1.0)
CASE = PROJECT_DIR / "tutorials/photo_test_m/single_k_test"
EMI = PROJECT_DIR / "build/emi"


def main():
    from tutorials.root_io import merge_scan

    output = CASE / "output"
    output.mkdir(parents=True, exist_ok=True)
    common = ["--settings", str(CASE / "settings.h"), "--photo-test"]
    points = []
    for index, asymmetry in enumerate(ASYMMETRIES):
        print(f"[{index + 1}/{len(ASYMMETRIES)}] s={asymmetry:+.2f}", flush=True)
        truth = output / f"s_{index:02d}_truth.root"
        fit = output / f"s_{index:02d}_truth_fit.root"
        subprocess.run(
            [str(EMI), "generate-mass-model", *common,
             "--k-asymmetry", str(asymmetry), "--polarization", "both",
             "--output", str(truth)],
            cwd=PROJECT_DIR, check=True,
        )
        subprocess.run(
            [str(EMI), "mass-fit", *common, "--polarization", "none",
             "--input", str(truth), "--output", str(fit)],
            cwd=PROJECT_DIR, check=True,
        )
        points.append({"truth_file": str(truth), "fit_file": str(fit),
                       "k_asymmetry": asymmetry})
    # Store scan coordinates alongside the results; analysis needs only scan.root.
    merge_scan(
        points, output / "scan.root",
        {"massMoments": "truth_file", "massFitResults": "fit_file",
         "massFitPredictions": "fit_file"},
        mass_truth_key="truth_file",
    )


if __name__ == "__main__":
    main()
