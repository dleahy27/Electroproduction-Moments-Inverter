#!/usr/bin/env python3

"""Generate PhotoTest mass-scan points and fit them with the single-k model."""

import csv
import subprocess
from pathlib import Path


# Settings: edit these values to define the scan and fit cost.
MASS_BINS = 20
MASS_MIN = 0.5
MASS_MAX = 4.0
SCALE_COUNT = 5
SCALE_MIN = 0.0
SCALE_MAX = 0.3
STARTS = 10_000
WORKERS = 8
SEED = 12_345
SETTINGS_FILE = "app/UserSettings.h"
SCAN_NAME = "photo_test_mass_scan"
INCLUDE_BACKGROUND = True
OVERWRITE = False
DRY_RUN = False


# Set up repository paths and evenly spaced scan values.
PROJECT_DIR = Path(__file__).resolve().parents[1]
EMI = PROJECT_DIR / "build" / "emi"
SETTINGS = PROJECT_DIR / SETTINGS_FILE
truth_dir = PROJECT_DIR / "InputFiles" / "Generated" / SCAN_NAME
fit_dir = PROJECT_DIR / "OutputFiles" / SCAN_NAME
truth_dir.mkdir(parents=True, exist_ok=True)
fit_dir.mkdir(parents=True, exist_ok=True)

mass_ratio = MASS_MAX / MASS_MIN
masses = [
    MASS_MIN * mass_ratio ** (index / (MASS_BINS - 1))
    for index in range(MASS_BINS)
]
scales = [
    SCALE_MIN + (SCALE_MAX - SCALE_MIN) * index / (SCALE_COUNT - 1)
    for index in range(SCALE_COUNT)
]

# Generate and fit every mass/scale pair, reusing completed files by default.
manifest_rows = []
total = len(masses) * len(scales)
point_number = 0
for scale_index, scale in enumerate(scales):
    for mass_index, mass in enumerate(masses):
        point_number += 1
        stem = f"scale_{scale_index:02d}_mass_{mass_index:02d}"
        truth_file = truth_dir / f"{stem}_truth.root"
        fit_file = fit_dir / f"{stem}_fit.root"
        run_seed = SEED + scale_index * len(masses) + mass_index

        manifest_rows.append({
            "scale_index": scale_index,
            "mass_index": mass_index,
            "mass_GeV": f"{mass:.17g}",
            "k_minus_scale": f"{scale:.17g}",
            "seed": run_seed,
            "truth_file": str(truth_file.relative_to(PROJECT_DIR)),
            "fit_file": str(fit_file.relative_to(PROJECT_DIR)),
        })
        print(
            f"[{point_number:3d}/{total}] mass={mass:.7g} GeV, "
            f"k-minus scale={scale:.5g}",
            flush=True,
        )

        generate_truth = OVERWRITE or not truth_file.exists()
        if generate_truth:
            command = [
                str(EMI), "generate-fixed",
                "--settings", str(SETTINGS),
                "--output", str(truth_file),
                "--mass", f"{mass:.17g}",
                "--k-minus-scale", f"{scale:.17g}",
                "--quiet",
            ]
            if not INCLUDE_BACKGROUND:
                command.append("--no-background")
            print("  " + " ".join(command), flush=True)
            if not DRY_RUN:
                subprocess.run(command, cwd=PROJECT_DIR, check=True)
        else:
            print(f"  using existing {truth_file.name}", flush=True)

        run_fit = OVERWRITE or generate_truth or not fit_file.exists()
        if run_fit:
            command = [
                str(EMI), "fit",
                "--settings", str(SETTINGS),
                "--input", str(truth_file),
                "--tree", "genMoments",
                "--output", str(fit_file),
                "--photo",
                "--polarization", "none",
                "--starts", str(STARTS),
                "--workers", str(WORKERS),
                "--seed", str(run_seed),
                "--no-hesse",
                "--quiet",
            ]
            print("  " + " ".join(command), flush=True)
            if not DRY_RUN:
                subprocess.run(command, cwd=PROJECT_DIR, check=True)
        else:
            print(f"  using existing {fit_file.name}", flush=True)

# Record the files and settings used by the analysis script.
manifest = fit_dir / "scan_points.csv"
with manifest.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=manifest_rows[0].keys())
    writer.writeheader()
    writer.writerows(manifest_rows)

print(f"\nWrote scan manifest to {manifest}")
print(
    "Analyse it with: conda run -n phdconda python "
    "AnalysisScripts/analyse-photo-test-mass-scan.py"
)
