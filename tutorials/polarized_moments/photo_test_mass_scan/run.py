"""Generate PhotoTest mass/asymmetry points and merge them for analysis.

settings.h owns the physics model. Constants below keep the original scan
and fit defaults; CLI options allow short runs or merging existing results.
"""

import argparse
import csv
import subprocess
import numpy as np
from tutorials.paths import PROJECT_DIR

# Settings: edit these values to define the scan and fit cost.
MASS_BINS = 15
MASS_MIN = 0.8
MASS_MAX = 1.8
ASYMMETRIES = (-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)
STARTS = 10_000
WORKERS = 10
SEED = 1441
SETTINGS_FILE = "tutorials/polarized_moments/photo_test_mass_scan/settings.h"
INCLUDE_BACKGROUND = False
# Set a relative moment-noise width, e.g. 0.05, or None for exact truth.
NOISE = None
# False resumes a partial scan. Set True only when the model or settings have
# changed and every existing truth/fit pair must be replaced.
OVERWRITE = True
DRY_RUN = False


CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/photo_test_mass_scan"
EMI = PROJECT_DIR / "build/emi"
SETTINGS = PROJECT_DIR / SETTINGS_FILE


def merge_completed_scan(points):
    """Merge available point pairs without shifting their manifest identifiers."""
    from tutorials.root_io import merge_scan

    completed = [
        {**point, "scan_case": case}
        for case, point in enumerate(points)
        if (PROJECT_DIR / point["truth_file"]).is_file()
        and (PROJECT_DIR / point["fit_file"]).is_file()
    ]
    if not completed:
        raise RuntimeError("The manifest contains no completed truth/fit pairs")
    merge_scan(
        completed,
        CASE_DIR / "output/root/scan.root",
        {"genMoments": "truth_file", "fitResults": "fit_file"},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--merge-only",
        action="store_true",
        help="merge existing manifest cases without fitting",
    )
    parser.add_argument("--starts", type=int, default=STARTS)
    parser.add_argument("--workers", type=int, default=WORKERS)
    parser.add_argument("--mass-bins", type=int, default=MASS_BINS)
    parser.add_argument("--asymmetries", type=float, nargs="+", default=ASYMMETRIES)
    parser.add_argument(
        "--resume", action="store_true", help="reuse completed truth/fit pairs"
    )
    parser.add_argument("--dry-run", action="store_true", default=DRY_RUN)
    options = parser.parse_args()
    truth_dir = CASE_DIR / "output/generated"
    fit_dir = CASE_DIR / "output/root"
    manifest = fit_dir / "scan_points.csv"
    if options.merge_only:
        with manifest.open(newline="") as stream:
            merge_completed_scan(list(csv.DictReader(stream)))
        return
    truth_dir.mkdir(parents=True, exist_ok=True)
    fit_dir.mkdir(parents=True, exist_ok=True)
    masses = np.linspace(MASS_MIN, MASS_MAX, options.mass_bins)
    overwrite = OVERWRITE and not options.resume
    # Generate and fit every mass/asymmetry pair, reusing completed files when resuming.
    manifest_rows = []
    total = len(masses) * len(options.asymmetries)
    point_number = 0
    for asymmetry_index, asymmetry in enumerate(options.asymmetries):
        for mass_index, mass in enumerate(masses):
            point_number += 1
            stem = f"asymmetry_{asymmetry_index:02d}_mass_{mass_index:02d}"
            truth_file = truth_dir / f"{stem}_truth.root"
            fit_file = fit_dir / f"{stem}_fit.root"
            # Reusing one seed holds the random-start sequence fixed between scan
            # points. Changes across mass or asymmetry then reflect the physics input
            # rather than a different draw of initial minimizer coordinates.
            run_seed = SEED

            manifest_rows.append(
                {
                    "asymmetry_index": asymmetry_index,
                    "mass_index": mass_index,
                    "mass_GeV": f"{mass:.17g}",
                    "k_asymmetry": f"{asymmetry:.17g}",
                    "fixed_moment_noise": ("" if NOISE is None else f"{NOISE:.17g}"),
                    "seed": run_seed,
                    "truth_file": str(truth_file.relative_to(PROJECT_DIR)),
                    "fit_file": str(fit_file.relative_to(PROJECT_DIR)),
                }
            )
            print(
                f"[{point_number:3d}/{total}] mass={mass:.7g} GeV, "
                f"k asymmetry s={asymmetry:.5g}",
                flush=True,
            )

            # Truth must be regenerated before its dependent fit.  The separate
            # booleans below also allow an interrupted scan to resume at its first
            # missing artifact when overwrite is false.
            generate_truth = overwrite or not truth_file.exists()
            if generate_truth:
                command = [
                    str(EMI),
                    "generate-fixed",
                    "--settings",
                    str(SETTINGS),
                    "--output",
                    str(truth_file),
                    "--mass",
                    f"{mass:.17g}",
                    "--k-asymmetry",
                    f"{asymmetry:.17g}",
                    "--quiet",
                ]
                if not INCLUDE_BACKGROUND:
                    command.append("--no-background")
                if NOISE is not None:
                    command.extend(["--noise", f"{NOISE:.17g}"])
                print("  " + " ".join(command), flush=True)
                if not options.dry_run:
                    subprocess.run(command, cwd=PROJECT_DIR, check=True)
            else:
                print(f"  using existing {truth_file.name}", flush=True)

            run_fit = overwrite or generate_truth or not fit_file.exists()
            if run_fit:
                command = [
                    str(EMI),
                    "fit",
                    "--settings",
                    str(SETTINGS),
                    "--input",
                    str(truth_file),
                    "--tree",
                    "genMoments",
                    "--output",
                    str(fit_file),
                    "--photo",
                    "--polarization",
                    "none",
                    "--starts",
                    str(options.starts),
                    "--workers",
                    str(options.workers),
                    "--seed",
                    str(run_seed),
                    "--no-hesse",
                    "--quiet",
                ]
                print("  " + " ".join(command), flush=True)
                if not options.dry_run:
                    subprocess.run(command, cwd=PROJECT_DIR, check=True)
            else:
                print(f"  using existing {fit_file.name}", flush=True)

    # Write the manifest only after commands have been considered for every point.
    # In dry-run mode it documents the proposed scan; otherwise it records the
    # generated file pairing consumed by analyse.py.
    manifest = fit_dir / "scan_points.csv"
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=manifest_rows[0].keys())
        writer.writeheader()
        writer.writerows(manifest_rows)

    if not options.dry_run:
        merge_completed_scan(manifest_rows)
    print(f"Wrote scan manifest to {manifest}")
    print(
        "Analyse it with: python -m tutorials.polarized_moments.photo_test_mass_scan.analyse"
    )


if __name__ == "__main__":
    main()
