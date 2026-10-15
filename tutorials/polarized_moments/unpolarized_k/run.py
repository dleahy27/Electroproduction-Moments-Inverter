#!/usr/bin/env python3

"""Generate polarized fixed points and fit their unpolarized moments.

Each repeat uses a deterministic seed derived from the base value.  Optional
suppression factors vary one truth component so the reduced model's response
can be measured by the companion analysis.
"""

import subprocess
from tutorials.paths import PROJECT_DIR

# Settings: edit these values before running the study.
SETTINGS_FILE = "tutorials/polarized_moments/unpolarized_k/settings.h"
EPSILON = None  # None uses the value in the settings file.
STARTS = 1_000
WORKERS = 8
SEED = 12_345
REPEATS = 200
SUPPRESSIONS = [None]  # For example: [0.2, 0.1, 0.05].
USE_HESSE = False


EMI = PROJECT_DIR / "build" / "emi"
SETTINGS = PROJECT_DIR / SETTINGS_FILE
CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/unpolarized_k"


def merge_completed_batch(input_dir, output_dir):
    """Record repetitions and merge their ROOT trees for the plotting pass."""
    import csv
    import re
    from tutorials.root_io import merge_scan

    pattern = re.compile(r"fixed_N(\d+)(?:_M(\d+))?_moments\.root")
    points = []
    for truth in sorted(input_dir.glob("fixed_N*_moments.root")):
        match = pattern.fullmatch(truth.name)
        fit = output_dir / truth.name.replace("_moments.root", "_unpolarized_fit.root")
        if match and fit.is_file():
            points.append(
                {
                    "truth_file": str(truth.relative_to(PROJECT_DIR)),
                    "fit_file": str(fit.relative_to(PROJECT_DIR)),
                    "repeat": int(match[1]),
                    "suppression": int(match[2] or 0),
                }
            )
    if not points:
        raise RuntimeError("No completed repetition pairs to merge")
    with (output_dir / "scan_points.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=points[0].keys())
        writer.writeheader()
        writer.writerows(points)
    merge_scan(
        points,
        output_dir / "scan.root",
        {"genMoments": "truth_file", "fitResults": "fit_file"},
    )


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--merge-only", action="store_true")
    parser.add_argument("--starts", type=int, default=STARTS)
    parser.add_argument("--workers", type=int, default=WORKERS)
    parser.add_argument("--repeats", type=int, default=REPEATS)
    options = parser.parse_args()
    # Name the batch after its repeat and suppression counts.
    batch = f"N{options.repeats}"
    if SUPPRESSIONS != [None]:
        batch += f"_M{len(SUPPRESSIONS)}"
    input_dir = CASE_DIR / "output/generated" / batch
    output_dir = CASE_DIR / "output/root" / batch
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    if options.merge_only:
        merge_completed_batch(input_dir, output_dir)
        return

    # Generate each configured fixed point, then fit its unpolarized projection.
    for suppression_index, suppression in enumerate(SUPPRESSIONS):
        for repeat in range(options.repeats):
            run_seed = SEED + suppression_index * options.repeats + repeat
            stem = f"fixed_N{repeat}"
            if suppression is not None:
                stem += f"_M{suppression_index}"
            input_file = input_dir / f"{stem}_moments.root"
            output_file = output_dir / f"{stem}_unpolarized_fit.root"

            generate_command = [
                str(EMI),
                "generate-fixed",
                "--settings",
                str(SETTINGS),
                "--output",
                str(input_file),
                "--seed",
                str(run_seed),
                "--quiet",
            ]
            if EPSILON is not None:
                generate_command.extend(["--epsilon", str(EPSILON)])
            if suppression is not None:
                generate_command.extend(["--suppression", str(suppression)])

            details = f"N={repeat}, seed={run_seed}"
            if suppression is not None:
                details += f", suppression={suppression}"
            print(f"Generating configured fixed point: {details}")
            subprocess.run(generate_command, cwd=PROJECT_DIR, check=True)

            fit_command = [
                str(EMI),
                "fit",
                "--settings",
                str(SETTINGS),
                "--input",
                str(input_file),
                "--tree",
                "genMoments",
                "--output",
                str(output_file),
                "--polarization",
                "none",
                "--starts",
                str(options.starts),
                "--workers",
                str(options.workers),
                "--seed",
                str(run_seed),
                "--quiet",
            ]
            if EPSILON is not None:
                fit_command.extend(["--epsilon", str(EPSILON)])
            if not USE_HESSE:
                fit_command.append("--no-hesse")

            print(f"Fitting its unpolarized moments: N={repeat}")
            subprocess.run(fit_command, cwd=PROJECT_DIR, check=True)

    print(f"Wrote generated files under {input_dir}")
    print(f"Wrote fitted files under {output_dir}")

    merge_completed_batch(input_dir, output_dir)


if __name__ == "__main__":
    main()
