#!/usr/bin/env python3

"""Generate and fit the double-polarized fixed-amplitude example.

Generation creates exact moments from a known two-k amplitude set.  Fitting
with the same polarization content then provides direct amplitude-level closure.
"""

import subprocess
from tutorials.paths import PROJECT_DIR

# Settings: edit these values before running the closure test.


def main():
    """Run the configured study when invoked as a Python module."""
    polarization = "both"
    EPSILON = 0.8
    STARTS = 20_000
    WORKERS = 10
    SEED = 0
    K_GAUGE = None  # For example, "b:T:1:0". It does not apply to "both".
    USE_HESSE = False

    EMI = PROJECT_DIR / "build" / "emi"
    CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/fixed_double"
    SETTINGS = CASE_DIR / "settings.h"
    GENERATED_DIR = CASE_DIR / "output/generated"
    OUTPUT_DIR = CASE_DIR / "output/root"
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Generate and then fit one closure sample for each polarization choice.
    input_file = GENERATED_DIR / f"fixed_{polarization}_moments.root"
    output_file = OUTPUT_DIR / f"fixed_{polarization}_fit.root"
    gauge = ["--k-gauge", K_GAUGE] if K_GAUGE and polarization != "both" else []

    print(f"Generating {polarization}-polarization closure input")
    subprocess.run(
        [
            str(EMI),
            "generate-example",
            "--settings",
            str(SETTINGS),
            "--epsilon",
            str(EPSILON),
            "--output",
            str(input_file),
            "--quiet",
        ],
        cwd=PROJECT_DIR,
        check=True,
    )

    print(f"Fitting {polarization}-polarization closure input")
    fit_command = [
        str(EMI),
        "fit",
        "--settings",
        str(SETTINGS),
        "--electro",
        "--input",
        str(input_file),
        "--tree",
        "genMoments",
        "--output",
        str(output_file),
        "--polarization",
        polarization,
        "--epsilon",
        str(EPSILON),
        "--starts",
        str(STARTS),
        "--workers",
        str(WORKERS),
        "--seed",
        str(SEED),
        "--quiet",
        *gauge,
    ]
    if not USE_HESSE:
        fit_command.append("--no-hesse")
    subprocess.run(fit_command, cwd=PROJECT_DIR, check=True)


if __name__ == "__main__":
    main()
