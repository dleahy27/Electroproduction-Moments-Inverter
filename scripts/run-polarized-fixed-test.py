#!/usr/bin/env python3

"""Generate and fit the three polarized fixed-amplitude examples."""

import subprocess
from pathlib import Path


# Settings: edit these values before running the closure test.
POLARIZATIONS = ["initial", "recoil", "both"]
EPSILON = 0.8
STARTS = 10_000
WORKERS = 1
SEED = 12_345
K_GAUGE = None  # For example, "b:T:1:0". It does not apply to "both".
USE_HESSE = False


PROJECT_DIR = Path(__file__).resolve().parents[1]
EMI = PROJECT_DIR / "build" / "emi"
GENERATED_DIR = PROJECT_DIR / "InputFiles" / "Generated"
OUTPUT_DIR = PROJECT_DIR / "OutputFiles"
GENERATED_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Generate and then fit one closure sample for each polarization choice.
for polarization in POLARIZATIONS:
    input_file = GENERATED_DIR / f"fixed_{polarization}_moments.root"
    output_file = OUTPUT_DIR / f"fixed_{polarization}_fit.root"
    gauge = ["--k-gauge", K_GAUGE] if K_GAUGE and polarization != "both" else []

    print(f"Generating {polarization}-polarization closure input")
    subprocess.run(
        [
            str(EMI), "generate-example",
            "--polarization", polarization,
            "--epsilon", str(EPSILON),
            "--output", str(input_file),
            "--quiet",
            *gauge,
        ],
        cwd=PROJECT_DIR,
        check=True,
    )

    print(f"Fitting {polarization}-polarization closure input")
    fit_command = [
        str(EMI), "fit",
        "--electro",
        "--input", str(input_file),
        "--tree", "genMoments",
        "--output", str(output_file),
        "--polarization", polarization,
        "--epsilon", str(EPSILON),
        "--starts", str(STARTS),
        "--workers", str(WORKERS),
        "--seed", str(SEED),
        "--quiet",
        *gauge,
    ]
    if not USE_HESSE:
        fit_command.append("--no-hesse")
    subprocess.run(fit_command, cwd=PROJECT_DIR, check=True)
