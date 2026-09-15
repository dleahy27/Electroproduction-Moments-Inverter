#!/usr/bin/env python3

"""Bootstrap one experimental moment bin."""

import subprocess
from pathlib import Path


# Settings: edit these values to choose the data and bootstrap run.
INPUT_FILE = "InputFiles/Experiment/e_rho_moments.root"
OUTPUT_FILE = "OutputFiles/e_rho_bootstrap_0.root"
TREE_NAME = "expMoments"
BIN = 0
EPSILON = 0.8
ELECTROPRODUCTION = True
TOYS = None  # None uses the value in app/UserSettings.h.
STARTS_PER_TOY = None  # None uses the value in app/UserSettings.h.
WORKERS = None  # None uses the value in app/UserSettings.h.


# Build the EMI command from the settings above.
PROJECT_DIR = Path(__file__).resolve().parents[1]
command = [
    str(PROJECT_DIR / "build" / "emi"),
    "bootstrap",
    "--input", str(PROJECT_DIR / INPUT_FILE),
    "--tree", TREE_NAME,
    "--bin", str(BIN),
    "--epsilon", str(EPSILON),
    "--output", str(PROJECT_DIR / OUTPUT_FILE),
]
if ELECTROPRODUCTION:
    command.append("--electro")
if TOYS is not None:
    command.extend(["--toys", str(TOYS)])
if STARTS_PER_TOY is not None:
    command.extend(["--starts-per-toy", str(STARTS_PER_TOY)])
if WORKERS is not None:
    command.extend(["--workers", str(WORKERS)])

# Run the bootstrap and stop if EMI reports an error.
subprocess.run(command, cwd=PROJECT_DIR, check=True)
