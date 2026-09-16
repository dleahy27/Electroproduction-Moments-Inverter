#!/usr/bin/env python3

"""Fit one experimental moment bin."""

import subprocess
from tutorials.electroproduction_paper.datasets import CHANNELS
from tutorials.paths import PROJECT_DIR


# Settings: edit these values to choose the data and fit.
CHANNEL = CHANNELS["e_rho"]
INPUT_FILE = CHANNEL.input_file
OUTPUT_FILE = "tutorials/electroproduction_paper/main/output/root/e_rho_fit_0.root"
TREE_NAME = "expMoments"
BIN = 0
EPSILON = CHANNEL.epsilon
ELECTROPRODUCTION = True
STARTS = None  # None uses the value in app/UserSettings.h.
WORKERS = None  # None uses the value in app/UserSettings.h.


# Build the EMI command from the settings above.
output_path = PROJECT_DIR / OUTPUT_FILE
output_path.parent.mkdir(parents=True, exist_ok=True)
command = [
    str(PROJECT_DIR / "build" / "emi"),
    "fit",
    "--input", str(INPUT_FILE),
    "--tree", TREE_NAME,
    "--bin", str(BIN),
    "--epsilon", str(EPSILON),
    "--output", str(output_path),
]
if ELECTROPRODUCTION:
    command.append("--electro")
if STARTS is not None:
    command.extend(["--starts", str(STARTS)])
if WORKERS is not None:
    command.extend(["--workers", str(WORKERS)])

# Run the fit and stop if EMI reports an error.
subprocess.run(command, cwd=PROJECT_DIR, check=True)
