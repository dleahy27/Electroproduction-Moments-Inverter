#!/usr/bin/env python3

"""Fit one experimental moment bin.

The constants below define one reproducible inversion. The local settings file
defines the P-wave model; run controls are passed as command-line overrides.
"""

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
STARTS = 10_000
WORKERS = 10
SEED = 0
USE_HESSE = True


# Build the EMI command from the settings above.


def main():
    """Run the configured study when invoked as a Python module."""
    output_path = PROJECT_DIR / OUTPUT_FILE
    output_path.parent.mkdir(parents=True, exist_ok=True)
    settings_file = (
        PROJECT_DIR / "tutorials/electroproduction_paper/main/experimental_settings.h"
    )
    command = [
        str(PROJECT_DIR / "build" / "emi"),
        "fit",
        "--settings",
        str(settings_file),
        "--input",
        str(INPUT_FILE),
        "--tree",
        TREE_NAME,
        "--bin",
        str(BIN),
        "--epsilon",
        str(EPSILON),
        "--output",
        str(output_path),
        "--starts",
        str(STARTS),
        "--workers",
        str(WORKERS),
        "--seed",
        str(SEED),
    ]
    if ELECTROPRODUCTION:
        command.append("--electro")
    if not USE_HESSE:
        command.append("--no-hesse")

    # Run the fit and stop if EMI reports an error.
    subprocess.run(command, cwd=PROJECT_DIR, check=True)


if __name__ == "__main__":
    main()
