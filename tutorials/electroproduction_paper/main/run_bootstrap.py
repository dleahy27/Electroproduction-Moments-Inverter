#!/usr/bin/env python3

"""Bootstrap one experimental moment bin.

This focused runner is the convenient place to tune toy count, fit starts, and
parallel workers before launching the complete channel sweep.
"""

import subprocess
from tutorials.electroproduction_paper.datasets import CHANNELS
from tutorials.paths import PROJECT_DIR

# Settings: edit these values to choose the data and bootstrap run.
CHANNEL = CHANNELS["e_rho"]
INPUT_FILE = CHANNEL.input_file
OUTPUT_FILE = (
    "tutorials/electroproduction_paper/main/output/root/e_rho_bootstrap_0.root"
)
TREE_NAME = "expMoments"
BIN = 0
EPSILON = CHANNEL.epsilon
ELECTROPRODUCTION = True
TOYS = 1_000
STARTS_PER_TOY = 1_000
WORKERS = 10
SEED = 0
USE_HESSE = False


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
        "bootstrap",
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
        "--toys",
        str(TOYS),
        "--starts-per-toy",
        str(STARTS_PER_TOY),
        "--workers",
        str(WORKERS),
        "--seed",
        str(SEED),
    ]
    if ELECTROPRODUCTION:
        command.append("--electro")
    if USE_HESSE:
        command.append("--hesse")

    # Run the bootstrap and stop if EMI reports an error.
    subprocess.run(command, cwd=PROJECT_DIR, check=True)


if __name__ == "__main__":
    main()
