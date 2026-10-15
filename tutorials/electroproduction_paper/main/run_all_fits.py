#!/usr/bin/env python3

"""Fit every configured experimental data bin.

One independent multi-start inversion is launched per channel and Q-squared
bin.  Output names follow the discovery convention used by the plot scripts.
"""

import subprocess
from tutorials.electroproduction_paper.datasets import CHANNELS, MAIN_CHANNEL_KEYS
from tutorials.paths import PROJECT_DIR

# Settings: edit the shared channel list in datasets.py to change this sweep.
TREE_NAME = "expMoments"
ELECTROPRODUCTION = True
STARTS = 10_000
WORKERS = 10
SEED = 0
USE_HESSE = True


EMI = PROJECT_DIR / "build" / "emi"
SETTINGS = (
    PROJECT_DIR / "tutorials/electroproduction_paper/main/experimental_settings.h"
)
OUTPUT_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/root"


def main():
    """Run the configured study when invoked as a Python module."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Run one fit for every bin, including each dataset's last bin.
    for name in MAIN_CHANNEL_KEYS:
        channel = CHANNELS[name]
        for bin_number in range(len(channel.q2)):
            command = [
                str(EMI),
                "fit",
                "--settings",
                str(SETTINGS),
                "--input",
                str(channel.input_file),
                "--tree",
                TREE_NAME,
                "--bin",
                str(bin_number),
                "--epsilon",
                str(channel.epsilon),
                "--output",
                str(OUTPUT_DIR / f"{name}_fit_{bin_number}.root"),
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
            subprocess.run(command, cwd=PROJECT_DIR, check=True)


if __name__ == "__main__":
    main()
