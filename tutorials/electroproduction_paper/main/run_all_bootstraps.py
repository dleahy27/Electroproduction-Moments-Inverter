#!/usr/bin/env python3

"""Bootstrap every configured experimental data bin.

The Cartesian product of configured channels and Q-squared bins defines the
batch.  Checked subprocess calls stop the sweep at the first failed analysis
instead of leaving a silently incomplete output collection.
"""

import subprocess
from tutorials.electroproduction_paper.datasets import CHANNELS, MAIN_CHANNEL_KEYS
from tutorials.paths import PROJECT_DIR

# Settings: edit the shared channel list in datasets.py to change this sweep.
TREE_NAME = "expMoments"
ELECTROPRODUCTION = True
TOYS = 1_000
STARTS_PER_TOY = 1_000
WORKERS = 10
SEED = 0
USE_HESSE = False


EMI = PROJECT_DIR / "build" / "emi"
SETTINGS = (
    PROJECT_DIR / "tutorials/electroproduction_paper/main/experimental_settings.h"
)
OUTPUT_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/root"


def main():
    """Run the configured study when invoked as a Python module."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Run one bootstrap for every bin, including each dataset's last bin.
    for name in MAIN_CHANNEL_KEYS:
        channel = CHANNELS[name]
        for bin_number in range(len(channel.q2)):
            command = [
                str(EMI),
                "bootstrap",
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
                str(OUTPUT_DIR / f"{name}_bootstrap_{bin_number}.root"),
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
            subprocess.run(command, cwd=PROJECT_DIR, check=True)


if __name__ == "__main__":
    main()
