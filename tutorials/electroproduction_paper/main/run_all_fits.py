#!/usr/bin/env python3

"""Fit every configured experimental data bin."""

import subprocess
from tutorials.electroproduction_paper.datasets import CHANNELS, MAIN_CHANNEL_KEYS
from tutorials.paths import PROJECT_DIR


# Settings: edit the shared channel list in datasets.py to change this sweep.
TREE_NAME = "expMoments"
ELECTROPRODUCTION = True
STARTS = None  # None uses the value in app/UserSettings.h.
WORKERS = None  # None uses the value in app/UserSettings.h.


EMI = PROJECT_DIR / "build" / "emi"
OUTPUT_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main/output/root"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Run one fit for every bin, including each dataset's last bin.
for name in MAIN_CHANNEL_KEYS:
    channel = CHANNELS[name]
    for bin_number in range(len(channel.q2)):
        command = [
            str(EMI), "fit",
            "--input", str(channel.input_file),
            "--tree", TREE_NAME,
            "--bin", str(bin_number),
            "--epsilon", str(channel.epsilon),
            "--output", str(OUTPUT_DIR / f"{name}_fit_{bin_number}.root"),
        ]
        if ELECTROPRODUCTION:
            command.append("--electro")
        if STARTS is not None:
            command.extend(["--starts", str(STARTS)])
        if WORKERS is not None:
            command.extend(["--workers", str(WORKERS)])
        subprocess.run(command, cwd=PROJECT_DIR, check=True)
