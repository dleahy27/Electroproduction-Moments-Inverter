#!/usr/bin/env python3

"""Fit every configured experimental data bin."""

import subprocess
from pathlib import Path


# Settings: each row is (file stem, last bin, epsilon).
DATASETS = [
    ("e_rho", 3, 0.8),
    ("mu_rho", 3, 0.9),
    ("e_omega", 2, 0.8),
    ("mu_omega", 2, 0.96),
]
TREE_NAME = "expMoments"
ELECTROPRODUCTION = True
STARTS = None  # None uses the value in app/UserSettings.h.
WORKERS = None  # None uses the value in app/UserSettings.h.


PROJECT_DIR = Path(__file__).resolve().parents[1]
EMI = PROJECT_DIR / "build" / "emi"

# Run one fit for every bin, including each dataset's last bin.
for name, last_bin, epsilon in DATASETS:
    for bin_number in range(last_bin + 1):
        command = [
            str(EMI), "fit",
            "--input", str(PROJECT_DIR / "InputFiles/Experiment" / f"{name}_moments.root"),
            "--tree", TREE_NAME,
            "--bin", str(bin_number),
            "--epsilon", str(epsilon),
            "--output", str(PROJECT_DIR / "OutputFiles" / f"{name}_fit_{bin_number}.root"),
        ]
        if ELECTROPRODUCTION:
            command.append("--electro")
        if STARTS is not None:
            command.extend(["--starts", str(STARTS)])
        if WORKERS is not None:
            command.extend(["--workers", str(WORKERS)])
        subprocess.run(command, cwd=PROJECT_DIR, check=True)
