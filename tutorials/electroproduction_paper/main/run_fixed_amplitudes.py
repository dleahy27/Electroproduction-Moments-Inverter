#!/usr/bin/env python3

"""Generate and fit the fixed electroproduction closure example.

The first command writes moments from known amplitudes; the second recovers an
equivalent unpolarized solution.  Running both steps together prevents a stale
truth file from being paired with changed generation settings.
"""

import subprocess

from tutorials.paths import PROJECT_DIR


# Settings: reduce the fit cost here for a quick check.
SETTINGS_FILE = "tutorials/electroproduction_paper/main/fixed_amplitudes_settings.h"
FIXED_TYPE = "large_L"
STARTS = 50_000
WORKERS = 10
SEED = 12_345


CASE_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main"
EMI = PROJECT_DIR / "build/emi"
truth_file = CASE_DIR / f"output/generated/fixed_{FIXED_TYPE}_moments.root"
fit_file = CASE_DIR / f"output/root/fixed_{FIXED_TYPE}_fit.root"
truth_file.parent.mkdir(parents=True, exist_ok=True)
fit_file.parent.mkdir(parents=True, exist_ok=True)

# Generate the deterministic single-k truth configured by the local header.
subprocess.run(
    [
        str(EMI), "generate-fixed",
        "--settings", str(PROJECT_DIR / SETTINGS_FILE),
        "--output", str(truth_file),
    ],
    cwd=PROJECT_DIR,
    check=True,
)

# Fit the generated unpolarized moments from near-truth magnitude distributions.
subprocess.run(
    [
        str(EMI), "fit",
        "--settings", str(PROJECT_DIR / SETTINGS_FILE),
        "--input", str(truth_file),
        "--tree", "genMoments",
        "--output", str(fit_file),
        "--electro",
        "--polarization", "none",
        "--starts", str(STARTS),
        "--workers", str(WORKERS),
        "--seed", str(SEED),
        "--no-hesse",
        "--quiet",
    ],
    cwd=PROJECT_DIR,
    check=True,
)
