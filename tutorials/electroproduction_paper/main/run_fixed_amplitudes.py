#!/usr/bin/env python3

"""Generate and fit the fixed photoproduction closure example."""

import subprocess

from tutorials.paths import PROJECT_DIR


# Settings: reduce the fit cost here for a quick check.
SETTINGS_FILE = "app/UserSettings.h"
STARTS = 10_000
WORKERS = 1
SEED = 12_345


CASE_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/main"
EMI = PROJECT_DIR / "build/emi"
truth_file = CASE_DIR / "output/generated/fixed_photomoments.root"
fit_file = CASE_DIR / "output/root/fixed_photomoments_unpolarized_fit.root"
truth_file.parent.mkdir(parents=True, exist_ok=True)
fit_file.parent.mkdir(parents=True, exist_ok=True)

# Generate the two-k truth configured by FixedMoments() in UserSettings.h.
subprocess.run(
    [
        str(EMI), "generate-fixed",
        "--settings", str(PROJECT_DIR / SETTINGS_FILE),
        "--output", str(truth_file),
        "--quiet",
    ],
    cwd=PROJECT_DIR,
    check=True,
)

# Fit the unpolarized projection of the generated moments.
subprocess.run(
    [
        str(EMI), "fit",
        "--settings", str(PROJECT_DIR / SETTINGS_FILE),
        "--input", str(truth_file),
        "--tree", "genMoments",
        "--output", str(fit_file),
        "--photo",
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
