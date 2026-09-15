#!/usr/bin/env python3

"""Generate polarized fixed points and fit their unpolarized moments."""

import subprocess
from pathlib import Path


# Settings: edit these values before running the study.
SETTINGS_FILE = "app/UserSettings.h"
EPSILON = None  # None uses the value in the settings file.
STARTS = 1_000
WORKERS = 8
SEED = 0
REPEATS = 200
SUPPRESSIONS = [None]  # For example: [0.2, 0.1, 0.05].
USE_HESSE = False


PROJECT_DIR = Path(__file__).resolve().parents[1]
EMI = PROJECT_DIR / "build" / "emi"
SETTINGS = PROJECT_DIR / SETTINGS_FILE

# Name the batch after its repeat and suppression counts.
batch = f"N{REPEATS}"
if SUPPRESSIONS != [None]:
    batch += f"_M{len(SUPPRESSIONS)}"
input_dir = PROJECT_DIR / "InputFiles" / "Generated" / "photo_k" / batch
output_dir = PROJECT_DIR / "OutputFiles" / "photo_k" / batch
input_dir.mkdir(parents=True, exist_ok=True)
output_dir.mkdir(parents=True, exist_ok=True)

# Generate each configured fixed point, then fit its unpolarized projection.
for suppression_index, suppression in enumerate(SUPPRESSIONS):
    for repeat in range(REPEATS):
        run_seed = SEED + suppression_index * REPEATS + repeat
        stem = f"fixed_N{repeat}"
        if suppression is not None:
            stem += f"_M{suppression_index}"
        input_file = input_dir / f"{stem}_moments.root"
        output_file = output_dir / f"{stem}_unpolarized_fit.root"

        generate_command = [
            str(EMI), "generate-fixed",
            "--settings", str(SETTINGS),
            "--output", str(input_file),
            "--seed", str(run_seed),
            "--quiet",
        ]
        if EPSILON is not None:
            generate_command.extend(["--epsilon", str(EPSILON)])
        if suppression is not None:
            generate_command.extend(["--suppression", str(suppression)])

        details = f"N={repeat}, seed={run_seed}"
        if suppression is not None:
            details += f", suppression={suppression}"
        print(f"Generating configured fixed point: {details}")
        subprocess.run(generate_command, cwd=PROJECT_DIR, check=True)

        fit_command = [
            str(EMI), "fit",
            "--settings", str(SETTINGS),
            "--input", str(input_file),
            "--tree", "genMoments",
            "--output", str(output_file),
            "--polarization", "none",
            "--starts", str(STARTS),
            "--workers", str(WORKERS),
            "--seed", str(run_seed),
            "--quiet",
        ]
        if EPSILON is not None:
            fit_command.extend(["--epsilon", str(EPSILON)])
        if not USE_HESSE:
            fit_command.append("--no-hesse")

        print(f"Fitting its unpolarized moments: N={repeat}")
        subprocess.run(fit_command, cwd=PROJECT_DIR, check=True)

print(f"Wrote generated files under {input_dir}")
print(f"Wrote fitted files under {output_dir}")
