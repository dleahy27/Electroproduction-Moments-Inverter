#!/usr/bin/env python3

"""Fit and bootstrap all GlueX photoproduction momentum-transfer bins.

Each input-tree entry is one momentum-transfer bin.  The same transverse model
is fitted independently in every bin, and deterministic output names allow the
companion analysis to recover the physical bin ordering.
"""

import subprocess

from tutorials.electroproduction_paper.datasets import INPUT_DIR
from tutorials.paths import PROJECT_DIR

# Settings: the GlueX input contains 18 bins numbered 0 through 17.
BINS = range(18)
STARTS = 10_000
WORKERS = 8
TOYS = 1_000
STARTS_PER_TOY = 1_000
RUN_FITS = True
RUN_BOOTSTRAPS = True


CASE_DIR = PROJECT_DIR / "tutorials/electroproduction_paper/gluex"
EMI = PROJECT_DIR / "build/emi"
SETTINGS = CASE_DIR / "settings.h"


def main():
    """Run the configured study when invoked as a Python module."""
    input_file = INPUT_DIR / "gluex_moments.root"
    output_dir = CASE_DIR / "output/root"
    output_dir.mkdir(parents=True, exist_ok=True)

    for bin_index in BINS:
        if RUN_FITS:
            subprocess.run(
                [
                    str(EMI),
                    "fit",
                    "--settings",
                    str(SETTINGS),
                    "--photo",
                    "--input",
                    str(input_file),
                    "--tree",
                    "expMoments",
                    "--bin",
                    str(bin_index),
                    "--output",
                    str(output_dir / f"gluex_fit_{bin_index}.root"),
                    "--starts",
                    str(STARTS),
                    "--workers",
                    str(WORKERS),
                ],
                cwd=PROJECT_DIR,
                check=True,
            )

        if RUN_BOOTSTRAPS:
            subprocess.run(
                [
                    str(EMI),
                    "bootstrap",
                    "--settings",
                    str(SETTINGS),
                    "--photo",
                    "--input",
                    str(input_file),
                    "--tree",
                    "expMoments",
                    "--bin",
                    str(bin_index),
                    "--output",
                    str(output_dir / f"gluex_bootstrap_{bin_index}.root"),
                    "--toys",
                    str(TOYS),
                    "--starts-per-toy",
                    str(STARTS_PER_TOY),
                    "--workers",
                    str(WORKERS),
                ],
                cwd=PROJECT_DIR,
                check=True,
            )


if __name__ == "__main__":
    main()
