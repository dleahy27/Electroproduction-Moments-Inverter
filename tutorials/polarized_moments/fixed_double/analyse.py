#!/usr/bin/env python3

"""Compare generated and fitted double-polarization amplitudes.

Both measured nucleon-spin indices remove the continuous spin-basis ambiguity,
allowing the resolved k=+1 and k=-1 complex amplitudes to be compared directly.
"""

from tutorials.analysis_utils import PROJECT_DIR
from tutorials.polarized_moments.fixed_analysis import analyse_fixed_closure

CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/fixed_double"
TRUTH_FILE = CASE_DIR / "output/generated/fixed_both_moments.root"
FIT_FILE = CASE_DIR / "output/root/fixed_both_fit.root"
OUTPUT_DIR = CASE_DIR / "output/plots"
SHOW_PLOT = False


def main():
    """Run the configured study when invoked as a Python module."""
    analyse_fixed_closure(
        TRUTH_FILE,
        FIT_FILE,
        OUTPUT_DIR,
        show_plot=SHOW_PLOT,
    )


if __name__ == "__main__":
    main()
