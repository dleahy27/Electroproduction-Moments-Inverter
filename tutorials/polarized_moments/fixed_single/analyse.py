#!/usr/bin/env python3

"""Compare the generated and fitted single-polarization examples.

Generation and fitting use the same explicit gauge, so the resolved k=+1 and
k=-1 complex amplitudes can be compared directly in each case.
"""

from tutorials.analysis_utils import PROJECT_DIR
from tutorials.polarized_moments.fixed_analysis import analyse_fixed_closure
from tutorials.polarized_moments.fixed_single.run import K_GAUGE

# Analyse both files produced by run.py so a stale selector cannot accidentally
# compare a newly generated case with an older fit from the other spin mode.
POLARIZATIONS = ("initial", "recoil")
CASE_DIR = PROJECT_DIR / "tutorials/polarized_moments/fixed_single"
SHOW_PLOT = False


def main():
    """Run the configured study when invoked as a Python module."""
    for polarization in POLARIZATIONS:
        analyse_fixed_closure(
            CASE_DIR / f"output/generated/fixed_{polarization}_moments.root",
            CASE_DIR / f"output/root/fixed_{polarization}_fit.root",
            CASE_DIR / "output/plots" / polarization,
            show_plot=SHOW_PLOT,
            polarization=polarization,
            k_gauge=K_GAUGE,
        )


if __name__ == "__main__":
    main()
