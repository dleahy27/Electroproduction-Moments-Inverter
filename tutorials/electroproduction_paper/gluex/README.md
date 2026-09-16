# GlueX photoproduction

This tutorial fits and bootstraps the 18 momentum-transfer bins in
`../input/gluex_moments.root`. The analysis plots the selected
S- and P-wave amplitudes in magnitude-phase and Argand coordinates and shows
their scaling with momentum transfer.

Edit the fit and bootstrap costs at the top of `run.py`. Either stage can be
disabled with `RUN_FITS` or `RUN_BOOTSTRAPS`.

```bash
python -m tutorials.electroproduction_paper.gluex.run
python -m tutorials.electroproduction_paper.gluex.analyse
```

ROOT results are written to `output/root`; figures are written to
`output/plots`. The analysis uses the momentum-transfer values declared near
the top of `analyse.py` in the same order as the input bins.

The fit is photoproduction, so no virtual-photon longitudinal contribution is
introduced. Magnitude/phase plots show numerical trends, while Argand plots
retain the correlated motion of the real and imaginary parts.
