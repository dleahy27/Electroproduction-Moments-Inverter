# PhotoTest mass-dependent studies

* [closure](closure/README.md): the one-wave Breit-Wigner closure. Its analysis
  checks mass, width, magnitude and moment predictions, and plots parameters,
  spectra, residuals and starts.
* [full_closure](full_closure/README.md): the four-resonance PhotoTest closure
  with explicit k=+/-1.
* [single_k_test](single_k_test/README.md): generate two-k truth and fit an
  effective k=0 model.
  This is a model-mismatch study. Effective couplings need not equal either
  injected sector, and general two-k truth need not have an exact k=0 solution.

Run from the repository root after `python build.py`:

```bash
python -m tutorials.photo_test_m.closure.run
python -m tutorials.photo_test_m.closure.analyse --check
python -m tutorials.photo_test_m.full_closure.run
python -m tutorials.photo_test_m.full_closure.analyse
python -m tutorials.photo_test_m.single_k_test.run
python -m tutorials.photo_test_m.single_k_test.analyse
```

Physics choices and defaults belong in each checked-in `settings.h`. The
runner never edits headers. EMI loads that header with ROOT's interpreter,
then applies the scan coordinates. Edit `ASYMMETRIES` at the top of
`single_k_test/run.py` to choose the scan. For a short run, reduce
`config.fit.starts` in its `settings.h`. The default single-k scan uses s=0
and s=1, 50,000 starts per asymmetry, 14 workers, and 20 mass points.

The first global fit start preserves every configured initial coordinate.
Later starts randomize free coupling magnitudes/phases within their bounds;
pole masses and widths always start from the header. The reference coupling
phase stays at its configured value. Each start uses a fresh Minuit2/Migrad
instance. `--max-calls` sets its budget; failed starts remain in the results
but are excluded from best-fit predictions. HESSE runs only after convergence.

Mass-generation scan overrides include `--mass-bins`, `--mass-min`,
`--mass-max`, and, with `--photo-test`, `--k-asymmetry`. `--polarization both`
selects explicit k sectors; `none` selects k=0. `--photo-test`
rebuilds the resonance table, so omit it when fitting a custom table.

All analyses use ROOT RDataFrames, NumPy, and Matplotlib. Use a Python
interpreter matching the Python version of the installed ROOT build.
All inputs are noiseless; positive assigned uncertainties define fit weights.
A small chi-square alone is insufficient evidence of parameter recovery.

Parameter markers show the best accepted fit. The bars span the minimum and
maximum values among the best 0.1% of converged fits, ranked by chi-square for
each scan case. The count is rounded up, with at least one fit retained. These
are ranges across minimizations, not statistical error estimates. The mass-fit
plots use the stored coupling phases; the mass-independent amplitude plots
use phases in [-pi, pi], so a range can cross that interval's boundary.
