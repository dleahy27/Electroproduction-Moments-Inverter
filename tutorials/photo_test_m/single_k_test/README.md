# Effective single-k PhotoTest study

The generator retains full k=+/-1 truth. The fit uses only its unpolarized
moments with k=0 amplitudes. General two-k truth has no exact single-k
representation: recovered effective couplings need not equal either sector.
At s=+1 the supplied unsplit k=0 couplings give an exact truth start.

Edit `ASYMMETRIES` in `run.py` (default: s=0 and s=1) and `settings.h`
for physics choices and job size (50,000 starts, 14 workers, 20 mass points).
Each minimization fits all mass points simultaneously. The runner generates
truth, fits each asymmetry, then merges the results:

```bash
python -m tutorials.photo_test_m.single_k_test.run
python -m tutorials.photo_test_m.single_k_test.analyse
```

For a short run, reduce `config.fit.starts` in `settings.h`.
The runner never rewrites settings headers.

`output/scan.root` contains the moments, fit results, predictions, truth
parameters, case identifiers, and asymmetries for the completed scan. Analysis
reads this file alone; no JSON manifest is needed. Individual case ROOT files
remain available for inspection. Older files outside the current scan are
excluded from the merged results.

Analysis selects the same best accepted finite minimum as the C++ writer and
plots intensity and residuals, pole masses and widths, and couplings. Parameter
bars span the minimum and maximum among the best 0.1% of converged fits, ranked
by chi-square. The marker stays at the best fit; the bars show the spread
across minimizations.
