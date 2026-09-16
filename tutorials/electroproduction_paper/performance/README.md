# Performance comparison

This analysis compares numerical and analytical gradient fits and, when the
external files are available, EMI and BruFit timing measurements.

Place the analytical and numerical ROOT fit files in `output/root` using the
names configured at the top of `analyse.py`. Timing CSV files belong directly
in `output/`. The BruFit location is also an editable setting.

Run:

```bash
python -m tutorials.electroproduction_paper.performance.analyse
```

Missing optional timing files simply omit that comparison. Available figures
are written to `output/plots`.

Use the same wave set, data bin, stopping tolerances, and number of starts for
meaningful timings. Analytical gradients reduce repeated finite-difference
evaluations; this study measures the complete fit cost rather than timing a
single objective call in isolation.
