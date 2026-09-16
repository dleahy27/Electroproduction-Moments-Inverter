# Bootstrap versus Hessian

This case compares bootstrap means and sample widths with the best-fit values
and propagated Hessian uncertainties for matching experimental bins. It plots
moments, magnitudes, phases, and `R` where the required branches are present.

First produce both result types with the runners in the main tutorial:

```bash
python -m tutorials.electroproduction_paper.main.run_all_fits
python -m tutorials.electroproduction_paper.main.run_all_bootstraps
```

Then run:

```bash
python -m tutorials.electroproduction_paper.bootstrap_vs_hessian.analyse
```

Matching files are discovered in `../main/output/root`. Figures are written to
`output/plots/<channel>/bin_<n>/`, with the `R` summary stored in the channel
directory.

Interpret agreement as a numerical diagnostic rather than an identity: the
Hessian linearizes the objective near one solution, while the bootstrap can
show skewness, wrapping phases, and migration between ambiguous minima. Both
inputs must use the same model, bin, and branch convention.
