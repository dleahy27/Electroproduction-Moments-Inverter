# Electroproduction paper analysis

These tutorials reproduce the electroproduction and comparison studies used
around the paper analysis. Published ROOT moment files are versioned in the
local `input/` directory. Central channel metadata (beam, meson, `Q^2`,
invariant mass, epsilon, and published `R`) is defined once in `datasets.py`.
Generated ROOT results and figures stay inside each case directory.

Build EMI before running a case:

```bash
python build.py --test
```

The cases are:

- [Main analysis](main/README.md): experimental Hessian and bootstrap fits,
  amplitude plots, moment closure, longitudinal-to-transverse ratio, and the
  fixed-amplitude closure example.
- [Bootstrap versus Hessian](bootstrap_vs_hessian/README.md): compare central
  values and uncertainties from the two error-estimation methods.
- [GlueX](gluex/README.md): photoproduction fits and bootstraps versus momentum
  transfer.
- [Performance](performance/README.md): analytical/numerical gradient and
  external BruFit timing comparisons.

Install the tutorial package as described in the parent README, or run scripts
as Python modules from the repository root:

```bash
python -m tutorials.electroproduction_paper.main.run_fit
```

The input tree is `expMoments`. Each entry is one kinematic bin and contains
moment central values plus their uncertainty/covariance information. The fit
and bootstrap runners preserve the source-bin index in their output names.
