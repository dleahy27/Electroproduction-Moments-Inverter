# EMI tutorials

Each case has a runner, an analysis script and a README describing its model.
The runner writes ROOT files; the analysis writes figures. Generated files
stay in the case's ignored `output/` directory.

## Setup

From the repository root:

```bash
conda env create -f environment.yml
conda activate emi
python build.py
python -m pip install -e ./tutorials
```

The runners use `build/emi` and locate their inputs and outputs relative to the
checkout. The editable Python installation lets you run them from another
working directory. C++ ROOT and PyROOT must come from the same installation;
the `emi` environment provides both.

## Running a case

Read the case README and check its settings before starting a large job.
Wave choices and generation models are in its settings header. Run counts,
worker counts, seeds and overwrite options are either in that header or near
the top of the runner.

Run the modules with `python -m`, for example:

```bash
python -m tutorials.polarized_moments.fixed_double.run
python -m tutorials.polarized_moments.fixed_double.analyse
```

A failed C++ command stops the runner. Importing a Python module does not start
a fit or create figures. Use a nonzero seed for repeatable results; seed zero
asks ROOT to choose a seed.

## Studies

- [Electroproduction paper](electroproduction_paper/README.md): published rho,
  omega, phi and GlueX inputs, fits, and Hessian/bootstrap comparisons.
- [Polarized moments](polarized_moments/README.md): spin-resolved closure tests
  and the PhotoTest mass/asymmetry scan.
- [Mass-dependent fits](photo_test_m/README.md): a single-resonance closure,
  the full PhotoTest closure, an effective single-k fit to two-k truth, and a
  Flatté/Breit–Wigner shape test.

The [configuration reference](configuration.md) covers settings, phase
conventions, input formats and worked examples. Each study keeps its own
physics discussion and numerical assumptions.

## Analysis and output files

The scan runners combine point trees into `scan.root`; `scan_case` identifies
the input point. The polarized-moment scans also keep CSV manifests and accept
`--merge-only` to rebuild the merged file without fitting again. PhotoTest
mass-fit scans store their asymmetries in the ROOT file itself. Individual
point files remain available for inspection and refitting.

[`root_io.py`](root_io.py) handles merging. Analyses read ROOT columns through
RDataFrame, calculate with NumPy and plot with Matplotlib. Shared tree readers,
fit selection and phase statistics are in
[`analysis_utils.py`](analysis_utils.py); published channel metadata is in
[`electroproduction_paper/datasets.py`](electroproduction_paper/datasets.py).

Generated ROOT files, data tables and figures are ignored by Git. The published
ROOT inputs already tracked under `electroproduction_paper/input/` remain part
of the repository. Keep new outputs in each case's `output/` directory.
