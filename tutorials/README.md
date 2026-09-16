# EMI tutorials

Each tutorial is a complete, reproducible physics study: a runner creates ROOT
results, an analysis module turns those results into figures, and a local
README records the model assumptions. Generated files remain below that
tutorial's ignored `output/` directory.

## Setup

Build and test EMI from the repository root:

```bash
python build.py --test
```

Install these Python modules in editable mode so they can be called from any
working directory:

```bash
python -m pip install -e ./tutorials
```

The package declares NumPy and Matplotlib. PyROOT is supplied by the ROOT
installation used for the C++ build and must be available in the selected
Python environment. Without installation, run modules from the repository root
so Python can resolve the `tutorials` package.

Editable installation is intentional: the runners locate the repository-local
`build/emi` executable and keep their results beside the tutorial sources. This
package is not intended as a standalone wheel detached from the checkout.

## Workflow

1. Read the case README and edit the uppercase settings near the top of its
   runner. In particular, check random starts, bootstrap samples, worker count,
   seed, and overwrite policy before a large job.
2. Run the case's `run` or `run_*` module. All commands use `subprocess.run`
   with `check=True`, so a failed C++ job stops the scan immediately.
3. Run the matching `analyse` module. It reads the produced ROOT trees and
   writes PDF figures below `output/plots`.

Use a nonzero seed when results must be reproducible. In EMI, seed zero asks
ROOT to initialize from a changing machine seed.

The shared [`analysis_utils.py`](analysis_utils.py) documents tree conversion,
best-minimum selection, circular phase statistics, and amplitude branch names.
Reaction-specific channel metadata is in
[`electroproduction_paper/datasets.py`](electroproduction_paper/datasets.py).

## Tutorial families

- [Electroproduction paper analysis](electroproduction_paper/README.md): fits
  published rho, omega, phi, and GlueX moment data and compares Hessian and
  bootstrap uncertainty estimates.
- [Polarized moments analysis](polarized_moments/README.md): closure tests for
  spin-resolved amplitudes, photon-helicity sectors, and a mass-dependent
  photoproduction model.

For example:

```bash
python -m tutorials.polarized_moments.fixed_double.run
python -m tutorials.polarized_moments.fixed_double.analyse
```

ROOT files can be large and are ignored by default. The published input files
under `electroproduction_paper/input/` are explicitly tracked; do not place
temporary outputs there.
