# Electroproduction Moments Inverter

This program reconstructs complex partial-wave amplitudes from measured
vector-meson moments. It can run ordinary fits or Gaussian bootstrap fits,
rebuild the supplied experimental moment tables, and generate synthetic/random moment
samples for closure tests.

The application is C++17 and ROOT based. CMake compiles the physics and fitting
code into a reusable library and links the small command-line program against
it.

## Requirements

- A C++17 compiler (GCC 9, Clang 10, or newer)
- CMake 3.18 or newer
- ROOT 6.26 or newer with Tree, MathMore, Minuit2, and multiprocessing support

Activate ROOT before configuring the project. For a binary ROOT installation:

```bash
source /path/to/root/bin/thisroot.sh
root-config --version
```

## Build

From the repository root:

```bash
./scripts/build.sh
```

This configures CMake and compiles a native executable at `build/emi`. The
`emi` executable is built from `app/main.C` and linked against the project's
`emi_core` library and ROOT; it is not a ROOT macro and does not need to be
launched through the `root` command.

Run it directly:

```bash
./build/emi --help
```

The supplied run scripts also use this repository-local executable, so no
installation or `PATH` changes are needed. To configure and build without the
helper script, use:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
alias emi = ./build/emi
```

### Optional installation

Install only if you want to run `emi` from any directory without writing its
path explicitly:

```bash
cmake --install build --prefix "$HOME/.local"
```

This places the executable at `$HOME/.local/bin/emi`. If that directory is not
already on `PATH`, add the following to your shell configuration. For Bash,
Zsh, and other Bourne-style shells:

```bash
export PATH="$HOME/.local/bin:$PATH"
emi --help
```

For Csh or Tcsh:

```csh
setenv PATH "${HOME}/.local/bin:${PATH}"
emi --help
```

In both cases, ROOT must still be active in the shell as described above.

## Choose the physics model

The routine choices are kept in `app/UserSettings.h`, next to the executable
source. Edit this file and rebuild. This avoids a long command line for settings
that normally stay fixed throughout a study.

The fit model is an explicit list of `(l,m)` waves:

```cpp
inline ModelConfig Model() {
  ModelConfig model;
  model
      .SetWaves({
          {1, -1},
          {1,  0},
          {1,  1},
      })
      .UseReflectivities(true, false) // positive only
      .EnforceLongitudinalParity(true);
  model.normalisationMoment = 2.0;
  return model;
}
```

`UseReflectivities(positive, negative)` accepts positive only `(true, false)`,
negative only `(false, true)`, or both `(true, true)`. At least one must be
enabled. Each selected wave must satisfy `l >= 0` and `|m| <= l`.

An explicit wave list is used instead of `lmax` and `mmax` as the max parameters resulted in null moments being created.
In particular, say you want to fit a purely D-wave process you do not want to limit time constructing S-wave and P-wave moments.
Moment models are also derived from the selected waves. Clebsch-Gordan-forbidden
moments have no terms and are not constructed. An input moment is used only when
the selected model can construct it and its quoted uncertainty is finite and
strictly positive. This is important for vector-meson-only tables containing
placeholder zero values with zero errors.

`GenerationModel()` in the same file controls the fixed and random generators.
`ConfigureFit()` and `ConfigureBootstrap()` contain the default run counts,
worker count, random seed, Hessian choice, and gradient choice:

```cpp
fit
    .SetStarts(10000)
    .SetWorkers(0)
    .SetSeed(0)
    .UseHesse(true)
    .UseNumericalGradients(false);
```

Parameter limits, minimizer tolerances, start-distribution widths, and step sizes
are implementation defaults. If one wants to change them, this will require going into the actual code itself to change.

## Run a fit

The supplied script fits electron-rho bin zero:

```bash
./scripts/run-fit.sh 0
```

It uses the model and run defaults in `app/UserSettings.h`. Environment
variables are available for quick, temporary overrides:

```bash
STARTS=100 WORKERS=2 ./scripts/run-fit.sh 0
```

For another input file:

```bash
emi fit \
  --input InputFiles/Experiment/e_rho_moments.root \
  --tree expMoments \
  --bin 0 \
  --epsilon 0.8 \
  --output OutputFiles/e_rho_fit_0.root
```

Paths, tree name, bin, and epsilon remain command-line values because they vary
between datasets. Run `emi --help` for optional runtime overrides.

Analytical gradients are used by default. Set
`.UseNumericalGradients(true)` in `UserSettings.h`, or pass
`--numerical-gradients` to use. An example script implementing a controlled comparison is provided:

```bash
./scripts/compare-gradients.sh 0 100
```

## Run a bootstrap

```bash
./scripts/run-bootstrap.sh 0
```

For a short installation check:

```bash
TOYS=2 STARTS_PER_TOY=2 WORKERS=1 ./scripts/run-bootstrap.sh 0
```

Each bootstrap independently samples from a Gaussian that is constructed from the published SDME values, with the width given by the quadrature of the published errors.
The standard minimization procedure outlined earlier is then ran for starts per worker and it retains the best valid minimum from its random starts. Hessian calculation is off
for bootstrap by default because the bootstrap distribution supplies the uncertainty.

The direct command is:

```bash
emi bootstrap \
  --input InputFiles/Experiment/e_rho_moments.root \
  --tree expMoments \
  --bin 0 \
  --epsilon 0.8 \
  --output OutputFiles/e_rho_bootstrap_0.root
```

## Experimental and synthetic inputs

The repository includes the source needed to regenerate both the fixed/random test data
and the experimental moments. Experimental tables are rebuilt from the SDME
values in `src/LeptoMoments.C` and `src/PhotoMoments.C`:

```bash
emi make-lepto --dataset e_rho
emi make-lepto --dataset mu_rho
emi make-lepto --dataset e_omega
emi make-lepto --dataset mu_omega
emi make-lepto --dataset e_phi
emi make-photo --dataset gluex
```

The default destination is `InputFiles/Experiment/<dataset>_moments.root`.
Use `--output FILE` to choose a different path.

Synthetic inputs use the waves in `GenerationModel()` and do not contain a
hard-coded list of amplitude names or moments:

```bash
emi generate-fixed \
  --epsilon 0.8 \
  --output InputFiles/Generated/fixed_test.root

emi generate-random \
  --events 100000 \
  --epsilon 0.8 \
  --seed 12345 \
  --output InputFiles/Generated/random_input_moments.root
```

Generated moment files use the same `RH04_*`/`RH_*` value and `_err` branch
convention as experimental input, so they can be passed straight back to `fit`.

## ROOT output

Single-run output contains a `fitResults` tree with one row per random start.
Bootstrap output contains a `PartialWaves` tree with one row per toy. Both store
only:

- fit identity and minimizer status (`fit_ok`/`valid`, `status`, `chi2`, etc.);
- fitted amplitude magnitudes and phases;
- the physical moments `H04` and `H1`, `H2`, `H3`, `H5`, `H6`, `H7`, `H8`
  that can be constructed from the selected waves;
- the longitudinal/transverse ratio `R` where applicable.

Single-run output also stores `err__*` parameter and physical-moment errors when
HESSE succeeds.

Raw input values always remain readable from the original experimental or
generated moment file.

## Code layout

```text
app/main.C                 command dispatch
app/UserSettings.h         editable model and run defaults
include/emi/Config.h       public configuration types and setters
include/emi/Runner.h       public library entry points
src/Input.C                ROOT input discovery and validation
src/Model.C                amplitudes and Clebsch-Gordan moment construction
src/Evaluation.C           moments, chi-square, and gradients
src/Context.C              fit context, normalisation, and Hessian products
src/Minimizer.C            Minuit2 setup
src/FitRunner.C            ordinary-fit orchestration
src/Bootstrap.C            bootstrap orchestration
src/LeptoMoments.C         electroproduction experimental tables
src/PhotoMoments.C         photoproduction experimental tables
src/FixedMoments.C         fixed synthetic input
src/RandomMoments.C        random synthetic inputs
```

`include/emi` is a normal public-header directory, not a second copy of the
program. `emi` is short for Electroproduction Moments Inverter and is also the
C++ namespace.
The library can be called directly from C++:

```cpp
#include "emi/Runner.h"

emi::ModelConfig model;
model.SetWaves({{1, -1}, {1, 0}, {1, 1}})
     .UseReflectivities(true, false);

emi::FitConfig fit;
fit.input = "InputFiles/Experiment/e_rho_moments.root";
fit.output = "OutputFiles/check.root";
fit.bin = 0;
fit.epsilon = 0.8;
fit.SetStarts(100).SetWorkers(1).SetSeed(123);

emi::RunFit(fit, model);
```

## Reproducibility

- Use a non-zero seed when comparing changes.
- Use one worker for the simplest analytical/numerical comparison.
- Multiple workers receive deterministic offsets from a non-zero base seed.
- Worker files are merged automatically and removed after a successful merge.
- Build products and generated ROOT output are ignored by Git.

## Analysis notebooks

`AnalysisScripts/` contains the existing PyROOT notebooks, these were used for my thesis work and for the plots within the paper ..... They are not required
to compile or run the inverter. Typical notebook dependencies are PyROOT, NumPy,
SciPy, and Matplotlib.

## License

See `LICENSE`.
