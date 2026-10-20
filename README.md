# Electroproduction Moments Inverter

EMI reconstructs complex partial-wave amplitudes from measured moments. It
supports binwise fits, Gaussian bootstrap fits, synthetic moment generation,
and simultaneous mass-dependent fits with shared Breit–Wigner or Flatté
resonances and complex backgrounds. It is a C++17 library with a ROOT-based
command-line program and Python analysis scripts.

## Build and setup

Use a C++17 compiler, CMake 3.18 or newer and ROOT 6.26 or newer with Tree,
MathMore, Minuit2 and multiprocessing support. The supplied environment also
includes PyROOT, NumPy and Matplotlib for the tutorials:

```bash
conda env create -f environment.yml
conda activate emi
python build.py
python -m pip install -e ./tutorials
./build/emi --help
```

For the pinned Linux x86-64-v3 builds, use
`conda create -n emi --file environment_linux_64.lock` instead of the YAML recipe.
With a separate ROOT installation, activate its `thisroot.sh` before building
and use a matching PyROOT/Python version. Avoid mixing ROOT installations.

The build helper creates `build/emi`; no installation is needed. The equivalent
CMake commands are:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

`python build.py --root-dir /path/to/root/cmake` selects a specific ROOT.
Use a fresh build directory when switching compilers or ROOT installations.
An optional `cmake --install build --prefix "$HOME/.local"` installs the command
in `~/.local/bin`. That installation still depends on this source checkout and
an active ROOT environment; rebuild after moving the checkout.

## Configure and run

Edit `app/UserSettings.h`, or pass `--settings FILE` to select another header.
Settings are loaded at runtime; editing them does not require recompilation.
Explicit command-line options override the settings. For example:

```bash
./build/emi fit --settings app/UserSettings.h \
  --input tutorials/electroproduction_paper/input/e_rho_moments.root \
  --tree expMoments --bin 0 --epsilon 0.8 --starts 100 --workers 2 \
  --output fit.root
```

The main commands are:

| Command | Purpose |
| --- | --- |
| `fit` | Fit one bin of measured moments |
| `bootstrap` | Fluctuate input moments and fit each sample |
| `mass-fit` | Fit resonance parameters across a mass grid |
| `generate-mass-model` | Generate moments from configured resonances |
| `generate-fixed`, `generate-example`, `generate-random` | Generate synthetic inputs |
| `make-lepto`, `make-photo` | Rebuild the supplied experimental moment tables |

Use `./build/emi --help` for options. Analytic gradients are the default;
`--numerical-gradients` enables a numerical cross-check. A nonzero seed makes
starts repeatable for a given worker count. Each worker owns its fit state and
part file; successful runs merge and remove the parts.

Ordinary fits write `fitResults`, bootstraps write `PartialWaves`, and mass fits
write `massFitResults`, `massFitPredictions` and optionally `massFitCurves`.
Results retain every ordinary/mass-fit start, including failures. Analyses select
accepted finite minima. HESSE enables covariance and propagated errors.
Generated files belong in a case's ignored `output/` directory.

## Examples and checks

The [tutorial index](tutorials/README.md) links runnable studies. A short
installation and mass-model check is:

```bash
python -m tutorials.photo_test_m.closure.run
python -m tutorials.photo_test_m.closure.analyse --check
python -m tutorials.photo_test_m.flatte.run
python -m tutorials.photo_test_m.flatte.analyse
```

See [mass-dependent studies](tutorials/photo_test_m/README.md) for the full
closure and effective single-k scan. Model equations, phase conventions and
worked examples live in the [configuration reference](tutorials/configuration.md)
and the individual tutorial documents.

## Code layout

- `include/emi/`: public configuration types and library entry points.
- `app/`: command parsing and runtime settings.
- `src/`: angular models, fitting, generation and ROOT output.
- `src/MassModels/`: resonance shapes, ordered parameter mapping and PhotoTest.
- `tutorials/`: runners, shared ROOT readers and plots.

The [source guide](src/README.md) describes the data flow and ownership rules.
See `LICENSE` for licensing.
