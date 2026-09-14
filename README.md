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

This configures CMake and compiles `build/emi`. It is a normal executable
linked against the project's `emi_core` library and ROOT; it does not need to
be launched through the `root` command.

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
alias emi='./build/emi'
```

### Optional installation

Install only if you want to run `emi` from any directory without writing its
path explicitly:

```bash
cmake --install build --prefix "$HOME/.local"
```

This places `emi` in `$HOME/.local/bin`. If that directory is not
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

The routine choices are kept in `app/UserSettings.h`. EMI loads this C++ file
at run time through ROOT, so edits take effect on the next invocation without
rebuilding the executable. Pass `--settings FILE` to use a different settings
file for a particular study.

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

Nucleon polarization is disabled by default. It can be selected in the model
configuration with

```cpp
model.UseNucleonPolarization(emi::NucleonPolarization::Initial);
```

The other polarized choices are `Recoil` and `Both`; `None` restores the
unpolarized model. The same choice can be made at run time for fits,
bootstraps, and synthetic generation:

```bash
emi fit ... --polarization initial
emi fit ... --polarization recoil
emi fit ... --polarization both
```

Polarized models retain separate nucleon non-flip (\(k=+1\)) and flip
(\(k=-1\)) amplitudes. Their branch suffix is ordered as `l_m_k`, for example
`a_T_1_0_1` and `a_T_1_0_m1`; any negative index is written with an `m`
prefix. With only one nucleon polarization measured, the unobserved
nucleon spin basis has an arbitrary common SO(2) mixing angle. Writing
\(c=\cos\theta\), \(s=\sin\theta\), and suppressing the wave labels, the
amplitudes transform as

\[
\begin{pmatrix}A'_+\\A'_-\end{pmatrix}_{\text{initial}}
=
\begin{pmatrix}c&-s\\s&c\end{pmatrix}
\begin{pmatrix}A_+\\A_-\end{pmatrix},
\qquad
\begin{pmatrix}A'_+\\A'_-\end{pmatrix}_{\text{recoil}}
=
\begin{pmatrix}c&\varepsilon s\\-\varepsilon s&c\end{pmatrix}
\begin{pmatrix}A_+\\A_-\end{pmatrix},
\]

where \(\varepsilon=+1\) for `a` reflectivity and \(-1\) for `b`. The same
angle acts on every wave because it is a change of the one unobserved nucleon
spin basis, not a wave-by-wave freedom. All single-polarized moments are
unchanged by this transformation, so the angle cannot be determined from the
fit. It must be fixed as a gauge convention.

The convention fixes the highest `(l,m)` natural transverse spin-non-flip wave
to be real and nonnegative,

\[
T^{(0)}_{l_{\max}m_{\max}}\in\mathbb{R}_{\geq0}.
\]

This removes the overall phase. For initial-only and recoil-only data, the
remaining SO(2) freedom is fixed with the highest `(l,m)` unnatural transverse
spin-non-flip wave,

\[
\operatorname{Im}U^{(0)}_{l_{\max}m_{\max}}=0.
\]

The second condition allows either sign of the real part. Internally that one
branch is therefore a signed real coordinate with its phase fixed to zero;
although its branch retains the amplitude-magnitude name, a negative stored
value means a phase of \(\pi\), not a negative physical magnitude. This avoids
introducing an extra discrete convention beyond `Im U = 0`. Both reference
waves use \(k=+1\), the spin non-flip component. A different SO(2) reference
wave can be selected at run time:

```bash
emi fit ... --polarization initial --k-gauge b:T:1:0
emi fit ... --polarization recoil --k-gauge b:L:1:1
```

The four fields are `reflectivity:orientation:l:m`. The requested wave must be
in the model; its spin non-flip component is used. The equivalent C++
configuration is

```cpp
model.UseNucleonPolarization(emi::NucleonPolarization::Initial)
     .SetKMixingGauge('b', 'T', 1, 0);
```

This selects the real reference amplitude that defines the SO(2) gauge; it does
not assign a measurable numerical value to \(\theta\). A numerical angle only has meaning
relative to another chosen spin basis. The polarized closure analysis below
therefore finds and reports the SO(2) rotation that best aligns a fitted result
with its generated reference. With both target and recoil polarization
measured, the SO(2) ambiguity is absent and only one overall reference phase is
fixed.

Polarized input moments use the tensor branch convention
`RH_<alpha>_<beta>_<delta>_<L>_<M>`, where `beta` is the initial-nucleon
Pauli index and `delta` is the recoil index. The inseparable transverse and
longitudinal response uses `RH04_<beta>_<delta>_<L>_<M>`. Initial-only data
have `delta=0`, recoil-only data have `beta=0`, and double-polarization data
may contain all values from zero to three. Every value branch has the usual
matching `_err` branch.

An unpolarized fit accepts either `RH_<alpha>_<L>_<M>` or the
`beta=delta=0` projection `RH_<alpha>_0_0_<L>_<M>`. The same fallback applies
to electroproduction `RH04` branches. A polarized generated file can therefore
be fitted directly with `--polarization none`, without duplicate alias branches
or a conversion step.

An explicit wave list is used instead of `lmax` and `mmax` (like in brufit) to avoid having to construct moments with placeholder zero values with zero errors.

`FixedMoments()` controls the user-defined fixed-moment generator described
below. `Fit()` and `Bootstrap()` contain the default run counts, worker count,
random seed, Hessian choice, and gradient choice:

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
  --electro \
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
python3 scripts/compare-gradients.py 0
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
  --electro \
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

The editable fixed-moment generator and the multi-event random generator are:

```bash
emi generate-fixed

emi generate-random \
  --events 100000 \
  --epsilon 0.8 \
  --seed 12345 \
  --output InputFiles/Generated/random_input_moments.root
```

`emi --generate-fixed` is accepted as an equivalent spelling.
For batch scripts, `generate-fixed` also accepts `--settings FILE`,
`--output FILE`, `--epsilon VALUE`, `--seed N`, `--mass VALUE`,
`--k-minus-scale VALUE`, `--no-background`, `--quiet`, and `--verbose`.
`--suppression VALUE` overrides the configured value only when
`FixedMoments()` selects `KPositiveDominant` or `KReflectivitySplit`.
The mass-model options apply only when `FixedMoments()` selects `PhotoTest`.
The wave set, polarization mode, generation mode, and amplitude lists remain
in the settings file.

Generated moment files use the same `RH04_*`/`RH_*` value and `_err` branch
convention as experimental input, so they can be passed straight back to `fit`.

### Editable fixed-moment generation

Edit `FixedMoments()` in `app/UserSettings.h` and run `emi generate-fixed`.
No rebuild is needed. All choices for one synthetic point are together there:

- output, photo/electroproduction, epsilon, and random seed;
- waves, reflectivities, and `None`, `Initial`, `Recoil`, or `Both` nucleon
  polarization;
- the generation mode;
- exact and individually random amplitudes in `Custom` mode.

The available modes are:

- `Custom`: combine exact and individually random amplitudes; everything else
  is zero;
- `AllRandom`: randomly populate every amplitude in the selected model;
- `KPositiveDominant`: suppress every `k=-1` amplitude;
- `KReflectivitySplit`: favour `k=+1` for `a` and `k=-1` for `b`;
- `PhotoTest`: evaluate the deterministic photoproduction S/P/D mass model
  described below.

The two dominant-K modes require `NucleonPolarization::Both`. For any
non-custom mode, leave both amplitude lists empty and use `randomMinimum`,
`randomMaximum`, and `suppression` to configure it. `AllRandom` with `Both`
polarization is the independent-random K baseline.

#### PhotoTest mass model

`PhotoTest` is an internal generation model under `src/MassModels`. It returns
a fresh complete complex-amplitude set from `Evaluate(double massGeV)` on every
call. Four constant-width Breit-Wigners provide the fixed S-, P-, and D-wave
contributions, and the two D-wave resonances are summed coherently before the
result is converted to magnitude and phase. Fixed complex backgrounds can be
enabled for the four S-wave sectors. `kMinusScale` multiplies each complete
`k=-1` amplitude, including its background.

Configure it in `FixedMoments()` or override its scan controls at runtime:

```bash
./build/emi generate-fixed --mass 1.306 --k-minus-scale 0.5
./build/emi generate-fixed --mass 1.306 --no-background
```

Run the standard two-K truth/single-K fit scan over 20 logarithmic mass points
from 0.5 to 12 GeV and five K-minus scales from 0 to 0.3 with:

```bash
python3 scripts/run-photo-test-mass-scan.py
python3 AnalysisScripts/analyse-photo-test-mass-scan.py
```

The runner writes a CSV manifest beside the fit outputs and resumes from
existing ROOT files unless `--overwrite` is supplied. Use `--starts` and
`--workers` to control the fit cost. The analysis uses PyROOT, NumPy, and
Matplotlib without pandas. For every fitted amplitude it writes absolute and
relative magnitude-error curves and shared-axis Argand panels under
`OutputFiles/photo_test_mass_scan/figures`.
The error reference is the incoherent two-K magnitude
`sqrt(|T+|^2 + |T-|^2)`, consistent with the bilinears used to generate the
moments. The Argand panels show both complex truth sectors separately and
rotate each truth reflectivity into the phase convention fixed by its
single-K fit.

Generation keeps both reflectivities and both independent K sectors in the
truth model. Before writing amplitudes, it rotates the complete amplitude set
by the phase of `a_T_2_2_1`, which is the phase reference used by the polarized
generation context. The rotation leaves all moments unchanged and writes
`aphi_T_2_2_1` as exactly zero. It then applies the existing common amplitude
normalisation so that `RH_0_0_0_0_0` equals the configured normalisation
moment. The ROOT file stores `mass_model`, `invariant_mass_GeV`,
`k_minus_scale`, `background_enabled`, the phase removed by the rotation, the
raw zeroth moment, and the common normalisation scale alongside the generated
truth amplitudes and moments.

In `Custom` mode an exact row is:

```cpp
// reflectivity, orientation, l, m, k, magnitude, phase [radians]
{'a', 'T', 1, 0, +1, 0.35, -0.40},
```

An individually random row is:

```cpp
// reflectivity, orientation, l, m, k, final minimum, final maximum
{'a', 'T', 1, 1, +1, 0.02, 0.10},
```

Its magnitude is sampled uniformly between the final `minimum` and `maximum`,
and its phase is sampled uniformly from \([-\pi,\pi]\). EMI rejects any
amplitude placed in both the fixed and random lists, including equivalent
longitudinal entries after parity is applied.

Use `k=0` for `NucleonPolarization::None` and `k=+1` or `k=-1` for a polarized
model. Photoproduction accepts transverse (`T`) amplitudes only. In `Custom`
mode, amplitudes not listed are exactly zero. Photoproduction models do not
construct or write longitudinal amplitude branches.

Do not list the normalization amplitude. EMI derives the
positive-reflectivity transverse amplitude at the highest selected `(l,m)`,
using `k=+1` for a polarized model and `k=0` otherwise. Both photo- and
electroproduction use

\[
w_{\mathrm{ref}}|A_{\mathrm{ref}}|^2 = \frac{N}{2}
 - \sum_{i\ne\mathrm{ref}} w_i |A_i|^2,
\]

where `N` is `model.normalisationMoment` (normally 2), and the weights are read
from the constructed \(H_0+\epsilon H_4\) model. For the usual complete wave
list this is the familiar sum of transverse intensities plus the
epsilon-weighted longitudinal intensities; the longitudinal part is absent in
photoproduction. Reading the weights from the model also keeps arbitrary
explicit wave lists normalized exactly. A fully random mode first samples all
amplitudes and then rescales them together, preserving their relative sizes.
The generator reports an error if custom amplitudes exceed the normalization,
or if an entry names a wave outside the model, uses an invalid `k`, or violates
a fixed reference phase.

### Polarized fixed-amplitude closure test

The polarized closure runner generates and fits matched copies of `Model()` for
initial, recoil, and double nucleon polarization:

```bash
STARTS=100 WORKERS=1 SEED=12345 \
  ./scripts/run-polarized-fixed-test.sh
```

Set `K_GAUGE=b:T:1:0` to use a particular single-polarized gauge reference.
`HESSE=1` enables the Hessian.
The runner uses the internal deterministic `generate-example` command so the
generated and fitted wave sets are the same. It does not depend on the current
contents of the user-editable `FixedMoments()` configuration.

Open `AnalysisScripts/polarized_fixed_amp_analysis.ipynb` and set its generated
and fitted ROOT filenames to inspect one of these cases. The notebook reads the
common polarized amplitude branches and plots the generated and fitted complex
amplitudes directly in the gauge imposed by generation and minimization. It
does not apply a phase rotation or a recoil-angle transformation.

### Generic polarized-to-unpolarized fixed-amplitude analysis

After generating a custom polarized point and fitting its unpolarized moments:

```bash
./build/emi generate-fixed
./build/emi fit
```

open `AnalysisScripts/unpolarized_k_dominance_R.ipynb`. The only user inputs
are the generated and fitted ROOT filenames near the top of the notebook. It
discovers the selected waves and common moments automatically, selects the
lowest-chi-square accepted fit, and uses no pandas. It produces:

- an Argand plot showing the two generated K amplitudes and the fitted
  unpolarized amplitude separately as phase context;
- fitted magnitude versus the quadrature of the generated K magnitudes;
- direct closure of the fitted moments against the generated
  `beta=delta=0` moments;
- a generated-versus-fitted `R` comparison when longitudinal waves are present.

The figures are written under
`AnalysisScripts/outputs/polarized_to_unpolarized_fixed/`. The notebook uses
PyROOT, NumPy, and Matplotlib, like the historic analysis notebooks.

### Configured polarized truth fitted as unpolarized

The historical K-dominance runner now uses the same `FixedMoments()` path as
ordinary fixed generation. It generates the configured point and fits only its
unpolarized `beta=delta=0` projection:

```bash
STARTS=10000 WORKERS=1 SEED=12345 \
  ./scripts/run-unpolarized-k-dominance-test.sh
```

In `Custom` mode, the exact and random amplitude rows in `FixedMoments()` are
used directly. `REPEATS=N` generates independent repetitions by incrementing
the seed. For either dominant-K mode, an optional list such as
`SUPPRESSIONS="0.2 0.1 0.05"` performs a suppression sweep; this option is
rejected for `Custom` and `AllRandom` because suppression has no meaning there.
`HESSE=1` enables Hessian calculation. Use `SETTINGS=path/to/settings.h` for a
separate study file. Ensure `Fit().photoproduction` describes the same process
as `FixedMoments().photoproduction`.

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
app/RuntimeSettings.C      runtime C++ settings loader
app/UserSettings.h         runtime fit and generation settings
include/emi/Config.h       public configuration types and setters
include/emi/Runner.h       public library entry points
src/Input.C                ROOT input discovery and validation
src/Model.C                amplitudes and Clebsch-Gordan moment construction
src/Evaluation.C           moments, chi-square, and gradients
src/Context.C              fit context, normalisation, and Hessian products
src/Minimizer.C            Minuit2 setup
src/FitRunner.C            ordinary-fit orchestration
src/Bootstrap.C            bootstrap orchestration
src/PolarizedModel.C       polarized tensor-moment construction
src/LeptoMoments.C         electroproduction experimental tables
src/PhotoMoments.C         photoproduction experimental tables
src/FixedMoments.C         fixed synthetic input
src/RandomMoments.C        random synthetic inputs
src/MassModels/MassModel.* internal shared mass-model value types and BW helper
src/MassModels/PhotoTest.* deterministic photoproduction S/P/D model
tests/MassModelsTest.C     PhotoTest formula and generation integration checks
scripts/run-polarized-fixed-test.sh
                           three polarized closure runs
AnalysisScripts/polarized_fixed_amp_analysis.ipynb
                           polarized fixed-point Argand closure
AnalysisScripts/unpolarized_k_dominance_R.ipynb
                           one polarized-truth/unpolarized-fit comparison
scripts/run-unpolarized-k-dominance-test.sh
                           configured fixed generation and unpolarized fit
scripts/run-photo-test-mass-scan.py
                           PhotoTest mass/scale generation and fitting scan
AnalysisScripts/analyse-photo-test-mass-scan.py
                           amplitude-error and Argand scan figures
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
