# Configuration and worked examples

This reference describes model settings, input conventions and the physics
behind the supplied examples. For installation, see the [main README](../README.md).

## Experimental mass-fit inputs

For experimental `mass-fit` inputs, the settings header can restore normalised
moments and select the fitted mass bins:

```cpp
auto& config = settings.massFit;
config.yieldBranch = "yield"; // Scalar unnormalised zeroth moment per entry.
config.normalisation = 2.0;   // Input zeroth-moment normalisation.
config.bins = {0, 2, 5};      // Optional original zero-based TTree entry indices.
config.minMassGeV = 1.1;      // Optional inclusive lower mass bound.
config.maxMassGeV = 1.5;      // Optional inclusive upper mass bound.
```

Setting `yieldBranch` multiplies every input moment and its error by
`yield / normalisation` before chi-square evaluation. Supply the unnormalised
zeroth moment in the same convention as the model; convert experimental event
yields to that convention first. Yield values and the denominator must be
positive and finite. Yield uncertainty and induced correlations are not
propagated. With no yield branch, moments/errors are read unchanged.

Empty `bins` and unset mass bounds select all entries. Index and mass selections
intersect, bounds include their endpoints, and duplicate indices do not repeat
observations. Either mass bound can be used alone. For a contiguous index range,
populate `bins` with a loop. Invalid indices/bounds and empty selections are
rejected. Moment discovery uses the first selected entry; predictions and `ndf`
use only selected entries. These controls apply only to `mass-fit`.


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

The Pauli indices follow the axial-vector reflection convention. Components
1 and 3 are the in-plane spin components and are odd under reflection in the
production plane; component 2 is normal to that plane and is even. This parity
assignment retains the interference between the `k=+1` and `k=-1` sectors,
which determines their relative phase. Only a phase common to every amplitude
remains unobservable. In a single-polarization fit, the additional real
unnatural reference described above fixes the SO(2) spin-basis angle; it is not
a second independently adjustable reflectivity phase.

The target--recoil replacement rules are listed in
[PolarizedModel.C](../src/PolarizedModel.C). Each of the 16 `(beta, delta)`
entries specifies the replacement helicities and their real or imaginary
prefactor.

An unpolarized fit accepts either `RH_<alpha>_<L>_<M>` or the
`beta=delta=0` projection `RH_<alpha>_0_0_<L>_<M>`. The same fallback applies
to electroproduction `RH04` branches. A polarized generated file can therefore
be fitted directly with `--polarization none`, without duplicate alias branches
or a conversion step.

EMI takes an explicit wave list. Only the moments needed for those waves
are constructed.

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

Parameter limits, minimizer tolerances, and step sizes are implementation
defaults. Magnitude start distributions can be set per amplitude in
`FitConfig::initialAmplitudes`:

```cpp
fit.initialAmplitudes = {
    // reflectivity, orientation, l, m, k, mean, standard deviation
    {'a', 'T', 0, 0, 0, 0.30, 0.05},
    {'b', 'T', 0, 0, 0, 0.15, 0.03},
    // ...one row for every independent magnitude...
};
```

Use `k=0` for an unpolarized model and `k=+1` or `k=-1` when nucleon
polarization is enabled. Omit the positive-reflectivity transverse amplitude
at the highest selected `(l,m)`: its magnitude is derived from the
normalization constraint. Fixed-zero parity amplitudes are omitted as well.
The configured Gaussian distributions are activated only when the list covers
every other independent magnitude. An empty or incomplete list keeps the
default broad random initialization for the entire start; phases remain
uniform on their allowed interval. Complete configured draws are rejected and
resampled when their weighted squared magnitudes lie outside the normalization
manifold, so their relative values are not silently rescaled.

## Run a fit

The supplied script fits electron-rho bin zero:

```bash
python -m tutorials.electroproduction_paper.main.run_fit
```

It uses the model and run defaults in `app/UserSettings.h`. Input paths, the
bin number, epsilon, and optional run overrides are variables at the top of
the script:

```bash
STARTS = 100
WORKERS = 2
```

For another input file:

```bash
emi fit \
  --electro \
  --input tutorials/electroproduction_paper/input/e_rho_moments.root \
  --tree expMoments \
  --bin 0 \
  --epsilon 0.8 \
  --output tutorials/electroproduction_paper/main/output/root/e_rho_fit_0.root
```

The underlying `emi` executable still accepts command-line values. Run
`emi --help` to see them.

Analytical gradients are used by default. Set
`.UseNumericalGradients(true)` in `UserSettings.h`, or pass
`--numerical-gradients` to use numerical derivatives. To compare their timing:

```bash
python -m tutorials.electroproduction_paper.performance.analyse
```

## Run a bootstrap

```bash
python -m tutorials.electroproduction_paper.main.run_bootstrap
```

For a short installation check, set these variables at the top of the script:

```bash
TOYS = 2
STARTS_PER_TOY = 2
WORKERS = 1
```

Each bootstrap sample draws SDMEs from Gaussians centred on the published
values, using the published errors added in quadrature as their widths. The
fitter keeps the best valid minimum from its random starts. HESSE is off by
default for bootstrap runs; uncertainties come from the sample distribution.

The direct command is:

```bash
emi bootstrap \
  --electro \
  --input tutorials/electroproduction_paper/input/e_rho_moments.root \
  --tree expMoments \
  --bin 0 \
  --epsilon 0.8 \
  --output tutorials/electroproduction_paper/main/output/root/e_rho_bootstrap_0.root
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

The default destination is
`tutorials/electroproduction_paper/input/<dataset>_moments.root`. Use
`--output FILE` to choose a different path. These published inputs are the only
ROOT files intentionally tracked by Git.

The editable fixed-moment generator and the multi-event random generator are:

```bash
emi generate-fixed

emi generate-random \
  --events 100000 \
  --epsilon 0.8 \
  --seed 12345 \
  --output random_input_moments.root
```

The deterministic closure-test point used by the polarized tutorials is also
available directly:

```bash
emi generate-example \
  --settings app/UserSettings.h \
  --polarization both \
  --output fixed_both_moments.root
```

Unlike `generate-fixed`, this command uses EMI's internal example amplitudes;
it still reads the wave set and polarization-related model choices from the
settings file.

`emi --generate-fixed` is accepted as an equivalent spelling.
For batch scripts, `generate-fixed` also accepts `--settings FILE`,
`--output FILE`, `--epsilon VALUE`, `--seed N`, `--mass VALUE`,
`--noise VALUE`, `--k-asymmetry VALUE`, `--no-background`, `--quiet`, and
`--verbose`.
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
- optional relative moment noise via `noise`;
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

Set `generation.noise` to a nonnegative relative width, for example `0.05`,
to smear every generated non-normalization moment with a Gaussian of width
`noise * abs(moment)`. The value is absent by default, leaving the generated
moments exact. The normalized `H04_0_0` moment stays exactly equal to the
configured normalization target. Nonzero smeared moments store their injected
width in the corresponding `_err` branch; use `--noise VALUE` to override the
setting for one generation run.

#### PhotoTest mass model

`PhotoTest` is an internal generation model under `src/MassModels`. It returns
a fresh complete complex-amplitude set from `Evaluate(double massGeV)` on every
call. Four constant-width Breit-Wigners provide the fixed S-, P-, and D-wave
contributions, and the two D-wave resonances are summed coherently before the
result is converted to magnitude and phase. Fixed complex backgrounds can be
enabled for the four S-wave sectors. The signed `kAsymmetry` control, `s`,
multiplies each complete `k=+1` amplitude by `(1+s)/2` and each complete
`k=-1` amplitude by `(1-s)/2`, including backgrounds. It accepts values in
`[-1, 1]`; the endpoints select one sector and zero splits the two equally.
Every `k=-1` amplitude receives a fixed phase offset for its `(l,m)` wave;
the hard-coded offsets are shared by both reflectivities and remain unchanged
throughout a scan.

Configure it in `FixedMoments()` or override its scan controls at runtime:

```bash
./build/emi generate-fixed --mass 1.306 --k-asymmetry 0.5
./build/emi generate-fixed --mass 1.306 --no-background
```

Run the standard two-K truth/single-K fit scan with:

```bash
python -m tutorials.polarized_moments.photo_test_mass_scan.run
python -m tutorials.polarized_moments.photo_test_mass_scan.analyse
```

The runner uses 20 linearly spaced masses from 0.5 to 2 GeV by default; edit
`MASS_MIN`, `MASS_MAX`, `MASS_BINS`, and `ASYMMETRIES` at the top of `run.py`.
It writes a CSV manifest beside the fit outputs and resumes existing ROOT files
by default. Set `OVERWRITE=True` only when old files are incompatible. Edit
`STARTS`, `WORKERS`, and the nonzero base `SEED` to control cost and
reproducibility.

For every wave, the analysis writes magnitude and phase curves. For every mass
it writes a 2-by-4 signed-error summary: the first two columns contain
`|A_fit|-A_quad`, the last two contain `|A_fit|^2-A_quad^2`, the top row shows
individual natural/unnatural waves, and the bottom row sums over the waves.
The summed panels autoscale independently. All products stay below the
tutorial's `output/` directory.
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
`k_asymmetry`, `background_enabled`, the phase removed by the rotation, the
unnormalised zeroth moment, and the common normalisation scale alongside the generated
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

The single- and double-polarization tutorials generate and fit matched model
copies:

```bash
python -m tutorials.polarized_moments.fixed_single.run
python -m tutorials.polarized_moments.fixed_double.run
```

Set `K_GAUGE = "b:T:1:0"` at the top of the script to use a particular
single-polarized gauge reference. Set `USE_HESSE = True` to enable the Hessian.
Each runner uses the internal deterministic `generate-example` command so the
generated and fitted wave sets are the same. They do not depend on the current
contents of the user-editable `FixedMoments()` configuration.

Run the corresponding analyses with:

```bash
python -m tutorials.polarized_moments.fixed_single.analyse
python -m tutorials.polarized_moments.fixed_double.analyse
```

The analysis reads the common polarized amplitude branches and plots the generated and fitted complex
amplitudes directly in the gauge imposed by generation and minimization. It
does not apply a phase rotation or a recoil-angle transformation.

### Generic polarized-to-unpolarized fixed-amplitude analysis

After generating a custom polarized point and fitting its unpolarized moments:

```bash
./build/emi generate-fixed
./build/emi fit
```

run the unpolarized-k tutorial. The user inputs are grouped near the top of
each module. The analysis
discovers the selected waves and common moments automatically, selects the
lowest-chi-square accepted fit, and uses no pandas. It produces:

- an Argand plot showing the two generated K amplitudes and the fitted
  unpolarized amplitude separately as phase context;
- fitted magnitude versus the quadrature of the generated K magnitudes;
- direct closure of the fitted moments against the generated
  `beta=delta=0` moments;
- a generated-versus-fitted `R` comparison when longitudinal waves are present.

The figures are written under
`tutorials/polarized_moments/unpolarized_k/output/plots/`. The script uses ROOT
RDataFrame, NumPy, and Matplotlib.

### Configured polarized truth fitted as unpolarized

The unpolarized-k runner uses the same `FixedMoments()` path as ordinary fixed
generation. It generates the configured point and fits only its unpolarized
`beta=delta=0` projection:

```bash
python -m tutorials.polarized_moments.unpolarized_k.run
python -m tutorials.polarized_moments.unpolarized_k.analyse
```

In `Custom` mode, the exact and random amplitude rows in `FixedMoments()` are
used directly. `REPEATS` generates independent repetitions by incrementing
the seed. For either dominant-K mode, a list such as
`SUPPRESSIONS = [0.2, 0.1, 0.05]` performs a suppression sweep; this option is
rejected for `Custom` and `AllRandom` because suppression has no meaning there.
`USE_HESSE = True` enables Hessian calculation. Edit `SETTINGS_FILE` to use a
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

Original input values always remain readable from the experimental or
generated moment file.
