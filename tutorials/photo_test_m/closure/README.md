# Single-Breit-Wigner mass-dependent closure

This tutorial performs one fit to moments from every mass bin. It differs
from a mass-independent inversion, where each bin receives its own unrelated
complex amplitudes: here the bins share pole masses, widths, and production
couplings. It is a closure test, so generated moments and fitted moments use
the same model.

Run from the repository root after building EMI:

```bash
python -m tutorials.photo_test_m.closure.run
python -m tutorials.photo_test_m.closure.analyse --check
```

`run.py` uses the editable [settings.h](settings.h) steering card. It
generates `output/moments.root`, fits it, and writes `output/fit.root`.

## Amplitude model

For each selected partial wave $w$, EMI evaluates

$$
A_w(m) = \sum_{r\in w} c_{wr}B_r(m)+b_w,\qquad
B_r(m)=\frac{M_r\Gamma_r}{M_r^2-m^2-iM_r\Gamma_r}.
$$

$m$ is the invariant mass in one tree entry. $M_r$ and $\Gamma_r$
are the constant-width Breit-Wigner pole mass and width. $c_{wr}$ is a
complex production coupling and $b_w$ is a complex constant background.
Contributions add before moments are formed, so their relative phases produce
physical interference.

`TruthModel()` generates a D wave with $M=1.320$ GeV and
$\Gamma=0.120$ GeV. `MassFit()` uses the same physics structure but begins
from displaced mass and coupling-magnitude coordinates. This closure has no background term.

| Config coordinate | Meaning |
| --- | --- |
| `poleMass` | $M_r$, in GeV |
| `width` | $\Gamma_r$, in GeV |
| `couplingMagnitude`, `couplingPhase` | $|c_{wr}|$, $\arg c_{wr}$ |

Every coordinate has a start, bounds, and a `fixed` flag. Each
`MassWave` must also occur in `Model()`; EMI rejects a mass model that
refers to an absent angular wave.

## From amplitudes to moments

The angular-moment formalism is the same one used by EMI's ordinary
mass-independent fit. At every mass, the BW model supplies complex amplitude
values to the sparse moment equations. Each response moment has the form

$$
H_k(m)=\sum_{ij}C_{ij}^{(k)}
\begin{cases}
\operatorname{Re}[A_i(m)A_j^*(m)]\\
\operatorname{Im}[A_i(m)A_j^*(m)],
\end{cases}
$$

where the known $C_{ij}^{(k)}$ contain angular momentum, reflectivity,
photon, and polarization factors. The mass model changes amplitudes; it does
not change those angular conventions.

The input ROOT tree has one scalar entry per mass bin. Its mass branch is
`mass_GeV` by default. Every measured moment branch, such as
`RH_0_0_0`, needs a positive matching uncertainty branch,
`RH_0_0_0_err`. This generator assigns constant absolute uncertainties of 0.001,
to set the chi-square weights; no noise is added.

## Simultaneous chi-square

For all mass points $n$ and observed moments $k$, the fitter minimizes

$$
\chi^2(\textbf p)=\sum_{n,k}
\left[
\frac{H_k^{\rm data}(m_n)-H_k^{\rm model}(m_n; \textbf p)}
{\sigma_{nk}}
\right]^2.
$$

$\textbf{p}$ contains every free pole and coupling coordinate. The same pole mass, width and couplings must describe every mass bin.

The current objective uses diagonal moment uncertainties. A real analysis
with correlated moments or mass bins needs a covariance-aware objective before
the chi-square and HESSE errors have their complete statistical meaning.

## Analytic derivatives and Minuit

EMI differentiates both the BW and the moment bilinears analytically. For
each fit coordinate $p$,

$$
\frac{\partial(A_iA_j^*)}{\partial p}=
\frac{\partial A_i}{\partial p}A_j^*
+A_i\left(\frac{\partial A_j}{\partial p}\right)^*.
$$

Taking the appropriate real or imaginary part gives
$\partial H_k/\partial p$, which is accumulated into the chi-square
gradient. Minuit2 uses this gradient to find the minimum. Set
`fit.useNumericalGradients = true` for a slower independent numerical
cross-check.

With `fit.hesse = true`, Minuit's HESSE calculation estimates the local
covariance at the selected minimum. The `err__<parameter>` output branches
are nominal local one-standard-deviation uncertainties; they rely on a
quadratic chi-square approximation and the input error model. Check
multi-resonance fits with toy experiments or bootstrap refits as well.

## Why one coupling phase is held fixed

Moments are unchanged if every amplitude receives the same phase:

$$
A_w(m)\longrightarrow e^{i\phi}A_w(m).
$$

All products $A_iA_j^*$ remain identical. The common phase is therefore
unmeasurable. Leaving it free would give Minuit an exact flat direction and a
singular covariance matrix.

EMI fixes the phase of the resonance on the highest-(l,m),
positive-reflectivity transverse wave. Here that is `D2_a_2_1_phase`, held at its configured initial
value. The reference phase defines the convention; its magnitude remains free.
In this one-wave test there are no other coupling phases. Set
`MassModelConfig::SetPhaseGauge(name, wave)` to override the
default; choose the coupling expected to have the largest magnitude.

The held phase is a convention and cannot be measured in this one-wave test.
It need not equal the generator phase, since neither changes the moments.
The fit writes `mass_fit_phase_gauge_coordinate` so that
the convention is recorded in the ROOT output.

## Output and visual closure check

`massFitResults` contains one row per start with Minuit status, total `chi2`, `ndf`,
fitted coordinates, and HESSE errors. `massFitPredictions` contains one row
per input mass and `model__<moment>` branches recalculated from the
best-fitting BW model.

`analyse.py` prints the recovered pole mass and width, then produces
paginated `output/plots/input_vs_model_*.pdf` figures. Every panel overlays
an input moment mass spectrum and its error bars with the recreated model.
Agreement checks the full chain: ROOT input, coherent BW amplitudes, angular
moments, minimization, and output serialization.

For data, change `MassFit().fit.input`, add the required resonances and
backgrounds, retain a physically motivated wave set, cover enough mass range
to constrain each line shape, inspect alternative minima, and use a
covariance-aware extension when correlated input uncertainties are available.
