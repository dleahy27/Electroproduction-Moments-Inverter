# Flatté conventions and published inputs

## Published input values

- f2(1270): BW mass **1.2754 GeV**, width **0.1866 GeV**, from the
  [PDG 2025 listing](https://pdg.lbl.gov/2025/listings/rpp2025-list-f2-1270.pdf).
  This is the state often described by its approximately 1275 MeV mass.
- f0(980): **m0 = 0.965 GeV**, **g_pi = 0.165**, **g_K/g_pi = 4.21**,
  from BES, M. Ablikim et al., Phys. Lett. B 607 (2005) 243,
  [Eq. (1) and Sec. 4](https://arxiv.org/pdf/hep-ex/0411001), the
  ABLIKIM 05 Flatté result cited in the
  [PDG f0 listing](https://pdg.lbl.gov/2025/listings/rpp2025-list-f0-980.pdf).
  There is no single PDG-average Flatté coupling set. We use this consistent
  measured set rather than combining the PDG mass estimate with unrelated
  decay parameters. Its m0 is not the peak position or a complex pole mass.

The BES denominator in Eq. (1) is m0²-s-i(g1 rho_pi+g2 rho_K), matching the
implemented convention. Thus the code's `pipi_g2` coefficient is **0.165 GeV²**
and `KK_g2` is **0.69465 GeV²**; do not square 0.165 or multiply it by m0.
The original paper labels g1 in mass units despite the denominator convention;
the values here are interpreted as its coefficients of rho in GeV².
For this simple example the two channels use charged daughter masses
m_pi = 0.13957039 GeV and m_K = 0.493677 GeV; it does not resolve neutral and
charged thresholds separately. Closed KK phase space is continued analytically.
The current BW is relativistic with constant width, without barrier factors or
a running breakup-momentum width.


## Configure another Flatté component

Select Flatté independently for any resonance contribution:

```cpp
MassResonanceConfig resonance;
resonance.name = "example_a0";
resonance.wave = {'a', 'T', 0, 0, 0};
resonance.poleMass = {0.980, 0.90, 1.08, false};
resonance.UseFlatte({
    {"eta_pi", 0.547862, 0.134977, {0.18, 0.01, 1.0, false}},
    {"KK", 0.493677, 0.493677, {0.40, 0.01, 2.0, false}},
});
config.model.resonances.push_back(resonance); // In MassFit().
// In MassGeneration(), add it to config.massModel.resonances instead.
```

Each channel is `{name, daughter1_GeV, daughter2_GeV, couplingSquared}`;
the last field uses the usual `{initial, lower, upper, fixed}` controls.
Daughter masses are fixed. At least two channels are required. Coupling-squared
bounds must be nonnegative. Contributions named alike share the mass and channel
couplings across wave projections, reflectivities and k sectors; their line-shape
choice and ordered channel definitions must agree. Production magnitude and phase
remain specific to each wave. BW and Flatté contributions may coexist coherently.
The `.width` member is unused for Flatté and produces no fitted coordinate.

The precise implemented convention, with s=m² and all masses in GeV, is

\[
F(m)=\frac{1\;\mathrm{GeV}^2}{m_0^2-s-i\sum_i g_i^2\rho_i(m)},\qquad
\rho_i(m)=\sqrt{1-\frac{(m_{i1}+m_{i2})^2}{s}}
            \sqrt{1-\frac{(m_{i1}-m_{i2})^2}{s}}=\frac{2q_i}{m}.
\]

The two square roots use the physical-sheet continuation: closed channels have
positive imaginary phase space between threshold and pseudothreshold, yielding
a real contribution to the denominator. They are retained below threshold.
The fitted channel parameters are g_i² in GeV², not g_i in GeV or m₀g_i; published
couplings in other conventions need conversion. The fixed numerator makes F
dimensionless; `couplingMagnitude` sets its production scale. `poleMass` denotes
the Flatté bare mass parameter here, not a solved complex pole position.
This is S-wave channel phase space; it includes no centrifugal-barrier factors.
For resonance and threshold conventions see the
[PDG resonance review](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-resonances.pdf).

The existing BW choice remains
`MassLineShape::ConstantWidthBreitWigner`, with
`m0*Gamma0/(m0*m0-m*m-i*m0*Gamma0)`. It has a relativistic denominator but
constant width: it contains no breakup-momentum dependence or barrier factors.
The hard-coded PhotoTest generator retains its existing BW prescription; to
select Flatté, use a custom mass-model table as above and omit `--photo-test`,
which replaces that table with the default BW model.

For the generic a0 example above, outputs contain `example_a0_mass`,
`example_a0_eta_pi_g2`, `example_a0_KK_g2` and per-wave production couplings. Channel couplings and mass can be fixed or fitted.
Analytic derivatives cover both the Flatté mass and channel couplings; the
`--numerical-gradients` option also works.
