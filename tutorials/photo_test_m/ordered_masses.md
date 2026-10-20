# Canonical resonance mass ordering

Mass-model generation and fitting assign resonance identities in ascending order
of configured starting mass, once during setup. Complete named components retain
their widths, decay channels, production couplings and metadata. Evaluation never
sorts fitted masses and uses no ordering penalty.

The default ordering sector is the partial-wave angular momentum `wave.l`.
For example, the two PhotoTest tensor components are ordered together; the
scalar and vector components are independent. BW and Flatté components may
belong to the same sector. A Flatté mass is its configured bare parameter m0,
not a peak location or a solved complex pole.

The schema does not specify every intrinsic quantum number. To distinguish
unrelated sectors with the same l, set explicit groups on their resonance records:

```cpp
resonance.massOrderingGroup = "tensor_2pp";
config.model.minimumMassGapGeV = 1.e-6; // In MassFit(), separation in GeV.
// In MassGeneration(), set config.massModel.minimumMassGapGeV instead.
```

A shared resonance is identified by its existing `name`: it has one mass and one
BW width or set of Flatté channel couplings across all helicities, orientations,
reflectivities and k sectors. Shared records must agree on l, group and shape
controls. An explicit group cannot mix different l values. Group identifiers use
letters, numbers and underscores. Equal starting masses are rejected with a clear
diagnostic; starting gaps below the configured positive minimum are also rejected.
The default gap is a numerical separation convention, not a physical spacing prior.
Single-resonance sectors retain the original mass coordinate.

## Bounds and fixed masses

Physical masses are represented by a base and positive successive gaps,
`m_i = m_base + sum(gaps before i)`. Original individual mass bounds and fixed
masses are preserved exactly, intersected with the ordering domain. There are
no independent rectangular bounds on the gaps that could silently reinterpret
those physical constraints.

MINUIT supports rectangular variable bounds. For arbitrary bounded masses,
independently fitting every gap would also require coupled constraints. The
implementation therefore fits the physical base mass and uses directly fitted
gaps when the preceding mass is fixed. Otherwise it uses a bounded gap fraction
`t_i` in [0,1], reconstructing the physical positive gap on evaluation.

Writing the original physical interval as [L_i, H_i], and a fixed mass as
[L_i,H_i]=[m_i,m_i], the feasible upper envelope is

```
U_i = min over k >= i of (H_k - (k-i)*minimum_gap)
m_base in [L_0, U_0]
a_i = max(L_i, m_previous + minimum_gap)
m_i = a_i + t_i*(U_i-a_i)
gap_i = m_i-m_previous
```

This covers the full allowed ordered domain, including overlapping or disjoint
mass windows and fixed masses at any rank. A fixed anchor restricts preceding
masses through the upper envelope and following masses through the positive gap.
Infeasible bounds are rejected. The map is continuous and piecewise differentiable
at changes of the active bound; the analytic Jacobian follows that map. As with
ordinary bounded MINUIT fits, covariance near an active constraint deserves care.

Mass and width starts remain the values specified by the steering card. The
existing random-start policy randomizes free production couplings after start
zero; each start uses the same ordered shape initialization. Bounds, fixed flags,
reference phases and amplitude normalization are preserved. Generation, analytic
and numerical fitting, and best-result loading use the same mapping. Repeated
mass fits, including external mass-toy loops, rebuild it identically. The existing
`RunBootstrap` is a bin-wise amplitude bootstrap; there is no separate built-in
mass-model bootstrap or mass-fit restart-file API to remap.

## Unbounded positive shape parameters

Call `resonance.poleMass.Unbounded()` or `resonance.width.Unbounded()` to remove
both configured fit limits while retaining the initial value and fixed flag.
Bounded configurations keep their previous exact-bound behavior.

A standalone unbounded positive shape uses `p = exp(q)`, with q registered as
a free MINUIT variable. An unbounded ordered group uses a free log
base and free log gaps, `m_base = exp(q_base)` and
`gap_i = minimum_gap + exp(q_i)` when no other physical mass bound is active.
Thus positivity and label ordering remain, without finite mass/width windows.
Mixed bounded/fixed groups continue to respect their exact feasible domain.

`positive__<physical_name>` stores free log coordinates for individual shapes;
`ordered__<group>__log_base_mass` and `ordered__<group>__gap_<name>_log` store
ordered log coordinates. The full analytic Jacobian includes the exponential
derivatives, so physical mass and width uncertainties still use `J V J^T`.
The physical output branch names remain unchanged. Analysis excludes both
`positive__` and `ordered__` diagnostics from physical parameter plots.

## ROOT output and uncertainty propagation

Existing physical branches such as `a2_1320_mass`, `a2_1700_mass` and
`err__a2_1700_mass` retain their names. Analysis helpers select these physical
branches and exclude the new `ordered__` diagnostics.

Additional scalar branches include:

- `ordered__l2__base_mass`: fitted base coordinate, when floating;
- `ordered__l2__gap_a2_1700_fraction`: fitted bounded gap coordinate, when needed;
- `ordered__l2__base_GeV` and `ordered__l2__gap_a2_1700_GeV`: reconstructed physical
  base and gap, with corresponding `err__` branches.

With successful HESSE covariance extraction, each row also stores:

- `fit_covariance`: full row-major covariance of MINUIT coordinates;
- `physical_covariance`: full row-major `J V_fit J^T`, including every mass/width/
  channel/production cross-correlation;
- `mass_jacobian`: full row-major physical-by-fit Jacobian;
- `covariance_ok`: whether those arrays were obtained for this start.

The `mass_fit_coordinate_names` and `mass_fit_physical_coordinate_names`
TObjStrings specify matrix row/column order, one name per line. The
`mass_fit_ordering_map` TObjString records `group resonance physical_index
previous_physical_index fit_index kind`, with kind 0=fixed, 1=base, 2=direct gap,
3=gap fraction, 4=free log base, 5=free log gap. Missing fit or preceding
indices are -1.
`mass_fit_minimum_gap_GeV` records the separation floor. Without HESSE, covariance
arrays are empty and uncertainty branches remain NaN. Fixed-mass covariance rows
are zero when a covariance is available. Original physical truth parameters are
retained in generated ROOT files alongside transformed diagnostics.

## Check a run

Run the closure analyses and the ROOT-output checker described in the
[mass-fit index](README.md). For example, after the single-k scan:

```bash
python -m tutorials.photo_test_m.verify \
  tutorials/photo_test_m/single_k_test/output/s_00_truth_fit.root \
  tutorials/photo_test_m/single_k_test/output/s_01_truth_fit.root \
  --starts 10000 --workers 10
```

The checker verifies saved start counts, reference phases, positive shapes,
ordered masses, dense predictions and the full covariance transformation when
available. Use `--require-covariance` for a run made with HESSE enabled. Failed
starts remain in the file and are counted separately.

Ordering assigns canonical resonance labels by mass. Other local minima, phase
freedoms and underdetermined production couplings can remain; it does not
establish a globally unique fit or physically unique couplings.
