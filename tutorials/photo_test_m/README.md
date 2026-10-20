# Mass-dependent fits

These studies generate exact moments and fit all mass points simultaneously.
The assigned moment errors set chi-square weights; they do not add noise.
Each case has a settings header, a runner and an analysis script.

| Case | What it checks |
| --- | --- |
| [closure](closure/README.md) | One BW, with displaced mass, width and magnitude starts |
| [full_closure](full_closure/README.md) | Full four-resonance PhotoTest model with explicit k=±1, starting at truth |
| [single_k_test](single_k_test/README.md) | Effective k=0 fits to two-k truth at s=0 and s=1 |
| [flatte](flatte/README.md) | S-wave f0(980) Flatté plus D-wave f2(1270) BW, starting away from truth |

From the repository root, after building EMI:

```bash
python -m tutorials.photo_test_m.closure.run
python -m tutorials.photo_test_m.closure.analyse --check
python -m tutorials.photo_test_m.full_closure.run
python -m tutorials.photo_test_m.full_closure.analyse --check
python -m tutorials.photo_test_m.flatte.run
python -m tutorials.photo_test_m.flatte.analyse
```

The single-k defaults are a larger study: 10,000 starts per asymmetry on 10
workers. For a short run that preserves the main output:

```bash
python -m tutorials.photo_test_m.single_k_test.run \
  --starts 20 --workers 10 --output /tmp/photo_single_check
python -m tutorials.photo_test_m.single_k_test.analyse --output /tmp/photo_single_check
```

Each runner leaves ROOT files and plots in its ignored `output/` directory.
`analysis.py` provides shared readers and parameter labels. `verify.py` checks
individual fit files without rerunning minimizations:

```bash
python -m tutorials.photo_test_m.verify \
  tutorials/photo_test_m/closure/output/fit.root --starts 8 --workers 2 \
  --require-covariance
```

It checks start/worker indices, fixed phase references, positive shapes, mass
ordering, finite predictions and the full covariance map when available. Closure
analyses also test injected parameter and moment recovery. The single-k study
can be a model mismatch, so its effective production couplings are not checked
against either generated sector.

## Model settings

Named resonances share shape parameters across their contributions; production
couplings remain specific to each wave. Start zero uses the supplied values.
Later starts randomize free production couplings while retaining configured
shape starts and phase references. `--max-calls` sets the per-start call budget.
HESSE is optional; failed starts remain available in the result tree.

For effective k=0, each populated reflectivity has an independent reference
phase. Explicit-k fits have one overall reference. Override a reference with
`MassModelConfig::SetPhaseGauge(name, wave)`; its `couplingPhase.initial` is
held fixed. `Unbounded()` removes finite shape windows while retaining physical
positivity. See [ordered masses](ordered_masses.md) for the mapping, constraints
and covariance output, and [Flatté model notes](flatte/model.md) for conventions.

`generate-mass-model --polarization none` supports configured k=0 models.
`MakePhotoTestMassModel(0.0, false)` supplies the unsplit PhotoTest table.
`--photo-test` rebuilds that table, so omit it when generating or fitting custom
resonances or preserving custom references. Generation accepts `--mass-bins`,
`--mass-min` and `--mass-max`; the PhotoTest scan also accepts `--k-asymmetry`.

Best-fit markers use accepted finite minima. Where parameter plots include
ranges, the case README states the selected fraction: the single-k plots use
5%, and the BW closures use 0.1%. These show variation across starts, not
statistical uncertainty. Moment overlays in the single-k and Flatté studies
show generated points and dense model curves without uncertainty bars.
