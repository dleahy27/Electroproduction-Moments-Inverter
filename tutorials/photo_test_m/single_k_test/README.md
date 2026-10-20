# Effective single-k PhotoTest study

The generator retains full k=+/-1 truth. The fit uses its unpolarized moments
with effective k=0 amplitudes. General two-k truth has no exact single-k
representation; effective couplings need not equal either injected sector.
At s=+1 the supplied unsplit k=0 model gives an exact truth start.

Edit `ASYMMETRIES` in `run.py` (default: s=0 and s=1) and `settings.h` for
physics and job size. Defaults are **10,000 minimizations per asymmetry,
10 workers, and 20 generated mass points**. Each minimization fits all mass
points simultaneously:

```bash
python -m tutorials.photo_test_m.single_k_test.run
python -m tutorials.photo_test_m.single_k_test.analyse
```

## Phase references and shape starts

The single-k fit explicitly fixes independent a and b reference phases:

```cpp
config.model.SetPhaseGauge("a2_1320", {'a', 'T', 2, 2, 0})
            .SetPhaseGauge("a2_1320", {'b', 'T', 2, 2, 0});
```

The corresponding `a2_1320_a_2_2_phase` and `a2_1320_b_2_2_phase` branches
stay at their configured initial value, zero here, for every random start.
The generated k=-1 phases are separate injected quantities; they need not
match the effective k=0 reference. Each named resonance still shares its mass
and width across projections and reflectivities.

The old resonance-specific mass/width windows have been removed:

```cpp
for (auto& resonance : config.model.resonances) {
    resonance.poleMass.Unbounded();
    resonance.width.Unbounded();
}
```

`Unbounded()` keeps the configured initial value and fixed flag and ignores the
lower/upper fields. Positive masses and widths use free log coordinates, with
no MINUIT upper/lower limits. Positive successive tensor mass gaps retain
canonical ordering. Positivity and ordering remain physical conventions; no
narrow mass or width windows remain. Initial masses and widths come from the
PhotoTest model, and each random start begins near those values. Subsequent
starts randomize free production couplings, preserving both phase references.
See [ordered masses](../ordered_masses.md) for mappings and covariance output.

Generation uses `--photo-test` for the two-k asymmetry scan. Fitting uses the
model from `settings.h` without `--photo-test`, preserving custom references,
starts and fixed flags. The assigned moment errors are chi-square weights for
this exact numerical test, not simulated statistical uncertainty.

## Outputs and plots

`output/scan.root` contains moments, fit results, injected truth, predictions,
scan coordinates and densely evaluated model curves. Analysis reads this file
alone. Individual ROOT files remain available for inspection.

The existing mass/width and coupling figures remain. Their bars show the range
among the best 5% of accepted starts, rather than measurement uncertainties.
The intensity figure now uses a smooth curve from the best-fit resonance
parameters, with residuals evaluated at the generated mass points.

`plots/moments/` contains one PDF per `(L,M)` per asymmetry:

- `M=0`: two panels, alpha=0 and alpha=1;
- `M>0`: four panels, alpha=0,1,2,3.

Generated values are markers; best-fit model predictions are continuous lines.
There are no moment error bars or uncertainty bands. `massFitCurves` stores
1,001 masses between the first and last fitted mass, evaluated directly from
the selected best resonance parameters. These lines are not interpolations
between the 20 generated points. `massFitPredictions` retains model values at
the original mass points. Set `config.curvePoints` in `MassFit()` to adjust
the dense-grid count.

For a separate short check without changing defaults:

```bash
python -m tutorials.photo_test_m.single_k_test.run --workers 10 --starts 20 --output /tmp/photo_single_check
python -m tutorials.photo_test_m.single_k_test.analyse --output /tmp/photo_single_check
```

`--settings`, `--asymmetries` and `--max-calls` override the steering card,
scan coordinates and per-start function-call budget. The runner never rewrites
settings headers. Re-run generation/fitting with the updated binary before
using the new analysis on older files that lack `massFitCurves`.
