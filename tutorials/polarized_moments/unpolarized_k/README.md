# Unpolarized k-sector test

This tutorial generates polarized amplitudes containing `k=+1` and `k=-1`
components, projects them into unpolarized moments, and fits those moments with
a single unpolarized amplitude per wave. The reference magnitude is

```text
A_quad = sqrt(|A_+|^2 + |A_-|^2).
```

Select `KPositiveDominant` or `KReflectivitySplit` in `FixedMoments()` inside
`app/UserSettings.h`. The suppression values and repeat count are controlled
at the top of `run.py`.

The default base seed is nonzero. Each suppression/repeat pair receives a
deterministic offset, so deleting and regenerating one point reproduces it.

```bash
python -m tutorials.polarized_moments.unpolarized_k.run
python -m tutorials.polarized_moments.unpolarized_k.analyse
```

The analysis reports mean generated and fitted magnitudes, plots the signed
error `|A_fit|-A_quad`, and uses up/down triangles for the two generated k
sectors. Set `REPEAT_TO_PLOT` and `SUPPRESSION_INDEX_TO_PLOT` to select the
example shown in the complex-amplitude plots.

Outputs are divided between `output/generated`, `output/root`, and
`output/plots`.
