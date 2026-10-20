# Flatté S-wave and Breit–Wigner D-wave test

This small closure study generates f0(980) in the S wave and f2(1270) in the
D wave, then fits their masses, decay parameters and production couplings.
The shapes use published PDG/BES values; production coefficients are toy values.
The [model notes](model.md) give the sources, equations and coupling convention.

Run from the repository root with a ROOT-enabled Python:

```bash
python -m tutorials.photo_test_m.flatte.run
python -m tutorials.photo_test_m.flatte.analyse
```

`settings.h` contains the model and defaults: 20 starts, 10 workers, both
reflectivities, effective k=0, no background or smearing. Starts are displaced
from the generated parameters. Each named resonance shares its shape across
projections and reflectivities; each sector has an independent phase reference.
Masses, the BW width and Flatté decay coefficients use positive free coordinates.
Override the run size with `run --starts N --workers N`.

The runner writes `output/moments.root` and `output/fit.root`. Analysis checks
shape recovery, reference phases and an independent intensity reconstruction,
then saves `output/verification_report.txt` and figures under `output/plots/`:

- masses, BW width, Flatté decay coefficients and production magnitudes/phases;
- total intensity and residuals, plus separate S/D contributions;
- `lineshapes.pdf`, showing shape moduli squared and phases;
- `moments/`, with one file per (L,M): two alpha panels for M=0, four otherwise.

The 82 generated points include the charged-kaon threshold. Model curves combine
2001 dense evaluations with original-grid predictions to retain that exact cusp.
Plots contain no uncertainty bars. Assigned errors only define chi-square weights.
The Flatté has no constant fitted width; the width plot contains only the BW.
