# Double-polarization fixed closure

This case retains both initial and recoil nucleon-polarization information. The
generated and fitted trees therefore contain explicit `k=+1` and `k=-1`
amplitudes without the continuous ambiguity present in a single-polarization
measurement.

Edit the numerical settings at the top of `run.py`, then execute:

```bash
python -m tutorials.polarized_moments.fixed_double.run
python -m tutorials.polarized_moments.fixed_double.analyse
```

The analysis selects the accepted minimum with the lowest chi square and plots
the generated and fitted complex amplitudes. Files are grouped beneath
`output/generated`, `output/root`, and `output/plots`.
