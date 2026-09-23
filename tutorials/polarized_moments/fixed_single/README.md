# Single-polarization fixed closure

This tutorial generates and fits separate initial-polarization and
recoil-polarization examples. In each case only one nucleon polarization is
observed, so the fit uses the gauge convention implemented by EMI to remove the
unobservable common spin-basis rotation.

Edit `EPSILON`, `STARTS`, `WORKERS`, `SEED`, or `K_GAUGE` at the top of
`run.py`, and edit the local `settings.h` to change the wave model, then run:

```bash
python -m tutorials.polarized_moments.fixed_single.run
```

The runner creates both the `initial` and `recoil` examples and passes the same
polarization and optional k-mixing gauge to generation and fitting. The analysis
processes both cases:

```bash
python -m tutorials.polarized_moments.fixed_single.analyse
```

Generated moments are stored in `output/generated`, fit trees in
`output/root`, and separate initial/recoil figures in `output/plots`.

For each case, the report shows the lowest accepted chi-square solution and the
solution closest to the known amplitudes among numerically equivalent minima.
The plots use the latter because this is a synthetic closure test. A very small
moment residual together with different amplitudes indicates an amplitude
ambiguity, rather than a failure of the fitted observables or gauge application.
