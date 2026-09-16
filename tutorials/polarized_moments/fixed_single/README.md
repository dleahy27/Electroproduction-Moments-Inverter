# Single-polarization fixed closure

This tutorial generates and fits separate initial-polarization and
recoil-polarization examples. In each case only one nucleon polarization is
observed, so the fit uses the gauge convention implemented by EMI to remove the
unobservable common spin-basis rotation.

Edit `EPSILON`, `STARTS`, `WORKERS`, `SEED`, or `K_GAUGE` at the top of
`run.py`, then run:

```bash
python -m tutorials.polarized_moments.fixed_single.run
```

The runner creates both the `initial` and `recoil` examples. To analyse one,
set `POLARIZATION` near the top of `analyse.py` to `"initial"` or `"recoil"`:

```bash
python -m tutorials.polarized_moments.fixed_single.analyse
```

Generated moments are stored in `output/generated`, fit trees in
`output/root`, and separate initial/recoil figures in `output/plots`.
