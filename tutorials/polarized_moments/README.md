# Polarized moments analysis

These tutorials exercise EMI with explicit nucleon spin-flip and non-flip
amplitudes. Each case keeps generated moments, fitted ROOT trees, and plots in
its own `output` directory.

Build EMI before running a case:

```bash
python build.py --test
```

The cases are:

- [Single polarization](fixed_single/README.md): initial-only and recoil-only
  fixed-amplitude closures, including the remaining spin-basis gauge choice.
- [Double polarization](fixed_double/README.md): simultaneous initial and
  recoil information with both nucleon-polarization indices resolved.
- [Unpolarized k test](unpolarized_k/README.md): fit the unpolarized projection
  of generated two-k amplitudes and study quadrature recovery.
- [PhotoTest mass scan](photo_test_mass_scan/README.md): deterministic S/P/D
  photoproduction amplitudes across invariant mass and k-minus scale.

Run scripts as modules from the repository root, for example:

```bash
python -m tutorials.polarized_moments.fixed_double.run
```

These are closure tests: moments are generated from known complex amplitudes
and then refitted. Agreement must be judged modulo the explicitly fixed phase
or spin-basis gauge. Each runner records its seed and keeps truth and fit ROOT
files separate beneath its local `output/` directory.
