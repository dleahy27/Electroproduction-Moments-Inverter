# PhotoTest mass and k-asymmetry scan

This tutorial studies a deterministic photoproduction model containing S-, P-,
and D-wave Breit-Wigner amplitudes. A signed asymmetry `s` splits each wave
between the two nucleon-helicity sectors: `k=+1` has factor `(1+s)/2` and
`k=-1` has factor `(1-s)/2`. Thus `s=+1` selects `k=+1`, `s=-1` selects
`k=-1`, and `s=0` gives an equal split. Every `k=-1` amplitude also has the
same fixed global phase offset, so the sectors are not coherent.

`run.py` scans linearly spaced invariant masses and the values in `ASYMMETRIES`. It
generates a two-k truth point and fits its unpolarized moments at every scan
point. Edit `MASS_BINS`, the mass range, `ASYMMETRIES`, number of starts, and workers
near the top of the file. The local `settings.h` owns the PhotoTest S/P/D model
and generation mode.

The nonzero base seed plus the asymmetry and mass indices uniquely determine every
fit seed recorded in the manifest.

```bash
python -m tutorials.polarized_moments.photo_test_mass_scan.run
python -m tutorials.polarized_moments.photo_test_mass_scan.analyse
```

`OVERWRITE=False` resumes a partial scan by reusing existing truth and fit
files. Set it to `True` only after changing the model or scan settings in a way
that invalidates those files.

The analysis writes one mass plot for every moment, amplitude magnitude, and
amplitude phase. It also writes a signed-error figure for every mass using a
2-by-4 layout. Columns 1-2 contain the linear differences
`|A_fit|-A_quad`; columns 3-4 contain squared-amplitude differences
`|A_fit|^2-A_quad^2`. The top row resolves individual waves into natural and
unnatural reflectivity, and the bottom row shows the signed sum over all waves.
Those sum axes autoscale independently so physically large variations remain
visible.

Colour identifies the k-sector asymmetry `s`; stars and dots identify generated/fitted
moments; up-triangles, down-triangles, and dots identify generated `k=+1`,
generated `k=-1`, and fitted unpolarized amplitudes. Small horizontal offsets
separate coincident mass points. `scan_points.csv` is the authoritative map
between masses, asymmetries, seeds, and ROOT files.

The directory layout is:

```text
output/generated/       generated truth ROOT files
output/root/            fitted ROOT files and scan manifest
output/plots/moments/   input/output moment comparisons
output/plots/amplitudes/
output/plots/signed_errors/
```
