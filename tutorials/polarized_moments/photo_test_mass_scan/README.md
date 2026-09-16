# PhotoTest mass and k-scale scan

This tutorial studies a deterministic photoproduction model containing S-, P-,
and D-wave Breit-Wigner amplitudes. The model specifies the `k=+1` amplitudes;
each `k=-1` complex amplitude is the corresponding positive-k amplitude times
`kMinusScale`.

`run.py` scans linearly spaced invariant masses and the values in `SCALES`. It
generates a two-k truth point and fits its unpolarized moments at every scan
point. Edit `MASS_BINS`, the mass range, `SCALES`, number of starts, and workers
near the top of the file.

The nonzero base seed plus the scale and mass indices uniquely determine every
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

Colour identifies `kMinusScale`; stars and dots identify generated/fitted
moments; up-triangles, down-triangles, and dots identify generated `k=+1`,
generated `k=-1`, and fitted unpolarized amplitudes. Small horizontal offsets
separate coincident mass points. `scan_points.csv` is the authoritative map
between masses, scales, seeds, and ROOT files.

The directory layout is:

```text
output/generated/       generated truth ROOT files
output/root/            fitted ROOT files and scan manifest
output/plots/moments/   input/output moment comparisons
output/plots/amplitudes/
output/plots/signed_errors/
```
