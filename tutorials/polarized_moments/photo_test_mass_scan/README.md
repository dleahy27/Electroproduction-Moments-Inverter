# PhotoTest mass and k-asymmetry scan

This tutorial studies a deterministic photoproduction model containing S-, P-,
and D-wave Breit-Wigner amplitudes. A signed asymmetry `s` splits each wave
between the two nucleon-helicity sectors: `k=+1` has factor `(1+s)/2` and
`k=-1` has factor `(1-s)/2`. Thus `s=+1` selects `k=+1`, `s=-1` selects
`k=-1`, and `s=0` gives an equal split. Every `k=-1` partial wave also has its
own fixed phase offset, shared by both reflectivities and unchanged for every
asymmetry setting.

`run.py` scans linearly spaced invariant masses and the values in
`ASYMMETRIES`. It generates a two-k truth point and fits its unpolarized
moments at every scan point. Edit `MASS_BINS`, the mass range,
`ASYMMETRIES`, number of starts, and workers near the top of the file. The
local `settings.h` owns the PhotoTest S/P/D model and generation mode.

Set `NOISE` in `run.py` to a relative generated-moment uncertainty, for
example `0.05`, or leave it as `None` for exact moment truth. The runner passes
the selected value through `--noise` and records it in the scan manifest.

The nonzero base seed is recorded in the manifest and reused for every fit, so
each scan point uses the same minimizer-start sequence.

```bash
python -m tutorials.polarized_moments.photo_test_mass_scan.run
python -m tutorials.polarized_moments.photo_test_mass_scan.analyse
python -m tutorials.polarized_moments.photo_test_mass_scan.analyse_unnormalised
python -m tutorials.polarized_moments.photo_test_mass_scan.analyse_ideal_cases
python -m tutorials.polarized_moments.photo_test_mass_scan.analyse_asymmetry_intensity_comparisons
```

`OVERWRITE=False` resumes a partial scan by reusing existing truth and fit
files. Set it to `True` only after changing the model or scan settings in a way
that invalidates those files.

`analyse.py` writes amplitudes and moments in the fit-normalized convention.
`analyse_unnormalised.py` restores the unnormalised model scale and writes the
corresponding amplitude and moment plots. Both analyses write a mass plot for
every moment, amplitude magnitude, and amplitude phase. They also write a
signed-error figure for every mass using a 2-by-4 layout. Columns 1-2 contain
the linear differences
`|A_fit|-A_quad`; columns 3-4 contain squared-amplitude differences
`|A_fit|^2-A_quad^2`. The top row resolves individual waves into natural and
unnatural reflectivity, and the bottom row shows the signed sum over all waves.
Those sum axes autoscale independently so physically large variations remain
visible. When Hessian correlations are present, it also writes one correlation
matrix for each completed mass and asymmetry point.

`analyse_ideal_cases.py` compares the unnormalised results at the ideal
settings `s=0` and `s=1`. For every wave it writes one magnitude figure and
one phase figure, with filled stars for the reconstructed `s=1` amplitude and
dots for `s=0`; both use the same colour. It also writes one
invariant-mass total/S/P/D unnormalised intensity plot and one two-panel
unnormalised moment plot for every common observable. The moment plots put
`s=0` on the left and `s=1` on the right, with stars for set values and dots
for reconstructed values. It also writes the mass-dependent sum of
`|A_fit(s=0)|^2-|A_fit(s=1)|^2` over all waves. The latter figure annotates
the same difference summed over every mass bin and includes a second panel
showing the ratio of the summed `s=0` and `s=1` fitted intensities.

`analyse_asymmetry_intensity_comparisons.py` extends this diagnostic to every
non-reference scan value of `s`. It writes one two-panel figure per value,
comparing its summed unnormalised fitted intensity with the `s=1` result by
difference and ratio.

Colour identifies the k-sector asymmetry `s`; stars and dots identify
generated and fitted moments; up-triangles, down-triangles, and dots identify
generated `k=+1`, generated `k=-1`, and fitted unpolarized amplitudes. Small
horizontal offsets separate coincident mass points. `scan_points.csv` is the
map between masses, asymmetries, seeds, and ROOT files.

Each generated PhotoTest ROOT tree also contains the model amplitudes before
the fit gauge rotation and the enforced `H_0000=2` normalization. The unnormalised ROOT branches use `raw_*` names: amplitude branches are `raw_a_*` and `raw_b_*`, and the tree
contains `raw_photo_intensity` and its S-, P-, and D-wave contributions
(`raw_photo_intensity_l0`,
`_l1`, and `_l2`). These intensities are sums of
amplitude moduli squared in the toy model's arbitrary units. They are useful
for displaying resonance structure, but are not physical differential cross
sections until flux, phase-space, and any desired detector factors are added.
The tree also records `raw_photo_H_0_0_0` and
`normalised_photo_H_0_0_0`. `analyse_unnormalised.py` derives the per-point
amplitude conversion factor
`sqrt(raw_photo_H_0_0_0 / normalised_photo_H_0_0_0)` and applies it to the
minimizer magnitudes and generated amplitudes. The same ratio without the
square root scales generated and fitted moments. The ideal-case analysis writes
`unnormalised_mass_intensity.pdf`, including the total and the three
orbital-wave contributions for both ideal asymmetries.

The directory layout is:

```text
output/generated/       generated truth ROOT files
output/root/            fitted ROOT files and scan manifest
output/plots/normalised/      fit-normalized analysis figures
output/plots/unnormalised/    unnormalised analysis figures
output/plots/ideal_cases_unnormalised/  s=0 and s=1 unnormalised comparisons
output/plots/asymmetry_intensity_comparisons/  each s compared with s=1
```

After fitting, the runner merges the point trees into `output/root/scan.root`.
All analyses load that file through `scan_data.py`. RDataFrame selects the
best 0.1% of accepted fits before reading their amplitudes and moments.
Correlations are plotted for the best fit.
Point files remain available for resuming and refitting. To merge an existing
scan without running new fits:

```bash
python -m tutorials.polarized_moments.photo_test_mass_scan.run --merge-only
```

Fitted coordinate/amplitude markers remain at the lowest-chi2 accepted fit.
Their bars span the minimum and maximum among the lowest-chi2 0.1% of
converged, finite fits, selected separately for each scan case (rounded up,
with at least one fit). The bars show the spread across minimizations, not statistical errors. Amplitude phases use the displayed principal
interval; a range crossing its boundary can therefore span almost 2*pi.
