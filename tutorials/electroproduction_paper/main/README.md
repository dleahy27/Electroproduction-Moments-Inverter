# Main electroproduction analysis

This directory contains the experimental fits and plots used for the main
electroproduction study. The default channels are electron and muon production
of rho and omega mesons. Input moment tables live in the parent `input/`
directory; all results from this tutorial are written locally.

## Directory contents

- `run_fit.py` and `run_bootstrap.py` run one editable data bin.
- `run_all_fits.py` and `run_all_bootstraps.py` run every configured channel.
- `run_fixed_amplitudes.py` generates and fits the fixed photoproduction
  closure point configured in `app/UserSettings.h`.
- `analyse_experimental_hessian.py` and
  `analyse_experimental_bootstrap.py` plot the fitted amplitudes.
- `analyse_hessian_moments.py` and `analyse_bootstrap_moments.py` compare fitted
  moments with the experimental inputs.
- `analyse_R_hessian.py` and `analyse_R_bootstrap.py` compare
  `R = sigma_L/sigma_T` with the published measurements.
- `analyse_fixed_amplitudes.py` compares the fixed truth with its unpolarized
  fit in amplitude and moment space.

## Run one bin

Edit the settings at the top of the two runners, then execute:

```bash
python -m tutorials.electroproduction_paper.main.run_fit
python -m tutorials.electroproduction_paper.main.run_bootstrap
```

## Run all configured experimental bins

```bash
python -m tutorials.electroproduction_paper.main.run_all_fits
python -m tutorials.electroproduction_paper.main.run_all_bootstraps
```

The fit cost comes from `app/UserSettings.h` unless the optional values near
the top of a runner are set explicitly.

Each experimental bin is fitted with many random starting points. The Hessian
analyses select the accepted entry with minimum chi square. Bootstrap analyses
use the sample distribution: ordinary means and standard deviations describe
magnitudes, while circular statistics describe phases across the `-pi/pi`
boundary.

## Analyse the experimental results

```bash
python -m tutorials.electroproduction_paper.main.analyse_experimental_hessian
python -m tutorials.electroproduction_paper.main.analyse_experimental_bootstrap
python -m tutorials.electroproduction_paper.main.analyse_hessian_moments
python -m tutorials.electroproduction_paper.main.analyse_bootstrap_moments
python -m tutorials.electroproduction_paper.main.analyse_R_hessian
python -m tutorials.electroproduction_paper.main.analyse_R_bootstrap
```

## Fixed-amplitude closure

Set `FixedMoments()` in `app/UserSettings.h`, then run:

```bash
python -m tutorials.electroproduction_paper.main.run_fixed_amplitudes
python -m tutorials.electroproduction_paper.main.analyse_fixed_amplitudes
```

The output layout is:

```text
output/generated/   generated fixed truth
output/root/        Hessian and bootstrap fit results
output/plots/       all figures, grouped by analysis
```

Delete a result or set the runner's overwrite option when changing the wave
model. Reusing files produced with a different `UserSettings.h` can otherwise
mix incompatible branch layouts.
