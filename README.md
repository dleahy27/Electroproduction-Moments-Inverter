# Electroproduction-Moments-Inverter — reconstructed README

## What this repository does

This project takes **published spin-density matrix elements (SDMEs)** for exclusive vector-meson production and turns them into a set of **moment observables** \(the `RH...` and `RH04...` quantities\), then **inverts** those moments back to complex production amplitudes by fitting a partial-wave model.

In practical terms, the repository supports three related tasks:

1. **Build ROOT input tables** of experimental moments from hard-coded SDME measurements.
2. **Fit amplitudes** to those moments with many randomized minimization starts.
3. **Study uncertainties and ambiguities** through bootstrap toys and synthetic “known-truth” datasets.

The code is written as a set of **ROOT macros** rather than as a compiled package with a build system.

---

## Physics background

### The underlying problem

For exclusive production of spin-1 vector mesons such as \(\rho^0\) and \(\omega\), the angular distribution of the decay products contains information about the helicity structure of the production amplitude. Experimental papers often report this information in terms of **SDMEs** rather than directly in terms of amplitudes.

This repository sits between those two descriptions:

- **Published SDMEs** \(\rightarrow\) converted into **moments**
- **Moments** \(\rightarrow\) fitted with a model written in terms of **complex amplitudes**

That inversion is useful because the amplitudes are much closer to the underlying reaction mechanism than the SDMEs are.

### What the moments are doing here

The fit is organized around moment labels like

- `RH04_L_M`
- `RH_alpha_L_M`
- internally reconstructed `H_alpha_L_M`

The code builds moment models from Clebsch–Gordan / Wigner-3j couplings and evaluates them from a chosen set of partial-wave amplitudes. In electroproduction mode, the code reconstructs a longitudinal/transverse ratio

\[
R = -\frac{H^4_{00}}{H^0_{00}}
\]

and uses it to rescale the raw `H` moments into the experimentally observed `RH` and `RH04` combinations.

### Amplitude basis used by the code

The parameter names tell you the basis directly:

- `a_*` and `b_*` = positive / negative reflectivity sectors
- `T` and `L` = transverse / longitudinal photon couplings
- `l,m` indices are encoded in names like `a_T_1_0`, `b_L_1_m1`
- phases are stored separately as `aphi_*`, `bphi_*`

With the default settings the code is effectively set up for **S- and P-wave spin-1 production**:
- \(l_{\max}=1\)
- \(m_{\max}=1\)

So the active wave content is basically the \(S\) and \(P\) sectors.

### Electroproduction vs photoproduction

The repository supports two physics modes:

- **Electroproduction / leptoproduction**
  - full set of moments including longitudinal and LT-interference structures
  - uses the `epsR4` parameter and the reconstructed \(R\)
- **Photoproduction**
  - turns on `photoProduction=true`
  - only keeps the lower moment sector (\(\alpha \le 3\))
  - all longitudinal amplitudes/phases are fixed to zero

---

## Repository structure

### Top level

- `MakeLeptoMoments.C`  
  Converts hard-coded **electro/leptoproduction SDMEs** into ROOT trees of experimental moments.

- `MakePhotoMoments.C`  
  Converts hard-coded **GlueX photoproduction SDMEs** into ROOT trees of experimental moments.

- `RunGivenMoments_Chi2Amps.C`  
  The **core fitter**. Reads an input ROOT moment table, builds the moment model, and runs many randomized minimizations.

- `RunGivenMoments_Chi2Amps_Bootstrap.C`  
  Wraps the fitter in a **toy/bootstrap procedure** by Gaussian-throwing the observed moments and refitting each toy.

- `GenerateMomentsFromFixedAmplitudes.C`  
  Creates a **synthetic input dataset** from a user-chosen set of amplitudes/phases. Useful for closure tests and ambiguity studies.

- `HERMES_analysis.ipynb`  
  A notebook in the repository root. It appears to be an older or lighter analysis notebook relative to the versions in `AnalysisScripts/`.

- `README.md`  
  Currently just a placeholder in the repository.

- `LICENSE`  
  GPL-3.0.

- `.gitignore`

### `AnalysisScripts/`

- `GlueX_analysis.ipynb`
- `HERMES_analysis.ipynb`
- `R_analysis.ipynb`
- `fixed_amp_analysis.ipynb`

These are post-processing / plotting notebooks, not part of the core fitting engine.

---

## Dependencies

### Required C++ / ROOT side

You need a ROOT installation that provides at least:

- `TFile`, `TTree`
- `TRandom3`
- `TBenchmark`
- `ROOT::Math::Minimizer`
- `Minuit2`
- `MathMore` / Wigner-3j support
- `ROOT::TProcessExecutor`
- `TFileMerger`

A working ROOT 6 installation with PyROOT available is the safest assumption.

### Required Python side for notebooks

The notebooks use:

- `numpy`
- `matplotlib`
- `ROOT` (PyROOT)
- `scipy` (`R_analysis.ipynb` uses `curve_fit`)

A minimal environment is something like:

```bash
python -m pip install numpy matplotlib scipy
```

with ROOT/PyROOT already available from your ROOT installation.

---

## Directory layout expected by the macros

The macros assume these folders already exist:

```bash
mkdir -p InputFiles/Experiment
mkdir -p InputFiles/Generated
mkdir -p OutputFiles
```

The code writes to those directories directly and does **not** create them for you.

- `MakeLeptoMoments.C` writes to `./InputFiles/Experiment/`
- `MakePhotoMoments.C` writes to `./InputFiles/Experiment/`
- `GenerateMomentsFromFixedAmplitudes.C` writes to `./InputFiles/Generated/`
- `RunGivenMoments_Chi2Amps.C` writes final fit results to `./OutputFiles/`

---

## End-to-end workflows

## 1. Electroproduction / leptoproduction workflow

### Step 1: build an experimental moment table

For electron \(\rho\):

```bash
root -l -q 'MakeLeptoMoments.C("e_rho")'
```

For muon \(\rho\):

```bash
root -l -q 'MakeLeptoMoments.C("mu_rho")'
```

For electron \(\omega\):

```bash
root -l -q 'MakeLeptoMoments.C("e_omega")'
```

For muon \(\omega\):

```bash
root -l -q 'MakeLeptoMoments.C("mu_omega")'
```

This creates a ROOT file in `InputFiles/Experiment/` containing one tree (default name: `expMoments`) with arrays of moment values and moment uncertainties vs. \(Q^2\).

### Step 2: fit one bin

Example: fit the first \(Q^2\) bin of HERMES \(\rho\)

```bash
root -l -q 'RunGivenMoments_Chi2Amps.C("InputFiles/Experiment/e_rho_moments.root","expMoments",0,"HERMES_results_0.root",1.0,false)'
```

Arguments are:

1. input ROOT file
2. tree name
3. bin index
4. output filename (written under `./OutputFiles/`)
5. `epsR4`
6. `photoProduction`

For electroproduction, use `photoProduction=false`.

### Step 3: loop over bins

For HERMES \(\rho\) (4 bins):

```bash
for i in 0 1 2 3; do
  root -l -q "RunGivenMoments_Chi2Amps.C(\"InputFiles/Experiment/e_rho_moments.root\",\"expMoments\",$i,\"HERMES_results_${i}.root\",1.0,false)"
done
```

Similarly:

- `mu_rho` has 4 bins
- `e_omega` has 3 bins
- `mu_omega` has 3 bins

### Step 4: bootstrap the uncertainties

```bash
root -l -q 'RunGivenMoments_Chi2Amps_Bootstrap.C("InputFiles/Experiment/e_rho_moments.root","expMoments",0,"HERMES_bootstrap_0.root",1.0,false)'
```

This generates toy datasets by fluctuating the observed moments and refitting each one.

---

## 2. Photoproduction workflow

### Step 1: build the GlueX moment table

```bash
root -l -q 'MakePhotoMoments.C("gluex")'
```

This writes a ROOT file under `InputFiles/Experiment/` with 18 \(-\bar t\) bins.

### Step 2: fit one GlueX bin

```bash
root -l -q 'RunGivenMoments_Chi2Amps.C("InputFiles/Experiment/gluex_moments.root","expMoments",0,"gluex_results_0.root",1.0,true)'
```

Important:
- for photoproduction, set the last argument to `true`

### Step 3: bootstrap one GlueX bin

```bash
root -l -q 'RunGivenMoments_Chi2Amps_Bootstrap.C("InputFiles/Experiment/gluex_moments.root","expMoments",0,"gluex_bootstrap_0.root",1.0,true)'
```

### Step 4: loop over all 18 bins

```bash
for i in $(seq 0 17); do
  root -l -q "RunGivenMoments_Chi2Amps.C(\"InputFiles/Experiment/gluex_moments.root\",\"expMoments\",$i,\"gluex_results_${i}.root\",1.0,true)"
done
```

and similarly for bootstrap.

---

## 3. Synthetic / closure-test workflow

This is the best way to test whether the inversion machinery can recover a known amplitude point.

### Step 1: generate fake moments from a fixed amplitude set

With defaults:

```bash
root -l -q 'GenerateMomentsFromFixedAmplitudes.C()'
```

This writes a generated ROOT file to `InputFiles/Generated/`.

### Step 2: fit the generated moments back

```bash
root -l -q 'RunGivenMoments_Chi2Amps.C("InputFiles/Generated/fixed_input_moments.root","genMoments",0,"fixedamps.root",1.0,false)'
```

### Step 3: inspect the output with the fixed-amplitude notebook

The `fixed_amp_analysis.ipynb` notebook is designed for this use case.

---

## File-by-file technical explanation

## `MakeLeptoMoments.C`

### Purpose

This macro hard-codes published SDME tables and converts them into the moment basis expected by the fitter.

### Datasets implemented in code

The file contains hard-coded datasets for:

- `e_rho`
- `mu_rho`
- `e_omega`
- `mu_omega`

with their associated \(Q^2\) bin centers and SDME values.

### Output

It writes a tree (default `expMoments`) containing:

- `Nbins`
- `Q2`
- all relevant electroproduction moment arrays:
  - `RH04_...`
  - `RH_1_...`
  - `RH_2_...`
  - `RH_3_...`
  - `RH_5_...`
  - `RH_6_...`
  - `RH_7_...`
  - `RH_8_...`
- corresponding `_err` branches for each moment
- metadata objects:
  - `dataset_key`
  - `dataset_title`
  - `NbinsMeta`

### Important implementation details

- Statistical and systematic uncertainties are combined in quadrature:
  \[
  \sigma = \sqrt{\sigma_{\rm stat}^2 + \sigma_{\rm syst}^2}
  \]
- Some \(L=1\) moments are filled with small placeholder values/errors via `SetMissingL1`.
- The SDME-to-moment conversion factors are encoded explicitly in `FillBin(...)`.

### Caveats

1. The comments at the top are stale:
   - they mention `MakeExperimentalMomentsTable_leptoproduction.C`
   - the actual file/function is `MakeLeptoMoments.C` / `MakeLeptoMoments(...)`

2. The comments/examples suggest that `"omega"` is a supported key, but the actual dataset dispatcher is written around `e_omega` and `mu_omega`. Use those explicit keys.

3. The fitter later assumes a fixed maximum of 18 bins internally. That matches GlueX, but it means you should only pass meaningful bin indices for the specific dataset:
   - HERMES/COMPASS \(\rho\): 0–3
   - \(\omega\): 0–2

---

## `MakePhotoMoments.C`

### Purpose

This macro creates the photoproduction input table, currently for GlueX \(\rho\) photoproduction.

### Dataset implemented

- `gluex`

with 18 mean-\(-t\) bins:
- 0.107, 0.121, 0.138, …, 0.940

### Output

The default output tree is `expMoments`, and it contains:

- `Nbins`
- `Q2`
- `mtbar`
- `t`
- photoproduction moments:
  - `RH_0_...`
  - `RH_1_...`
  - `RH_2_...`
  - `RH_3_...`
- corresponding `_err` branches
- metadata:
  - `dataset_key`
  - `dataset_title`
  - `NbinsMeta`

### Why does it store `Q2`, `mtbar`, and `t` all together?

For GlueX photoproduction there is no \(Q^2\)-bin structure like in electroproduction, so the file reuses the same mean-\(-t\) value in all three arrays. This appears to be for downstream convenience / compatibility with plotting code that expects a `Q2` branch.

### Important implementation details

- `RH_0_0_0` is hard-wired to `2.0`
- several formally absent moments are forced to zero with a small default error
- one bin (`i == 15`) is manually overridden with hard-coded numbers after the standard conversion

That bin-specific override is worth knowing if you are auditing exact numerical agreement.

---

## `RunGivenMoments_Chi2Amps.C`

### Purpose

This is the central inversion engine.

It:

1. reads an experimental or generated moment table,
2. builds the allowed amplitude parameters,
3. constructs the moment model \(H^\alpha_{LM}\),
4. runs many random minimization starts,
5. writes all fit attempts to a ROOT tree.

### Main internal components

#### `FitConfig`
Controls:
- wave truncation (`lmax`, `mmax`)
- electro vs photo mode
- number of starts
- minimizer settings
- MCMC prescan options
- random seed
- input file / tree / bin
- `epsR4`

#### `BuildObservedMoments(...)`
Reads the observed moments from the input ROOT file for one bin.

#### `BuildAmplitudePhaseParameters(...)`
Creates the full list of magnitude and phase parameters:
- `a_T_*`, `a_L_*`, `b_T_*`, `b_L_*`
- `aphi_*`, `bphi_*`

It also fixes several phases/amplitudes to remove trivial gauge/normalization freedoms and to impose the current model choices.

#### `BuildMomentModels(...)`
Builds the full model of the reconstructed moments `H_alpha_L_M` using Clebsch–Gordan coefficients and cached phase-pair trigonometric factors.

#### `Chi2Function`
Evaluates the fit objective.

#### `RunGivenMoments_Chi2Amps(...)`
The top-level user macro:
- splits starts across cores with `ROOT::TProcessExecutor`
- launches worker fits
- merges `.part_*.root` files into one final output

### How the parameterization works

The code uses a **simplex / logit parameterization** for the free magnitude variables so that the total magnitude norm is constrained. In other words, the fit is not just independently varying every magnitude between 0 and 1; it is enforcing a normalized amplitude set.

That is physically sensible because the moments are sensitive mostly to relative magnitudes/phases and an overall normalization freedom must be fixed somewhere.

### Output

The main fit macro writes a ROOT file under `OutputFiles/` with a tree named:

- `fitResults`

Each entry corresponds to **one random start / one local minimization attempt**.

The tree contains:

- `log_val`
- all amplitude parameters
- all reconstructed model moments `H_alpha_L_M`
- the observed-like `RH...` / `RH04...` values reconstructed from the fitted amplitudes

### Important caveats

1. **The objective is not a conventional error-weighted chi-square.**  
   The code reads uncertainties into `moment.sigma`, but in the main fit objective it minimizes plain sums of squared residuals:
   \[
   \chi^2_{\text{code}} = \sum_i (m_i^{\rm obs} - m_i^{\rm model})^2
   \]
   rather than \(\sum_i [(m_i^{\rm obs} - m_i^{\rm model})/\sigma_i]^2\).

2. **Bin numbering is zero-based.**  
   Pass `bin=0` for the first bin, `bin=1` for the second, etc.

3. **The internal bin upper bound is hard-coded to 18.**  
   This is compatible with GlueX, but for smaller electroproduction datasets you should only use the bins that actually exist.

4. **Photoproduction mode is special.**  
   You must set `photoProduction=true`; otherwise the code will try to reconstruct electroproduction scaling quantities that do not belong in the pure photo case.

---

## `RunGivenMoments_Chi2Amps_Bootstrap.C`

### Purpose

This file quantifies uncertainties by throwing toy moment sets and refitting them.

### What it does

For each toy:

1. clone the nominal observed moments,
2. Gaussian-throw each moment using its stored uncertainty,
3. run the same fit machinery,
4. keep the best-fit solution for that toy.

### Output

The output tree is:

- `PartialWaves`

Each row is one toy fit result.

Branches include:
- `toy`
- `log_val`
- the toy-thrown observed moments
- fitted amplitudes
- reconstructed model moments

### Parallelization

The code parallelizes over workers, writes `.part_*.root` files, merges them, then deletes the temporary parts.

### Default settings inside the macro

At present the wrapper hard-codes:

- `nToys = 1000`
- `nStartsPerToy = 1000`
- `nCores = 10`

If you need other values, you will likely want to edit the macro or call the lower-level setup helper.

### Caveat

Some notebooks still look for a tree called `toyFitResults` for part of the GlueX workflow, but the current bootstrap macro writes `PartialWaves`.

---

## `GenerateMomentsFromFixedAmplitudes.C`

### Purpose

This file creates a known-truth synthetic dataset.

It is useful for:

- closure tests,
- ambiguity studies,
- checking whether two apparently different amplitude solutions generate the same observable moments.

### What you edit

The main editable part is `FillUserAmplitudes(...)`, where the code explicitly sets values like:

- `a_T_1_1`
- `a_T_1_0`
- `a_T_1_m1`
- `a_L_1_1`
- etc.
- and the corresponding phases.

After that, the code normalizes the amplitudes, computes all model moments, reconstructs `RH` / `RH04`, and writes them to a ROOT tree.

### Output

The output tree is:

- `genMoments`

and includes:
- `R`
- all parameter branches
- full `RH04_*_*`
- the exact observed subset used by the fitter
- corresponding `_err` branches
- `Q2`

### Best use case

Use this file when you want to answer questions like:

- “If I generate data from this amplitude point, does the fitter recover it?”
- “Does flipping these two phases leave the observed moments unchanged?”
- “Are two minima physically ambiguous or just numerically different?”

### Caveat

The `Q2vals` argument is declared as a `std::vector<double>`, but the branch creation only names a single `Q2` branch. In practice, treat this as a **single-value input** unless you refactor the branch-writing logic.

---

## Analysis notebooks

## `AnalysisScripts/HERMES_analysis.ipynb`

### Purpose

Post-processes the electroproduction fits for:

- `e_rho`
- `mu_rho`
- `e_omega`
- `mu_omega`

### What it does

- reads `fitResults` trees from files like
  - `HERMES_results_i.root`
  - `MuRho_results_i.root`
  - `ElOmega_results_i.root`
  - `Omega_results_i.root`
- reads bootstrap outputs from files like
  - `HERMES_bootstrap_i.root`
  - etc.
- converts amplitude magnitudes/phases into complex-plane points
- makes “elliptical” uncertainty plots
- studies amplitude scaling vs. \(Q^2\)

### Important note

The notebook uses **absolute file paths from the original author’s machine**. You will need to edit those paths.

---

## `AnalysisScripts/GlueX_analysis.ipynb`

### Purpose

Post-processes GlueX photoproduction fits and compares the extracted amplitudes with another fit framework.

### What it does

- reads `gluex_results_i.root`
- reads `gluex_bootstrap_i.root`
- plots complex amplitudes vs. \(-\bar t\)
- compares this project’s output against an external BruFit-style result file

### Important note

This notebook also contains hard-coded absolute file paths and a tree-name assumption that appears older than the current bootstrap macro.

---

## `AnalysisScripts/R_analysis.ipynb`

### Purpose

Studies the ratio \(R\) as a function of \(Q^2\) for the electroproduction channels.

### What it does

- reads bootstrap outputs,
- extracts means / spreads,
- compares against externally quoted \(R\) values,
- fits simple VMD-style parameterizations using `scipy.optimize.curve_fit`.

This notebook is more physics-summary oriented than inversion-engine oriented.

---

## `AnalysisScripts/fixed_amp_analysis.ipynb`

### Purpose

This is the natural follow-up notebook for `GenerateMomentsFromFixedAmplitudes.C`.

### What it does

- loads a closure-test output file like `fixedamps.root`
- compares recovered amplitudes to the known truth
- visualizes the result in magnitude/phase space and in the complex plane

This is probably the best notebook to start with if you want to understand solution ambiguities.

---

## Recommended practical usage order

If you are new to the repository, I would use it in this order:

1. **Create the directory structure**
2. **Run `GenerateMomentsFromFixedAmplitudes.C`**
3. **Fit that generated file with `RunGivenMoments_Chi2Amps.C`**
4. **Inspect with `fixed_amp_analysis.ipynb`**
5. **Then move to a real dataset**
   - `MakeLeptoMoments.C` for electroproduction
   - `MakePhotoMoments.C` for GlueX photoproduction
6. **Run the bootstrap wrapper**
7. **Use the analysis notebooks only after editing paths/tree names**

That order lets you debug the inversion with a known truth before touching real data.

---

## Suggested naming convention for outputs

The notebooks already expect a naming pattern close to:

### Electroproduction
- `HERMES_results_0.root`, ..., `HERMES_results_3.root`
- `HERMES_bootstrap_0.root`, ..., `HERMES_bootstrap_3.root`
- `MuRho_results_0.root`, ...
- `MuRho_bootstrap_0.root`, ...
- `ElOmega_results_0.root`, ...
- `Omega_results_0.root`, ...

### Photoproduction
- `gluex_results_0.root`, ..., `gluex_results_17.root`
- `gluex_bootstrap_0.root`, ..., `gluex_bootstrap_17.root`

If you follow that convention, the plotting notebooks will need fewer edits.

---

## Common pitfalls

1. **Forgetting to create the output directories**
2. **Using 1-based rather than 0-based bin numbering**
3. **Passing `photoProduction=false` for GlueX**
4. **Using notebook paths from the author’s machine without editing them**
5. **Assuming `"omega"` is a working dataset key in `MakeLeptoMoments.C`**
6. **Assuming the main fit uses experimental errors in the chi-square weights**
7. **Trying to interpret every local minimum as unique physics rather than as a discrete ambiguity**
8. **Using more bins than the dataset actually contains**

---

## Short “how do I run everything?” checklist

### Electroproduction example
```bash
mkdir -p InputFiles/Experiment InputFiles/Generated OutputFiles

root -l -q 'MakeLeptoMoments.C("e_rho")'

for i in 0 1 2 3; do
  root -l -q "RunGivenMoments_Chi2Amps.C(\"InputFiles/Experiment/e_rho_moments.root\",\"expMoments\",$i,\"HERMES_results_${i}.root\",1.0,false)"
  root -l -q "RunGivenMoments_Chi2Amps_Bootstrap.C(\"InputFiles/Experiment/e_rho_moments.root\",\"expMoments\",$i,\"HERMES_bootstrap_${i}.root\",1.0,false)"
done
```

### Photoproduction example
```bash
root -l -q 'MakePhotoMoments.C("gluex")'

for i in $(seq 0 17); do
  root -l -q "RunGivenMoments_Chi2Amps.C(\"InputFiles/Experiment/gluex_moments.root\",\"expMoments\",$i,\"gluex_results_${i}.root\",1.0,true)"
  root -l -q "RunGivenMoments_Chi2Amps_Bootstrap.C(\"InputFiles/Experiment/gluex_moments.root\",\"expMoments\",$i,\"gluex_bootstrap_${i}.root\",1.0,true)"
done
```

### Closure-test example
```bash
root -l -q 'GenerateMomentsFromFixedAmplitudes.C()'
root -l -q 'RunGivenMoments_Chi2Amps.C("InputFiles/Generated/fixed_input_moments.root","genMoments",0,"fixedamps.root",1.0,false)'
```

---

## Final assessment

This repository is not a polished package yet; it is much closer to a **research working codebase**. But the core logic is clear and scientifically meaningful:

- translate SDMEs into moments,
- encode the moment formalism in a reflectivity amplitude basis,
- fit amplitudes with many randomized starts,
- propagate uncertainties with toy Monte Carlo,
- analyze discrete ambiguities with notebooks and synthetic tests.

If you wanted to improve it next, the highest-impact upgrades would be:

1. switch the fit objective to a truly error-weighted chi-square,
2. clean up dataset-key handling in `MakeLeptoMoments.C`,
3. remove hard-coded notebook paths,
4. move the macros into headers/sources or at least a cleaner reusable structure,
5. add a real README and a reproducible environment file.

