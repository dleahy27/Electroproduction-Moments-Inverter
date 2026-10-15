# Source guide

The public C++ interface is in `include/emi/Config.h` and `Runner.h`. The CLI
in `app/main.C` loads a settings header through `RuntimeSettings.C`, applies
explicit CLI overrides, and calls the library. The source files below `src/` implement the library operations.

## Ordinary amplitude fits

1. `Input.C` discovers moment branches and reads the requested bin. Published
   inputs use array branches in entry zero; generated fixed points use scalars.
2. `Context.C` connects branch names, the complete amplitude vector, and Minuit's
   independent coordinates. It validates configured start distributions and
   prepares model/observation indices once per worker.
3. `Model.C` and `PolarizedModel.C` build sparse angular bilinears, sharing the
   Clebsch–Gordan and amplitude-naming conventions. Explicit Pauli indices apply
   to target/recoil polarization; unpolarized moments use the reduced form.
4. `Evaluation.C` reconstructs constrained amplitudes, evaluates moments and
   chi-square, and applies the normalization chain rule to analytic gradients.
5. `Minimizer.C` creates fresh Minuit2 state per start. `Uncertainties.C` propagates
   its covariance to full amplitudes, moments, correlations, and R.
6. `FitRunner.C` saves every start; `Bootstrap.C` fluctuates moment values and
   saves one selected finite minimum per toy. Each runner applies its own
   fit-selection rules.

`Detail.h` declares these shared internal operations and value types. Sparse
terms, packed quantum-number keys, and cached phase pairs avoid repeated angular
algebra/string searches in the minimizer. They are prepared once; the evaluation
loops work with ordinary indexed vectors.

## Generation and mass fits

`FixedMoments.C` creates configured or example amplitude points and converts them
to the same measured moments consumed by fits. `RandomMoments.C` creates an
ensemble. `LeptoMoments.C` and `PhotoMoments.C` convert the published SDME tables;
the numerical tables and conversion factors are in those source files.
`MassModels/` contains deterministic line shapes and the PhotoTest truth model.

`MassFit.C` reuses the angular bilinears and replaces binwise amplitudes with
coherent Breit–Wigner/background functions across a mass grid. It builds shared
pole coordinates, chooses the phase reference, evaluates the mass chi-square and
gradient, and writes predictions from the best accepted global start. Generation,
gradients, and saved predictions share the H0 + epsilon H4 response contraction.
The first global start keeps the supplied coordinates; later starts randomize
free couplings using the existing worker seed sequence.

## ROOT ownership and workers

ROOT branches keep addresses, so their backing vectors must be allocated before
branch creation and must retain that storage until `Fill()` finishes. Assign
values in place when an expression could otherwise move a temporary vector.
Each worker owns its context, minimizer, and ROOT file. `WorkerOutput.h` supplies
only part-file naming and merging; runners retain scheduling, seed arithmetic,
selection, and cleanup order. Parts remain available if final processing fails.

`Progress.h` uses a shared lock-free atomic because a normal counter would be
copied when ROOT forks. Only the parent prints progress. The counter publishes
no physics data; relaxed atomic ordering is sufficient for this display.
