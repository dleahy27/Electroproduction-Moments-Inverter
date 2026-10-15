# Full PhotoTest closure

This case generates and fits the same four-resonance PhotoTest model, including
both reflectivities and the k=+1 and k=-1 sectors. The mass grid has 21 points
from 0.8 to 1.8 GeV. Moments are noiseless; assigned uncertainties of 0.001 set
the chi-square weights.

After building EMI, run:

```bash
python -m tutorials.photo_test_m.full_closure.run
python -m tutorials.photo_test_m.full_closure.analyse --check
```

The defaults in [settings.h](settings.h) use one start at the generated
parameters, with HESSE disabled. This checks the full model and output path;
it does not test recovery from random initial parameters. Increase the number
of starts for that study.

The runner writes `output/moments.root` and `output/fit.root`. Analysis checks
the fitted parameters and moments and writes figures to `output/plots/`.
