# flatppl-physics-neutrinos

[Function reference](docs/reference.md)

Three-flavour vacuum oscillations, coherent constant and layered matter
propagation, spherical-shell paths and neutrino mass observables in FlatPPL.
Library version **0.1**.

```flatppl
flatppl_compat = "0.1"
nu = load_module("https://raw.githubusercontent.com/flatppl-community/flatppl-physics-neutrinos/main/main.flatppl")
p = nu.parameters(0.59, 0.15, 0.78, 1.2, 7.5e-5, 2.5e-3)
P = nu.matter_probabilities(0.6, 295, p, 1.4, false)
appearance = P[1, 2] % nu_mu -> nu_e
```

Energy uses GeV, baseline km, mass splittings eV² and electron density mol/cm³.
Probability rows are final flavour and columns initial flavour, ordered
electron, muon, tau. All parameters are explicit. The library supplies no
experiment data, density profile, flux model or priors.

The [appearance example](examples/playground.flatppl) fits one mixing angle
from synthetic counts with an explicit prior, observation kernel and posterior.

Run the independent complex-rotation, Taylor-propagation and analytic checks
with Python 3.10+ and Node 24+:

```sh
FLATPPL_JS_DIR=/path/to/flatppl-js python tests/check_math.py
```

The [runtime fixture](tests/runtime.flatppl) and `runtime_values` in the checker
also support StableHLO value, gradient and batching checks.
