# flatppl-physics-neutrinos

[Function reference](docs/reference.md) · [Live example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-physics-neutrinos%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)

Three-flavour neutrino oscillations, coherent matter propagation, spherical-shell paths, and mass observables for FlatPPL. Library version **0.1**.

```flatppl
flatppl_compat = "0.1"
nu = load_module("https://raw.githubusercontent.com/flatppl-community/flatppl-physics-neutrinos/608fa0076fe22a9dd7634e35a5b1365b12af2403/main.flatppl")
p = nu.parameters(0.59, 0.15, 0.78, 1.2, 7.5e-5, 2.5e-3)
P = nu.matter_probabilities(0.6, 295, p, 1.4, false)
appearance = P[1, 2] % nu_mu -> nu_e
```

Energy uses GeV, baseline km, mass splittings eV², and electron density mol/cm³. Matrix rows index final flavour. Columns index initial flavour. Flavour order is electron, muon, tau.

The [appearance example](examples/playground.flatppl) infers one mixing angle from synthetic counts with an explicit prior and observation kernel.

See the reference for definitions, domains, numerical limits, and validation commands.
