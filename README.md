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

Propagation returns $P_{\beta\alpha}=|S_{\beta\alpha}|^2$. At constant density, $S=\exp(-iCLH/E)$, with $C=2.533865358078198$ and mass-squared Hamiltonian $H$ in eV². Layers compose coherently as $S=S_n\cdots S_1$ in source-to-detector order.

Energy uses GeV, baseline km, mass splittings eV², and electron density mol/cm³. Probabilities have shape `(3, 3)`, with final-flavour rows and initial-flavour columns, ordered electron, muon, tau. Amplitudes use `record(re, im)` with the same matrix shape.

The [appearance example](examples/playground.flatppl) infers one mixing angle from synthetic counts with an explicit prior and observation kernel.

See the reference for definitions, domains, numerical limits, and validation commands.
