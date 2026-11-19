# Neutrino function reference

Library version **0.1** · [Package overview](../README.md)

## Conventions

Flavours are ordered electron, muon, tau. A mixing matrix has flavour rows and
mass-eigenstate columns. An amplitude or probability matrix has **final-flavour
rows and initial-flavour columns**: `P[1, 2]` means muon to electron.

| Quantity | Unit or convention |
| --- | --- |
| `energy` | GeV, strictly positive |
| `baseline`, layer `lengths` | km, nonnegative |
| Mixing angles and phases | radians |
| `dm21` | $m_2^2-m_1^2>0$, eV² |
| `dm31` | $m_3^2-m_1^2$, eV²; positive normal ordering, negative inverted ordering |
| Absolute masses | eV, nonnegative |
| `electron_density` | mol of electrons per cm³, nonnegative |
| `antineutrino` | Boolean; `false` selects neutrinos |

Complex matrices use `record(re = R, im = I)`, with two real $3\times3$
matrices. This representation preserves the real tensor interface in both
engines. Propagation keeps amplitudes coherent across layers and applies the
absolute square only at the end. No clipping, renormalization, phase averaging,
absorption or decoherence is applied.

| Function | Purpose |
| --- | --- |
| [`parameters`](#parameters) | Collect six oscillation parameters |
| [`pmns`](#pmns) | PDG mixing matrix |
| [`vacuum_amplitude`](#vacuum_amplitude) | Coherent vacuum amplitude |
| [`vacuum_probabilities`](#vacuum_probabilities) | Vacuum transition probabilities |
| [`matter_hamiltonian`](#matter_hamiltonian) | Standard matter mass-squared Hamiltonian |
| [`matter_amplitude`](#matter_amplitude) | Constant-density amplitude |
| [`matter_probabilities`](#matter_probabilities) | Constant-density probabilities |
| [`layered_amplitude`](#layered_amplitude) | Ordered coherent layer propagation |
| [`layered_probabilities`](#layered_probabilities) | Layered transition probabilities |
| [`masses`](#masses) | Three masses for either ordering |
| [`beta_mass`](#beta_mass) | Kinematic beta-decay effective mass |
| [`majorana_mass`](#majorana_mass) | Light-Majorana effective mass |
| [`atmospheric_baseline`](#atmospheric_baseline) | Surface-detector source distance |
| [`shell_path`](#shell_path) | Full chord through supplied spherical shells |

## Mixing and vacuum propagation

### `parameters`

`parameters(theta12, theta13, theta23, delta, dm21, dm31)` returns a record
with those six fields. It chooses no prior or best-fit value. Mixing angles
lie in $[0,\pi/2]$, while `delta` is periodic modulo $2\pi$. For normal
ordering require `dm31 > dm21`; for inverted ordering require `dm31 < 0`.
The constructor collects values; these domains are caller preconditions.

### `pmns`

`pmns(theta12, theta13, theta23, delta)` returns the real and imaginary parts
of $U=R_{23}U_{13}(\delta)R_{12}$. With $c_{ij}=\cos\theta_{ij}$ and
$s_{ij}=\sin\theta_{ij}$:

$$
U=\begin{pmatrix}
c_{12}c_{13}&s_{12}c_{13}&s_{13}e^{-i\delta}\\
-s_{12}c_{23}-c_{12}s_{23}s_{13}e^{i\delta}&c_{12}c_{23}-s_{12}s_{23}s_{13}e^{i\delta}&s_{23}c_{13}\\
s_{12}s_{23}-c_{12}c_{23}s_{13}e^{i\delta}&-c_{12}s_{23}-s_{12}c_{23}s_{13}e^{i\delta}&c_{23}c_{13}
\end{pmatrix}.
$$

Majorana phases are omitted here because they cancel from oscillation probabilities.

### `vacuum_amplitude`

`vacuum_amplitude(energy, baseline, p, antineutrino)` takes the parameter
record `p` and returns a complex amplitude record. Define
$d=(0,\Delta m_{21}^2,\Delta m_{31}^2)$ and $C=2.5338653580781976$:

$$S_{\beta\alpha}=\sum_{j=1}^3 U_{\beta j}
 e^{-iCd_jL/E}U^*_{\alpha j}.$$

The common $m_1^2$ phase is removed. Antineutrinos use $U^*$ in place of $U$.
At zero baseline the amplitude is the identity. This function uses the
closed-form spectral expression, without a matrix exponential.

### `vacuum_probabilities`

`vacuum_probabilities(energy, baseline, p, antineutrino)` returns the real
$3\times3$ matrix $P_{\beta\alpha}=|S_{\beta\alpha}|^2$ from
`vacuum_amplitude`. Each column sums to one up to floating-point error.

## Matter propagation

### `matter_hamiltonian`

`matter_hamiltonian(energy, p, electron_density, antineutrino)` returns the
complex **mass-squared** Hamiltonian $H$ in eV²:

```math
H=U\,\mathrm{diag}(0,\Delta m_{21}^2,\Delta m_{31}^2)U^\dagger
 +\mathrm{diag}(a,0,0),\qquad a=2E V_{CC}.
```

The energy in the expression for $a$ is converted to eV.
$V_{CC}=\sqrt{2}G_FN_A n_e$, using
$G_F=8.961877245622253\times10^{-38}$ eV cm³ and
$N_A=6.02214076\times10^{23}$ mol⁻¹. Numerically,
$a=1.52649324494995\times10^{-4}E n_e$ eV² for the documented input units.
For antineutrinos conjugate $U$ and reverse the sign of $a$.

Electron density is not bulk mass density. For ordinary matter the usual
conversion is approximately $n_e=\rho Y_e$ mol/cm³ when $\rho$ is in g/cm³
and $Y_e$ is electrons per nucleon. Supply that conversion explicitly.

### `matter_amplitude`

`matter_amplitude(energy, baseline, p, electron_density, antineutrino)` returns
$S=\exp(-iCLH/E)$ as a complex amplitude record, using `matter_hamiltonian`.
It evaluates the equivalent real $6\times6$ exponential:

```math
\exp\left[\frac{CL}{E}
\begin{pmatrix}\Im H&\Re H\\-\Re H&\Im H\end{pmatrix}\right].
```

The upper-left and lower-left blocks are $\Re S$ and $\Im S$. At zero
density this equals `vacuum_amplitude`, including its phase convention.

### `matter_probabilities`

`matter_probabilities(energy, baseline, p, electron_density, antineutrino)`
returns the absolute square of `matter_amplitude`, entry by entry.

### `layered_amplitude`

`layered_amplitude(energy, p, lengths, electron_densities, antineutrino)`
takes **nonempty**, aligned vectors in source-to-detector traversal order.
It returns the complex amplitude

$$S=S_nS_{n-1}\cdots S_1,$$

where $S_j$ uses `lengths[j]` and `electron_densities[j]`. Zero-length and
zero-density layers are allowed. The library computes the vacuum Hamiltonian
once and uses `scan` to compose the layer amplitudes. Reversing unequal layers
can change the result.

### `layered_probabilities`

`layered_probabilities(energy, p, lengths, electron_densities, antineutrino)`
returns the absolute square of `layered_amplitude`. Multiplying the individual
layer probability matrices would discard interference and gives a different model.

## Mass observables

### `masses`

`masses(lightest, dm21, dm31)` returns `[m1, m2, m3]`. Define
$q=\max(0,-\Delta m_{31}^2)$. Then

$$m_j=\sqrt{m_{\mathrm{lightest}}^2+q+d_j}.$$

This assigns `lightest` to $m_1$ in normal ordering and $m_3$ in inverted
ordering. Changing the ordering is a model choice, not an interpolation through
`dm31 = 0`. Use `sum(masses(...))` for the sum of masses.

### `beta_mass`

`beta_mass(m, theta12, theta13)` takes the three masses and returns

$$m_\beta=\sqrt{\sum_j |U_{ej}|^2m_j^2}.$$

This kinematic effective mass is independent of the CP and Majorana phases.

### `majorana_mass`

`majorana_mass(m, theta12, theta13, delta, alpha21, alpha31)` uses
$U_{\mathrm{PDG}}\mathrm{diag}(1,e^{i\alpha_{21}/2},e^{i\alpha_{31}/2})$
and returns

$$m_{\beta\beta}=\left|c_{13}^2c_{12}^2m_1
 +c_{13}^2s_{12}^2m_2e^{i\alpha_{21}}
 +s_{13}^2m_3e^{i(\alpha_{31}-2\delta)}\right|.$$

It is the effective mass for light-Majorana exchange, not a complete
double-beta-decay rate. The absolute value is nondifferentiable at exact cancellation.

## Spherical geometry

### `atmospheric_baseline`

`atmospheric_baseline(cos_zenith, radius, height)` assumes a detector at the
sphere's surface and production at radius `radius + height`. Distances use
one consistent unit, with `radius > 0`, `height >= 0` and
`-1 <= cos_zenith <= 1`. Down-going vertical means `cos_zenith = +1`:

$$L=\sqrt{(R+h)^2-R^2(1-c_z^2)}-Rc_z.$$

The down-going branch uses an equivalent rationalized expression to avoid
subtracting nearly equal distances. Underground detector depth and a
production-height distribution are separate model inputs.

### `shell_path`

`shell_path(impact, radii, electron_densities)` returns a table with columns
`length` and `electron_density`. `radii` is a nonempty vector of strictly
increasing positive outer radii. Its aligned densities describe shells from
the centre outward. Require `0 <= impact <= radii[last]`.

Let $R_0=0$ and

$$h_j=\sqrt{\max(R_j^2-b^2,0)}-\sqrt{\max(R_{j-1}^2-b^2,0)}.$$

The returned lengths are $(h_n,\ldots,h_1,h_1,\ldots,h_n)$, with densities
in the same order. Their sum is the full chord length
$2\sqrt{R_n^2-b^2}$. Zero intersections remain as zero-length layers.
Derivatives are well-defined away from shell tangencies; the square-root
boundary itself is not differentiable.

```flatppl
path = nu.shell_path(3000, [3480, 6371], [5.5, 2.4])
P = nu.layered_probabilities(2.0, p, path.length, path.electron_density, false)
```

This is a **full through-sphere chord**, not an automatic atmospheric Earth
path. Use it only for up-going Earth crossings. Down-going trajectories have
no Earth chord. A vacuum atmosphere segment can be supplied explicitly before
the shell segments. No Earth density model or production profile is bundled.

## Numerical scope and sources

Vacuum and standard coherent matter propagation are included. Sterile
states, nonstandard interactions, absorption, resolution smearing, fluxes and
detector responses are separate extensions or analysis inputs. Choose a
precision that resolves the oscillation phases of the supplied $L/E$ range.

Conventions follow the [PDG neutrino review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-neutrino-mixing.pdf).
The explicit physics/experiment separation and physical conversion constants
also follow [Newtrinos.jl](https://github.com/philippeller/Newtrinos.jl), described
by [Eller and Schultheiß](https://arxiv.org/abs/2608.12940). The library uses its
own FlatPPL implementation; its checker composes complex rotation matrices
and a complex Taylor propagator independently.

## Examples and validation

[Run the appearance example in the playground](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-physics-neutrinos%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post).
It fits one mixing angle from synthetic counts with an explicit prior,
observation kernel and posterior, and pins
[library revision 608fa00](https://github.com/flatppl-community/flatppl-physics-neutrinos/tree/608fa0076fe22a9dd7634e35a5b1365b12af2403).

Run the independent complex-rotation, Taylor-propagation and analytic checks
with Python 3.10+ and Node 24+:

```sh
FLATPPL_JS_DIR=/path/to/flatppl-js python tests/check_math.py
```

The [runtime fixture](../tests/runtime.flatppl) and `runtime_values` in the checker
also support StableHLO value, gradient and batching checks.
