# Cross section of the 2D model: from the image I(Qx,Qy) to sigma

How the `DirectLoad2D` model turns the tabulated intensity `I(Qx,Qy)` — what
SasView exports: a differential cross section per unit solid angle over
detector-plane coordinates — into the total cross section `σ` (barn/atom) that
transport needs, for every incident direction and energy. Extracted from the
specification, section "Cross sections → Anisotropic case: from the 2D image
to the cross section"
file format on [Data format](data-format), usage in [Tutorial: 2D models](tutorial-2d).
A byproduct of the `σ` computation is used by the rejection method to sample
the scattered neutron momentum `k_f`, following the rejection method of
X.-X. Cai *et al.*,
[Rejection-based sampling of inelastic neutron scattering, J. Comput. Phys. (2019)](https://www.sciencedirect.com/science/article/pii/S0021999118307885).

## Symbols


| symbol           | meaning                                                                                                        |
| ---------------- | -------------------------------------------------------------------------------------------------------------- |
| `E`              | neutron kinetic energy                                                                                         |
| `k`              | wavenumber:`k = sqrt(E / (2.07214 meV·Aa²))`, i.e. `E[meV] = 2.07214 · k²` for `k` in 1/Aa                 |
| `k_i, k_f`       | incident / final wave vectors,`k_i = k·k̂_i`, `k_f = k·k̂_f`                                               |
| `k̂_i, k̂_f`   | incident / final unit directions; elasticity:`|k_f| = |k_i| = k`                                               |
| `a_x, a_y, a_z`  | components of the incident unit direction,`k̂_i = (a_x, a_y, a_z)`, in the material frame                     |
| `Q`              | scattering vector,`Q = k_f − k_i`; magnitude `Q = 2k·sin(θ/2)`, `0 ≤ Q ≤ 2k`                              |
| `Qx, Qy`         | in-plane components`Q·x̂` and `Q·ŷ` — the two coordinates of the 2D table                                 |
| `x̂, ŷ, ẑ`    | material-frame axes; the table represents`I` on the plane `Q·ẑ = 0` (the `q_z = 0` convention)               |
| `θ, ψ`         | outcome polar angle (from`ẑ`) and azimuth: `k̂_f = (sinθ·cosψ, sinθ·sinψ, cosθ)`                      |
| `I_table(Qx,Qy)` | the tabulated intensity, barn/(atom sr);`I = 0` outside the grid                                               |
| `D`              | disc of in-plane points reachable by elastic scattering at this`k̂_i`, `E` (Step 2)                           |
| `ρ'`            | distance from the disc centre:`ρ'² = (Qx + k·a_x)² + (Qy + k·a_y)²`                                      |
| `α, φ`         | polar parameters of the disc (Step 4); for the design beam`α` is the outcome polar angle and `φ` its azimuth |
| `σ(k̂_i, E)`   | total cross section, barn/atom — the quantity transport uses                                                  |

## The algorithm in five steps

### Step 1 — definition: integrate the image over the sphere

An outcome whose scattering vector has in-plane components `(Qx,Qy)`
contributes `I_table(Qx,Qy)·dΩ` to `dσ` when its direction covers solid angle
`dΩ`. The total cross section integrates this over all elastic outcomes:

```
σ(k̂_i, E) = ∮ I( Q(k̂_f)·x̂ , Q(k̂_f)·ŷ ) dΩ ,      Q = k_f − k_i ,
```

an integral over the unit sphere of `k̂_f`; elasticity `|k̂_f| = 1` links `Q`
to `k̂_f` at fixed `k̂_i`, `E`.

### Step 2 — elasticity folds the sphere onto the plane

Writing `k̂_f = k̂_i + Q/k` and imposing `\|k̂_f\| = 1`: every in-plane point of
the disc

```
D = { (Qx,Qy) :  (Qx + k·a_x)² + (Qy + k·a_y)² ≤ k² }
```

(radius `k`, centred at `(−k·a_x, −k·a_y)`) is reached by exactly **two**
elastic outcomes, the two branches of the elastic sphere:

```
k̂_f·ẑ = ± sqrt(k² − ρ'²) / k .
```

Each plane pixel of `D` is thus the image of two equal-weight sphere patches,
and the model *defines* `I` on the backward branch by the same plane value
(the `q_z = 0` mirror extension — exact for structures mirror-symmetric about
the `xy` plane, and exact in the closed-form limits of Step 5 for any
structure).

### Step 3 — Jacobian: the plane integral

The central projection along `ẑ` from the unit sphere to in-plane coordinates
has area element `dA = k²·|k̂_f·ẑ|·dΩ` per branch (the shift by
`−k·(a_x,a_y)` preserves area elements). Inverting, `dΩ = dA / (k²·|k̂_f·ẑ|)`,
and summing the two branches:

```
σ(k̂_i, E) = 2 · ∬_D  I_table(Qx,Qy) / ( k² · sqrt(k² − ρ'²) )  dQx dQy
```

This is the definition implemented by the plugin: it uses nothing but the
table, the kinematics, and the disc. The integrand diverges (integrably) as
`ρ' → k`, where `|k̂_f·ẑ| → 0` (the horizon).

### Step 4 — the smooth angular form (what is actually evaluated)

Parametrise the disc by

```
(Qx, Qy) = (−k·a_x, −k·a_y) + k·sinα·(cosφ, sinφ) ,      α ∈ [0, π] :
```

then `dA = k²·sinα·cosα·dα·dφ` and `sqrt(k² − ρ'²) = k·cosα`; the range
`α ∈ [0, π/2]` covers the forward branch and `α ∈ [π/2, π]` traverses the
disc a second time (the backward branch), giving

```
σ(k̂_i, E) = ∫₀^π ∫₀^{2π}  I_table( Qx(α,φ), Qy(α,φ) ) · sinα  dφ  dα .
```

No singularity anywhere. In the material frame `α` is exactly the outcome
polar angle and `φ` its azimuth, so this integral *is* `∮ I dΩ` itself; the
reference harness evaluates it by Gauss–Legendre quadrature in `α`.

### Step 5 — closed-form limits (the validation anchors)

- **Constant table `I ≡ I₀`:** `σ = 4π·I₀`, independent of `k` and of beam
  direction.
- **`E → 0` (`k → 0`):** the disc shrinks onto the origin, so
  `σ → 4π·I(0)` — for any table, any beam direction.
- **`E → 0` expansion:** with `g = ∇I(0)`, `H = ∇²I(0)`, `a_⊥ = (a_x, a_y)`,

```
σ = 4π·I₀ − 4π·k·a_⊥·g + k²·( 2π·a_⊥ᵀ·H·a_⊥ + (2π/3)·tr H ) + O(k³) .
```

The linear term is the response of `σ` to beam tilt, generated by the
disc-centre shift `−k·a_⊥`. For a linear ramp `I = I₀ + g·Q` the series
terminates, so `σ = 4π·(I₀ − k·a_⊥·g)` holds exactly at any `k`.

Both branches are required: a forward-only convention would give `2π·I₀`
instead of `4π·I₀` and break the `E → 0` isotropic limit. The cost of the
mirrored backward branch — up to a factor-2 inflation of the *total* `σ` for
narrow patterns at `2k ≫` pattern width — is a documented convention effect,
derived in the specification section "The backward branch and the relation to
the 1D model" and summarised in [Tutorial: 2D models](tutorial-2d)
(Conventions).

## Result

`σ` depends on the incident **direction**, through the disc centre `−k·a_⊥`,
i.e. through the angle between the beam and the material pattern — like Bragg
reflection in a single crystal. This `σ(k̂_i, E)` is what the plugin's
`crossSection` returns and what transport uses for collision rates. How
outcome directions are *sampled* from the same image (per-energy CDFs, exact
and rejection-free on axis) is specified in the specification section
"Cross-section evaluation and sampling" and described hands-on in
[Tutorial: 2D models](tutorial-2d).
