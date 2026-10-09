# Tutorial: 2D anisotropic SANS with DirectLoad2D

The `DirectLoad2D` model of the SasCSNS plugin transports neutrons through a
**tabulated, anisotropic** 2D intensity `I(Qx,Qy)`: the detector-plane image of
an oriented-particle SANS model. As with the 1D [DirectLoad](tutorial-1d)
model, the plugin develops no new physics models; the pattern is computed
upstream and the plugin transports neutrons through it. The specification puts
it this way:

> **SasView computes physics; the plugin transports it.** A single new model type, `DirectLoad2D`, ingests a generic tabulated 2D intensity `I(qx,qy)`; everything shape-specific happens upstream, in SasView, and in the converter tool `sasview2ncmat.py`.
>
> — specification, section "Scope and design principle" ([doc/anisotropic_directload2d.pdf](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/doc/anisotropic_directload2d.pdf))

This tutorial walks through the full 2D workflow — exporting `I(Qx,Qy)` from
SasView and converting it (`sasview2ncmat.py --2d`), loading and sampling it
through the NCrystal Python API, cross-checking the plugin cross section
against independent quadrature, and transport-integrated attenuation — and
closes with the model's coordinate **Conventions**, quoted verbatim from the
specification. Data producers should read that section before deciding whether
the 2D model is appropriate for their use case.

## Prerequisites

- The SasCSNS plugin installed and importable, per [Installation](installation)
  (`pip install .` in the repo; `python` must see the plugin).
- Steps 1 and 3 need `numpy` and `scipy` (step 1 only via the chain-script
  shortcut suggested there — the converter itself needs only the standard
  library); `sasmodels` is optional (one guarded cross-check is skipped when
  absent). No network access and no local
  SasView installation are needed — the example scripts emulate the exports.
- The NCMAT data format itself is described on [Data format](data-format).

## Step 1: export I(Qx,Qy) from SasView and convert

From SasView, save the 2D dataset as ASCII. The converter
[script/sasview2ncmat.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/sasview2ncmat.py)
reads rows of `Qx Qy I` (header lines start with `#`; extra columns are
ignored), with `Qx`,`Qy` in 1/Aa. No local SasView is needed to follow
along: running `script/example_sasview_chain_2d.py` first writes the
emulated export `cylinder2d_sasview_export.dat` into
`$SASVIEW_CHAIN_2D_WORKDIR` (default `/tmp/sasview_chain_2d`) — use it in
place of `sasview_2d_export.dat` below. Then run:

```sh
python3 script/sasview2ncmat.py sasview_2d_export.dat --2d \
        -o cylinder2d.ncmat --material sio2 --density 2.2 \
        --phi 0.4 --volume 4.19e6
```

The flags shown are the converter's standard set; the complete reference
table lives on [Data format, Converter
CLI](data-format#converter-cli-sasview2ncmatpy). In short: `--2d` selects
the `DirectLoad2D` model (without it the 1D pipeline runs), `--phi` and
exactly one of `--radius`/`--volume` are required since the unit conversion
needs `V_p`; `--solvent` is recorded as a `# solvent: NAME` comment in both modes (see
[Data format, Pitfalls](data-format#pitfalls)).

### The unit conversion is mandatory and automatic

The converter multiplies every `I` value by `scale = (phi / V_p) / n_d`
(barn/object -> barn/atom) — deliberately with **no flag** — recording the
applied scale in the output file header; all values stored in the
`DirectLoad2D` section are barn/(atom sr). The normative description (the
rationale, the `1 angstrom^2 = 1e8 barn` trap and the invariant `Sigma = n_d * sigma`) is on
[Data format, Units](data-format#units-barnatom-not-barnobject).

### Data requirements for the 2D export

The rules are normative on [Data format](data-format): the converter's
per-file checks and warnings — complete rectangular grid (hard error, fill
masked pixels first), at least 2 strictly ascending points per axis,
bilinear resampling of non-uniform axes onto a uniform grid (warning),
non-finite `I` a hard error, negative `I` clamped to 0 (warning), a warning
when either axis has fewer than 8 points, the `Q = 0` start and the
`--emax` coverage gate — are enumerated under
[Sanity checks performed](data-format#sanity-checks-performed);
the Q-range rule — with `E[meV] = 2.07214·k²` (`k` in 1/Aa), each half-side
of the table must cover `[-2k, 2k]` for the highest energy of interest,
otherwise cross sections are silently biased — is specified under
[Q-range and coverage requirements](data-format#q-range-and-coverage-requirements).

Producer guidance, not checked by the `--2d` converter: keep the Q spacing at
or below `π/(10·R)` — at least 10 points per form-factor oscillation period
`π/R`. That resolution rule is enforced (as a warning) only by the 1D
pipeline, and only when `--radius` is given.

## Step 2: load and sample with the Python API

```python
import NCrystal as NC

scat = NC.createScatter('cylinder2d.ncmat')

beam = (0.0, 0.0, 1.0)            # unit beam direction in the material frame
ekin = 2.07214e-3                 # eV  (E = 2.07214 meV  ->  k = 1 Aa^-1)

sigma = float(scat.crossSection(ekin, beam))          # barn/atom

ek_out, dirs = scat.sampleScatter(ekin, beam, repeat=200_000)
# ek_out equals ekin (scattering is elastic, exactly);
# dirs holds the outcome unit directions (x, y, z components in the
# material frame as separate arrays).
```

The bare relative path works because step 1 wrote `cylinder2d.ncmat` into
the current working directory (the repo root): NCrystal resolves relative
paths against the cwd by default, so no `addCustomSearchDirectory` call is
needed here — that call ([Tutorial: 1D](tutorial-1d), step 4) is only for
files living outside the cwd.

Like all NCrystal (3.x and 4.x) Python calls, `crossSection` and `sampleScatter` take
`ekin` in **eV** — pass `E_meV * 1e-3`.

Sampling behaviour, as stated in the README model table and the specification
(section "Cross-section evaluation and sampling"): on-axis sampling is exact
and rejection-free; tilted beams are handled via dilated-proposal rejection.
In implementation terms (constants in
[include/NCDirectLoad2D.hh](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/include/NCDirectLoad2D.hh);
defaults, subject to change): beams within 5 mrad of the +z material axis are
drawn from per-energy CDFs built from the image — exactly on axis this is a
single table lookup with no rejection, and the residual tilt inside the cone
is handled by a small envelope mask with near-unity acceptance; beams tilted
by up to
60 mrad are proposed from a dilated copy of the table and accepted with a
pointwise mask, with a measured acceptance of about 0.68 at 50 mrad for a
Gaussian table; beyond that, a uniform-sphere rejection fallback runs — slow
but correct. The design energy window is 0.1–100 meV; queries outside the
design range are clamped and flagged with a one-time warning. The whole plugin
model is gated by the `sans` request parameter (`sans=0` disables it).

Accuracy expectations, as measured by the validation harness (see
[Validation](validation)): on-axis cross sections are exact at the energy
nodes; off-node and tilted values carry grid interpolation errors at the 1e-3
to 1e-2 level for smooth tables (median 8.2e-4, worst 1.0e-2 over 150
randomised checks).

## Step 3: the full 2D chain

[script/example_sasview_chain_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/example_sasview_chain_2d.py)
(validation rung 6, see [Validation](validation)) is the complete runnable version of steps 1–2: it emulates
a SasView 2D export of a cylinder (R = 25 Aa, L = 30 Aa, axis tilted 30° from
+z, `phi = 0.13`) on a 161x161 grid, converts it with `sasview2ncmat.py --2d`,
reads the file back to verify `I_file = I_object · (phi/V_p)/n_d`, and compares
the plugin cross section against an independent midpoint quadrature over the
same table:

```python
# excerpt (abridged) from script/example_sasview_chain_2d.py
for e_mev, beam, tol in [(0.3, (0, 0, 1), 5e-3),
                         (0.7, (0, 0, 1), 5e-3),
                         (3.0, (0, 0, 1), 3e-2),   # disc truncated by table edge
                         (0.7, (math.sin(0.03), 0, math.cos(0.03)), 2e-2),
                         (0.7, (math.sin(0.06), 0, math.cos(0.06)), 2e-2)]:
    ref = sigma_quad(k_of_e(e_mev), beam)   # independent quadrature of the file table
    got = float(scat.crossSection(e_mev * 1e-3, beam))
```

The script also checks a far-out-of-range beam (90 degrees: the clamped result
must stay finite and positive), that the sampled outcomes reproduce the
elongated image (`<Qx^2>/<Qy^2>` must differ from 1 by more than 5 %), and that
sampling is exactly elastic. If `sasmodels` is installed, the emulated export
is cross-checked against the real cylinder kernel: shapes must agree to machine
precision and the ratio must be one constant — the unit convention that the
converter absorbs.

One lesson worth keeping: the script's reference quadrature uses a dense
**midpoint** rule, not Gauss-Legendre — on ring-structured tables the
alpha-integrand varies by orders of magnitude inside the forward boundary
layer, and a GL-64 rule was measured to overestimate sigma by about 12x while
the midpoint rule converges within a few hundred nodes. Always
convergence-check reference quadratures themselves.

The companion figure script
[script/fig_chain_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/fig_chain_2d.py)
renders `doc/fig_chain_2d.pdf` (input image vs. density of 2M sampled
outcomes) and auto-runs the chain script if its NCMAT file is missing.

## Step 4: attenuation (transport-integrated, validation rung 8)

[script/example_transmission_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/example_transmission_2d.py)
exercises, as validation rung 8 of the [validation ladder](validation), the
absolute cross-section scale the way a transport code sees it:
analog Monte Carlo through a slab, with free flight lengths drawn from
`Exp(1/Sigma)`:

```python
# excerpt from script/example_transmission_2d.py
info = NCrystal.createInfo(path)
n_inv_a3 = float(info.getNumberDensity())   # atoms/Aa^3, as transport sees it
sigma_tot = sig_s + sig_a                    # barn/atom
# unit cancellation: 1 barn = 1e-24 cm^2 and 1 Aa^3 = 1e-24 cm^3, so
# n * sigma IS the macroscopic cross section in cm^-1, with no factor:
sigma_inv = n_inv_a3 * sigma_tot             # 1/cm
s = -np.log(rng.random(n_neutron)) / sigma_inv
uncollided = int(np.sum(s > length_cm))      # must follow exp(-Sigma*L)
```

The material is a narrow radial Gaussian pattern (`s_Q = 0.1` 1/Aa, amplitude
40 barn/(atom sr)) at `E = 2.07214` meV; `1e6` neutrons per run. All four
checks pass: the uncollided fraction reproduces Beer-Lambert at
L = 0.25/1/3 mfp (pulls 0.1/0.9/1.7 sigma), the absorbed fraction of
interactions matches `Sigma_a/Sigma_tot` (+1.9 sigma), a 30 mrad tilted beam
attenuates with the tilted `sigma` (+1.7 sigma; `sigma` rises by 1.0 % at
30 mrad for this pattern — anisotropic attenuation survives transport), and
the 1D-vs-2D convention comparison of check 4 finds
`sigma_2D/sigma_1D = 2.026`, inside the accepted [1.8, 2.3] band. That last
row operationalises the convention effect described in the Conventions section
below.

## Conventions: the q_z = 0 convention and the backward branch

The specification is the authoritative statement of the model's coordinate
conventions; the key passages are reproduced below verbatim. In the
blockquotes, text in square brackets is **editorial**: it expands the
specification's internal cross-references (to its own sections and equation
labels, which do not resolve on this wiki) and is not part of the quoted
words. The convention checks belong to validation rung 7
(`script/validate_1d_vs_2d.py`, see [Validation](validation)).

**The convention, stated plainly.** A SasView 2D grid stores only the
*transverse* components of the scattering vector: a detector pixel labelled
`(qx, qy)` means the model-frame vector `R^-1 · [qx, qy, 0]^T` — the
longitudinal entry is set to zero before the orientation rotation `R` is
applied. The `DirectLoad2D` table is exactly that transverse pair,
`I(Qx, Qy)`, with `Qx`, `Qy` the in-plane components of `Q = kf − ki` in the
material frame, evaluated identically for forward and backward outcomes.

This is verifiable, not folklore: the specification (section "The
DirectLoad2D data format", paragraph "Provenance of the q_z=0 convention
(receipts)") reads the convention directly out of the sasmodels sources
(v1.1.0) at three independent levels:

> **Spec, section "The DirectLoad2D data format", paragraph "Provenance of the q_z=0 convention (receipts)":** "(i) *Manual*: the orientation guide derives the model-frame scattering vector from detector coordinates as `[qa,qb,qc]^T = R^-1 · [qx,qy,0]^T` (`docs-source/guide/orientation/orientation.rst`) — the longitudinal entry is set to zero before the rotation. (ii) *Kernel source*: the compiled C kernel implements exactly that equation; the helper `qac_apply` in `kernel_iq.c` computes `δq_c = R31·qx + R32·qy` and `q_ab² = qx² + qy² − δq_c²`, with no `q_z` input anywhere in the call chain. (iii) *Data model*: the 2D data container `Data2D` (`data.py`) carries only `qx`/`qy` arrays [...]" ([...] marks elided text.)

> **Same paragraph:** "As an end-to-end empirical check, the emulated export of rung 6 reproduces a real sasmodels kernel at 200 detector pixels to `4e-16` in shape."

### Accuracy of the convention (inherited from SasView)

The neglected longitudinal component is `Q_z = −2k·sin²(θ/2)`, so
`|Q_z|/Q_⊥ = tan(θ/2)`: a second-order-in-θ error, negligible for typical
SANS (θ ≲ 3°, particles up to a few hundred angstrom) and *not introduced by the
plugin* — it is SasView's own 2D convention. Quoted verbatim from the
specification (same section, paragraph "Accuracy of the q_z=0 convention
(inherited from SasView)"):

> **Spec, section "The DirectLoad2D data format", paragraph "Accuracy of the q_z=0 convention (inherited from SasView)":** "Quantitatively, `|Q|/Q_⊥ = 1/cos(θ/2)`, so the wave-vector magnitude is misstated by only `θ²/8`: `3.8e-5` at `θ = 1°`, `3.4e-4` at `3°`, `3.8e-3` at `10°`. In *intensity* the error carries a second growth direction, toward the pattern edge: for a feature of relative width `s_Q` (Gaussian), `I(Q_⊥)/I(|Q|) = exp[Q_⊥²·tan²(θ/2)/(2·s_Q²)]`, i.e. it is the feature sharpness in units of `Q_⊥·tan(θ/2)` that decides. For a sharp feature `Q_⊥ = 3·s_Q` at the table edge:"

| θ | 1° | 3° | 5° | 10° | 20° | 30° |
|---|---|---|---|---|---|---|
| `\|Q\|/Q_⊥ − 1` | 3.8e-5 | 3.4e-4 | 9.5e-4 | 3.8e-3 | 1.5e-2 | 3.5e-2 |
| edge intensity error | 0.03 % | 0.31 % | 0.86 % | 3.5 % | 15 % | 38 % |

So the forward-cone intuition holds to within a few % up to `θ ~ 10°` even at
the sharpest table edge, and is essentially exact at `3°`; beyond `~20°` a
table generated for the true geometry (or the 1D description) is required.
Two caveats bound this statement. First, it is about *detector-visible
forward outcomes*: at the lowest design energy (`E = 0.1` meV,
`k = 0.22` 1/Aa) even the table edge lies beyond `90°`, where no forward
intuition applies — there the two conventions happen to agree anyway, because
every outcome samples `I ≈ I(0)`. Second, the error that is genuinely *not*
small is not a curvature correction at all but the mirror branch of the
backward-branch section below, which inflates the *total* `σ` by up to a
factor two regardless of how small `θ` is in the forward cone — a convention
artifact tracing to the dropped longitudinal component, not a curvature
effect. Whenever it matters, fall back to a 1D description or supply tables
generated for the true geometry.

### The backward branch, and what the convention does to attenuation

> **README.md:** "Rung 7 documents the one known convention effect: the 2D model's `q_z = 0` mirror extension inflates the *total* cross section of narrow patterns by up to a factor two (attenuation only; detector-visible forward scattering is unaffected)."

> **Spec, section "The DirectLoad2D data format", subsection "What the data means":** "A SasView 2D export tabulates `I(qx,qy)` on the detector plane *perpendicular to the dataset's beam axis*, fixed here to the material ẑ: `qx` along material x̂, `qy` along material ŷ. sasmodels evaluates oriented models at `q_z=0`, i.e. the table represents `I` restricted to the plane `Q·ẑ = 0` — which is precisely the coordinate pair `(Qx,Qy)` [see the specification's section "Kinematics and notation"]. The plugin evaluates for every outcome `I(kf) = I_table(Qx, Qy)` [spec eq. "eval2d"], i.e. *the same two numbers whether the outcome sits on the forward or the backward branch* of [spec eq. "branches"] — the mirror extension of [Step 2 in the specification's section "Cross sections"]. For the design beam along ẑ, `(Qx,Qy)` coincide with the detector coordinates (same sign convention `Q = kf − ki`)."

> **Spec, section "The backward branch and the relation to the 1D model", paragraph "The two conventions":** "The 1D `DirectLoad` model evaluates the table at the *full* elastic scattering-vector magnitude, `Q = kf − ki`, `Q = 2k·sin(θ/2)`, and integrates it over the sphere once: `σ_1D(k) = (2π/k²)·∫₀^{2k} I(Q)·Q dQ`, `μ = 1 − Q²/(2k²)` [spec eq. "sigma1d"]. The 2D model evaluates the table at the *transverse* components only (detector coordinates, `q_z=0` convention); since `Q_⊥` does not distinguish θ from π−θ, the two elastic branches of [spec eq. "branches"] both contribute: `σ_2D(k) = ∫ I_table(Q_⊥) dΩ = 2π·∫₀^π I_table(k·sin α)·sin α dα` [spec eq. "sigma2d"]."

> **Spec, section "The backward branch and the relation to the 1D model", paragraph "Where they differ — the mirror doubling":** "For a three-dimensional particle whose pattern core lies at `Q ≪ 2k`, a neutron scattered towards θ≈π has `|Q| ≈ 2k`, so the physical amplitude is `I(|Q|)≈0`: [spec eq. "sigma1d"] correctly gives it negligible weight. The 2D model instead sees `Q_⊥≈0` and assigns this backward outcome the *forward* intensity `I(0)` — the mirror extension of Step 2 in action. [...] This is not an implementation error on either side — both models reproduce their own conventions to 10⁻⁵ — but a consequence of the `q_z=0` convention: the 1D model computes the *physical* total cross section (validated at 0.14 % against NCrystal's builtin `@CUSTOM_HARDSPHERESANS`), while the 2D model's total cross section includes the mirrored backward lobe." (Editorial note: the elided sentence derives the doubling, [spec eq. "double"]: `σ_2D ≃ 2·σ_1D` for `Q_pattern ≪ 2k`.)

> **Spec, section "The backward branch and the relation to the 1D model", paragraph "Operational impact and options":** "Detector-visible forward scattering — intensities, anisotropy, shape of the SANS signal — is unaffected. What is inflated is the collision rate, i.e. simulated attenuation: for narrow patterns at `2k ≫ Q_pattern` the DirectLoad2D material attenuates up to twice as strongly as the physical material. Whether this matters, and which remedy to adopt if so, is a model-design decision that belongs to the data producer; the candidates are (a) keep the present convention (exact for E→0 and for ẑ-invariant structures; documented inflation otherwise), (b) restrict the mirrored branch by the elastic horizon weight `|kf·ẑ|`, which restores [spec eq. "sigma1d"] for radial patterns but breaks the `4π·I(0)` anchor, or (c) extend the data format to carry a `Q_z` dependence (a full 3D table), removing the ambiguity at the cost of format complexity and upstream data volume. The validation harness now pins the behaviour of the current choice either way."

> **Spec, section "Validation", subsection "Rung 8: transport-integrated attenuation (pass)":** "And the doubling row of [the specification's section "The backward branch and the relation to the 1D model"] becomes operational: for the *same* narrow-Gaussian particles and beam, the 1D-convention material transmits 57.9 % through one mfp of the 2D-convention material, against the 2D material's 36.8 % — up to a factor-two difference in attenuation from the convention alone, which is the number the data producer needs in order to judge whether the `q_z=0` convention is acceptable for their use case."

Details: see the specification
[doc/anisotropic_directload2d.pdf](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/doc/anisotropic_directload2d.pdf),
sections "The DirectLoad2D data format" (convention, receipts, accuracy) and
"The backward branch and the relation to the 1D model" (derivation and
options).

## Complete runnable versions

- 2D SasView chain (rung 6):
  [script/example_sasview_chain_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/example_sasview_chain_2d.py)
- Transport-integrated attenuation (rung 8):
  [script/example_transmission_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/example_transmission_2d.py)

Run them from the repo root with the plugin environment active; both print
per-check PASS/FAIL lines and exit nonzero on failure. For the 1D counterpart
and the validation ladder, see [Tutorial: 1D](tutorial-1d) and
[Validation](validation); the NCMAT section grammar is on
[Data format](data-format); return to [Home](home) for the overview.
