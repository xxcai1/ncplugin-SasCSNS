# Validation

The SasCSNS models ship with their validation: a special-case harness — the **validation ladder,
8 main rungs plus the half-rungs 4b/4c, all passing**. The doctrine, the special cases and every measured number are part of the
specification
([doc/anisotropic_directload2d.pdf](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/doc/anisotropic_directload2d.pdf),
section "Validation"); this page summarises the philosophy, gives the command and expected result
for every rung, and says what to do when a check fails. The rungs that drive the full chain are
walked through on the tutorial pages: rung 6 (1D) in
[Tutorial: 1D, "Run the complete example"](tutorial-1d#run-the-complete-example); rung 6 (2D) and
the rung 7-8 material in [Tutorial: 2D](tutorial-2d), steps 3-4 and Conventions.

## Philosophy

The specification (subsection "Doctrine") states the principle in one line: **model development
without designed validation is worthless — the validation is part of the model, delivered together
with it.** The method: find SPECIAL CASES whose parameter sets simplify the physics into something
computable in closed form, then run the model through its *completely ordinary numerical path* — it
still reads its numerical table, still interpolates, still integrates — and compare against the
hand-derived value. The model never special-cases anything; the reference always comes from
outside. Two pillars are validated separately: **Pillar A — cross section, in absolute terms**
(`sigma` decides collision probabilities; its scale is physical and is checked against closed forms
and independent quadrature), and **Pillar B — sampler, in shape only** (a sampled distribution is a
normalised density; its absolute scale is meaningless, so only its shape is tested, by binned chi^2
against the exact target in outcome-angle coordinates).

The same physical case is evaluated by four independent paths and compared pairwise: (1) closed
form, (2) continuum quadrature (no table involved), (3) table quadrature and samplers — the python
harness, a reference twin of the C++ model, (4) the plugin itself, with the same assertions
re-pointed. 1 vs 2 validates the mathematics, 2 vs 3 quantifies table discretisation, 3 vs 4
validates the implementation. Two further elements complete the design: **cross-model consistency**
(rung 7 — the two SANS models share no code or tables, so each partially checks the other) and
**regression anchors** (rung 5 — pinned values locking the released behaviour). Data files are
validated at conversion time too: the converter's gate (uniformity, `I >= 0`, monotone axes,
coverage `Qmax >= 2k(Emax)` — the error that historically cost a silent factor 2-14) is described
on [Data format](data-format), with a truncation control case in rung 6 (1D).

## The ladder at a glance

Run from the repo root with the NCrystal+plugin environment active (`python` must see the plugin,
i.e. `pip install .` first — see [Installation](installation)):

```bash
python script/validate_directload2d.py    # rungs 2-3: python model vs closed forms (29 checks)
python script/validate_plugin_2d.py       # rung 4: C++ plugin vs reference (32 checks)
python script/validate_hard_sphere_equiv.py  # rung 4b: constant table = hard sphere (15 checks)
python script/validate_physical_closedforms.py  # rung 4c: Guinier + Zimm/OZ closed forms (30 checks)
python script/validate_1d_vs_2d.py        # rung 7: cross-model consistency + conventions (7 checks)
python script/example_sasview_chain.py    # rung 6 (1D): full SasView chain
python script/example_sasview_chain_2d.py # rung 6 (2D): full SasView chain + figure data
python script/example_transmission_2d.py  # rung 8: transport-integrated attenuation
ncrystal-pluginmanager --test SasCSNS     # rung 5: regression anchors
python script/fig_chain_2d.py             # render doc/fig_chain_2d.pdf
python script/fig_conventions.py          # render doc/fig_conventions.pdf
```

| rung | what is tested | command | expected result (numbers as measured in the specification) |
|---|---|---|---|
| 1 | closed forms, hand-derived for special parameter sets | (expected values inside `script/validate_directload2d.py`) | S1 `sigma = 4*pi*I0 = 37.6991` and the S8 linear ramp: machine precision, rel. err. `2.2e-14` |
| 2-3 | continuum quadrature and the python reference twin | `python script/validate_directload2d.py` | 29 checks, 0 failures; table discretisation `1.8e-5` (S5, 401^2); tilted-beam acceptance `0.676` at 50 mrad |
| 4 | the C++ plugin, same assertions re-pointed | `python script/validate_plugin_2d.py` | 32 checks, 0 failures; worst `sigma` error `1.0e-2`, median `8.2e-4` (S9 fuzz, 150 checks); all shape `p >= 0.04` |
| 4b | the hard-sphere equivalence, plugin only (no python twin) | `python script/validate_hard_sphere_equiv.py` | 15 checks, 0 failures; `sigma = 4*pi*I0` machine-exact (`1e-13`) on axis and in the `E->0` clamp, `1e-4` tilted rows; `1D == 2D`; outcomes uniform on the sphere |
| 4c | real physics with closed-form sigma (Guinier, Zimm/OZ), both models | `python script/validate_physical_closedforms.py` | 30 checks, 0 failures; 2D worst `1.2e-3`, 1D worst `4e-6` (eq:isoxs to quadrature accuracy); mirror ratio from the closed forms to `2e-4` |
| 5 | regression anchors pinned in `src/NCTestPlugin.cc` | `ncrystal-pluginmanager --test SasCSNS` | ends `All tests of plugin were successful!`; anchor tolerances `1e-4` / `2e-4` |
| 6 | the full SasView chain, 1D and 2D | `python script/example_sasview_chain.py` and `python script/example_sasview_chain_2d.py` | 1D: worst `0.14%` against the builtin HardSphereSANS over 0.5-20 meV; 2D: worst `1.2e-2` against an independent quadrature |
| 7 | cross-model consistency + conventions | `python script/validate_1d_vs_2d.py` | 7 checks, 0 failures; mirror ratio `2.026` (`E = 2.07` meV) and `2.014` (`E = 5` meV), band `[1.8, 2.3]` |
| 8 | transport-integrated attenuation (Beer-Lambert) | `python script/example_transmission_2d.py` | 4 checks (`1e6` neutrons each), all pass; pulls within `1.9 sigma` |
| — | figures for the specification | `python script/fig_chain_2d.py`, `python script/fig_conventions.py` | renders `doc/fig_chain_2d.pdf` and `doc/fig_conventions.pdf` |

The three `validate_*.py` harnesses print one row per check and a final `ALL PASS` /
`N FAILURE(S)` line; the rung 4b/4c and chain scripts print per-check PASS/FAIL rows and
`ALL CHECKS PASS` / `N check(s) FAILED` to the same effect. All exit nonzero on any failure —
so the ladder can run unattended in CI.

## Rungs 1-3: closed forms, quadrature, and the python twin

[script/validate_directload2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/validate_directload2d.py)
is pure python (no plugin needed) and re-derives the model twice: once by continuum quadrature of
analytic references, once as a bilinear table twin of the C++ code. Its special cases, and what
each locks down:

- **S1, constant table**: `sigma = 4*pi*I0`, outcomes uniform on the sphere — the overall
  normalisation and the full-sphere (two-branch) convention; a one-sided model would give
  `2*pi*I0` and is excluded, at machine precision (`2.2e-14`).
- **S2, `E -> 0`**: `sigma -> 4*pi*I(0)` for any table and any beam (`5.1e-5` on a tilted beam) —
  the elastic-kinematics prefactor `1/(k^2 |kf . z|)` and the anchor shared with the 1D model.
- **S4, single pixel**: outcomes cluster at the two elastic branch points (measured fraction
  `1.000`); the two branches and their Jacobian weighting.
- **S5, Gaussian**: no closed form, so three *independent* evaluations must agree: continuum vs
  table `1.8e-5` (the discretisation), table 401^2 vs 801^2 `1.4e-5`, sampling-grid vs quadrature
  `4.8e-4`; tilted-beam `sigma` via table-rotation invariance.
- **S6, Airy rings**: `I = [2*J1(q*R)/(q*R)]^2` — sharp ring structure through both samplers
  (continuum vs table `5.1e-5`), psi uniformity.
- **S7, conventions**: a blob at `+Qx` scatters to `+Qx` (sign convention `Q = kf - ki`);
  `sigma(+z) = sigma(-z)` to below `1e-12`.
- **S8, `E -> 0` expansion**: exact for a linear ramp at any `k` — the linear tilt coefficient to
  machine precision (`2.2e-14`, on-axis and at 30/-30/200 mrad), the curvature coefficient to
  `1.3e-4`.

Total: **29 checks, 0 failures**. The harness paid for itself: it caught a factor-2 bug in the
reference quadrature itself (S1's machine-precision anchor), excluded the one-sided convention,
showed that plane-coordinate sampling cannot be exact at the elastic horizon (the replacement was
designed before any C++ existed), and exposed a Gauss-Legendre reference rule that overestimated a
ring-structured table by about 12x. References must be validated too.

## Rung 4: the C++ plugin

[script/validate_plugin_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/validate_plugin_2d.py)
exports the same special cases to NCMAT files, loads them through the real plugin factory and
re-points the identical assertions at the plugin's `crossSection` / `sampleScatter`: **32 checks,
0 failures**. The absolute scale is set by the S9 fuzz envelope — 30 random smooth tables times 5
random kinematics each: median error `8.2e-4`, worst `1.0e-2` over 150 checks, exactly the
documented grid-interpolation envelope — and every shape test passes with `p >= 0.04`. This rung
caught the grid-resolution failure modes (nearest-node shape `p ~ 1e-14`; ~1% `sigma` curvature
errors from a too-coarse energy grid), fixed before release.

## Rungs 4b and 4c: closed forms straight through the plugin

Two half-rungs validate the integration method itself, with no python twin in the loop.

**Rung 4b — the hard-sphere equivalence**
([script/validate_hard_sphere_equiv.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/validate_hard_sphere_equiv.py)):
a constant table `I(Qx,Qy) = I0` *is* an isotropic hard sphere — `sigma = 4*pi*I0` exactly, for
every `k` and every beam direction. 15 checks: `sigma` machine-exact (`1e-13`) on axis at four
energies and in the `E->0` clamp, `1e-4` on the tilted rows (the s-grid quadrature floor),
`1D == 2D`, the table-edge truncation row, and outcomes uniform on the sphere (chi2 `p = 0.48/0.69`
at the 500k-event convention).

**Rung 4c — physical closed forms**
([script/validate_physical_closedforms.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/validate_physical_closedforms.py)):
two real SANS physics models are closed form for BOTH DirectLoad models at all `k` — the Guinier
law (`erfi` / elementary forms) and Zimm/Ornstein-Zernike (`atanh` / `ln` forms); the formulas are
in the specification, section "Rung 4c". 30 checks, after the closed forms are themselves
machine-checked against independent quadratures: 2D worst `1.2e-3`, 1D worst `4e-6`, `E->0` node
clamp `1.4e-6`, and the mirror ratio `sigma_2D/sigma_1D` (1.11-1.26 for these wide patterns,
continuously connecting the constant-table value 1 to the narrow-pattern limit 2) reproduced from
the closed forms to `2e-4`. Rung 4c also pins a measured model detail: a 1D file violating the
coverage rule (`2k > qmax`) serves a legacy `sigma ~ E^-1/2` lookup extrapolation above the
boundary energy, not the `I=0` truncation — documented in the specification, section "Behaviour
outside the table".

## Rung 5: regression anchors

`ncrystal-pluginmanager --test SasCSNS` pins fixed-value anchors into
[src/NCTestPlugin.cc](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/src/NCTestPlugin.cc):
the constant-table `sigma = 4*pi*I0` on axis, at 30 mrad and at a perpendicular beam (tolerances
`1e-4` / `2e-4`); exact elasticity and unit outcome norms; the `sans=0` bypass; and a peaked table
with pinned `sigma` values at two energies plus a forward-scattering assertion at `k >>` table
extent (either elastic branch, `|mu| > 0.9`). A successful run ends with `All tests of plugin were
successful!` (full expected output on [Installation](installation)). These anchors make later
refactors comparable against the released behaviour.

## Rung 6: the full SasView chain

- **1D** —
  [script/example_sasview_chain.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/example_sasview_chain.py)
  (walkthrough: [Tutorial: 1D](tutorial-1d)): emulated SasView export, conversion, and the
  macroscopic cross section `Sigma = n_d * sigma` compared against NCrystal's builtin
  `@CUSTOM_HARDSPHERESANS` for the same physical system at six energies over 0.5-20 meV: worst
  agreement **0.14%**; the same script runs the truncated-table control (biases of factors 2 up
  to 14).
- **2D** —
  [script/example_sasview_chain_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/example_sasview_chain_2d.py)
  (walkthrough: [Tutorial: 2D](tutorial-2d), step 3): emulated 2D export of a tilted cylinder,
  conversion with `--2d`, then plugin `sigma` against an independent midpoint quadrature of the
  same file table: worst relative error `1.2e-2`; sampled outcomes reproduce the elongated image
  (`<Qx^2>/<Qy^2> = 0.47`); where `sasmodels` is installed, the emulation is verified against the
  real kernel (shapes to `4e-16`).

## Rung 7: cross-model consistency (the finding)

[script/validate_1d_vs_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/validate_1d_vs_2d.py)
checks the two independent SANS models against each other in the small-angle regime (narrow
Gaussian, `s_Q = 0.1` 1/Aa, 40 pixels per `s_Q`), where their conventions must agree on the forward
branch. Seven checks, all passing: each model's `sigma` against its own continuum (`9.4e-6` on both
1D rows; `2.7e-3` twice for 2D), the mirror ratio in the accepted `[1.8, 2.3]` band
(`2.026` at `E = 2.07` meV, `2.014` at `E = 5` meV), and design-cone shape agreement (max relative
bin deviation `0.039` — statistical noise — at `2e5` events). The doubling rows are the campaign's
headline *finding*, not a failure (see below); the shape row had earlier exposed an angular-grid
aliasing and driven its fix.

## Rung 8: transport-integrated attenuation

[script/example_transmission_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/example_transmission_2d.py)
exercises Pillar A the way a transport code sees it: analog Monte Carlo through a slab, free paths
drawn from `Exp(1/Sigma)`, the uncollided fraction must reproduce Beer-Lambert. Four checks at
`1e6` neutrons each, all passing: uncollided fraction at `L = 0.25/1/3` mfp (pulls
`0.1/0.9/1.7 sigma`), absorbed fraction (`+1.9 sigma`), a 30 mrad tilted beam (`+1.7 sigma`;
`sigma` rises by `1.0%` at 30 mrad for this pattern), and the mirror doubling (`2.026`, in
`[1.8, 2.3]`). This pins the unit chain barn/atom -> cm^-1 end to end — `n [atoms/Aa^3] *
sigma [barn/atom]` *is* `Sigma` in 1/cm, no factor at all.

## Figures for the specification

- [script/fig_chain_2d.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/fig_chain_2d.py)
  renders `doc/fig_chain_2d.pdf` — rung 6 at the shape level: the emulated input image
  `I(Qx,Qy)` (tilted cylinder) next to the density of 2M outcomes sampled from the C++ plugin, on
  a common colour scale with the input contours overlaid; auto-runs the 2D chain script when its
  NCMAT is missing.
- [script/fig_conventions.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/fig_conventions.py)
  renders `doc/fig_conventions.pdf` — the picture of the rung 7 finding: left, the sampled
  `p(theta)` of the same narrow-Gaussian particles in both models, with the 2D model's backward
  lobe carrying mirrored forward intensity (the `sigma` doubling); right, the Beer-Lambert
  transmission of a 1D- vs a 2D-convention material — the operational consequence.

Both need `matplotlib` (and `scipy`); they render the specification's figures — the numbers live
in the harness scripts above.

## The known convention effect: q_z = 0 mirror extension

Rung 7 documents the one known convention effect of the 2D model. Quoted verbatim from the
specification, section "The backward branch and the relation to the 1D model" (`[...]` marks elided text):

> "The 2D model instead sees `Q_perp ~ 0` and assigns this backward outcome the *forward* intensity
> `I(0)` — the mirror extension of Step 2 in action. [...] and `sigma_2D ~= 2 * sigma_1D`
> (`Q_pattern << 2k`). [...] This is not an implementation error on either side — both models reproduce
> their own conventions to `10^-5` — but a consequence of the `q_z = 0` convention: the 1D model
> computes the *physical* total cross section (validated at `0.14%` against NCrystal's builtin
> `@CUSTOM_HARDSPHERESANS`), while the 2D model's total cross section includes the mirrored
> backward lobe."

> "Detector-visible forward scattering — intensities, anisotropy, shape of the SANS signal — is
> unaffected. What is inflated is the collision rate, i.e. simulated attenuation: for narrow
> patterns at `2k >> Q_pattern` the DirectLoad2D material attenuates up to twice as strongly as
> the physical material."

So: the *total* cross section of narrow patterns can be inflated by up to a factor two (rung 7
measures `2.026` / `2.014`; rung 8 makes it operational — through one mfp, the 1D-convention
material transmits `57.9%` where the 2D-convention material transmits `36.8%`), while
detector-visible forward scattering is unaffected. Whether that matters is a data-producer
decision; the full derivation and the options are in the specification section quoted above and
summarised on [Tutorial: 2D](tutorial-2d) (Conventions).

## If a check fails

1. **Check the environment first.** A stale or missing plugin build is the most common cause:
   re-run `pip install .`, confirm `ncrystal-pluginmanager --test SasCSNS` still passes, and make
   sure the `python` running the scripts is the one the plugin was installed into (see
   [Installation](installation), troubleshooting).
2. **Compare the failing row with the specification** — section "Validation" of the
   [PDF](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/doc/anisotropic_directload2d.pdf)
   reproduces every check table with expected value and tolerance.
3. **Read the failure kind.** Pillar A rows (cross sections vs closed forms or quadrature) are
   deterministic: a failure is a real regression. Pillar B rows (shape p-values, binomial pulls)
   carry documented noise margins (`p > 0.001`, `< 0.08`, `4 sigma`); a marginal breach should
   reproduce on re-run before it is called one. The python harnesses seed their own random draws,
   so two runs of the same script are directly comparable.
4. **Report with the printed row** — it carries the check name, expected value, measured value and
   error; quote it (with the NCrystal and plugin versions) when filing an issue.

Back to [Home](home).
