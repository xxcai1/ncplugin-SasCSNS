# Data format: the `@CUSTOM_SASCSNS` NCMAT section

This page is the reference for the `@CUSTOM_SASCSNS` custom section, by which an NCMAT file selects one of the SasCSNS SANS models. Worked usage lives in [Tutorial: 1D models](tutorial-1d) and [Tutorial: 2D models](tutorial-2d); how the models are verified is in [Validation](validation).

| model type | input | character |
|---|---|---|
| `DirectLoad` | tabulated 1D `I(Q)` | isotropic; full elastic sphere, `Q = 2k sin(theta/2)`; physical total cross section (README: validated against the NCrystal builtin HardSphereSANS at 0.14%) |
| `DirectLoad2D` | tabulated 2D `I(Qx,Qy)` | anisotropic pattern in detector coordinates; elastic |
| `HardSphere` | analytic `radius`, contrast | analytic reference model |

The plugin develops no new physics models: SasView (or any other calculator) evaluates the form factor and exports a tabulated intensity; the plugin transports neutrons through it numerically. Most files are produced by the converter [script/sasview2ncmat.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/sasview2ncmat.py), which always writes `DirectLoad` (1D) or `DirectLoad2D` (2D). `HardSphere` files are hand-written; the shipped example is [data/sascsns_sio2_spheres.ncmat](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/data/sascsns_sio2_spheres.ncmat). Two shipped `DirectLoad2D` examples carry real physics with a closed-form total cross section (doc section "Rung 4c"): [data/sascsns_guinier2d.ncmat](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/data/sascsns_guinier2d.ncmat) (Guinier, Rg = 1.0/Aa) and [data/sascsns_oz2d.ncmat](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/data/sascsns_oz2d.ncmat) (Zimm/OZ, xi = 0.7/Aa), each loading as `plugins::SasCSNS/<name>.ncmat` after install.

## Where the section goes in the file

`@CUSTOM_SASCSNS` goes after the standard sections, as in every example below. Only the **first** `@CUSTOM_SASCSNS` block of a file is read. Loading requires the plugin to be installed ([Installation](installation)) and honours two gates implemented in [src/NCPluginFactory.cc](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/src/NCPluginFactory.cc) (factory priority 999):

- The config parameter `sans` must not be disabled; pass `sans=0` to fall back to the standard NCrystal models.
- The file must be single-phase: files with `@OTHERPHASES` are declined here and served by NCrystal's builtin multi-phase SANS factories.

Dispatch order: (1) any line whose first token is `DirectLoad2D` -> anisotropic `DirectLoad2D`; (2) otherwise the first keyword line in file order decides — first token `DirectLoad` -> isotropic `DirectLoad`, `HardSphere` -> analytic `HardSphere` (the picker loop in `src/NCSansModelPicker.cc` scans the section lines in order); (3) anything else -> `BadInput` with the message `" @CUSTOM_SASCSNS with undefined load method"` (the leading space is part of the code literal in `src/NCSansModelPicker.cc`). Quick load check from Python (all usage snippets on this page use the Python API):

```python
import NCrystal
proc = NCrystal.createLoadedMaterial('plugins::SasCSNS/sascsns_sio2_spheres.ncmat')
```

The `plugins::SasCSNS/...` path resolves to the data file bundled with the installed plugin (the same smoke test as in [Installation](installation); a repo-relative path like `data/...` would only resolve with the process running from the repo root). `createLoadedMaterial` is also exported under its convenience alias `NCrystal.load`.

## Type `DirectLoad` — tabulated 1D I(Q)

Two keyword lines, located by their first token:

| line | content | unit |
|---|---|---|
| `Q <q1> <q2> ...` | scattering-vector magnitudes | 1/angstrom |
| `I <i1> <i2> ...` | intensity per atom per steradian (after conversion, see [Units](#units-barnatom-not-barnobject)) | barn/(atom sr) |

Parsing notes ([src/NCSansModelPicker.cc](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/src/NCSansModelPicker.cc)):

- Both lines are mandatory (missing one raises `findCustomLineIter can not find parameter <keyword>`); a non-numeric token raises `Invalid values specified in the @CUSTOM_SASCSNS section Q` / `... section I`.
- **Keep `Q` and `I` exactly the same length**: the parser loops over the `Q` tokens while indexing the `I` line, so an `I` line **shorter** than `Q` surfaces as an out-of-range exception, not a friendly error message, while a **longer** `I` line is silently truncated (extra values are never read).
- No minimum point count in the parser (the converter requires >= 2); other lines in the section (e.g. `solvent`) are ignored here.

Minimal example (constructed in exactly the format the converter emits; the `Q`/`I` values are placeholders, not physics data — produce real files with `sasview2ncmat.py`):

```text
NCMAT v6
@DENSITY
  2.2 g_per_cm3
@DYNINFO
  element Si
  fraction 0.3333333
  type freegas
@DYNINFO
  element O
  fraction 0.6666667
  type freegas
@CUSTOM_SASCSNS
  DirectLoad
  Q 0.001 0.002 0.004 0.008
  I 1e-03 8e-04 3.2e-04 5.6e-05
```

## Type `HardSphere` — analytic reference model

| line | content | unit |
|---|---|---|
| `radius <R>` | sphere radius; exactly one parameter | angstrom |
| `solvent <cfgstring> <fraction>` | optional; exactly two parameters; `fraction` = solvent volume fraction | dimensionless |

The contrast is not in the section: SLD and number density come from the rest of the file (`@DENSITY` + `@DYNINFO`, or a crystal structure with atom positions). I(Q) is generated internally on a fixed `Q = logspace(-6, 10, 1000)` grid (1/angstrom) from the solid-sphere form factor `P = 3(sin(qR) - qR cos(qR))/(R^3 q^3)`, with `I = V^2 sld^2 P^2 / (V n_d)`, `V = 4/3 pi R^3`; points with `q < 1e-5` are replaced by the `q -> 0` limit.

The shipped example [data/sascsns_sio2_spheres.ncmat](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/data/sascsns_sio2_spheres.ncmat) (10 nm-radius amorphous silica spheres), comments trimmed:

```text
NCMAT v6
@DENSITY
  2.2 g_per_cm3
@DYNINFO
  element Si
  fraction 0.3333333
  type freegas
@DYNINFO
  element O
  fraction 0.6666667
  type freegas
@CUSTOM_SASCSNS
  HardSphere
  radius 100
```

Caveats: a file with neither `@DENSITY`+`@DYNINFO` nor atom positions never reaches the plugin — NCrystal itself rejects it with `BadInput: missing information about both unit cell and material dynamics` (`@DYNINFO` without `@DENSITY`, and a `@CELL` without atom positions, are likewise rejected before the plugin runs). The wrong-parameter-count error messages are mislabelled in the code — "radius in the @CUSTOM_SASCSNS radius field should prove one parameter" and "radius in the @CUSTOM_SASCSNS solvent field should prove two parameters" (the latter checks the *solvent* line) — quote them as-is when debugging.

## Type `DirectLoad2D` — tabulated anisotropic I(Qx,Qy)

The grammar, quoted verbatim from the specification ([doc/anisotropic_directload2d.tex](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/doc/anisotropic_directload2d.tex), section "The DirectLoad2D data format", subsection "NCMAT grammar"):

> ```text
> @CUSTOM_SASCSNS
>   DirectLoad2D
>   Qx <x_1> <x_2> ... <x_Nx>     # Nx ascending q_x values [1/Aa]; may be negative
>   Qy <y_1> <y_2> ... <y_Ny>     # Ny ascending q_y values [1/Aa]
>   I  <I_1> <I_2> ...            # Nx*Ny values [barn/(atom sr)], ROW-MAJOR:
>   <continuation lines allowed>  # I[qy_index*Nx + qx_index] (qy outer loop)
> ```

Hard constraints enforced by [src/NCDirectLoad2D.cc](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/src/NCDirectLoad2D.cc) (violations are `BadInput` parse errors, message text in brackets):

- `Qx`, `Qy` and `I` lines are all mandatory ["DirectLoad2D section must contain Qx, Qy and I lines"].
- Each axis needs >= 2 values ["... axis needs at least 2 values"], strictly ascending ["... axis must be strictly ascending"], and uniformly spaced within relative tolerance `1e-6 * |d0|` ["... axis must be uniformly spaced"] — uniformity is what enables the O(1) bilinear lookup, and why the converter writes axes with 12 significant digits (fewer digits can shift consecutive pixel centres right onto the tolerance on fine grids).
- I values may start on the `I` line and wrap freely over following lines; the block ends at the first line whose first token is `Qx`, `Qy`, `I` or `DirectLoad2D`. The total must be exactly `Nx*Ny` ["DirectLoad2D I section must contain exactly N values (Nx x Ny), got M"], every value >= 0 ["DirectLoad2D I values must be non-negative"], and not all zero ["DirectLoad2D I values are all zero"].
- Values are stored row-major with `qy` as the outer axis, evaluated bilinearly, and **I = 0 outside the tabulated grid**.

Minimal example (2x2 grid — the smallest the parser accepts; physically useless, the converter warns if either axis has fewer than 8 points; `I` values are placeholders). Row-major order means row 1 is `(I(Qx1,Qy1), I(Qx2,Qy1))`:

```text
NCMAT v7
@DENSITY
  2.2 g_per_cm3
@DYNINFO
  element Si
  fraction 0.3333333
  type freegas
@DYNINFO
  element O
  fraction 0.6666667
  type freegas
@CUSTOM_SASCSNS
  DirectLoad2D
  Qx 0 0.01
  Qy 0 0.01
  I 1 0.5
    0.5 0.1
```

(All four values on the single `I` line is equally valid — wrapping is free. The converter emits NCMAT **v7** for 2D files, v6 for 1D, so DirectLoad2D files need a v7-capable NCrystal.)

### The q_z = 0 convention of the table coordinates

What the convention assumes and how accurate it is, quoted verbatim from the specification (doc/anisotropic_directload2d.tex, section "The DirectLoad2D data format", subsection "What the data means"; `[...]` marks elided cross-references):

> sasmodels evaluates oriented models at `q_z=0`, i.e. the table represents `I` restricted to the plane `Q . z_hat = 0` `[...]`. The plugin evaluates for every outcome `I(k_f) = I_table(Qx, Qy)`, i.e. *the same two numbers whether the outcome sits on the forward or the backward branch* `[...]`.

> The neglected component is `Q_z = -2k sin^2(theta/2)`, i.e. `|Q_z|/Q_perp = tan(theta/2)`: a second-order-in-theta error, negligible for typical SANS (`theta <~ 3 degrees`, particles up to a few hundred angstrom) and *not introduced by the plugin* — it is SasView's own 2D convention. `[...]` Whenever it matters, the user should fall back to a 1D description or supply tables generated for the true geometry.

Details: see the specification, section "The backward branch and the relation to the 1D model". The convention's provenance (the sasmodels `R^-1 [qx, qy, 0]^T` receipts) and its quantitative accuracy table are covered on [Tutorial: 2D models](tutorial-2d) (Conventions).

## Units: barn/atom, not barn/object

SasView reports intensity per **object** (particle), in barn/object, while the values in this section must be per **atom** (barn/(atom sr)). The converter applies the conversion **always and automatically — there is deliberately no flag for it**:

```text
scale    = (phi / V_p) / n_d          # n_objects / n_atoms
n_d      = density * AVOGADRO * atoms_batch / (mass_batch * 1.0e24)   # 1/angstrom^3
I_ncmat  = I_sasview * scale
```

Here `phi` is the particle volume fraction (`--phi`, mandatory, in (0,1)), `V_p` the particle volume (`--radius` in angstrom with `V_p = 4/3 pi R^3`, or `--volume` in angstrom^3), and `n_d` the atomic number density from `--material` (presets: sio2, al2o3, fe3o4, air) and `--density` (g/cm3). The applied scale is recorded in the output file header. The conversion is easy to get wrong by hand: a form-factor scale `V_p^2 dRho^2` with `dRho` in 1/angstrom^2 already has dimension angstrom^2, and **1 angstrom^2 = 1e8 barn**.

When comparing files (different `phi`, different batching), compare the macroscopic cross section `Sigma = n_d * sigma`, never the per-atom `sigma`: `Sigma` is the invariant quantity transport actually uses.

## Q-range and coverage requirements

Kinematic relations behind all checks here: `E[meV] = 2.07214 * k^2` (k in 1/angstrom), elastic scattering, largest reachable momentum transfer `Q = 2k`.

- **1D coverage rule:** the `I(Q)` table must span `[0, 2k(E_max)]` for the highest neutron energy `E_max` of the use case. The converter prints `File validation: Qmax = ... 1/Aa -> fully covers neutron energies up to E = ... meV (Q = 2k rule)` and, given `--emax`, **aborts** if `Qmax < 2k(Emax)`. A truncated table makes the plugin extrapolate outside it — historically the silent factor 2-14 error — so extend the table or lower `--emax`. A table starting above 5% of `Qmax` triggers a WARNING (the forward SANS intensity peaks at `Q = 0`).
- **2D coverage rule:** each half-side must cover `2k` (checked as `qmax = min(qx[-1], qy[-1])`; the disc of radius `2k` must fit inside the table). Export grids should start at `Q = 0`: without a `Q = 0` pixel the converter warns that the forward SANS intensity is amputated and interpolated from the first pixel inward.
- **Resolution rule:** use at least 10 points per form-factor oscillation period (period `pi/R`), i.e. spacing `<= pi/(10 R)`; the 1D converter checks this whenever `--radius` is given and warns if the minimum spacing is too coarse (it never runs in `--2d` mode).
- Outside its grid DirectLoad2D uses `I = 0`, so a truncated table stays internally consistent (truncated angles contribute nothing to sigma and are never sampled) — but the missing forward intensity is real physics loss. At load time the plugin warns if the table half-width does not even cover `2k(0.1 meV)`; runtime queries outside the calibrated energy window [0.1, 100] meV get boundary-node values plus a one-time warning.

## Converter CLI: `sasview2ncmat.py`

Most `@CUSTOM_SASCSNS` files are produced by [script/sasview2ncmat.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/sasview2ncmat.py). This table is the single reference for its command line (the tutorials link here instead of repeating it):

| flag | meaning | default |
|---|---|---|
| `sasview_file` | input SasView ASCII export (positional); 1D rows `Q I [dI] [dQ]`, 2D rows `Qx Qy I` | required |
| `--2d` | select the 2D `I(Qx,Qy)` path -> `DirectLoad2D` model; without it the 1D `DirectLoad` pipeline runs (1D is the default — there is no `--1d` flag) | off |
| `-o`, `--output` | output NCMAT file name | `sas_model.ncmat` |
| `--material` | material preset: `sio2`, `al2o3`, `fe3o4`, `air` | `sio2` |
| `--density` | particle material density [g/cm3] | `2.2` |
| `--phi` | particle volume fraction in the sample; must satisfy 0 < phi < 1 | required |
| `--radius` | sphere radius [Angstrom], `V_p = 4/3 pi R^3`; mutually exclusive with `--volume` — exactly one of the two is required, since the unit conversion needs `V_p` | — |
| `--volume` | particle volume [Angstrom^3], for non-spherical particles | — |
| `--emax` | highest neutron energy [meV] the file must cover; aborts if the table does not cover it (see [Q-range and coverage requirements](#q-range-and-coverage-requirements)) | off |
| `--solvent` | solvent name recorded for reference (the `DirectLoad` I(Q) already includes the solvent contrast); recorded as a trailing comment line `  # solvent: NAME`, which the NCMAT parser strips before the section is tokenised, so combining it with `--2d` is safe | none |

## Sanity checks performed

By the converter [script/sasview2ncmat.py](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/sasview2ncmat.py) (1D is the default mode — there is no `--1d` flag, only `--2d`, see [Converter CLI](#converter-cli-sasview2ncmatpy)):

- 1D: at least 2 (Q,I) points (hard error); drops `Q<=0`, `I<=0` or NaN points (warning; `+inf` passes this filter); sorts ascending in Q and drops Q duplicates within 1e-12 (warning).
- 2D: at least 4 (Qx,Qy,I) rows; rows must form a complete rectangular grid (hard error — fill masked pixels first); non-finite I is a hard error; negative I is clamped to 0 (warning); axes need >= 2 strictly ascending points and are resampled bilinearly onto a uniform grid if non-uniform (warning); warns if either axis has fewer than 8 points (the check is `min(len(Qx), len(Qy)) < 8`, so a 6x200 grid triggers it too), and warns if the grid does not start at `Q = 0` (the forward SANS intensity is then amputated — see [Q-range and coverage requirements](#q-range-and-coverage-requirements)).
- Both modes enforce `--emax` (hard error) and print the applied scale. The `File validation: Qmax = ...` coverage line quoted above is the 1D form; `--2d` prints its own coverage summary (`File validation: table Nx x Ny over [...]x[...] 1/Aa; fully covers neutron energies up to E = ... meV (disc radius 2k must fit inside the table)`).

Plugin-side: the per-type `BadInput` constraint lists above, plus one-time runtime warnings (DirectLoad2D energy-window and beam-tilt clamp notices, low-energy table half-width check).

## Pitfalls

- `--solvent` is harmless in every mode: the converter records the name as a true comment line `  # solvent: NAME` (verified by round-trip through NCrystal in both modes). Hand-written live `solvent` lines remain dangerous: the `DirectLoad2D` parser does **not** stop the I block on `solvent`, so it would swallow such a line as an extra I value and reject the file. `HardSphere` expects the two-parameter form `solvent <cfgstring> <volume_fraction>`.
- `HardSphere` needs `@DENSITY`+`@DYNINFO` (or a crystal structure); files without them are rejected outright by NCrystal before the plugin is consulted.

## See also

[Home](home) | [Installation](installation) | [Tutorial: 1D models](tutorial-1d) | [Tutorial: 2D models](tutorial-2d) | [Validation](validation)
