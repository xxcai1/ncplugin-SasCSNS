# NCrystal plugin SasCSNS

**SasCSNS** is an NCrystal plugin providing general-purpose **small-angle
neutron scattering (SANS) models** for particle transport.
Models are declared in a `@CUSTOM_SASCSNS` section of an
ordinary NCMAT file; the material is then used like any other NCrystal
material through the normal NCrystal Python API (`import NCrystal`, then e.g.
`NCrystal.createInfo` / `NCrystal.createScatter` / `NCrystal.createLoadedMaterial`
(the latter also available as `NCrystal.load`)).

> **Full specification (data format, mathematics, sampling, validation):**
> [doc/anisotropic_directload2d.pdf](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/doc/anisotropic_directload2d.pdf)

## The DirectLoad philosophy

The specification, section "Scope and design principle", states the design in
one line:

> **SasView computes physics; the plugin transports it.** A single new model
> type, `DirectLoad2D`, ingests a generic tabulated 2D intensity
> I(q_x,q_y); everything shape-specific happens upstream, in SasView, and
> in the converter tool `sasview2ncmat.py`.

That is, the plugin develops **no new physics models**: SasView (or any other
calculator) evaluates the form factor and exports a tabulated intensity, and
the plugin transports neutrons through that table numerically.

## Model types (`@CUSTOM_SASCSNS` NCMAT section)

| type | input | notes |
|---|---|---|
| `DirectLoad` | tabulated 1D `I(Q)` | full elastic sphere, `Q = 2k sin(theta/2)`; physical total cross section (validated against builtin HardSphereSANS at 0.14%) |
| `DirectLoad2D` | tabulated 2D `I(Qx,Qy)` | anisotropic patterns in detector coordinates (`q_z = 0` convention); on-axis sampling exact and rejection-free; tilted beams via dilated-proposal rejection; see [doc/anisotropic_directload2d.pdf](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/doc/anisotropic_directload2d.pdf) for the full specification |
| `HardSphere` | analytic `radius`, contrast | analytic reference model |

The `q_z = 0` convention of the `DirectLoad2D` row carries one known
convention effect, quantified in the specification, section "The backward
branch and the relation to the 1D model", and summarised in the
[README](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/README.md);
the convention's provenance (the sasmodels `R^-1 [qx, qy, 0]^T` receipts) and
its quantitative accuracy are covered on
[Tutorial: 2D models](tutorial-2d) (Conventions).

## Quickstart

Run from the repo root with the NCrystal+plugin environment active (`python`
must see the plugin, i.e. `pip install .` first):

```bash
pip install .                                              # 1. install the plugin
python script/sasview2ncmat.py EXPORT.dat --2d \
    --material sio2 --density 2.2 --phi 0.4 --radius 100 \
    -o sans_model.ncmat                                    # 2. SasView export -> NCMAT
python -c "import NCrystal; print(NCrystal.createInfo('sans_model.ncmat'))"   # 3. load (smoke test)
```

1. **Install** the plugin: [Installation](installation). The converter is the
   repo script `script/sasview2ncmat.py`, run directly in step 2.
2. **Convert** a SasView ASCII export (1D by default; `--2d` selects the 2D
   `I(Qx,Qy)` path) into an NCMAT
   file: converter behavior and the `@CUSTOM_SASCSNS` grammar are described in
   [Data format](data-format); end-to-end walkthroughs in
   [Tutorial: 1D models](tutorial-1d) and [Tutorial: 2D models](tutorial-2d).
3. **Load**: the resulting NCMAT file loads through the NCrystal Python
   API exactly like any other material.

## Units and data requirements

- SasView intensities are per **object** (`barn/object`), NCrystal needs per
  **atom** (`barn/(atom sr)`); the converter `script/sasview2ncmat.py` applies
  this conversion **always and automatically**, with no scale option
  (normative copy: [Data format, Units](data-format#units-barnatom-not-barnobject)).
- The tabulated intensity must cover the scattering vectors reached at the
  highest simulated energy, i.e. up to `2k(E_max)`. The converter aborts with
  a hard error if `--emax` is given and the table does not cover it (Qmax <
  `2k(E_max)`); otherwise it reports the highest fully covered energy. In the
  2D path (`--2d`), non-uniform axes are resampled onto uniform grids and
  negative `I` values are clamped to 0, each with a warning; in the 1D path,
  points with negative `I` are dropped (with a warning) and a non-uniform 1D
  `Q` grid is accepted as-is. For `DirectLoad2D`, the
  uniform-grid, `I >= 0` and strictly-ascending-axis rules are hard parse
  errors when the NCMAT file is loaded; the 1D `DirectLoad` parser enforces
  only the monotone axis (via NCrystal's `PointwiseDist`, which requires
  increasing `Q` values).

## Status

The models are validated by a special-case harness — the **validation ladder,
8 main rungs plus the half-rungs 4b/4c, all passing**: closed forms, continuum
quadrature, a python twin, the C++ plugin (plus the hard-sphere-equivalence and
physical-closed-form half-rungs), regression anchors, the full SasView chain,
cross-model consistency, and transport-integrated attenuation (Beer-Lambert).
The doctrine and every measured number are on the [Validation](validation) page.

## Wiki pages

| Page | Content |
|---|---|
| [Installation](installation) | Building and installing the plugin; checking that NCrystal sees it. |
| [Tutorial: 1D models](tutorial-1d) | End-to-end walkthrough: SasView 1D export, `sasview2ncmat.py` (1D is the default; `--2d` selects 2D), `DirectLoad` transport via the NCrystal Python API. |
| [Tutorial: 2D models](tutorial-2d) | End-to-end walkthrough: 2D export, `sasview2ncmat.py --2d`, anisotropic `DirectLoad2D` transport. |
| [Data format](data-format) | The `@CUSTOM_SASCSNS` grammar for `DirectLoad` and `DirectLoad2D`, the stored units, and what the converter does and checks. |
| [Validation](validation) | The eight-rung validation ladder: doctrine, the special-case harness, and commands. |
