# NCrystal plugin SasCSNS

A [NCrystal](https://github.com/mctools/ncrystal) plugin providing
general-purpose **small-angle neutron scattering (SANS) models** for
particle transport. See the wiki page for more information.
https://github.com/xxcai1/ncplugin-SasCSNS/wiki

## Model types (`@CUSTOM_SASCSNS` NCMAT section)

| type | input | notes |
|---|---|---|
| `DirectLoad` | tabulated 1D `I(Q)` | full elastic sphere, `Q = 2k sin(theta/2)`; physical total cross section (validated against builtin HardSphereSANS at 0.14%) |
| `DirectLoad2D` | tabulated 2D `I(Qx,Qy)` | anisotropic patterns in detector coordinates (`q_z = 0` convention); on-axis sampling exact and rejection-free; tilted beams via dilated-proposal rejection; see `doc/anisotropic_directload2d.pdf` for the full specification |
| `HardSphere` | analytic R, contrast | analytic reference model |

Following the *DirectLoad* philosophy the plugin develops **no new
physics models**: SasView (or any other calculator) evaluates the form
factor and exports a tabulated intensity; the plugin transports
neutrons through it numerically.

## Conversion tool

`script/sasview2ncmat.py` converts SasView 1D (`--1d`) and 2D (`--2d`)
ASCII exports to NCMAT files, applying the mandatory
`barn/object -> barn/atom` unit conversion always and automatically.

## Validation

The models are validated by a special-case harness ("validation
ladder", 8 rungs, all passing); the design doctrine and every measured
number are documented in
[`doc/anisotropic_directload2d.pdf`](doc/anisotropic_directload2d.pdf).
Run from the repo root with the NCrystal+plugin environment active
(`python` must see the plugin, i.e. `pip install .` first):

```
python script/validate_directload2d.py    # rungs 2-3: python model vs closed forms (29 checks)
python script/validate_plugin_2d.py       # rung 4: C++ plugin vs reference (32 checks)
python script/validate_1d_vs_2d.py        # rung 7: cross-model consistency + conventions (7 checks)
python script/example_sasview_chain.py    # rung 6 (1D): full SasView chain
python script/example_sasview_chain_2d.py # rung 6 (2D): full SasView chain + figure data
python script/example_transmission_2d.py  # rung 8: transport-integrated attenuation
ncrystal-pluginmanager --test SasCSNS     # rung 5: regression anchors
python script/fig_chain_2d.py             # render doc/fig_chain_2d.pdf
python script/fig_conventions.py          # render doc/fig_conventions.pdf
```

Rung 7 documents the one known convention effect: the 2D model's
`q_z = 0` mirror extension inflates the *total* cross section of
narrow patterns by up to a factor two (attenuation only;
detector-visible forward scattering is unaffected). Data producers
should read `doc/anisotropic_directload2d.pdf`, section "The backward
branch and the relation to the 1D model".
