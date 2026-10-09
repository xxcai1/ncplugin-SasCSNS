# Tutorial: 1D DirectLoad — from SasView to a NCrystal material

This page walks through the complete 1D chain: a SasView `I(Q)` export is
converted by `script/sasview2ncmat.py` into an NCMAT file with an
`@CUSTOM_SASCSNS` `DirectLoad` section, which is then used for neutron
transport through NCrystal's Python API. Following the *DirectLoad*
philosophy the plugin develops **no new physics models**: SasView (or any
other calculator) evaluates the form factor and exports a tabulated
intensity; the plugin transports neutrons through it numerically. The
physical system used throughout is the chain example's: SiO2 spheres of
R = 100 Angstrom radius at 40 vol% in vacuum (dry nanopowder), material
density 2.2 g/cm3, probed at 0.5-20 meV. The complete runnable version of
everything below is
[`script/example_sasview_chain.py`](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/script/example_sasview_chain.py);
it emulates the SasView export itself, so it needs neither a SasView
installation nor network access. For environment setup see
[Installation](installation); the NCMAT section syntax is specified in
[Data format](data-format); the anisotropic counterpart is
[Tutorial: 2D DirectLoad2D](tutorial-2d).

## Step 1: export a 1D I(Q) table from SasView

In SasView, evaluate a sphere model and use *Save Data* (ASCII export).
The converter `script/sasview2ncmat.py` expects:

- one row per Q value with columns `Q I [dI] [dQ]`. Only the first two
  tokens (Q and I) are read; extra columns (`dI`, `dQ`) are ignored.
  Header/metadata lines start with `#`; commas are tolerated as
  separators; lines that do not parse as numbers (column titles etc.) are
  silently skipped. At least 2 `(Q, I)` points are required.
- Q and I must be positive: points with `Q <= 0`, `I <= 0` or NaN
  values are dropped with a warning; the table is sorted ascending in Q,
  and duplicated Q values (within 1e-12) are dropped with a warning.
- the I column is interpreted as **barn per object** (per particle) —
  SasView's per-particle intensity for `scale = 1`, `background = 0`.

You do not need SasView installed to exercise the chain: the example
script `script/example_sasview_chain.py` (needs numpy) emulates exactly
such an export and writes the files `sasview_export_<tag>.dat` into
`$SASVIEW_CHAIN_WORKDIR`
(default `/tmp/sasview_chain`) — run it first and use any of them in place of
`sasview_iq.dat` in step 2. Its write loop (excerpt):

```python
i_obj = V_p**2 * delta_rho**2 * BARN_PER_AA2 * sphere_form_factor(qgrid)
...
    fh.write('# Q (1/A)    I (barn/object)   dI   dQ\n')
    for qi, ii in zip(qgrid, i_obj):
        fh.write(f'{qi:.8g}\t{ii:.8g}\t0\t0\n')
```

> **Callout — the intensity is in Angstrom^2, not barn.** If you generate
> the export yourself: the per-object forward amplitude `V_p**2 *
> delta_rho**2`, with the SLD contrast `delta_rho` in 1/Angstrom^2, has
> dimension **Angstrom^2**. Since 1 Angstrom^2 = 1e8 barn, the exported
> intensity must include a factor `1e8` (the example defines
> `BARN_PER_AA2 = 1.0e8`). Omitting it is the classic
> 8-orders-of-magnitude error.

## Step 2: convert with sasview2ncmat.py

From the repo root:

```bash
python3 script/sasview2ncmat.py sasview_iq.dat -o sas_model.ncmat \
        --material sio2 --density 2.2 --phi 0.4 --radius 100.0
```

1D is the default mode — there is no `--1d` flag (the `--2d` flag selects
the anisotropic 2D path, see [Tutorial: 2D DirectLoad2D](tutorial-2d)).
The complete flag reference is the table on
[Data format, Converter CLI](data-format#converter-cli-sasview2ncmatpy); in
short: `--radius`/`--volume` are mutually exclusive, and one of them (plus
`--phi`) is mandatory because the unit conversion needs the particle volume.
Note that `--solvent` (optional) is recorded as a comment line
`# solvent: NAME` in the generated file, in both 1D and 2D modes; the
pitfall concerns hand-written live `solvent` section lines instead (see
[Data format, Pitfalls](data-format#pitfalls)).

### What the automatic unit conversion does

The converter multiplies every I value by
`scale = n_objects / n_atoms = (phi / V_p) / n_d` — **always applied, there
is deliberately no flag for it** (the normative copy, including the
`1 Angstrom^2 = 1e8 barn` trap and the `Sigma = n_d * sigma` invariant, is
on [Data format, Units](data-format#units-barnatom-not-barnobject); here
`n_d = 0.066152 atoms/Aa^3` for the `sio2` preset at 2.2 g/cm3). For the
system above (R = 100 Aa): `V_p = 4.18879e6 Angstrom^3`,
`n_obj/n_d = 1.4435377e-06`, so the forward intensity `I(0) = 2.1176e10
barn/object` becomes about `3.06e4 barn/atom` in the NCMAT file.

A real converter run on the emulated export prints (paths shortened; the
two `File validation:` lines summarize the data requirements of step 3 for
your actual file):

```text
File validation: Qmax = 6.6 1/Aa -> fully covers neutron energies up to E = 22.57 meV (Q = 2k rule)
File validation: min spacing 4.55e-06 1/Aa resolves form-factor oscillations of R = 100 Aa (need <= 0.00314) OK
Wrote sas_model.ncmat with 2120 (Q,I) points from sasview_iq.dat
Unit conversion (always applied): scale = (phi/V_p)/n_d = (0.4/4.18879e+06)/0.066152 = 1.4435377e-06 [barn/object -> barn/atom]
Load it with NCrystal, e.g.:  nctool sas_model.ncmat -x 1e-5:1.0 -p
```

### What the converter writes

The generated NCMAT file (same run, abridged — the header carries one
`Q=... I=...` comment line per table point):

```text
NCMAT v6
# Generated by sasview2ncmat.py from a SasView-exported I(Q) curve.
# Scattering model: SasCSNS plugin, DirectLoad mode (vacuum solvent)
#
# Unit conversion applied (SasView barn/object -> NCrystal barn/atom):
#   I_ncmat = I_sasview * 1.4435377e-06
#   scale = (phi / V_p) / n_d with V_p = 4/3*pi*R^3, R = 100 Aa (sphere)
#
#   Q=0.0001       I=30568.1
#   ...            (one comment line per table point)
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
  Q 0.0001 0.00010455089 0.00010930888 ... 6.6
  I 30568.104 30568.048 30567.986 ... 1.3489498e-06
```

Keep the Q and I lists the **same length**: the parser loops over the `Q`
tokens while indexing the `I` line, so an `I` list **shorter** than `Q`
surfaces as an out-of-range exception rather than a clear error message,
while an `I` list **longer** than `Q` is silently truncated (the extra
values are never read).

## Step 3: data requirements for the I(Q) table

The general rules — the `[0, 2k(E_max)]` coverage rule with the `--emax`
gate, the `pi/(10 R)` resolution rule, and starting the table at (or near)
`Q = 0` — are normative on [Data format, Q-range and coverage
requirements](data-format#q-range-and-coverage-requirements). For the chain
example's system they mean concretely:

- Kinematics `E[meV] = 2.07214 * k^2` (k in 1/Angstrom) makes the table
  reach `2*sqrt(E_max/2.07214)` 1/Angstrom — e.g. 6.2 1/Aa for the
  E_max = 20 meV top of the 0.5-20 meV range probed below.
- Violating the coverage rule is costly: over that range, the chain
  example's truncation control (table cut at Qmax = 0.45 1/Aa where
  6.2 1/Aa is needed) shows cross-section biases of factors 2 up to 14.
- R = 100 Aa means the Q spacing must be <= pi/(10*R) = 0.00314 1/Aa
  (oscillation period pi/R = 0.0314 1/Aa); coarser tables trigger a
  converter WARNING (`too coarse to resolve form-factor oscillations of
  R = ...`).

## Step 4: load and transport in NCrystal (Python API)

```python
import NCrystal
workdir = '/path/to/dir/holding/sas_model.ncmat'     # the repo-root cwd of step 2 works
NCrystal.addCustomSearchDirectory(workdir)           # needed only for files outside the cwd
scat = NCrystal.createScatter('sas_model.ncmat')     # DirectLoad SANS model

e = 5.0e-3                                           # E = 5 meV, in eV!
sigma = float(scat.crossSection(e, (0, 0, 1)))       # [barn/atom], beam +z
ek, dirs = scat.sampleScatter(e, (0, 0, 1), repeat=200_000)  # outcomes
```

- Bare relative paths resolve against the current working directory by
  default, so from step 2's repo-root cwd the plain
  `NCrystal.createScatter('sas_model.ncmat')` needs no search directory at
  all (this is what [Tutorial: 2D](tutorial-2d) step 2 relies on);
  `addCustomSearchDirectory` matters when the file lives elsewhere — e.g.
  in the chain example's `$SASVIEW_CHAIN_WORKDIR` (default
  `/tmp/sasview_chain`), which its own script registers as `WORK`.
- The NCrystal Python API (3.x and 4.x) takes the energy argument
  (`ekin`) in **eV**, not meV — every script in the repo passes
  `E_meV * 1e-3` (5 meV becomes `5.0e-3`).
- Scattering is elastic: the sampled outcome energies `ek` equal the
  incident energy. The model is isotropic and consumes the table directly
  (no fitting). A quick command-line check: `nctool sas_model.ncmat
  -x 1e-5:1.0 -p` (`-x/--xrange` takes a single colon-syntax argument).
  Note: the hint line printed by sasview2ncmat.py uses this same colon
  syntax; the two-argument form (`-x 1e-5 1.0`) is rejected by nctool
  with "unrecognized arguments: 1.0" — use the colon-syntax command
  above.
- The plugin is gated by the `sans` config parameter: active by default,
  `sans=0` disables the SANS models, leaving NCrystal's standard
  scattering. It serves single-phase materials only — files with
  `@OTHERPHASES` are declined and handled by NCrystal's builtin SANS
  factories.
- The resulting material is **SANS-only**: the plugin model replaces the
  standard scattering of the carrier in this file (step 5 shows why that
  matters in comparisons).

## Step 5: the trust anchor — 0.14% against the builtin model

Why believe the numbers? The full chain was cross-checked against
NCrystal's official builtin `@CUSTOM_HARDSPHERESANS` model for exactly the
system above (same R, phi, density): the SANS macroscopic cross sections
agree to **0.14%** at worst over 0.5-20 meV. Two normalisation conventions
must be respected in such a comparison:

1. The builtin reference divides by the number density of the void-diluted
   mixture (declared via `@OTHERPHASES` with a vacuum phase), so its
   per-atom values are a factor 1/(1-phi) larger than those of the
   pure-SiO2 DirectLoad file, where phi is folded into the conversion
   scale instead. The invariant comparison quantity is therefore the
   **macroscopic cross section Sigma = n_d * sigma** — never raw sigma
   per atom.
2. The builtin factory adds the standard free-gas scattering on top of the
   SANS process, while the DirectLoad file is SANS-only. The standard part
   is measured via `sans=0` and subtracted.

The chain example implements exactly this (excerpt):

```python
NCrystal.addCustomSearchDirectory(WORK)
dl = NCrystal.createScatter(ncmat_full)               # SANS only (plugin model)
of = NCrystal.createScatter(builtin)                  # std + SANS (builtin model)
of_std = NCrystal.createScatter(builtin + ';sans=0')  # std part only
n_file = n_atoms                    # number density of the DirectLoad file
n_mix = (1.0 - PHI) * n_atoms       # number density of the void-diluted
                                    # mixture the builtin model divides by
...
    s_chain, s_of = sig_dl * n_file, sig_of * n_mix
    s_of_sans = s_of - sig_std * n_mix
    ratio = s_chain / s_of_sans
```

Unit bookkeeping for Sigma, once and for all: 1 barn = 1e-24 cm^2 and
1 Angstrom^3 = 1e-24 cm^3 cancel exactly, so
`n_d [atoms/Aa^3] * sigma [barn/atom]` **is** Sigma in 1/cm — no extra
factor. A 0.14% agreement at the Sigma level validates the entire chain:
export units, the barn/object -> barn/atom conversion, and the plugin
model itself. The same script also runs the truncated-Q control test of
step 3. See [Validation](validation) for the complete validation ladder.

## Run the complete example

```bash
python script/example_sasview_chain.py
```

executes all five steps in one go (validation rung 6, see
[Validation](validation)): emulated SasView export, unit
conversion, builtin-model reference, the Sigma comparison table at six
energies (0.5-20 meV) with the worst-agreement percentage, and the
truncated-Q control test. Output files land in `$SASVIEW_CHAIN_WORKDIR`
(default `/tmp/sasview_chain`). See also
[Data format](data-format) and [Tutorial: 2D DirectLoad2D](tutorial-2d).
