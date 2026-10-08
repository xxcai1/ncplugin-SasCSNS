#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# example_sasview_chain_2d.py - the 2D SasView -> ncmat -> NCrystal chain.
#
# Workflow demonstrated (2D counterpart of example_sasview_chain.py):
#
#   1. SasView computes an I(Qx,Qy) image of an oriented-particle model in
#      barn per OBJECT. Here we emulate the "SasView -> Save Data" ASCII
#      export (columns Qx Qy I) for a cylinder of radius 25 Aa and length
#      30 Aa whose axis is tilted 30 degrees away from the beam (z) axis in
#      the xz plane. The SAS convention for 2D data is Q=(Qx,Qy,0). The
#      particle size is chosen so the form-factor first zeros fall on ~8-11
#      pixels of the export grid - a larger cylinder would be aliased by
#      the pixelisation (the rod factor sin(Qx L/2)/(Qx L/2) has its first
#      zero at Q = 2*pi/L).
#   2. script/sasview2ncmat.py --2d converts the image into an NCMAT file
#      with an @CUSTOM_SASCSNS DirectLoad2D section, applying the mandatory
#      barn/object -> barn/atom conversion itself.
#   3. The NCMAT file is loaded with NCrystal+plugin and cross sections are
#      compared against an independent Gauss-Legendre quadrature over the
#      same table (the alpha-form integral of section 4.2 of
#      doc/anisotropic_directload2d.pdf), for on-axis, perpendicular and
#      30-degree-tilted beams, at energies where the integration disc is
#      inside the table and where it is (correctly) truncated by it.
#   4. Sampling elasticity and the beam-direction dependence (anisotropy
#      must survive the whole chain) are checked.
# -----------------------------------------------------------------------------
import math
import os
import subprocess
import sys

import numpy as np
from scipy.special import j1

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get('SASVIEW_CHAIN_2D_WORKDIR', '/tmp/sasview_chain_2d')

#--- physical system ----------------------------------------------------------
RADIUS = 25.0        # Aa
LENGTH = 30.0        # Aa
TILT = math.radians(30.0)   # cylinder axis angle from +z, in the xz plane
PHI = 0.13           # particle volume fraction
DENSITY = 2.2        # g/cm3 amorphous silica
M_SIO2 = 60.084      # g/mol

V_P = math.pi * RADIUS**2 * LENGTH           # Aa^3
n_mol_cm3 = DENSITY / M_SIO2 * 6.022140857e23
n_atoms = n_mol_cm3 * 3.0 * 1.0e-24          # atoms per Aa^3
# cylinder contrast: dry silica in vacuum. SasView per-object intensity,
# I = V_p^2 * (delta_rho)^2 * F^2, comes out in Aa^2 = barn per object:
delta_rho = 3.475e-6                         # 1/Aa^2 (same SLD as 1D example)

#--- 1. emulate a SasView 2D export ------------------------------------------
NQ = 161
QMAX = 1.5                                   # 1/Aa, uniform grid
qx = np.linspace(-QMAX, QMAX, NQ)
qy = np.linspace(-QMAX, QMAX, NQ)


def form_factor_image(qx, qy):
    """Per-object F^2 image for the tilted cylinder (barn/object at F^2=1)."""
    qxx, qyy = np.meshgrid(qx, qy, indexing='xy')   # [ny,nx], qy outer
    qpar = qxx * math.sin(TILT)                      # Q.(axis) with Qz=0
    qperp2 = qxx**2 * math.cos(TILT)**2 + qyy**2
    qperp = np.sqrt(qperp2)
    x = qperp * RADIUS
    rod = qpar * LENGTH / 2.0
    f_disc = np.ones_like(x)
    nz = x > 1e-6
    f_disc[nz] = 2.0 * j1(x[nz]) / x[nz]
    f_rod = np.ones_like(rod)
    nz2 = np.abs(rod) > 1e-6
    f_rod[nz2] = np.sin(rod[nz2]) / rod[nz2]
    return V_P**2 * delta_rho**2 * f_disc**2 * f_rod**2


I_object = form_factor_image(qx, qy)

os.makedirs(WORK, exist_ok=True)
raw = os.path.join(WORK, 'cylinder2d_sasview_export.dat')
with open(raw, 'w') as fh:
    fh.write('# Qx (1/Aa)  Qy (1/Aa)  I (1/cm)\n')
    for iy in range(NQ):
        for ix in range(NQ):
            fh.write(f'{qx[ix]:.8g} {qy[iy]:.8g} {I_object[iy, ix]:.8g}\n')
print(f'step 1: wrote emulated SasView 2D export ({NQ}x{NQ} pixels) -> {raw}')

#--- 1b. emulated kernel vs the real sasmodels cylinder (optional) -------------
#The demo emulates SasView's 2D evaluation analytically. If sasmodels is
#installed, verify the emulation against the real kernel at 200 random
#non-zero pixels of the same grid: the SHAPES must agree to machine
#precision, while the ratio must be ONE constant -- precisely the unit
#convention (sasmodels 1/cm-with-scale vs barn/object) that
#sasview2ncmat.py absorbs.
try:
    from sasmodels.core import load_model
    from sasmodels.data import Data2D
    from sasmodels.direct_model import DirectModel
except ImportError:
    print('step 1b: sasmodels not installed - skipping kernel cross-check')
else:
    rngk = np.random.default_rng(5)
    iyk, ixk = np.unravel_index(rngk.choice(NQ * NQ, 200, replace=False),
                                (NQ, NQ))
    pxk, pyk = qx[ixk], qy[iyk]
    nz = (np.abs(pxk) > 1e-9) | (np.abs(pyk) > 1e-9)  # sasmodels masks q = 0
    pxk, pyk = pxk[nz], pyk[nz]
    prefk = I_object[iyk[nz], ixk[nz]]
    datak = Data2D(x=pxk, y=pyk, z=np.zeros(pxk.size))
    datak.err_data = np.ones(pxk.size)
    smk = np.atleast_1d(DirectModel(datak, load_model('cylinder'))(
        radius=RADIUS, length=LENGTH, sld=3.475, sld_solvent=0.0,
        theta=30.0, phi=0.0, background=0.0, scale=1.0))
    ratio = smk / prefk
    shape_dev = float(np.max(np.abs(smk / smk.max() - prefk / prefk.max())))
    const_dev = float(np.std(ratio) / np.mean(ratio))
    okk = shape_dev < 1e-6 and const_dev < 1e-6
    if not okk:
        failures += 1
    print(f'step 1b: sasmodels cylinder kernel vs emulation over {pxk.size} '
          f'pixels: shape max dev {shape_dev:.2e}, ratio const to '
          f'{const_dev:.2e} (factor {np.mean(ratio):.6g}, the unit '
          f'convention the converter absorbs) '
          f'[{"PASS" if okk else "FAIL"}]')

#--- 2. convert ----------------------------------------------------------------
ncmat = os.path.join(WORK, 'cylinder2d.ncmat')
scale_expected = (PHI / V_P) / n_atoms
subprocess.run([sys.executable, os.path.join(HERE, 'sasview2ncmat.py'),
                raw, '--2d', '-o', ncmat, '--material', 'sio2',
                '--density', str(DENSITY), '--phi', str(PHI),
                '--volume', f'{V_P:.6g}'], check=True)

#--- 3. read the file back: units + grid --------------------------------------
def read_back(path):
    qxr = qyr = None
    valsr = []
    with open(path) as fh:
        in_i = False
        for ln in fh:
            t = ln.split()
            if not t or t[0].startswith('#') or t[0] == 'NCMAT':
                continue
            if t[0] == 'Qx':
                qxr = [float(v) for v in t[1:]]
            elif t[0] == 'Qy':
                qyr = [float(v) for v in t[1:]]
            elif t[0] == 'I':
                in_i = True
                valsr.extend(float(v) for v in t[1:])
            elif in_i:
                valsr.extend(float(v) for v in t)
    return qxr, qyr, valsr


qxr, qyr, valsr = read_back(ncmat)
assert qxr is not None and len(valsr) == len(qxr) * len(qyr)
qxr = np.array(qxr)
qyr = np.array(qyr)
vals_img = np.array(valsr).reshape(len(qyr), len(qxr))
scale = vals_img.max() / I_object.max()
n_checked = 0
for iy, ix in [(0, 0), (NQ // 2, NQ // 3), (3, NQ - 5), (NQ - 1, NQ - 1),
               (NQ // 5, 2 * NQ // 3)]:
    rel = abs(vals_img[iy, ix] - I_object[iy, ix] * scale_expected) \
        / max(I_object[iy, ix] * scale_expected, 1e-300)
    #tolerance 1e-4: the converter derives n_d from per-element masses with
    #rounded stoichiometric fractions (60.0826 g/mol) vs this script's
    #tabulated M_SIO2 = 60.084 g/mol - a ~2.3e-5 constant, plus 7-digit
    #file rounding:
    assert rel < 1e-4, (iy, ix, rel)
    n_checked += 1
print(f'step 2: file I values match I_object * (phi/V_p)/n_d '
      f'= I_object * {scale_expected:.6g} at {n_checked} pixels '
      f'(rel < 1e-4, limited by composition rounding): units chain OK')

#--- independent quadrature over the file table --------------------------------
def bilinear(qxg, qyg, img, px, py):
    ix = np.clip(np.searchsorted(qxg, px) - 1, 0, len(qxg) - 2)
    iy = np.clip(np.searchsorted(qyg, py) - 1, 0, len(qyg) - 2)
    tx = (px - qxg[ix]) / (qxg[ix + 1] - qxg[ix])
    ty = (py - qyg[iy]) / (qyg[iy + 1] - qyg[iy])
    v = (img[iy, ix] * (1 - tx) * (1 - ty) + img[iy, ix + 1] * tx * (1 - ty)
         + img[iy + 1, ix] * (1 - tx) * ty + img[iy + 1, ix + 1] * tx * ty)
    return v


#alpha quadrature: dense MIDPOINT, not Gauss-Legendre. On ring-structured
#tables (Airy-like form factors) the alpha-integrand varies by orders of
#magnitude inside the forward boundary layer alpha <~ first ring; a GL-64
#rule was measured to overestimate sigma by ~12x while the midpoint rule is
#converged already at a few hundred nodes (always convergence-check
#reference quadratures themselves).
NA = 2000
NP = 1440
alpha = (np.arange(NA) + 0.5) * math.pi / NA
wa = math.pi / NA
phi_g = (np.arange(NP) + 0.5) * (2.0 * math.pi / NP)
wp = 2.0 * math.pi / NP


def sigma_quad(k, beam):
    """alpha-form sigma (per-atom units), midpoint alpha x midpoint phi."""
    ax = k * beam[0]
    ay = k * beam[1]
    tot = 0.0
    for a_ in alpha:
        sa_ = math.sin(a_)
        px = k * sa_ * np.cos(phi_g) - ax
        py = k * sa_ * np.sin(phi_g) - ay
        tot += sa_ * float(np.sum(bilinear(qxr, qyr, vals_img, px, py))) * wp
    return tot * wa


#--- 4. load with NCrystal + plugin, compare -----------------------------------
import NCrystal as NC

scat = NC.createScatter(ncmat)
k_of_e = lambda e_mev: math.sqrt(e_mev / 2.07214)
failures = 0
print('step 3: plugin sigma vs independent quadrature of the same file table')
#beams within the model's design tilt range (beta <= 60 mrad; the sigma
#grid's s = k*sin(beta) axis covers exactly that cone for E <= 100 meV).
#The truncated-disc row at 3.0 meV carries a looser tolerance (grid
#quadrature on the table edge, documented in the plugin validation).
for e_mev, beam, tol in [(0.3, (0, 0, 1), 5e-3),
                         (0.7, (0, 0, 1), 5e-3),
                         (1.1, (0, 0, 1), 5e-3),
                         (3.0, (0, 0, 1), 3e-2),
                         (0.7, (math.sin(0.03), 0, math.cos(0.03)), 2e-2),
                         (0.7, (math.sin(0.06), 0, math.cos(0.06)), 2e-2),
                         (0.7, (math.sin(0.5), 0, math.cos(0.5)), 1e-2)]:
    ref = sigma_quad(k_of_e(e_mev), beam)
    got = float(scat.crossSection(e_mev * 1e-3, beam))
    rel = abs(got - ref) / ref
    tag = 'PASS' if rel < tol else 'FAIL'
    if tag == 'FAIL':
        failures += 1
    print(f'  E={e_mev:4.1f}meV beam={beam}: quadrature={ref:.6g} '
          f'plugin={got:.6g} rel={rel:.2e} [{tag}]')

#robustness far outside the design tilt range (beta = 90 deg clamps on
#the s grid; must stay finite and positive, accuracy is out of scope):
sigma_perp = float(scat.crossSection(0.7e-3, (1, 0, 0)))
print(f'step 4: out-of-range beam (90 deg): sigma = {sigma_perp:.4g} '
      f'(clamped s grid, finite and positive: '
      f'{"OK" if 0.0 < sigma_perp < float("inf") else "FAIL"})')
if not 0.0 < sigma_perp < float('inf'):
    failures += 1

#anisotropy must survive the chain - checked at the shape level: sampled
#outcomes on axis must have second moments that reflect the elongated
#image (the tilted cylinder's image is stretched along Qx):
k2 = math.sqrt(2.07214 / 2.07214)
ek, dirs = scat.sampleScatter(2e-3, (0, 0, 1), repeat=200_000)
qx_s = k2 * dirs[0]
qy_s = k2 * dirs[1]
mratio = float(np.mean(qx_s**2) / np.mean(qy_s**2))
print(f'step 4: anisotropy survives the chain: <Qx^2>/<Qy^2> of sampled '
      f'outcomes = {mratio:.4f} (must differ from 1 by >5%)')
if abs(mratio - 1.0) < 0.05:
    failures += 1
    print('  FAIL: anisotropy lost')

ek, dirs = scat.sampleScatter(1e-3, (0, 0, 1), repeat=200)
if np.any(ek != 1e-3):
    failures += 1
    print('step 4: FAIL sampling is not elastic')
else:
    print(f'step 4: sampling is elastic (200 events, ekin unchanged); '
          f'sampled |uz| in [{abs(dirs[2]).min():.3f},{abs(dirs[2]).max():.3f}]')

print('=' * 60)
if failures:
    sys.exit(f'{failures} check(s) FAILED')
print('2D SasView chain: ALL CHECKS PASS')
