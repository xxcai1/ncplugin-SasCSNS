#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig_chain_2d.py - figure for doc/anisotropic_directload2d.pdf (rung 6):
# the anisotropic I(Qx,Qy) input image of the 2D chain demo (tilted cylinder)
# next to the (Qx,Qy) distribution of outcomes sampled by the C++ plugin from
# it. Linear colour scale (SasView-style), same contours in both panels: the
# sampler reproduces
# the anisotropic pattern -- shape-level proof that it survives the whole
# SasView -> NCMAT -> NCrystal chain.
# -----------------------------------------------------------------------------
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get('SASVIEW_CHAIN_2D_WORKDIR', '/tmp/sasview_chain_2d')
NCMAT = os.path.join(WORK, 'cylinder2d.ncmat')
OUT = os.path.join(HERE, '..', 'doc', 'fig_chain_2d.pdf')

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

if not os.path.exists(NCMAT):
    subprocess.run([sys.executable, os.path.join(HERE, 'example_sasview_chain_2d.py')],
                   check=True)

import NCrystal  # noqa: E402

#--- regenerate the input image (same constants as the chain demo) ------------
RADIUS, LENGTH = 25.0, 30.0
TILT = math.radians(30.0)
V_P = math.pi * RADIUS**2 * LENGTH
DELTA_RHO = 3.475e-6
from scipy.special import j1  # noqa: E402

qx1 = qy1 = None
with open(NCMAT) as fh:
    for ln in fh:
        t = ln.split()
        if t and t[0] == 'Qx':
            qx1 = np.array([float(v) for v in t[1:]])
        elif t and t[0] == 'Qy':
            qy1 = np.array([float(v) for v in t[1:]])
        elif t and t[0] == 'I':
            break
X, Y = np.meshgrid(qx1, qy1)          # [ny,nx], qy outer
qpar = X * math.sin(TILT)
qperp = np.sqrt(X**2 * math.cos(TILT)**2 + Y**2)
x, rod = qperp * RADIUS, qpar * LENGTH / 2.0
f_disc = np.where(x > 1e-6, np.sinc(x / math.pi), 1.0)      # 2*J1(x)/x
f_rod = np.where(np.abs(rod) > 1e-6, np.sinc(rod / math.pi), 1.0)
I_img = V_P**2 * DELTA_RHO**2 * f_disc**2 * f_rod**2

#--- sample outcomes and map them back to (Qx,Qy) -----------------------------
sc = NCrystal.createScatter(NCMAT + ';dcutoff=0')
E_MEV = 2.07214                        # k = 1 /Aa
K = math.sqrt(E_MEV / 2.07214)
_, dirs = sc.sampleScatter(E_MEV * 1e-3, (0, 0, 1), repeat=2_000_000)
ux, uy, _ = dirs
QX, QY = K * ux, K * uy

#--- figure: both panels as probability densities on one colour scale --------
#input normalised to a probability density over the (Qx,Qy) plane:
dA = (qx1[1] - qx1[0]) * (qy1[1] - qy1[0])
P_in = I_img / I_img.sum() / dA
#histogram with bins CENTRED on the pixel centres:
dx, dy = qx1[1] - qx1[0], qy1[1] - qy1[0]
xe = np.append(qx1 - dx / 2, qx1[-1] + dx / 2)
ye = np.append(qy1 - dy / 2, qy1[-1] + dy / 2)
Hh, _, _ = np.histogram2d(QX, QY, bins=[xe, ye])
P_out = Hh.T / Hh.sum() / dA          # qy outer to match panel (a)

WIN = 0.5                    # display window around the pattern /Aa
mask = P_in > 0
Pk = P_in[mask].max()
norm = None                        # linear colour scale, SasView-style
lv = [Pk * f for f in (0.02, 0.1, 0.3, 0.6, 0.9)]

fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.35), sharey=True)
fig.subplots_adjust(left=0.07, right=0.86, top=0.86, bottom=0.14, wspace=0.06)

ax = axes[0]
im = ax.pcolormesh(qx1, qy1, np.where(mask, P_in, np.nan),
                   norm=norm, cmap='viridis', shading='nearest')
ax.contour(qx1, qy1, np.where(mask, P_in, np.nan), levels=lv,
           colors='k', linewidths=0.6)
ax.set_xlim(-WIN, WIN)
ax.set_ylim(-WIN, WIN)
ax.set_title(r'input density $I(Q_x,Q_y)/\!\int I$ (cylinder at $30^\circ$)',
             fontsize=10)
ax.set_xlabel(r'$Q_x$ [$\mathrm{\AA}^{-1}$]')
ax.set_ylabel(r'$Q_y$ [$\mathrm{\AA}^{-1}$]')

ax = axes[1]
mout = P_out > 0
ax.pcolormesh(qx1, qy1, np.where(mout, P_out, np.nan),
              norm=norm, cmap='viridis', shading='nearest')
#overlay the INPUT contours (white) to make pattern agreement visible
ax.contour(qx1, qy1, np.where(mask, P_in, np.nan), levels=lv,
           colors='k', linewidths=0.6)
ax.set_xlim(-WIN, WIN)
ax.set_ylim(-WIN, WIN)
ax.set_title(r'sampled outcome density, 2M events', fontsize=10)
ax.set_xlabel(r'$Q_x$ [$\mathrm{\AA}^{-1}$]')

cbar = fig.colorbar(im, ax=axes, shrink=0.92, pad=0.02)
cbar.set_label(r'probability density per $\mathrm{\AA}^2$', fontsize=9)
fig.savefig(OUT, bbox_inches='tight')
try:
    fig.savefig('/tmp/fig_chain_2d_preview.png', dpi=110, bbox_inches='tight')
except Exception:
    pass
print(f'wrote {OUT}')
