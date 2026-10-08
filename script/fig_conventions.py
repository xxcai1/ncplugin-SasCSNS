#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig_conventions.py - figure for doc/anisotropic_directload2d.pdf (rung 7 /
# the backward-branch section): the one picture of the q_z=0 convention
# finding. Same narrow-Gaussian particles, same beam, both plugin models:
#   left  : sampled p(theta) -- the 2D model's backward lobe (theta ~ pi)
#           carries the mirrored forward intensity where the physical
#           I(|Q|~2k) ~ 0 (1D model); this is the sigma doubling.
#   right : Beer-Lambert transmission through a slab -- the operational
#           consequence: same particles, up to 2x different attenuation.
# -----------------------------------------------------------------------------
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get('PLUGIN2D_WORKDIR', '/tmp/plugin_2d')
OUT = os.path.join(HERE, '..', 'doc', 'fig_conventions.pdf')

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, HERE)
from example_transmission_2d import transport, write_ncmat_1d, write_ncmat_2d  # noqa: E402

import NCrystal  # noqa: E402

E_MEV = 2.07214                        # k = 1 /Aa
S_Q = 0.1
AMP = 40.0

q1d = np.linspace(0.0, 3.2, 6401)
p1 = write_ncmat_1d('figconv1d', q1d, AMP * np.exp(-q1d ** 2 / (2 * S_Q ** 2)))
g2 = np.linspace(-3.2, 3.2, 1601)
X, Y = np.meshgrid(g2, g2)
p2 = write_ncmat_2d('figconv2d', g2, g2,
                    AMP * np.exp(-(X ** 2 + Y ** 2) / (2 * S_Q ** 2)))

sc1 = NCrystal.createScatter(p1)
sc2 = NCrystal.createScatter(p2)
ekin = E_MEV * 1e-3
n = 2_000_000
_, d1 = sc1.sampleScatter(ekin, (0, 0, 1), repeat=n)
_, d2 = sc2.sampleScatter(ekin, (0, 0, 1), repeat=n)
th1 = np.arccos(np.clip(np.asarray(d1[2]), -1, 1))
th2 = np.arccos(np.clip(np.asarray(d2[2]), -1, 1))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.4))

#--- left: p(theta) over the full sphere, log y --------------------------------
nb = 120
edges = np.linspace(0, math.pi, nb + 1)
h1, _ = np.histogram(th1, bins=edges, density=True)
h2, _ = np.histogram(th2, bins=edges, density=True)
mid = 0.5 * (edges[1:] + edges[:-1])
dth = edges[1] - edges[0]
floor = 1.0 / (n * dth)
ax1.plot(mid, np.degrees(mid) * 0 + np.sin(mid) / 2, color='0.75', lw=1.2,
         ls=':', label=r'uniform sphere $\sin\theta/2$')
ax1.plot(mid, np.where(h1 > 0, h1, np.nan), lw=1.6,
         label=rf'1D DirectLoad: $\sigma_\mathrm{{1D}}='
               rf'{float(sc1.crossSection(ekin, (0, 0, 1))):.2f}$ barn')
ax1.plot(mid, np.where(h2 > 0, h2, np.nan), lw=1.6,
         label=rf'2D DirectLoad2D: $\sigma_\mathrm{{2D}}='
               rf'{float(sc2.crossSection(ekin, (0, 0, 1))):.2f}$ barn')
ax1.set_yscale('log')
ax1.set_ylim(2e-5, 40)
ax1.set_xlabel(r'outcome polar angle $\theta$ (rad)')
ax1.set_ylabel(r'sampled density  $p(\theta)$')
ax1.set_title(r'narrow Gaussian $s_Q=0.1\,\mathrm{\AA}^{-1}$, $E=2.07$ meV, '
              r'$k=1\,\mathrm{\AA}^{-1}$', fontsize=9)
ax1.legend(fontsize=8, frameon=False, loc='lower left',
           bbox_to_anchor=(0.02, 0.02))
ax1.annotate('backward lobe:\nmirrored forward\nintensity (factor 2)',
             xy=(3.02, 5.0), xytext=(1.62, 6.0), fontsize=8,
             arrowprops=dict(arrowstyle='->', lw=0.9))
ax1.annotate('physical:\n$I(|Q|\\!\\approx\\!2k)\\approx 0$',
             xy=(3.08, 2.5e-4), xytext=(2.05, 3e-3), fontsize=8,
             arrowprops=dict(arrowstyle='->', lw=0.9))

#--- right: transmission vs thickness ------------------------------------------
r2 = transport(p2, ekin, (0, 0, 1), 1.0, 1000)
r1 = transport(p1, ekin, (0, 0, 1), 1.0, 1000)
L = np.linspace(0, 3.2 * r2['mfp'], 300)
T2 = np.exp(-r2['sigma_inv'] * L)
T1 = np.exp(-r1['sigma_inv'] * L)
ax2.plot(L, T1, lw=1.6, label='1D DirectLoad material')
ax2.plot(L, T2, lw=1.6, label='2D DirectLoad2D material')
# MC points at 0.25/1/3 mfp (1e6 neutrons, binomial 4-sigma whiskers tiny)
rng = np.random.default_rng(7)
for tag, r, c in (('1D', r1, 'C0'), ('2D', r2, 'C1')):
    for fac in (0.25, 1.0, 3.0):
        Lc = fac * r2['mfp']
        rr = transport({'1D': p1, '2D': p2}[tag], ekin, (0, 0, 1), Lc, 1_000_000)
        p_mc = rr['uncollided'] / rr['n']
        p_err = math.sqrt(p_mc * (1 - p_mc) / rr['n'])
        ax2.errorbar(Lc, p_mc, yerr=p_err, fmt='o', color=c, ms=3.5, capsize=2)
ax2.set_xlabel(r'slab thickness $L$ (cm)')
ax2.set_ylabel('uncollided fraction')
ax2.set_title(r'$\Sigma=n(\sigma_\mathrm{sc}+\sigma_\mathrm{ab})$: '
              rf'1D ${r1["sigma_inv"]:.3f}\,$cm$^{{-1}}$, '
              rf'2D ${r2["sigma_inv"]:.3f}\,$cm$^{{-1}}$ (same particles!)',
              fontsize=9)
ax2.legend(fontsize=8, frameon=False, loc='upper right')
ax2.set_xlim(0, L[-1])
ax2.set_ylim(0, 1.02)

fig.tight_layout()
fig.savefig(OUT)
print(f'wrote {OUT}')
