#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig_sigma_closedforms.py - figure for doc/anisotropic_directload2d.pdf
# (rung 4c, sweep view): the closed-form total cross sections sigma_2D(k) of
# the Guinier and Zimm/OZ tables (lines) against the C++ plugin evaluated at
# 14 log-spaced on-axis energies (dots), linear scales; below each, the
# relative deviation with the rung's 1.5e-3 tolerance band. The rise of the
# deviation toward k = 1.5 is the bilinear table's truncation error as the
# sampling disc (radius k) approaches the table edge (qmax = 1.6/Aa).
# Same model constants as script/validate_physical_closedforms.py.
# -----------------------------------------------------------------------------
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get('PLUGIN2D_WORKDIR', '/tmp/plugin_2d')
OUT = os.path.join(HERE, '..', 'doc', 'fig_sigma_closedforms.pdf')

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from scipy.special import erfi  # noqa: E402

I0, RG, XI = 2.0, 1.0, 0.7
HALF, E_OVER_K2 = 1.6, 2.07214
TOL = 1.5e-3                       # the rung 4c pass band


def sigma2d_guinier(k):
    a = k ** 2 * RG ** 2 / 3.0
    if a < 1e-12:
        return 4.0 * math.pi * I0
    return 4.0 * math.pi * I0 * math.sqrt(math.pi) / 2.0 \
        * math.exp(-a) * float(erfi(math.sqrt(a))) / math.sqrt(a)


def sigma2d_oz(k):
    a = k ** 2 * XI ** 2
    if a < 1e-12:
        return 4.0 * math.pi * I0
    return 4.0 * math.pi * I0 * \
        math.atanh(math.sqrt(a / (1.0 + a))) / math.sqrt(a * (1.0 + a))


def write_2d(name, i_func, n=501):
    g = np.linspace(-HALF, HALF, n)
    X, Y = np.meshgrid(g, g)
    flat = i_func(np.hypot(X, Y)).ravel()      # row-major, qy outer
    path = os.path.join(WORK, f'{name}.ncmat')
    with open(path, 'w') as fh:
        fh.write('NCMAT v7\n@DENSITY\n  2.2 g_per_cm3\n@DYNINFO\n  element Si\n'
                 '  fraction 1\n  type freegas\n@CUSTOM_SASCSNS\n'
                 '  DirectLoad2D\n  Qx ')
        fh.write(' '.join(f'{v:.10g}' for v in g) + '\n  Qy ')
        fh.write(' '.join(f'{v:.10g}' for v in g) + '\n  I ')
        for j, v in enumerate(flat):
            fh.write(f'{v:.8g}')
            fh.write('\n  ' if (j + 1) % 10 == 0 and j + 1 < flat.size else ' ')
        fh.write('\n')
    return path


def i_guinier(q):
    return I0 * np.exp(-q ** 2 * RG ** 2 / 3.0)


def i_oz(q):
    return I0 / (1.0 + q ** 2 * XI ** 2)


os.makedirs(WORK, exist_ok=True)
KS = np.geomspace(0.25, 1.5, 14)               # E in [0.13, 4.7] meV, in range
KD = np.linspace(0.25, 1.5, 300)               # dense line for the closed form

import NCrystal  # noqa: E402

MODELS = (('Guinier  ($R_g=1.0\\,\\mathrm{{\\AA}}^{{-1}}$)',
           'guinier2d', i_guinier, sigma2d_guinier),
          ('Zimm/OZ  ($\\xi=0.7\\,\\mathrm{{\\AA}}^{{-1}}$)',
           'oz2d', i_oz, sigma2d_oz))

fig, axes = plt.subplots(2, 2, figsize=(9.6, 5.6), sharex='col',
                         gridspec_kw=dict(height_ratios=[2.2, 1.05]))
fig.subplots_adjust(left=0.075, right=0.985, top=0.92, bottom=0.115,
                    wspace=0.18, hspace=0.10)

worst_all = 0.0
for col, (title, name, i_f, s_cf) in enumerate(MODELS):
    path = write_2d(name, i_f)
    sc = NCrystal.createScatter(path)
    cf = np.array([s_cf(k) for k in KS])
    pl = np.array([float(sc.crossSection(k * k * E_OVER_K2 * 1e-3, (0, 0, 1)))
                   for k in KS])
    dev = (pl - cf) / cf
    worst_all = max(worst_all, float(np.abs(dev).max()))
    for k, s in zip(KS, pl):
        print(f'{name}  k={k:.4f}  closed={s_cf(k):11.6f}  '
              f'plugin={s:11.6f}  rel.dev={abs(s / s_cf(k) - 1.0):.2e}')
    kw = dev.max(), dev.min(), KS[int(np.argmax(np.abs(dev)))]

    ax = axes[0, col]
    ax.plot(KD, [s_cf(k) for k in KD], '-', lw=1.8, color=f'C{col}',
            label='closed form')
    ax.plot(KS, pl, 'o', ms=4.6, mfc='white', mec=f'C{col}', mew=1.4,
            label='C++ plugin', ls='none')
    ax.set_ylabel(r'$\sigma_{2D}$  [barn/atom]')
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc='best', framealpha=0.9)

    ax = axes[1, col]
    ax.axhspan(-TOL * 1e3, TOL * 1e3, color='0.88', zorder=0)
    ax.axhline(0.0, color='0.35', lw=0.8)
    ax.plot(KS, dev * 1e3, 'o-', ms=3.6, lw=1.1, color=f'C{col}')
    ax.annotate(f'worst {max(abs(kw[0]), abs(kw[1])) * 1e3:.2f}' + r'$\times10^{-3}$'
                f' at $k={kw[2]:.2f}$', xy=(0.97, 0.08),
                xycoords='axes fraction', ha='right', fontsize=8)
    ax.set_ylabel(r'$10^3\,\dfrac{\sigma_\mathrm{plugin}-\sigma_\mathrm{cf}}{\sigma_\mathrm{cf}}$')
    ax.set_xlabel(r'$k$  [$\mathrm{\AA}^{-1}$]')

for ax in axes[0]:
    ax.tick_params(labelbottom=False)

fail = int(worst_all > TOL)
print(f'worst relative deviation over the sweep: {worst_all:.2e} '
      f'(band +-{TOL:.1e})  [{fail} gated check failed]')
fig.savefig(OUT, bbox_inches='tight')
try:
    fig.savefig('/tmp/fig_sigma_closedforms_preview.png', dpi=110,
                bbox_inches='tight')
except Exception:
    pass
print(f'wrote {OUT}')
sys.exit(1 if fail else 0)
