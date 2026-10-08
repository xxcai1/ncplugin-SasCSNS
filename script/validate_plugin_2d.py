#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# validate_plugin_2d.py - validation ladder rung 4: the C++ DirectLoad2D
#                         plugin against the validated python reference twin.
#
# The SAME special cases as validate_directload2d.py (S1 constant table, S2
# blob, S5 Gaussian) are exported to NCMAT files with @CUSTOM_SASCSNS
# DirectLoad2D sections, loaded through the real plugin factory, and pushed
# through the plugin's crossSection()/sampleScatter():
#   Pillar A: absolute cross sections vs the python quadrature (agreement
#             at the level of the sigma-grid quadrature + interpolation).
#   Pillar B: sampled outcome shapes vs the exact normalised targets
#             (binned chi2 in outcome-angle coordinates, as in rung 3).
# -----------------------------------------------------------------------------
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_directload2d import (  # noqa: E402
    Table2D, chi2_angular, k_of_e, sigma_table, make_grid)

import NCrystal  # noqa: E402

WORK = os.environ.get('PLUGIN2D_WORKDIR', '/tmp/plugin_2d')
E_OVER_K2 = 2.07214


def write_ncmat(name, qx, qy, vals, per_line=10):
    """NCMAT file with a DirectLoad2D section (Si freegas carrier; the plugin
    factory replaces the standard models, so the material is SANS-only)."""
    path = os.path.join(WORK, f'{name}.ncmat')
    qx = np.asarray(qx); qy = np.asarray(qy)
    flat = np.asarray(vals, float).ravel()  # row-major, qy outer
    with open(path, 'w') as fh:
        fh.write('NCMAT v7\n@DENSITY\n  2.2 g_per_cm3\n')
        fh.write('@DYNINFO\n  element Si\n  fraction 1\n  type freegas\n')
        fh.write('@CUSTOM_SASCSNS\n  DirectLoad2D\n  Qx ')
        fh.write(' '.join(f'{v:.10g}' for v in qx) + '\n')
        fh.write('  Qy ' + ' '.join(f'{v:.10g}' for v in qy) + '\n  I ')
        for i, v in enumerate(flat):
            fh.write(f'{v:.8g}')
            fh.write('\n  ' if (i + 1) % per_line == 0 and i + 1 < flat.size else ' ')
        fh.write('\n')
    return path


def outcomes(sc, ekin_mev, direction, n):
    _, dirs = sc.sampleScatter(ekin_mev * 1e-3, direction, repeat=n)
    return np.stack(dirs, axis=1)


def main():
    os.makedirs(WORK, exist_ok=True)
    rows = []
    np.seterr(all='ignore')

    def row(name, exp, got, kind, tol):
        rows.append((name, exp, got, kind, tol))
        if kind == 'rel':
            err = abs(got / exp - 1.0)
            ok = err < tol
            met = f'{err:.2e}'
        elif kind == 'p':
            ok = got >= exp and got == got
            met = f'{got:.4f}' if got == got else 'nan'
        else:
            err = abs(got - exp)
            ok = err <= tol
            met = f'{err:.2e}'
        print(f'{name:<66}{exp:>13.6g}{got:>13.6g}{met:>11}  '
              + ('PASS' if ok else 'FAIL'), flush=True)

    # ---------------- S1: constant table ----------------
    I0 = 3.0
    g = make_grid(2.2, 401)
    vals = np.full((401, 401), I0)
    p1 = write_ncmat('const', g, g, vals)
    sc1 = NCrystal.createScatter(p1)
    print(f'--- S1 constant table ({p1})', flush=True)
    tab1 = Table2D(g, g, vals, 'const')
    k1 = k_of_e(2.07214)
    for e_mev, d, tag, ref, tol in [
            (1.0, (0, 0, 1), 'E=1meV beam z', 4 * math.pi * I0, 1e-3),
            (25.0, (0, 0, 1), 'E=25meV beam z (disc>table: truncation)',
             sigma_table(tab1, (0, 0, 1), k_of_e(25.0)), 0.02),
            (100.0, (0, 0, 1), 'E=100meV beam z (disc>table: truncation)',
             sigma_table(tab1, (0, 0, 1), k_of_e(100.0)), 0.05),
            (2.07214, (math.sin(0.03), 0, math.cos(0.03)),
             'E=2.07meV beta=30mrad', 4 * math.pi * I0, 1e-3),
            (2.07214, (0, 0, -1), 'E=2.07meV beam -z', 4 * math.pi * I0, 1e-3),
            (2.07214, (1, 0, 0), 'E=2.07meV beam x (any tilt, s<=k)', 4 * math.pi * I0, 1e-3)]:
        row(f'S1 plugin sigma {tag}',
            ref, float(sc1.crossSection(e_mev * 1e-3, d)), 'rel', tol)
    kf = outcomes(sc1, 2.07214, (0, 0, 1), 200_000)
    row('S1 plugin elasticity max||kf|-1|', 0.0,
        float(np.abs(np.linalg.norm(kf, axis=1) - 1).max()), 'near', 1e-12)
    z = np.clip(kf[:, 2], -1, 1)
    cnt, _ = np.histogram(z, bins=20, range=(-1, 1))
    from scipy.stats import chi2 as chi2dist
    exp = np.full(20, z.size / 20)
    c2 = float(((cnt - exp) ** 2 / exp).sum())
    row('S1 plugin outcomes uniform on sphere (chi2 p)', 0.01,
        chi2dist.sf(c2, 19), 'p', None)

    # ---------------- S2: blob at low energy ----------------
    n = 401
    g = make_grid(0.8, n)
    X, Y = np.meshgrid(g, g)
    vals = 1.0 + np.exp(-((X - 0.25) ** 2 + (Y + 0.3) ** 2) / (2 * 0.12 ** 2))
    tab = Table2D(g, g, vals, 'blob')
    p2 = write_ncmat('blob', g, g, vals)
    sc2 = NCrystal.createScatter(p2)
    print(f'--- S2 blob, E->0 behaviour ({p2})', flush=True)
    ki_t = np.array([0.6, -0.5, 1.1]); ki_t /= np.linalg.norm(ki_t)
    e_mev = 0.1  # plugin's lowest energy node
    row('S2 plugin sigma(0.1meV, tilted) vs python quadrature',
        sigma_table(tab, tuple(ki_t), k_of_e(e_mev)),
        float(sc2.crossSection(e_mev * 1e-3, tuple(ki_t))), 'rel', 5e-3)
    row('S2 plugin sigma(0.1meV, tilted) vs 4*pi*I(0) (curvature ~3e-3)',
        4 * math.pi * tab.q0(),
        float(sc2.crossSection(e_mev * 1e-3, tuple(ki_t))), 'rel', 1e-2)

    # ---------------- S5: Gaussian ----------------
    s_q = 0.15
    g = make_grid(2.0, 401)
    X, Y = np.meshgrid(g, g)
    vals = np.exp(-(X ** 2 + Y ** 2) / (2 * s_q ** 2))
    tab = Table2D(g, g, vals, 'gauss')
    p3 = write_ncmat('gauss', g, g, vals)
    sc3 = NCrystal.createScatter(p3)
    print(f'--- S5 Gaussian ({p3})', flush=True)
    for e_mev, tol in [(1.0, 5e-3), (2.07214, 5e-3), (5.0, 5e-3),
                       (20.0, 0.02), (55.0, 0.05)]:
        row(f'S5 plugin sigma({e_mev}meV, z) vs python quadrature'
            + ('' if tol < 0.01 else ' [disc>table]'),
            sigma_table(tab, (0, 0, 1), k_of_e(e_mev)),
            float(sc3.crossSection(e_mev * 1e-3, (0, 0, 1))), 'rel', tol)
    beta = 0.03
    ki_t = (math.sin(beta), 0.0, math.cos(beta))
    #off-node in both s and lnE -> combined (s,E) grid interpolation error,
    #documented; on-axis rows carry the E-curvature alone (~3e-3).
    row('S5 plugin sigma(2.07meV, 30mrad) vs python quadrature',
        sigma_table(tab, ki_t, 1.0),
        float(sc3.crossSection(2.07214e-3, ki_t)), 'rel', 1e-2)
    row('S5 plugin sigma(2.07meV, -z) vs python quadrature',
        sigma_table(tab, (0, 0, -1), 1.0),
        float(sc3.crossSection(2.07214e-3, (0, 0, -1))), 'rel', 5e-3)

    kf = outcomes(sc3, 2.07214, (0, 0, 1), 200_000)
    row('S5 plugin on-axis shape chi2 p (200k)', 1e-3,
        chi2_angular(tab, (0, 0, 1), 1.0, kf, nbins=24)[1], 'p', None)
    kf = outcomes(sc3, 2.07214, ki_t, 150_000)
    row('S5 plugin tilted 30mrad shape chi2 p (150k)', 1e-3,
        chi2_angular(tab, ki_t, 1.0, kf, nbins=24)[1], 'p', None)
    ki_t1 = (math.sin(0.002), 0.0, math.cos(0.002))
    kf = outcomes(sc3, 2.07214, ki_t1, 120_000)
    row('S5 plugin tier-1 tilt 2mrad shape chi2 p (120k)', 1e-3,
        chi2_angular(tab, ki_t1, 1.0, kf, nbins=24)[1], 'p', None)
    kf = outcomes(sc3, 3.3, (0, 0, 1), 120_000)  # between energy nodes
    row('S5 plugin between-nodes E=3.3meV shape chi2 p (120k)', 1e-3,
        chi2_angular(tab, (0, 0, 1), k_of_e(3.3), kf, nbins=24)[1], 'p', None)

    print('\n' + '=' * 100)
    nfail = sum(1 for name, exp, got, kind, tol in rows
                if (kind == 'rel' and not abs(got / exp - 1.0) < tol)
                or (kind == 'p' and not (got >= exp and got == got))
                or (kind == 'near' and not abs(got - exp) <= tol))
    print('ALL PASS' if nfail == 0 else f'{nfail} FAILURE(S)')
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
