#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# validate_physical_closedforms.py - rung 4c: REAL physics models with
#                                    closed-form total cross sections.
#
# The synthetic special cases (constant/ramp/Gaussian tables) test the
# machinery; this rung tests it on real physics. Both standard analytic SANS
# models are closed form for BOTH DirectLoad models at ALL k (not only in
# limits):
#
#  2D (plane reduction, exact: rho' = k sin(theta)):
#      sigma_2D = 4*pi * int_0^{pi/2} I(k sin t) sin t dt
#    Guinier  I = I0 exp(-q^2 Rg^2/3):
#      sigma = 4*pi*I0 * sqrt(pi)/2 * exp(-a) * erfi(sqrt(a)) / sqrt(a),
#      a = k^2 Rg^2 / 3
#    Zimm/OZ  I = I0 / (1 + q^2 xi^2):
#      sigma = 4*pi*I0 * atanh(sqrt(a/(1+a))) / sqrt(a (1+a)),  a = k^2 xi^2
#
#  1D (true 3D reduction eq:isoxs, Q = 2k sin(theta/2), Jacobian exact):
#      sigma_1D = (2*pi/k^2) * int_0^{min(2k,qmax)} Q I(Q) dQ
#    Guinier:  sigma = (2*pi/k^2) * I0 * (3/(2 Rg^2)) * (1-exp(-Q*^2 Rg^2/3))
#    Zimm/OZ:  sigma = (2*pi/k^2) * I0 * ln(1 + xi^2 Q*^2) / (2 xi^2)
#
# Both models -> 4*pi*I0 as k -> 0 (the S2/Guinier-theorem anchor); the
# constant table of rung 4b is their a->0 plateau. NOTE the two models'
# sigma DIFFER for the same I values (the rung-7 mirror doubling): at
# k=0.5/Aa the ratio is ~1.11 for these wide patterns. The closed forms are
# themselves machine-checked here against independent quadratures before use.
#
# Checks (Guinier Rg=1.0/Aa, OZ xi=0.7/Aa, I0=2.0, tables +-1.6/Aa, so the
# 2D sampling disc (radius k <= 1.5) and the 1D interval 2k (k <= 0.8) are
# fully covered):
#   1. 2D plugin sigma vs closed form, on axis, k = 0.25/0.5/1.0/1.5,
#   2. same at tilts 30 mrad (phi=0) and 60 mrad (phi=17deg): radially
#      symmetric physics must be direction-independent,
#   3. E -> 0 clamp: the plugin serves the lowest energy node, i.e. the
#      closed form at k(0.1 meV),
#   4. 1D plugin sigma vs its closed form (2k <= qmax, the coverage rule),
#   5. mirror ratio sigma_2D/sigma_1D: plugin ratio == ratio of the two
#      closed forms (rung-7 doubling, quantified for wide patterns),
#   6. finding, documented not silent: for a 1D file violating the coverage
#      rule (2k > qmax) the plugin serves its sigma ~ E^-1/2 lookup
#      extrapolation above the boundary energy E(qmax) -- NOT the I=0
#      truncation of doc sec:edges; the row pins the measured behaviour.
# -----------------------------------------------------------------------------
import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import erfi

import NCrystal  # noqa: E402

WORK = os.environ.get('PLUGIN2D_WORKDIR', '/tmp/plugin_2d')
I0 = 2.0
RG, XI = 1.0, 0.7
HALF = 1.6
E_OVER_K2 = 2.07214
failures = 0


def i_guinier(q):
    return I0 * np.exp(-q ** 2 * RG ** 2 / 3.0)


def i_oz(q):
    return I0 / (1.0 + q ** 2 * XI ** 2)


#--- 2D closed forms (plane backbone) -----------------------------------------
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


#--- 1D closed forms (true 3D reduction; Q* = min(2k, qmax)) ------------------
def sigma1d_guinier(k):
    qs = min(2.0 * k, HALF)
    return 2.0 * math.pi / k ** 2 * I0 * (3.0 / (2.0 * RG ** 2)) \
        * (1.0 - math.exp(-qs ** 2 * RG ** 2 / 3.0))


def sigma1d_oz(k):
    qs = min(2.0 * k, HALF)
    return 2.0 * math.pi / k ** 2 * I0 \
        * math.log(1.0 + qs ** 2 * XI ** 2) / (2.0 * XI ** 2)


def backbone2d(i_func, k):
    v, _ = quad(lambda t: math.sin(t) * float(i_func(k * math.sin(t))),
                0.0, math.pi / 2.0, epsabs=1e-13, epsrel=1e-13, limit=400)
    return 4.0 * math.pi * v


def backbone1d(i_func, k):
    v, _ = quad(lambda q: q * float(i_func(q)), 0.0, min(2.0 * k, HALF),
                epsabs=1e-13, epsrel=1e-13, limit=400)
    return 2.0 * math.pi / k ** 2 * v


def write_2d(name, i_func, n=501):
    g = np.linspace(-HALF, HALF, n)
    X, Y = np.meshgrid(g, g)
    vals = i_func(np.hypot(X, Y))
    path = os.path.join(WORK, f'{name}.ncmat')
    flat = vals.ravel()
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


def write_1d(name, i_func, n=1001):
    q = np.linspace(0.0, HALF, n)
    path = os.path.join(WORK, f'{name}.ncmat')
    with open(path, 'w') as fh:
        fh.write('NCMAT v7\n@DENSITY\n  2.2 g_per_cm3\n@DYNINFO\n  element Si\n'
                 '  fraction 1\n  type freegas\n@CUSTOM_SASCSNS\n  DirectLoad\n'
                 '  Q ')
        fh.write(' '.join(f'{v:.10g}' for v in q) + '\n  I ')
        fh.write(' '.join(f'{float(i_func(qq)):.8g}' for qq in q) + '\n')
    return path


def check(name, got, exp, tol):
    global failures
    dev = abs(got - exp) / abs(exp)
    ok = dev <= tol
    failures += 0 if ok else 1
    print(f'  {name:<62s} got {got:12.6f}  exp {exp:12.6f}  '
          f'dev {dev:.2e}  tol {tol:.0e}  [{"PASS" if ok else "FAIL"}]')
    return dev


def main():
    global failures
    os.makedirs(WORK, exist_ok=True)
    print('closed forms vs independent quadratures (self-check)')
    worst_f = 0.0
    for k in (0.25, 0.5, 1.0, 1.5):
        for nm, s2, s1, i_f in (('Guinier', sigma2d_guinier(k),
                                 sigma1d_guinier(k), i_guinier),
                                ('OZ     ', sigma2d_oz(k),
                                 sigma1d_oz(k), i_oz)):
            d2 = abs(s2 / backbone2d(i_f, k) - 1.0)
            d1 = abs(s1 / backbone1d(i_f, k) - 1.0)
            worst_f = max(worst_f, d2, d1)
            print(f'  {nm} k={k:4.2f}: 2D {s2:10.6f} vs quad {backbone2d(i_f, k):10.6f} '
                  f'(dev {d2:.1e});  1D {s1:10.6f} vs quad {backbone1d(i_f, k):10.6f} '
                  f'(dev {d1:.1e})')
    if worst_f > 1e-11:
        sys.exit(f'closed-form self-check failed: {worst_f}')

    p2g, p2z = write_2d('guinier2d', i_guinier), write_2d('oz2d', i_oz)
    p1g, p1z = write_1d('guinier1d', i_guinier), write_1d('oz1d', i_oz)
    sc = {p: NCrystal.createScatter(p) for p in (p2g, p2z, p1g, p1z)}
    ks = (0.25, 0.5, 1.0, 1.5)
    ks1d = (0.25, 0.5, 0.75, 0.8)
    worst = 0.0
    for tag, p2, p1, s2cf, s1cf in (('Guinier', p2g, p1g, sigma2d_guinier,
                                     sigma1d_guinier),
                                    ('OZ     ', p2z, p1z, sigma2d_oz,
                                     sigma1d_oz)):
        print(f'{tag}: plugin sigma vs closed forms (real-physics special case)')
        for k in ks:
            e_mev = k * k * E_OVER_K2
            worst = max(worst, check(
                f'{tag} 2D k={k:.2f} (E={e_mev:.3g} meV) on-axis',
                float(sc[p2].crossSection(e_mev * 1e-3, (0, 0, 1))),
                s2cf(k), 1.5e-3))
        for beta, phi in ((0.030, 0.0), (0.060, math.radians(17.0))):
            b = (math.sin(beta) * math.cos(phi), math.sin(beta) * math.sin(phi),
                 math.cos(beta))
            worst = max(worst, check(
                f'{tag} 2D k=1.0 tilt {math.degrees(beta):.0f} mrad '
                f'toward phi={math.degrees(phi):.0f}',
                float(sc[p2].crossSection(1e-3 * E_OVER_K2, b)),
                s2cf(1.0), 5e-3))
        k_low = math.sqrt(0.1 / E_OVER_K2)   # lowest calibrated energy node
        worst = max(worst, check(
            f'{tag} 2D E->0 clamp serves lowest node (k={k_low:.4f})',
            float(sc[p2].crossSection(1e-6, (0, 0, 1))),
            s2cf(k_low), 1e-3))
        for k in ks1d:
            assert 2.0 * k <= HALF            # coverage rule respected
            e_mev = k * k * E_OVER_K2
            worst = max(worst, check(f'{tag} 1D k={k:.2f} (2k<=qmax)',
                                     float(sc[p1].crossSection(e_mev * 1e-3,
                                                               (0, 0, 1))),
                                     s1cf(k), 5e-5))
        for k in (0.5, 0.75):
            e_mev = k * k * E_OVER_K2
            r_pl = float(sc[p2].crossSection(e_mev * 1e-3, (0, 0, 1))) \
                / float(sc[p1].crossSection(e_mev * 1e-3, (0, 0, 1)))
            worst = max(worst, check(f'{tag} mirror ratio s2D/s1D k={k:.2f}',
                                     r_pl, s2cf(k) / s1cf(k), 3e-3))
        # finding: 1D beyond the coverage rule (2k > qmax) serves the
        # E^-1/2 lookup extrapolation above E(qmax) -- pinned here so any
        # change of that legacy behaviour breaks loudly instead of silently.
        k_b = 0.8                             # last k with 2k <= qmax
        e_b = k_b * k_b * E_OVER_K2
        for k in (1.0, 1.2):
            e_mev = k * k * E_OVER_K2
            extrap = float(sc[p1].crossSection(e_b * 1e-3, (0, 0, 1))) \
                * math.sqrt(e_b / e_mev)
            worst = max(worst, check(
                f'{tag} 1D beyond coverage k={k:.2f} '
                '(serves E^-1/2 extrapolation)',
                float(sc[p1].crossSection(e_mev * 1e-3, (0, 0, 1))),
                extrap, 2e-3))
    print('=' * 78)
    print(f'worst relative deviation: {worst:.2e}')
    if failures:
        sys.exit(f'{failures} check(s) FAILED')
    print('rung 4c (physical closed forms through the plugin): ALL CHECKS PASS')


if __name__ == '__main__':
    sys.exit(main())
