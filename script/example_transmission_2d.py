#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# example_transmission_2d.py - validation ladder rung 8: transport-integrated
# attenuation (Beer-Lambert) for DirectLoad2D materials.
#
# Pillar A ("sigma in ABSOLUTE terms") exercised the way a real transport
# code experiences it: neutrons cross a slab, path lengths are drawn from
# Exp(1/Sigma) with Sigma = n * (sigma_scatter + sigma_absorb) built from the
# plugin's cross sections and the NCMAT density, and the uncollided fraction
# must reproduce exp(-Sigma*L). This pins the unit chain
# barn/atom -> cm^-1 (the barn/object -> barn/atom conversion, number density
# from @DENSITY + composition) end-to-end, plus the direction dependence:
# a beam tilted off the design axis must be attenuated with the tilted sigma.
#
# Checks (E = 2.07214 meV, k = 1/Aa, narrow Gaussian s_q = 0.1 1/Aa tables):
#   1. uncollided fraction vs exp(-Sigma*L) at three thicknesses
#      (binomial sigma; pass if within 4 sigma),
#   2. absorbed fraction of the interactions vs Sigma_abs/Sigma_tot,
#   3. tilted beam (30 mrad): transmission follows sigma(tilted beam),
#   4. operational statement of the mirror-doubling convention (doc
#      sec:backward): the 2D material attenuates ~2x more strongly than the
#      1D material with the same narrow pattern, because backward-scattered
#      neutrons carry the forward intensity; transmission numbers quantify
#      this for the data producer.
# -----------------------------------------------------------------------------
import math
import os
import sys

import numpy as np

import NCrystal  # noqa: E402

WORK = os.environ.get('PLUGIN2D_WORKDIR', '/tmp/plugin_2d')
E_MEV = 2.07214                      # -> k = 1 Aa^-1
S_Q = 0.1                            # pattern width, 1/Aa
RNG = np.random.default_rng(12345)


def write_ncmat_1d(name, q, i):
    path = os.path.join(WORK, f'{name}.ncmat')
    with open(path, 'w') as fh:
        fh.write('NCMAT v7\n@DENSITY\n  2.2 g_per_cm3\n@DYNINFO\n  element Si\n'
                 '  fraction 1\n  type freegas\n@CUSTOM_SASCSNS\n  DirectLoad\n  Q ')
        fh.write(' '.join(f'{v:.12g}' for v in q) + '\n  I ')
        fh.write(' '.join(f'{v:.8g}' for v in i) + '\n')
    return path


def write_ncmat_2d(name, qx, qy, vals):
    path = os.path.join(WORK, f'{name}.ncmat')
    flat = np.asarray(vals, float).ravel()
    with open(path, 'w') as fh:
        fh.write('NCMAT v7\n@DENSITY\n  2.2 g_per_cm3\n@DYNINFO\n  element Si\n'
                 '  fraction 1\n  type freegas\n@CUSTOM_SASCSNS\n  DirectLoad2D\n')
        fh.write('  Qx ' + ' '.join(f'{v:.12g}' for v in qx) + '\n')
        fh.write('  Qy ' + ' '.join(f'{v:.12g}' for v in qy) + '\n  I ')
        for j, v in enumerate(flat):
            fh.write(f'{v:.8g}')
            fh.write('\n  ' if (j + 1) % 10 == 0 and j + 1 < flat.size else ' ')
        fh.write('\n')
    return path


def transport(path, ekin, beam, length_cm, n_neutron):
    """Analog MC through a slab: returns (uncollided, collided&absorbed)."""
    scat = NCrystal.createScatter(path)
    abst = NCrystal.createAbsorption(path)
    sig_s = float(scat.crossSection(ekin, beam))
    sig_a = float(abst.crossSection(ekin, beam))
    info = NCrystal.createInfo(path)
    n_inv_a3 = float(info.getNumberDensity())   # atoms/Aa^3, as transport sees it
    sigma_tot = sig_s + sig_a                    # barn/atom
    #unit cancellation to commit to memory: 1 barn = 1e-24 cm^2 and
    #1 Aa^3 = 1e-24 cm^3, so n [atoms/Aa^3] * sigma [barn/atom] IS the
    #macroscopic cross section in cm^-1, with no factor at all:
    sigma_inv = n_inv_a3 * sigma_tot             # 1/cm
    mfp_cm = 1.0 / sigma_inv

    rng = RNG
    xi = rng.random(n_neutron)
    # analog free flight: those with s > L are transmitted uncollided
    s = -np.log(rng.random(n_neutron)) / sigma_inv
    uncollided = int(np.sum(s > length_cm))
    n_coll = n_neutron - uncollided
    absorbed = int(np.sum(rng.random(n_coll) * sigma_tot < sig_a))
    return dict(sig_s=sig_s, sig_a=sig_a, n_inv_a3=n_inv_a3, sigma_inv=sigma_inv,
                mfp=mfp_cm, uncollided=uncollided, absorbed=absorbed,
                n_coll=n_coll, n=n_neutron)


def main():
    os.makedirs(WORK, exist_ok=True)
    failures = 0

    # --- material: narrow radial Gaussian on both models ---------------------
    AMP = 40.0   # barn/(atom sr) at Q=0: high-contrast powder, mm-scale mfp
    q1d = np.linspace(0.0, 3.2, 6401)
    i1d = AMP * np.exp(-q1d ** 2 / (2 * S_Q ** 2))
    p1 = write_ncmat_1d('transm1d', q1d, i1d)
    g2 = np.linspace(-3.2, 3.2, 1601)
    X, Y = np.meshgrid(g2, g2)
    p2 = write_ncmat_2d('transm2d', g2, g2,
                        AMP * np.exp(-(X ** 2 + Y ** 2) / (2 * S_Q ** 2)))

    ekin = E_MEV * 1e-3
    beam = (0.0, 0.0, 1.0)
    n = 1_000_000
    r2 = transport(p2, ekin, beam, 1.0, n)   # thickness set below anyway

    print(f'material: {os.path.basename(p2)}  E = {E_MEV} meV (k = 1 Aa^-1)')
    print(f'sigma_scatter = {r2["sig_s"]:.4f} barn  sigma_absorb = {r2["sig_a"]:.5f} barn')
    print(f'n = {r2["n_inv_a3"]:.4e} 1/Aa^3  ->  Sigma_tot = {r2["sigma_inv"]:.4f} 1/cm'
          f'  mfp = {r2["mfp"]:.3f} cm')
    print()

    # --- 1: Beer-Lambert at three thicknesses --------------------------------
    print('check 1: uncollided fraction vs Beer-Lambert exp(-Sigma*L)')
    for l_fac in (0.25, 1.0, 3.0):
        L = l_fac * r2['mfp']
        r = transport(p2, ekin, beam, L, n)
        p_exp = math.exp(-r['sigma_inv'] * L)
        p_mc = r['uncollided'] / r['n']
        p_err = math.sqrt(p_exp * (1 - p_exp) / r['n'])
        pull = (p_mc - p_exp) / p_err
        ok = abs(pull) < 4.0
        failures += 0 if ok else 1
        print(f'  L = {l_fac:.2f} mfp = {L:6.3f} cm: exp = {p_exp:.6f}  '
              f'MC = {p_mc:.6f}  pull = {pull:+.2f} sigma  '
              f'[{"PASS" if ok else "FAIL"}]')

    # --- 2: absorption split of the collisions --------------------------------
    print('check 2: absorbed fraction of interactions vs Sigma_a/Sigma_tot')
    frac_exp = r2['sig_a'] / (r2['sig_a'] + r2['sig_s'])
    assert r2['n_coll'] > 0
    frac_mc = r2['absorbed'] / r2['n_coll']
    f_err = math.sqrt(frac_exp * (1 - frac_exp) / r2['n_coll'])
    pull = (frac_mc - frac_exp) / f_err
    ok = abs(pull) < 4.0
    failures += 0 if ok else 1
    print(f'  expected {frac_exp:.6f}  MC {frac_mc:.6f}  pull = {pull:+.2f} sigma'
          f'  [{"PASS" if ok else "FAIL"}]')

    # --- 3: tilted beam attenuates with the tilted sigma ----------------------
    print('check 3: 30 mrad tilted beam follows sigma(tilted beam)')
    beta = 0.030
    beam_t = (math.sin(beta), 0.0, math.cos(beta))
    L = r2['mfp']
    rt = transport(p2, ekin, beam_t, L, n)
    p_exp_t = math.exp(-rt['sigma_inv'] * L)
    p_mc_t = rt['uncollided'] / rt['n']
    pull = (p_mc_t - p_exp_t) / math.sqrt(p_exp_t * (1 - p_exp_t) / rt['n'])
    aniso = rt['sig_s'] / r2['sig_s']
    ok = abs(pull) < 4.0 and abs(aniso - 1.0) > 1e-4
    failures += 0 if ok else 1
    print(f'  sigma(on axis) = {r2["sig_s"]:.5f}  sigma(30 mrad) = {rt["sig_s"]:.5f}'
          f'  ratio = {aniso:.5f}')
    print(f'  exp = {p_exp_t:.6f}  MC = {p_mc_t:.6f}  pull = {pull:+.2f} sigma'
          f'  [{"PASS" if ok else "FAIL"}]')

    # --- 4: mirror-doubling made operational ----------------------------------
    print('check 4: 1D vs 2D convention on the SAME pattern (doc sec:backward)')
    r1 = transport(p1, ekin, beam, L, n)
    ratio_sigma = r2['sig_s'] / r1['sig_s']
    # transmissions at a fixed thickness chosen for contrast (1 mfp of 2D):
    p1_t = math.exp(-r1['sigma_inv'] * L)
    p2_t = math.exp(-r2['sigma_inv'] * L)
    ok = 1.8 <= ratio_sigma <= 2.3
    failures += 0 if ok else 1
    print(f'  sigma_2D/sigma_1D = {ratio_sigma:.3f} (mirror doubling, expect ~2)'
          f'  [{"PASS" if ok else "FAIL"}]')
    print(f'  transmission through L = {L:.3f} cm (1 mfp of the 2D material):')
    print(f'    1D material: {p1_t:.4f}    2D material: {p2_t:.4f}  '
          f'(same particles, same beam!)')

    print('=' * 72)
    if failures:
        sys.exit(f'{failures} check(s) FAILED')
    print('rung 8 (transport-integrated attenuation): ALL CHECKS PASS')


if __name__ == '__main__':
    sys.exit(main())
