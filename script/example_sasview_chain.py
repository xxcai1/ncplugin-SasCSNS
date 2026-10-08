#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# example_sasview_chain.py - the SasView -> ncmat -> NCrystal chain, validated
#                           against NCrystal's builtin @CUSTOM_HARDSPHERESANS.
#
# Workflow demonstrated:
#
#   1. SasView computes I(Q) of a sphere model in barn per OBJECT (particle).
#      Here we emulate what "SasView -> Save Data" would export for a sphere
#      model (scale=1, background=0), so this example runs without SasView
#      installed. If you have SasView, simply replace this step with its
#      exported ASCII file.
#   2. script/sasview2ncmat.py converts the exported I(Q) curve into an NCMAT
#      file with an @CUSTOM_SASCSNS DirectLoad section. The --scale flag
#      converts the SasView "barn per object" unit into the NCrystal "barn per
#      atom" unit: I_ncmat = I_sasview * n_objects/n_atoms.
#   3. The resulting NCMAT file is loaded with NCrystal and its cross section
#      is compared against the official builtin HARDSPHERESANS model for the
#      same physical system (SiO2 spheres of 100 Angstrom radius, 40% volume
#      fraction in vacuum, i.e. dry nanopowder). The physics must agree, and
#      does, but note the two normalisation conventions:
#       - the builtin model divides by the number density of the whole
#         (void-diluted) mixture, so its per-atom values are 1/(1-phi) larger
#         than those of our pure-SiO2 DirectLoad file (where the particle
#         volume fraction is folded into the scale factor). The macroscopic
#         cross section, Sigma = n_d * sigma, is the invariant to compare.
#       - the builtin factory adds the standard free-gas scattering on top of
#         the SANS process, while the DirectLoad file is SANS-only. We
#         therefore measure the standard part via ";sans=0" and subtract it.
#   4. Finally, the importance of a sufficiently large Q range is shown by
#      repeating the chain with the table truncated at Qmax = 0.45/Aa: the
#      probed neutrons have Q = 2k up to 6.2/Aa (20 meV), so the plugin must
#      extrapolate beyond the table edge and the cross sections come out
#      badly biased. Always cover Q in [0, 2k(E_max)] with the table!
# -----------------------------------------------------------------------------
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get('SASVIEW_CHAIN_WORKDIR', '/tmp/sasview_chain')

#--- physical system: SiO2 spheres, R = 100 A, 40 vol% in vacuum -------------
RADIUS = 100.0                 # Angstrom
PHI = 0.4                      # particle volume fraction (<0.5)
DENSITY = 2.2                  # g/cm3, amorphous silica
M_SIO2 = 60.084                # g/mol
B_SI, B_O = 4.1491e-5, 5.803e-5  # coherent b in Angstrom (4.1491 fm, 5.803 fm)

n_mol_cm3 = DENSITY / M_SIO2 * 6.022140857e23         # molecules per cm3
n_atoms = n_mol_cm3 * 3.0 * 1.0e-24                   # atoms per A^3
sld_sio2 = (n_atoms / 3) * (B_SI + 2 * B_O)           # 1/A^2 (~3.475e-6)
sld_void = 0.0
delta_rho = sld_sio2 - sld_void
V_p = 4.0 / 3.0 * math.pi * RADIUS**3                 # sphere volume in A^3
BARN_PER_AA2 = 1.0e8                                  # 1 barn = 1e-8 Aa^2

ENERGIES_MEV = (0.5, 1.0, 2.0, 5.0, 9.0, 20.0)
SCALE = (PHI / V_p) / n_atoms    # barn/object -> barn/atom (n_objects/n_atoms)


def sphere_form_factor(q):
    """SasView 'sphere' model P(Q) (scale=1, background=0)."""
    x = q * RADIUS
    with np.errstate(divide='ignore', invalid='ignore'):
        f = 3.0 * (np.sin(x) - x * np.cos(x)) / np.where(x == 0, 1.0, x**3)
        f[x == 0] = 1.0  # limit: F(0)=1
    return f**2


def run_chain(qgrid, tag):
    """Emulate a SasView export on the given Q grid, convert it, return path."""
    i_obj = V_p**2 * delta_rho**2 * BARN_PER_AA2 * sphere_form_factor(qgrid)
    dat = os.path.join(WORK, f'sasview_export_{tag}.dat')
    with open(dat, 'w') as fh:
        fh.write(f'# SasView emulated export: sphere model, R={RADIUS:g} A\n')
        fh.write('# Q (1/A)    I (barn/object)   dI   dQ\n')
        for qi, ii in zip(qgrid, i_obj):
            fh.write(f'{qi:.8g}\t{ii:.8g}\t0\t0\n')
    ncmat = os.path.join(WORK, f'sas_sphere_{tag}.ncmat')
    subprocess.run([sys.executable, os.path.join(HERE, 'sasview2ncmat.py'),
                    dat, '-o', ncmat, '--material', 'sio2',
                    '--density', str(DENSITY), '--radius', str(RADIUS),
                    '--scale', f'{SCALE:.8g}'], check=True,
                   stdout=subprocess.DEVNULL)
    return ncmat


def write_builtin_reference():
    void = os.path.join(WORK, 'void.ncmat')
    with open(void, 'w') as fh:
        fh.write('NCMAT v6\n@DENSITY\n  1e-12 g_per_cm3\n@DYNINFO\n'
                 '  element He\n  fraction 1\n  type freegas\n')
    builtin = os.path.join(WORK, 'builtin_hardsphere.ncmat')
    with open(builtin, 'w') as fh:
        fh.write(f"""NCMAT v6
# SiO2 with {PHI:.0%} volume of vacuum spheres, official builtin model
@DENSITY
  {DENSITY} g_per_cm3
@DYNINFO
  element Si
  fraction 0.3333333
  type freegas
@DYNINFO
  element O
  fraction 0.6666667
  type freegas
@OTHERPHASES
  {PHI} {void}
@CUSTOM_HARDSPHERESANS
  {RADIUS} #sphere radius in angstrom.
""")
    return builtin


def main():
    os.makedirs(WORK, exist_ok=True)

    #Q grids for the I(Q) table: logarithmic in the Guinier regime, then
    #linear with ~10 points per oscillation period of P(Q) (pi/R~0.031/Aa).
    #The full grid reaches past Q=2k of the highest probed energy
    #(2k(20meV) = 6.2/Aa), the truncated grid demonstrates step 4:
    q_full = np.concatenate([np.logspace(-4.0, -1.7, 120),
                             np.linspace(0.0205, 6.6, 2000)])
    q_trunc = np.concatenate([np.logspace(-4.0, -1.7, 120),
                              np.linspace(0.0205, 0.45, 220)])

    print('[1/4] emulating SasView sphere-model exports (barn per object)...')
    ncmat_full = run_chain(q_full, 'full')
    print(f'      I(Q->0) after --scale: {SCALE * V_p**2 * delta_rho**2 * BARN_PER_AA2:.6g}'
          ' barn/atom (n_objects/n_atoms conversion applied)')
    print(f'      -> {ncmat_full}')

    print('[2/4] writing official builtin-model reference material...')
    builtin = write_builtin_reference()
    print(f'      -> {builtin}')

    print('[3/4] comparing cross sections (macroscopic Sigma = n_d * sigma)...')
    import NCrystal
    NCrystal.addCustomSearchDirectory(WORK)
    dl = NCrystal.createScatter(ncmat_full)          # SANS only (plugin model)
    of = NCrystal.createScatter(builtin)             # std + SANS (builtin model)
    of_std = NCrystal.createScatter(builtin + ';sans=0')  # std part only
    n_file = n_atoms                    # number density of the DirectLoad file
    n_mix = (1.0 - PHI) * n_atoms       # number density of the void-diluted
                                        # mixture the builtin model divides by

    print('\n  E(meV)  Sigma_chain  Sigma_builtin  std_part(sans=0)  '
          'Sigma_builtin_SANS  chain/builtin_SANS')
    worst = 0.0
    for e_mev in ENERGIES_MEV:
        e = e_mev * 1e-3
        sig = lambda sc: float(sc.crossSection(e, (0, 0, 1)))
        sig_dl, sig_of, sig_std = sig(dl), sig(of), sig(of_std)
        s_chain, s_of = sig_dl * n_file, sig_of * n_mix
        s_of_sans = s_of - sig_std * n_mix
        ratio = s_chain / s_of_sans
        worst = max(worst, abs(ratio - 1.0))
        print(f'{e_mev:9.1f} {s_chain:12.4f} {s_of:13.4f} {sig_std * n_mix:15.4f} '
              f'{s_of_sans:18.4f} {ratio:19.4f}')
    print(f'\n   -> SANS macroscopic cross sections agree to {worst:.2%} (the '
          'residual is linear-interpolation error of the I(Q) table).'
          if worst < 0.02 else
          f'\n   -> MISMATCH of {worst:.2%}: investigate!')
    print('   NB: per-atom values differ by 1/(1-phi) by construction - the '
          'invariant is Sigma = n_d * sigma.')

    print('[4/4] control test: same chain but Q table truncated at 0.45/Aa...')
    ncmat_trunc = run_chain(q_trunc, 'trunc')
    dl_trunc = NCrystal.createScatter(ncmat_trunc)
    print('\n  E(meV)  2k(1/Aa)  Sigma_truncQ  Sigma_builtin_SANS   ratio')
    for e_mev in ENERGIES_MEV:
        e = e_mev * 1e-3
        k = math.sqrt(e_mev / 2.07214)  # E[meV] = 2.07214 * k[Aa^-1]^2
        sig = lambda sc: float(sc.crossSection(e, (0, 0, 1)))
        s_trunc = sig(dl_trunc) * n_file
        s_of_sans = (sig(of) - sig(of_std)) * n_mix
        print(f'{e_mev:9.1f} {2*k:9.3f} {s_trunc:13.4f} {s_of_sans:18.4f} '
              f'{s_trunc/s_of_sans:11.4f}')
    print('   -> with Qmax below 2k the plugin extrapolates beyond the table'
          ' edge and the cross section is biased. Cover [0,2k(Emax)]!')


if __name__ == '__main__':
    main()
