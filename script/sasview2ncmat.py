#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# sasview2ncmat.py - convert a SasView-exported I(Q) curve (1D) or I(Qx,Qy)
#                    image (2D, --2d) into an NCMAT file that uses the
#                    DirectLoad / DirectLoad2D mode of the SasCSNS NCrystal
#                    plugin (@CUSTOM_SASCSNS section).
#
# SasView (https://www.sasview.org) can export 1D fitting data as plain ASCII
# files with columns like "Q (1/A)  I (1/cm)  dI  dQ" (header lines start
# with "#"). This script converts such a file into NCMAT v6 data:
#
#   python3 sasview2ncmat.py sasview_iq.dat -o sas_model.ncmat \
#           --material sio2 --density 2.2 --phi 0.4 --radius 100.0
#
# Unit conversion (ALWAYS applied - there is deliberately no flag for it):
#   SasView reports intensity in barn per OBJECT (particle), while NCrystal
#   expects barn per ATOM. The converter therefore multiplies every I value
#   by scale = n_objects / n_atoms = (phi / V_p) / n_d, where phi is the
#   particle volume fraction, V_p the particle volume and n_d the atomic
#   number density of the file material (computed from --material and
#   --density). The scale used is recorded in the output file header.
#
# 2D mode (--2d): the input must be an ASCII file with rows "Qx Qy I"
#   (SasView 2D export: columns Qx, Qy in 1/Aa and I in 1/cm; header lines
#   start with '#'). The rows must form a complete rectangular grid with
#   ascending axes; non-uniform axes are resampled onto uniform grids (the
#   DirectLoad2D model requires uniform grids). Qx and Qy are the in-plane
#   components of the scattering vector for a beam along the +z axis of the
#   material frame (see doc/anisotropic_directload2d.pdf, section 5.1).
#
# The resulting NCMAT file requires the SasCSNS plugin:
#   nctool sas_model.ncmat
# -----------------------------------------------------------------------------
import argparse
import math
import sys

AVOGADRO = 6.02214076e23  # 1/mol


def read_sasview_iq(filename):
    """Read a SasView-exported ASCII I(Q) file, return sorted [(Q, I)] list."""
    points = []
    with open(filename) as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.replace(',', ' ').split()
            if len(parts) < 2:
                continue
            try:
                q, i = float(parts[0]), float(parts[1])
            except ValueError:
                continue  # column-title or metadata line
            points.append((q, i))
    if len(points) < 2:
        sys.exit(f'ERROR: found fewer than 2 (Q,I) points in {filename}')
    # Sort ascending in Q, drop NaN/negative-I points and duplicated Q values:
    bad = [(q, i) for q, i in points if not (q > 0.0 and i > 0.0)]
    if bad:
        print(f'WARNING: dropped {len(bad)} of {len(points)} data points with '
              f'Q<=0, I<=0 or non-finite values (first: Q={bad[0][0]!r} I={bad[0][1]!r})')
    points = sorted((q, i) for q, i in points if q > 0.0 and i > 0.0)
    deduped = []
    n_dup = 0
    for q, i in points:
        if deduped and abs(q - deduped[-1][0]) < 1.0e-12:
            n_dup += 1
            continue
        deduped.append((q, i))
    if n_dup:
        print(f'WARNING: dropped {n_dup} duplicated Q values')
    return deduped


def read_sasview_iqxy(filename):
    """Read a SasView-exported ASCII I(Qx,Qy) file.

    Returns (qx_list, qy_list, vals) with vals[qy_index][qx_index] (row-major,
    qy outer). Rows must form a complete rectangular grid with ascending
    axes; non-uniform axes are resampled onto uniform grids.
    """
    rows = []
    with open(filename) as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.replace(',', ' ').split()
            if len(parts) < 3:
                continue
            try:
                qx, qy, i = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue  # column-title or metadata line
            rows.append((qx, qy, i))
    if len(rows) < 4:
        sys.exit(f'ERROR: found fewer than 4 (Qx,Qy,I) rows in {filename}')
    qxs = sorted({r[0] for r in rows})
    qys = sorted({r[1] for r in rows})
    if len(rows) != len(qxs) * len(qys):
        sys.exit(f'ERROR: rows do not form a complete rectangular grid '
                 f'({len(rows)} rows vs {len(qxs)}x{len(qys)} axes points). '
                 f'Masked or missing pixels must be filled before conversion.')
    lookup = {(r[0], r[1]): r[2] for r in rows}
    n_bad = sum(1 for v in lookup.values()
                if v != v or v in (float('inf'), float('-inf')))
    if n_bad:
        sys.exit(f'ERROR: {n_bad} non-finite I values in {filename} '
                 f'(masked pixels must be filled before conversion)')
    vals = [[lookup[(qx, qy)] for qx in qxs] for qy in qys]
    n_neg = sum(1 for row in vals for v in row if v < 0.0)
    if n_neg:
        print(f'WARNING: clamped {n_neg} negative I values to 0 '
              f'(noise at high Q; the model requires I >= 0)')
        vals = [[max(v, 0.0) for v in row] for row in vals]
    qx_u = _resample_axis_uniform(qxs, 'Qx')
    qy_u = _resample_axis_uniform(qys, 'Qy')
    if qx_u is not None or qy_u is not None:
        print('WARNING: non-uniform axis detected, resampling onto a uniform '
              'grid (DirectLoad2D requires uniform axes)')
        vals = _bilinear_resample(qxs, qys, vals,
                                  qx_u if qx_u is not None else qxs,
                                  qy_u if qy_u is not None else qys)
    return (qx_u if qx_u is not None else qxs,
            qy_u if qy_u is not None else qys, vals)


def _resample_axis_uniform(ax, name):
    """Return a uniform version of an ascending axis, or None if uniform."""
    n = len(ax)
    if n < 2:
        sys.exit(f'ERROR: {name} axis needs at least 2 points')
    d0 = ax[1] - ax[0]
    if d0 <= 0.0:
        sys.exit(f'ERROR: {name} axis must be strictly ascending')
    if all(abs((ax[i] - ax[i - 1]) - d0) <= 1.0e-6 * abs(d0) for i in range(2, n)):
        return None
    step = (ax[-1] - ax[0]) / (n - 1)
    return [ax[0] + i * step for i in range(n)]


def _interp1(x0, x1, v0, v1, x):
    if x1 == x0:
        return v0
    t = (x - x0) / (x1 - x0)
    return v0 + t * (v1 - v0)


def _bilinear_resample(qxs, qys, vals, qx_u, qy_u):
    """Resample values defined on (qxs,qys) onto (qx_u,qy_u) (bilinear).

    The input data is treated as the bilinear interpolation of its grid;
    resampling along each axis in turn evaluates that interpolant.
    """
    import bisect

    def axis_weights(axis, target):
        # segment index and fraction for each target coordinate:
        out = []
        for x in target:
            j = min(max(bisect.bisect_right(axis, x) - 1, 0), len(axis) - 2)
            out.append((j, (x - axis[j]) / (axis[j + 1] - axis[j])))
        return out

    wx = axis_weights(qxs, qx_u)
    wy = axis_weights(qys, qy_u)
    # first interpolate along x for each input y row:
    inter = []
    for row in vals:
        newrow = []
        for (j, t) in wx:
            v0, v1 = row[j], row[j + 1]
            newrow.append(v0 + t * (v1 - v0))
        inter.append(newrow)
    # then along y:
    out = []
    for (i, t) in wy:
        r0, r1 = inter[i], inter[i + 1]
        out.append([v0 + t * (v1 - v0) for v0, v1 in zip(r0, r1)])
    return out


def write_ncmat_2d(out, qx, qy, vals, material, density, solvent, scale):
    composition = resolve_material(material)
    lines = ['NCMAT v7',
             '#',
             '# Generated by sasview2ncmat.py --2d from a SasView-exported '
             'I(Qx,Qy) image.',
             '# Scattering model: SasCSNS plugin, DirectLoad2D mode '
             '(tabulated anisotropic 2D SANS).',
             '# Convention: Qx,Qy are the in-plane components of the '
             'scattering vector for a beam',
             '# along the +z axis of the material frame; elastic scattering; '
             'I=0 outside the table.',
             '#',
             '# Unit conversion applied (SasView barn/object -> NCrystal barn/atom):',
             f'#   I_ncmat = I_sasview * {scale:.8g}',
             '#   scale = (phi / V_p) / n_d',
             '#']
    lines += ['',
              '@DENSITY',
              f'  {density} g_per_cm3']
    for element, frac in composition:
        lines += ['@DYNINFO',
                  f'  element {element}',
                  f'  fraction {frac}',
                  '  type freegas']
    #12 significant digits on the axes: the DirectLoad2D parser enforces
    #uniform spacing with a 1e-6*|d0| tolerance, which %.8g can violate on
    #finer grids (rounding shifts consecutive pixel centres by ~1e-8 -- right
    #at the tolerance for 400-point axes):
    lines += ['@CUSTOM_SASCSNS',
              '  DirectLoad2D',
              '  Qx ' + ' '.join(f'{q:.12g}' for q in qx),
              '  Qy ' + ' '.join(f'{q:.12g}' for q in qy)]
    # I values, row major with qy as the outer axis; start on the I line:
    flat = ['%.7g' % (v * scale) for row in vals for v in row]
    per_line = 10
    lines.append('  I ' + ' '.join(flat[:per_line]))
    for b in range(per_line, len(flat), per_line):
        lines.append('    ' + ' '.join(flat[b:b + per_line]))
    if solvent:
        lines.append(f'  solvent {solvent}')
    lines.append('')
    with open(out, 'w') as fh:
        fh.write('\n'.join(lines))


MATERIAL_PRESETS = {
    # name: (list of (element, fraction))
    'sio2': [('Si', 0.3333333), ('O', 0.6666667)],
    'al2o3': [('Al', 0.4), ('O', 0.6)],
    'fe3o4': [('Fe', 0.4285714), ('O', 0.5714286)],
    'air': [('N', 0.76), ('O', 0.24)],
}

ATOMIC_MASS = {  # g/mol (standard atomic weights)
    'Si': 28.085, 'O': 15.999, 'Al': 26.9815, 'Fe': 55.845, 'N': 14.007,
}


def number_density(material, density):
    """Atomic number density [1/Aa^3] of a material preset at given density."""
    composition = resolve_material(material)
    mass_batch = sum(frac * ATOMIC_MASS[el] for el, frac in composition)
    atoms_batch = sum(frac for _, frac in composition)
    # n_d = rho * N_A * atoms_batch / M_batch, converted cm^-3 -> Aa^-3:
    return density * AVOGADRO * atoms_batch / (mass_batch * 1.0e24)


def resolve_material(material):
    try:
        return MATERIAL_PRESETS[material.lower()]
    except KeyError:
        sys.exit(f'ERROR: unknown material {material!r} '
                 f'(known: {", ".join(sorted(MATERIAL_PRESETS))})')


def write_ncmat(out, points, material, density, radius, solvent, scale):
    composition = resolve_material(material)
    lines = ['NCMAT v6',
             '#',
             '# Generated by sasview2ncmat.py from a SasView-exported I(Q) curve.',
             '# Scattering model: SasCSNS plugin, DirectLoad mode'
             + (f' with solvent {solvent}' if solvent else ' (vacuum solvent)'),
             '#',
             f'# Unit conversion applied (SasView barn/object -> NCrystal barn/atom):',
             f'#   I_ncmat = I_sasview * {scale:.8g}',
             f'#   scale = (phi / V_p) / n_d'
             + (f' with V_p = 4/3*pi*R^3, R = {radius:g} Aa (sphere)' if radius else ''),
             '#']
    for q, i in points:
        lines.append(f'#   Q={q:<12.6g} I={i * scale:<12.6g}')
    lines += ['',
              f'@DENSITY',
              f'  {density} g_per_cm3']
    for element, frac in composition:
        lines += ['@DYNINFO',
                  f'  element {element}',
                  f'  fraction {frac}',
                  '  type freegas']
    #%.12g on Q: the DirectLoad parser's 1e-6*|d0| uniformity tolerance
    #makes %.6g unsafe for fine grids (see write_ncmat_2d):
    lines += ['@CUSTOM_SASCSNS',
              '  DirectLoad',
              '  Q ' + ' '.join(f'{q:.12g}' for q, _ in points),
              '  I ' + ' '.join(f'{i * scale:.8g}' for _, i in points)]
    if solvent:
        lines.append(f'  solvent {solvent}')
    lines.append('')
    with open(out, 'w') as fh:
        fh.write('\n'.join(lines))


def main():
    ap = argparse.ArgumentParser(
        description='Convert a SasView-exported I(Q) curve into an NCMAT file '
                    'for the SasCSNS NCrystal plugin (DirectLoad mode). The '
                    'barn/object -> barn/atom unit conversion is always applied, '
                    'computed from the mandatory physical parameters below.')
    ap.add_argument('sasview_file', help='SasView ASCII data file (1D: Q I '
                    '[dI] [dQ] columns; --2d: Qx Qy I rows)')
    ap.add_argument('--2d', dest='mode2d', action='store_true',
                    help='input is a 2D I(Qx,Qy) image -> DirectLoad2D model '
                         '(anisotropic); see header of this script')
    ap.add_argument('-o', '--output', default='sas_model.ncmat',
                    help='output NCMAT file name (default: %(default)s)')
    ap.add_argument('--material', default='sio2',
                    help='particle material preset: ' + ', '.join(sorted(MATERIAL_PRESETS))
                         + ' (default: %(default)s)')
    ap.add_argument('--density', type=float, default=2.2,
                    help='particle material density in g/cm3 (default: %(default)s)')
    ap.add_argument('--phi', type=float, required=True,
                    help='particle volume fraction in the sample (mandatory)')
    vgroup = ap.add_mutually_exclusive_group(required=True)
    vgroup.add_argument('--radius', type=float,
                        help='particle radius in Angstrom, for spherical particles '
                             '(V_p = 4/3*pi*R^3)')
    vgroup.add_argument('--volume', type=float,
                        help='particle volume in Angstrom^3, for non-spherical particles')
    ap.add_argument('--emax', type=float, default=None,
                    help='highest neutron energy [meV] the file must cover. '
                         'If given and Qmax < 2k(Emax), abort with an error '
                         '(a truncated Q table silently biases cross '
                         'sections; see example_sasview_chain.py step 4)')
    ap.add_argument('--solvent', default='',
                    help='solvent material name, stored as a comment for '
                         'reference (DirectLoad I(Q) already includes solvent '
                         'contrast; default: none)')
    args = ap.parse_args()

    if not 0.0 < args.phi < 1.0:
        sys.exit(f'ERROR: --phi must be in (0,1), got {args.phi}')

    n_d = number_density(args.material, args.density)
    if args.radius is not None:
        if args.radius <= 0.0:
            sys.exit(f'ERROR: --radius must be positive, got {args.radius}')
        v_p = 4.0 / 3.0 * math.pi * args.radius**3
    else:
        if args.volume <= 0.0:
            sys.exit(f'ERROR: --volume must be positive, got {args.volume}')
        v_p = args.volume
    scale = (args.phi / v_p) / n_d  # barn/object -> barn/atom

    if args.mode2d:
        main_2d(args)
        return

    points = read_sasview_iq(args.sasview_file)

    # --- data-file validation (a truncated or coarse table biases results
    #     silently - always report what the file can and cannot describe) ---
    qmax = points[-1][0]
    kmax = 0.5 * qmax                      # largest k with 2k <= Qmax
    e_cover = 2.07214 * kmax**2            # [meV] highest fully covered energy
    print(f'File validation: Qmax = {qmax:.6g} 1/Aa -> fully covers neutron '
          f'energies up to E = {e_cover:.4g} meV (Q = 2k rule)')
    if args.emax is not None and qmax < 2.0 * math.sqrt(args.emax / 2.07214):
        sys.exit(f'ERROR: Qmax = {qmax:.6g} 1/Aa does not cover --emax '
                 f'{args.emax:g} meV (needs Q >= {2.0 * math.sqrt(args.emax / 2.07214):.6g} '
                 f'1/Aa). Extend the table or lower --emax.')
    if len(points) > 1:
        dq = min(b[0] - a[0] for a, b in zip(points, points[1:]))
        if args.radius:
            dq_needed = math.pi / (10.0 * args.radius)  # >=10 pts per oscillation
            if dq < dq_needed:
                print(f'File validation: min spacing {dq:.3g} 1/Aa resolves '
                      f'form-factor oscillations of R = {args.radius:g} Aa '
                      f'(need <= {dq_needed:.3g}) OK')
            else:
                print(f'WARNING: min Q spacing {dq:.3g} 1/Aa is too coarse to '
                      f'resolve form-factor oscillations of R = {args.radius:g} Aa '
                      f'(want >= 10 points per pi/R period, i.e. spacing <= '
                      f'{dq_needed:.3g} 1/Aa)')
        if points[0][0] > 0.05 * qmax:
            print(f'WARNING: table starts at Q = {points[0][0]:.3g} 1/Aa; the '
                  f'forward region below it is not tabulated (SANS intensity '
                  f'peaks at Q = 0)')

    write_ncmat(args.output, points, args.material, args.density,
                args.radius, args.solvent, scale)
    print(f'Wrote {args.output} with {len(points)} (Q,I) points from '
          f'{args.sasview_file}')
    print(f'Unit conversion (always applied): scale = (phi/V_p)/n_d = '
          f'({args.phi:g}/{v_p:.6g})/{n_d:.6g} = {scale:.8g} '
          f'[barn/object -> barn/atom]')
    print('Load it with NCrystal, e.g.:  nctool '
          f'{args.output} -x 1e-5 1.0 -p')


def main_2d(args):
    """The --2d code path: I(Qx,Qy) image -> DirectLoad2D NCMAT file."""
    if not 0.0 < args.phi < 1.0:
        sys.exit(f'ERROR: --phi must be in (0,1), got {args.phi}')
    if args.radius is not None and args.volume is not None:
        sys.exit('ERROR: --radius and --volume are mutually exclusive')
    n_d = number_density(args.material, args.density)
    if args.radius is not None:
        if args.radius <= 0.0:
            sys.exit(f'ERROR: --radius must be positive, got {args.radius}')
        v_p = 4.0 / 3.0 * math.pi * args.radius**3
    elif args.volume is not None:
        if args.volume <= 0.0:
            sys.exit(f'ERROR: --volume must be positive, got {args.volume}')
        v_p = args.volume
    else:
        sys.exit('ERROR: --2d mode needs --volume (or --radius) so the '
                 'barn/object -> barn/atom conversion can be computed')
    scale = (args.phi / v_p) / n_d

    qx, qy, vals = read_sasview_iqxy(args.sasview_file)
    qmax = min(qx[-1], qy[-1])   # square tables: each half-side must cover 2k
    kmax = 0.5 * qmax
    e_cover = 2.07214 * kmax**2
    print(f'File validation: table {len(qx)}x{len(qy)} over '
          f'[{qx[0]:.6g},{qx[-1]:.6g}]x[{qy[0]:.6g},{qy[-1]:.6g}] 1/Aa; '
          f'fully covers neutron energies up to E = {e_cover:.4g} meV '
          f'(disc radius 2k must fit inside the table)')
    if args.emax is not None and qmax < 2.0 * math.sqrt(args.emax / 2.07214):
        sys.exit(f'ERROR: table half-width {qmax:.6g} 1/Aa does not cover '
                 f'--emax {args.emax:g} meV (needs >= '
                 f'{2.0 * math.sqrt(args.emax / 2.07214):.6g} 1/Aa). Extend '
                 f'the table or lower --emax.')
    if min(len(qx), len(qy)) < 8:
        print(f'WARNING: only {len(qx)}x{len(qy)} pixels; anisotropic '
              f'structure is poorly resolved on such a coarse grid')
    if qx[0] > 0.0 or qy[0] > 0.0:
        print(f'WARNING: table starts at Qx={qx[0]:.6g}, Qy={qy[0]:.6g} '
              f'-- it does not contain Q=0, so the forward SANS intensity '
              f'(dominant at low neutron energies) is amputated; the model '
              f'interpolates from the first pixel inward. Export grids '
              f'should start at Q=0.')
    write_ncmat_2d(args.output, qx, qy, vals, args.material, args.density,
                   args.solvent, scale)
    print(f'Wrote {args.output} with {len(qx) * len(qy)} I(Qx,Qy) values '
          f'from {args.sasview_file}')
    print(f'Unit conversion (always applied): scale = (phi/V_p)/n_d = '
          f'({args.phi:g}/{v_p:.6g})/{n_d:.6g} = {scale:.8g} '
          f'[barn/object -> barn/atom]')
    print('Load it with NCrystal, e.g.:  nctool '
          f'{args.output} -x 1e-5 1.0 -p')


if __name__ == '__main__':
    main()
