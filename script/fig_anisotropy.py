#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig_anisotropy.py - the anisotropy of the DirectLoad2D model at a glance
# (four panels, doc/fig_anisotropy.pdf):
#   (a) anisotropic input table I(Qx,Qy): elliptical Gaussian, axis ratio 4:1,
#       major axis at 30 deg in the detector plane,
#   (b) (Qx,Qy) density of 2M outcomes sampled by the C++ plugin from it --
#       the sampler reproduces the anisotropic pattern (shape level),
#   (c) the anisotropy a detector ring sees: d sigma / d psi at theta = 10 deg
#       (ring |Q_perp| = k sin theta) -- sampled points against the exact
#       prediction I(Q(psi)), modulation depth ~300x here,
#   (d) counterpoint: the INTEGRATED sigma is direction-independent --
#       relative change along the two principal axes over the 60 mrad design
#       cone, plugin (dots) vs an INDEPENDENT midpoint quadrature (lines);
#       both sit inside the tabulated model's ~1% numerical floor.
# -----------------------------------------------------------------------------
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get('ANISO_FIG_WORKDIR', '/tmp/aniso_fig')
OUT = os.path.join(HERE, '..', 'doc', 'fig_anisotropy.pdf')

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

#--- material: tilted 4:1 elliptical Gaussian ---------------------------------
AMP = 40.0            # barn/(atom sr) at the peak
S_MAJ, S_MIN = 0.20, 0.05   # 1/Aa, axis ratio 4:1
TAU = math.radians(30.0)    # major-axis azimuth in the plane
EXT, NPTS = 1.6, 801        # table: +-1.6 1/Aa on an 801^2 grid (pitch 0.004)
E_MEV = 2.07214             # k = 1 /Aa
K = 1.0

os.makedirs(WORK, exist_ok=True)
g = np.linspace(-EXT, EXT, NPTS)
PITCH = g[1] - g[0]
X, Y = np.meshgrid(g, g)                 # [ny,nx]: x inner (j,i) = (qy,qx)
U = X * math.cos(TAU) + Y * math.sin(TAU)
V = -X * math.sin(TAU) + Y * math.cos(TAU)
I_TAB = AMP * np.exp(-0.5 * ((U / S_MAJ) ** 2 + (V / S_MIN) ** 2))

ncmat = os.path.join(WORK, 'aniso4to1.ncmat')
flat = I_TAB.ravel()                     # row-major, qy outer (grammar: I[qy*Nx+qx])
with open(ncmat, 'w') as fh:
    fh.write('NCMAT v7\n@DENSITY\n  2.2 g_per_cm3\n@DYNINFO\n  element Si\n'
             '  fraction 1\n  type freegas\n@CUSTOM_SASCSNS\n  DirectLoad2D\n')
    fh.write('  Qx ' + ' '.join(f'{v:.12g}' for v in g) + '\n')
    fh.write('  Qy ' + ' '.join(f'{v:.12g}' for v in g) + '\n  I ')
    for j, v in enumerate(flat):
        fh.write(f'{v:.8g}')
        fh.write('\n  ' if (j + 1) % 10 == 0 and j + 1 < flat.size else ' ')
    fh.write('\n')

import NCrystal  # noqa: E402
sc = NCrystal.createScatter(ncmat)
EKIN = E_MEV * 1e-3

#--- panel (b): sample outcomes, map to (Qx,Qy) --------------------------------
N_EVT = 2_000_000
_, dirs = sc.sampleScatter(EKIN, (0.0, 0.0, 1.0), repeat=N_EVT)
QXS, QYS = K * dirs[0], K * dirs[1]

#--- independent reference: midpoint quadrature over the plane disc ------------
def bilinear(qx, qy):
    """table values at (qx,qy); 0 outside the grid"""
    fx = (qx - g[0]) / PITCH
    fy = (qy - g[0]) / PITCH
    ok = (fx >= 0) & (fx <= NPTS - 1.0001) & (fy >= 0) & (fy <= NPTS - 1.0001)
    fx = np.clip(fx, 0.0, NPTS - 1.0001)
    fy = np.clip(fy, 0.0, NPTS - 1.0001)
    i0 = fx.astype(np.int64)
    j0 = fy.astype(np.int64)
    dx, dy = fx - i0, fy - j0
    v = (I_TAB[j0, i0] * (1 - dx) * (1 - dy) + I_TAB[j0, i0 + 1] * dx * (1 - dy)
         + I_TAB[j0 + 1, i0] * (1 - dx) * dy + I_TAB[j0 + 1, i0 + 1] * dx * dy)
    return np.where(ok, v, 0.0)

NQ = 501                                  # pitch 0.004 in n-space = table pitch
_nc = -1.0 + (np.arange(NQ) + 0.5) * 2.0 / NQ
NX, NY = np.meshgrid(_nc, _nc)
MASK = NX ** 2 + NY ** 2 < 1.0
NXm, NYm = NX[MASK], NY[MASK]
DA = (2.0 / NQ) ** 2
W = 2.0 / np.sqrt(1.0 - NXm ** 2 - NYm ** 2)   # two elastic branches

def sigma_quad(beta, phi):
    """sigma_scatter [barn/atom] for beam tilted by beta toward azimuth phi,
    from scratch:  sigma = 2*double-int I(Q)/|kf.z| dA  over the plane disc"""
    ax, ay = K * math.sin(beta) * math.cos(phi), K * math.sin(beta) * math.sin(phi)
    q = bilinear(K * NXm - ax, K * NYm - ay)
    return float(np.sum(q * W) * DA)

def beam(beta, phi):
    return (math.sin(beta) * math.cos(phi), math.sin(beta) * math.sin(phi),
            math.cos(beta))

#--- self-checks (gated; printed) ----------------------------------------------
s_on_pl = float(sc.crossSection(EKIN, (0.0, 0.0, 1.0)))
s_on_qd = sigma_quad(0.0, 0.0)
dev_on = abs(s_on_qd - s_on_pl) / s_on_pl
PHI_MAJ, PHI_MIN = TAU, TAU + math.pi / 2.0
BETAS = np.arange(0.0, 60.0001, 2.5) * 1e-3   # mrad -> rad
devs = []
for phi in (PHI_MAJ, PHI_MIN):
    for b in BETAS:
        sp = float(sc.crossSection(EKIN, beam(b, phi)))
        sq = sigma_quad(b, phi)
        devs.append(abs(sq - sp) / sp)
dev_max = max(devs)
failures = int(dev_on > 0.015) + int(dev_max > 0.015)

#--- figure --------------------------------------------------------------------
WIN = 0.6                     # display window: +-3 s_maj
dA = PITCH * PITCH
P_IN = I_TAB / I_TAB.sum() / dA
xe = np.append(g - PITCH / 2, g[-1] + PITCH / 2)
Hh, _, _ = np.histogram2d(QXS, QYS, bins=[xe, xe])
P_OUT = Hh.T / Hh.sum() / dA
inwin = (np.abs(g) <= WIN)
Pk = P_IN[np.ix_(inwin, inwin)].max()
norm = None                       # linear colour scale, SasView-style
lv = [Pk * f for f in (0.02, 0.1, 0.3, 0.6, 0.9)]

fig = plt.figure(figsize=(10.2, 8.6))
gs = GridSpec(2, 2, figure=fig, left=0.07, right=0.90, top=0.94, bottom=0.07,
              wspace=0.22, hspace=0.28, width_ratios=[1, 1])

# (a) input table
ax = fig.add_subplot(gs[0, 0])
im = ax.pcolormesh(g, g, P_IN, norm=norm, cmap='viridis',
                   shading='nearest', rasterized=True)
ax.contour(g, g, P_IN, levels=lv, colors='k', linewidths=0.6)
ax.annotate('', xy=(S_MAJ * 1.6 * math.cos(TAU), S_MAJ * 1.6 * math.sin(TAU)),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color='w', lw=1.2))

ax.set_xlim(-WIN, WIN)
ax.set_ylim(-WIN, WIN)
ax.set_title(r'(a) input $I(Q_x,Q_y)$: 4:1 Gaussian at $30^\circ$', fontsize=10)
ax.set_xlabel(r'$Q_x$ [$\mathrm{\AA}^{-1}$]')
ax.set_ylabel(r'$Q_y$ [$\mathrm{\AA}^{-1}$]')

# (b) sampled outcomes
ax = fig.add_subplot(gs[0, 1])
mout = P_OUT > 0
im2 = ax.pcolormesh(g, g, np.where(mout, P_OUT, np.nan), norm=norm,
                    cmap='viridis', shading='nearest', rasterized=True)
ax.contour(g, g, P_IN, levels=lv, colors='k', linewidths=0.6)
ax.set_xlim(-WIN, WIN)
ax.set_ylim(-WIN, WIN)
ax.set_title(r'(b) sampled outcome density, 2M events', fontsize=10)
ax.set_xlabel(r'$Q_x$ [$\mathrm{\AA}^{-1}$]')
ax.set_ylabel(r'$Q_y$ [$\mathrm{\AA}^{-1}$]')

# (c) azimuthal profile on a detector ring: d sigma / d psi at theta = 10 deg
ax = fig.add_subplot(gs[1, 0])
TH_RING = math.radians(10.0)
R0 = K * math.sin(TH_RING)
DR = 0.010
inring = np.abs(np.hypot(QXS, QYS) - R0) < DR / 2.0
psi = np.mod(np.arctan2(QYS[inring], QXS[inring]), 2 * math.pi)
NPB = 60
hh, edges = np.histogram(psi, bins=NPB, range=(0.0, 2 * math.pi))
ctr = 0.5 * (edges[1:] + edges[:-1])
dens = hh / (hh.sum() * (edges[1] - edges[0]))   # unit-integral density
ax.semilogy(np.degrees(ctr), dens, 'o', ms=3, color='0.35',
            label=r'sampled, 2M events')
ps_rad = np.linspace(0.0, 2.0 * math.pi, 721)
I_ring = bilinear(R0 * np.cos(ps_rad), R0 * np.sin(ps_rad))
C = 1.0 / (I_ring.sum() * (ps_rad[1] - ps_rad[0]))   # unit integral in psi
ax.semilogy(np.degrees(ps_rad), I_ring * C, '-', lw=1.6, color='C1',
            label=r'table prediction $I(Q(\psi))$')
modsens = I_ring.max() / I_ring.min()
ax.annotate(f'modulation depth {modsens:.0f}' + r'$\times$',
            xy=(0.97, 0.25), xycoords='axes fraction', ha='right', fontsize=8)
ax.set_xlim(0, 360)
ax.set_xticks([0, 90, 180, 270, 360])
ax.set_xlabel(r'$\psi$ [deg] (azimuth of $\vec Q_\perp$, $\psi=0$ along $Q_x$)')
ax.set_ylabel(r'$\mathrm{d}\sigma/\mathrm{d}\psi$  [arb.]')
ax.set_title(r'(c) anisotropy on the ring $\theta=10^\circ$'
             r' ($|\vec Q_\perp|=%.2f\,\mathrm{\AA}^{-1}$)' % R0, fontsize=10)
ax.legend(fontsize=7.5, loc='upper center', framealpha=0.9)

# (d) the integrated sigma is direction-independent: relative change + floor
ax = fig.add_subplot(gs[1, 1])
s0 = {phi: sigma_quad(0.0, phi) for phi in (PHI_MAJ, PHI_MIN)}
for phi, lab in ((PHI_MAJ, 'major axis ($30^\circ$)'),
                 (PHI_MIN, 'minor axis ($120^\circ$)')):
    sq = np.array([sigma_quad(b, phi) for b in BETAS])
    sp = np.array([float(sc.crossSection(EKIN, beam(b, phi))) for b in BETAS])
    ln, = ax.plot(BETAS * 1e3, (sq / s0[phi] - 1) * 1e3, '-', lw=1.6,
                  label=f'quadrature: {lab}')
    ax.plot(BETAS * 1e3, (sp / s0[phi] - 1) * 1e3, 'o', ms=3.5, mfc='none',
            color=ln.get_color(), label=f'plugin: {lab}')
floor = 15.0
ax.axhspan(-floor, floor, color='0.85', zorder=0)
ax.annotate('tabulated-model numerical floor (rung 6)'
            f' $\\approx\\pm{floor:.0f}\\times10^{{-3}}$',
            xy=(0.03, 0.78), xycoords='axes fraction', fontsize=8)
ax.set_xlabel('beam tilt magnitude $\\beta$ [mrad]')
ax.set_ylabel(r'$(\sigma(\beta)-\sigma(0))/\sigma(0)$  [$10^{-3}$]')
ax.set_ylim(-20, 25)
ax.set_title('(d) integrated $\\sigma$: direction-independent', fontsize=10)
ax.legend(fontsize=7.5, loc='lower right', framealpha=0.9)

cbar = fig.colorbar(im, ax=[fig.axes[0], fig.axes[1]], shrink=0.90,
                    pad=0.025, location='right')
cbar.set_label(r'probability density per $\mathrm{\AA}^2$', fontsize=8)
fig.savefig(OUT, bbox_inches='tight', dpi=220)
try:
    fig.savefig('/tmp/fig_anisotropy_preview.png', dpi=110, bbox_inches='tight')
except Exception:
    pass

# moment-based shape check along the pattern axes (printed, not drawn)
u_s = QXS * math.cos(TAU) + QYS * math.sin(TAU)
v_s = -QXS * math.sin(TAU) + QYS * math.cos(TAU)
ratio = np.mean(u_s ** 2) / np.mean(v_s ** 2)
print(f'on-axis sigma: plugin {s_on_pl:.6f}  quadrature {s_on_qd:.6f}  '
      f'dev {dev_on:.2e}')
print(f'sigma(tilt) curve: max |quadrature-plugin|/plugin = {dev_max:.2e} '
      f'[{"PASS" if dev_max <= 0.015 else "FAIL"}]')
print(f'sampled 2nd-moment ratio along pattern axes: {ratio:.2f} '
      f'(table: {(S_MAJ / S_MIN) ** 2:.1f})')
print(f'ring psi modulation depth (theta=10 deg): {modsens:.0f}x')
dm = (sigma_quad(60e-3, PHI_MAJ) / s0[PHI_MAJ] - 1) * 100
dn = (sigma_quad(60e-3, PHI_MIN) / s0[PHI_MIN] - 1) * 100
print(f'sigma(60 mrad) relative change: {dm:+.2f}% major, '
      f'{dn:+.2f}% minor -- below the {dev_max:.1e} numerical floor')
print(f'wrote {OUT}  [{failures} gated check(s) failed]')
