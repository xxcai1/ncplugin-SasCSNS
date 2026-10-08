#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# validate_directload2d.py - validation harness for the DirectLoad2D anisotropic
#                            SANS model (pure python; no plugin required).
#
# VALIDATION DOCTRINE (X.-X. Cai): find SPECIAL CASES whose parameter sets
# simplify the physics into something computable in closed form, then run the
# model through its completely ordinary numerical path (read the table, do the
# numerics) and compare. The model never special-cases anything; the reference
# always comes from outside.
#
# Two pillars:
#   Pillar A (cross section, ABSOLUTE): sigma against closed-form special
#           cases and against independent quadrature of continuum references.
#   Pillar B (sampler, SHAPE ONLY): sampled outcome distributions compared to
#           normalised reference densities (chi2); absolute scale is
#           meaningless for a distribution and is not tested here.
#
# Ladder of independent evaluations of the same physical case:
#   1. closed form       (special parameter sets; hand-derived in the doc)
#   2. continuum quad.   (analytic cylinder/Gaussian; no table involved)
#   3. table quad.       (bilinear table; exact python twin of the C++ model)
#   4. [later] plugin    (NCrystal C++ process; same assertions re-pointed)
# 1 vs 2 validates the math, 2 vs 3 quantifies table discretisation,
# 3 vs 4 validates the implementation.
#
# Model definition (doc/anisotropic_directload2d.tex):
#   - I(Qx,Qy) [barn/(atom sr)] on a uniform ascending grid in the material
#     frame; bilinear; I = 0 outside; Qx = Q.x-hat, Qy = Q.y-hat with
#     Q = kf - ki (wave-vector transfer, k = sqrt(E/2.07214meV) = 2pi/lambda).
#   - elastic; outcomes live on the FULL elastic sphere. Each plane point
#     (Qx,Qy) inside the disc has TWO elastic preimages (kf.z = +-sqrt)/k
#     with equal solid angle dA/(k^2*|kf.z|): the model extends the qz=0
#     dataset to the backward branch by mirror symmetry (exact for
#     z-symmetric structures and in both closed-form limits:
#     constant table -> 4*pi*I0, E->0 -> 4*pi*I(0)).
#   - sigma(ki,E) = double integral over the plane disc
#     D = {(Qx+ka_x)^2+(Qy+ka_y)^2 <= k^2} of I/(k^2*|kf.z|) dA. With
#     rho = k*sin(alpha) this is exactly
#         sigma = int_0^{2pi} int_0^{pi/2} I(disc point) sin(alpha) da dphi,
#     smooth (no singularity), the reference quadrature below.
#   - sampler (doc section 7.3, v3): per-energy CDF over the OUTCOME ANGLES
#     (theta,psi) of the design beam z-hat:
#         f(theta,psi) = I_table(k sin(th)cos(ps), k sin(th)sin(ps)) * sin(th),
#     so elasticity and the hemisphere are built into the coordinates
#     (kf = (sin th cos ps, sin th sin ps, cos th)) and no rejection is needed.
#     Tilted beams (neutrons already deflected before the SANS collision):
#     propose from the DILATED table CDF (max-filter of I by the maximal
#     in-plane beam shift k*sin(beta_max)) and accept with the pointwise mask
#         xi <= I(Q(kf-ki)) / I_dilated(Q(kf-zhat))   (<= 1 by construction),
#     exact for any tilt <= beta_max; far-off beams fall back to uniform-
#     hemisphere rejection (slow but correct).
# -----------------------------------------------------------------------------
import argparse
import math
import sys

import numpy as np

E_OVER_K2 = 2.07214  # E[meV] = E_OVER_K2 * k^2 [1/Aa^2]


def k_of_e(e_mev):
    return math.sqrt(e_mev / E_OVER_K2)


# =============================================================================
# The model, as pure-python reference (ladder item 3)
# =============================================================================
class Table2D:
    """Uniform-grid I(Qx,Qy) in barn/(atom*sr), bilinear, 0 outside."""

    def __init__(self, qx, qy, vals, name='table'):
        self.qx = np.asarray(qx, float)
        self.qy = np.asarray(qy, float)
        self.vals = np.asarray(vals, float).reshape(self.qy.size, self.qx.size)
        for nm, arr in (('Qx', self.qx), ('Qy', self.qy)):
            d = np.diff(arr)
            if arr.size < 2 or np.any(d <= 0):
                raise ValueError(f'{nm} must be strictly ascending, >=2 points')
            if not np.allclose(d, d[0], rtol=1e-6):
                raise ValueError(f'{nm} grid must be uniformly spaced')
        self.name = name
        self.dx = float(self.qx[1] - self.qx[0])
        self.dy = float(self.qy[1] - self.qy[0])

    def eval(self, x, y):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        inside = ((x >= self.qx[0]) & (x <= self.qx[-1])
                  & (y >= self.qy[0]) & (y <= self.qy[-1]))
        fx = np.clip((x - self.qx[0]) / self.dx, 0, self.qx.size - 1 - 1e-12)
        fy = np.clip((y - self.qy[0]) / self.dy, 0, self.qy.size - 1 - 1e-12)
        i0 = fx.astype(int).clip(0, self.qx.size - 2)
        j0 = fy.astype(int).clip(0, self.qy.size - 2)
        tx = fx - i0
        ty = fy - j0
        v = (self.vals[j0, i0] * (1 - tx) * (1 - ty)
             + self.vals[j0, i0 + 1] * tx * (1 - ty)
             + self.vals[j0 + 1, i0] * (1 - tx) * ty
             + self.vals[j0 + 1, i0 + 1] * tx * ty)
        return np.where(inside, v, 0.0)

    def q0(self):
        return float(self.vals[self.qy.size // 2, self.qx.size // 2])

    def dilated(self, radius):
        """Max-filter of I over all grid cells within `radius` (for the tilted
        beam proposal). Grid-based; exact within the table resolution."""
        if radius <= 0:
            return self
        nx = int(radius / self.dx) + 1
        ny = int(radius / self.dy) + 1
        out = self.vals.copy()
        for j in range(-ny, ny + 1):
            for i in range(-nx, nx + 1):
                if i * i * self.dx * self.dx + j * j * self.dy * self.dy > radius**2:
                    continue
                sh = np.full_like(out, 0.0)
                js = slice(max(0, j), self.qy.size + min(0, j))
                jd = slice(max(0, -j), self.qy.size + min(0, -j))
                iss = slice(max(0, i), self.qx.size + min(0, i))
                idd = slice(max(0, -i), self.qx.size + min(0, -i))
                sh[jd, idd] = self.vals[js, iss]
                out = np.maximum(out, sh)
        t = Table2D(self.qx, self.qy, out, self.name + ' dilated')
        return t


def sigma_table(tab, ki, k, nalpha=400, nphi=720):
    """sigma by Gauss-Legendre quadrature of the smooth alpha-form (reference).
    Independent of the samplers and of the sigma used by them."""
    ax, ay, _ = ki
    cx, cy = -k * ax, -k * ay
    xg, wg = np.polynomial.legendre.leggauss(nalpha)
    alpha = 0.5 * math.pi * (xg + 1.0)           # [0, pi]: full sphere
    wa = 0.5 * math.pi * wg
    phi = (np.arange(nphi) + 0.5) * (2 * math.pi / nphi)
    A, P = np.meshgrid(alpha, phi, indexing='ij')
    Qx = cx + k * np.sin(A) * np.cos(P)
    Qy = cy + k * np.sin(A) * np.sin(P)
    I = tab.eval(Qx.ravel(), Qy.ravel()).reshape(A.shape)
    return float(np.sum(I * np.sin(A) * wa[:, None]) * (2 * math.pi / nphi))


def sigma_cont_generic(func, ki, k, nalpha=600, nphi=1080):
    """Continuum sigma of an analytic I(Qx,Qy) via the same alpha-quadrature."""
    ax, ay, _ = ki
    cx, cy = -k * ax, -k * ay
    xg, wg = np.polynomial.legendre.leggauss(nalpha)
    alpha = 0.5 * math.pi * (xg + 1.0)           # [0, pi]: full sphere
    wa = 0.5 * math.pi * wg
    phi = (np.arange(nphi) + 0.5) * (2 * math.pi / nphi)
    A, P = np.meshgrid(alpha, phi, indexing='ij')
    Qx = cx + k * np.sin(A) * np.cos(P)
    Qy = cy + k * np.sin(A) * np.sin(P)
    I = func(Qx.ravel(), Qy.ravel()).reshape(A.shape)
    return float(np.sum(I * np.sin(A) * wa[:, None]) * (2 * math.pi / nphi))


def plane_target_weight(tab, ki, k, Qx, Qy):
    """Unnormalised exact target density over plane points: 2*I/|kf.z|
    (two elastic branches preimage each plane point)."""
    ax, ay, _ = ki
    rho2 = (Qx + k * ax) ** 2 + (Qy + k * ay) ** 2
    with np.errstate(invalid='ignore', divide='ignore'):
        kfz = np.sqrt(np.clip(k * k - rho2, 0.0, None)) / k
        w = np.where(kfz > 0, 2 * tab.eval(Qx, Qy) / np.maximum(kfz, 1e-300), 0.0)
    return w


# =============================================================================
# Samplers (Pillar B)
# =============================================================================
class AngularCDF:
    """Per-energy sampling table over the outcome angles (theta,psi) of the
    design beam z-hat. masses = I_eff(k*sin(th)cos(ps), k*sin(th)sin(ps))*sin(th)
    on a uniform grid; sigma_grid = total mass (same discretisation => exact
    xs/sampler consistency). dilate_r > 0 builds the tilted-beam proposal."""

    def __init__(self, tab, k, nth=192, npsi=384, dilate_r=0.0):
        src = tab.dilated(dilate_r) if dilate_r > 0 else tab
        th = (np.arange(nth) + 0.5) * (math.pi / nth)
        ps = (np.arange(npsi) + 0.5) * (2 * math.pi / npsi)
        T, P = np.meshgrid(th, ps, indexing='ij')
        Qx = k * np.sin(T) * np.cos(P)
        Qy = k * np.sin(T) * np.sin(P)
        f = src.eval(Qx.ravel(), Qy.ravel()).reshape(T.shape) * np.sin(T)
        self.dth = math.pi / nth
        self.dps = 2 * math.pi / npsi
        self.k = k
        mass = f * self.dth * self.dps
        self.cum = np.cumsum(mass.ravel())
        self.total = float(self.cum[-1])
        if self.total <= 0:
            raise ValueError(f'empty angular CDF at k={k}')
        self.cum /= self.total
        self.sigma_grid = self.total        # sigma on this discretisation
        self.nth, self.npsi = nth, npsi

    def sample_kf(self, rng, n):
        u = rng.random(n)
        idx = np.searchsorted(self.cum, u)
        np.clip(idx, 0, self.cum.size - 1, out=idx)
        jth, jps = np.unravel_index(idx, (self.nth, self.npsi))
        th = (jth + rng.random(n)) * self.dth
        ps = (jps + rng.random(n)) * self.dps
        st = np.sin(th)
        return np.stack([st * np.cos(ps), st * np.sin(ps), np.cos(th)], axis=1)


def sample_baseline(tab, ki, k, rng, n, oversample=16):
    """Obviously-correct rejection sampler: uniform on the full sphere,
    accept with prob I(plane projection)/I_max. Independent of all model code
    (acceptance can be small for peaked tables; tests use it at moderate
    acceptance or with generous draws)."""
    imax = float(tab.vals.max())
    out = np.empty((n, 3))
    got = 0
    while got < n:
        m = max(int((n - got) * oversample * 1.5), 4096)
        z = 2.0 * rng.random(m) - 1.0
        phi = 2 * math.pi * rng.random(m)
        r = np.sqrt(1 - z * z)
        kf = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)
        Qx = k * (kf[:, 0] - ki[0])
        Qy = k * (kf[:, 1] - ki[1])
        acc = tab.eval(Qx, Qy) / imax
        sel = rng.random(m) < acc
        take = kf[sel][: n - got]
        out[got:got + take.shape[0]] = take
        got += take.shape[0]
        if got < n and m > 200 * (n - got) * max(1.0 / max(float(tab.vals.max()), 1e-300), 1):
            pass  # keep looping; hard failure handled by caller timeout
    return out


def sample_tilted(tab, tab_dil, cdf_dil, ki, k, rng, n):
    """Tilted-beam fast path: propose kf from the dilated per-E CDF, accept with
    xi <= I(Q(kf-ki)) / I_dil(Q(kf-zhat))  (<= 1 by the dilation argument:
    the target support is the table dilated by the maximal in-plane beam
    shift, which is exactly where the proposal lives).
    Returns (outcomes, acceptance)."""
    out = np.empty((n, 3))
    got = 0
    draws = 0
    ax, ay, _ = ki
    while got < n:
        m = max(int((n - got) * 1.4), 2048)
        kf = cdf_dil.sample_kf(rng, m)
        Qx = k * (kf[:, 0] - ax)
        Qy = k * (kf[:, 1] - ay)
        I_t = tab.eval(Qx, Qy)
        Qx0 = k * kf[:, 0]
        Qy0 = k * kf[:, 1]
        I_p = tab_dil.eval(Qx0, Qy0)
        acc = np.where(I_p > 0, I_t / np.maximum(I_p, 1e-300), 0.0)
        sel = rng.random(m) < acc
        take = kf[sel][: n - got]
        out[got:got + take.shape[0]] = take
        got += take.shape[0]
        draws += m
        if draws > 2000 * n + 100000:
            raise RuntimeError('tilted acceptance collapsed')
    return out, got / draws


# =============================================================================
# Shape statistics (Pillar B)
# =============================================================================
def chi2_angular(tab, ki, k, kf, nbins=30, nquad=24):
    """Binned chi2 of outcome DIRECTIONS against the exact normalised target
    over the upper hemisphere, in the outcome's own angles (theta from z-hat,
    phi). In these coordinates the target density
        w(theta,phi) = I_table(plane components of k*(kf_hat-ki_hat)) * sin(theta)
    is smooth everywhere (no 1/|kf.z| horizon singularity), so plain
    Gauss-Legendre per bin converges fast. Works for tilted beams too."""
    th = np.arctan2(np.hypot(kf[:, 0], kf[:, 1]), kf[:, 2])
    ps = np.arctan2(kf[:, 1], kf[:, 0])
    thb = np.linspace(0.0, math.pi, nbins + 1)
    psb = np.linspace(-math.pi, math.pi, nbins + 1)
    xg, xw = np.polynomial.legendre.leggauss(nquad)
    yg, yw = np.polynomial.legendre.leggauss(nquad)
    h = np.zeros((nbins, nbins))
    for i in range(nbins):
        gt = thb[i] + (xg + 1) / 2 * (thb[i + 1] - thb[i])
        wt = xw * (thb[i + 1] - thb[i]) / 2
        for j in range(nbins):
            gp = psb[j] + (yg + 1) / 2 * (psb[j + 1] - psb[j])
            wp = yw * (psb[j + 1] - psb[j]) / 2
            T, P = np.meshgrid(gt, gp, indexing='ij')
            kfx = np.sin(T) * np.cos(P)
            kfy = np.sin(T) * np.sin(P)
            Qx = k * (kfx - ki[0])
            Qy = k * (kfy - ki[1])
            w = tab.eval(Qx.ravel(), Qy.ravel()).reshape(T.shape) * np.sin(T)
            h[i, j] = (w * wt[:, None] * wp[None, :]).sum()
    tot = h.sum()
    if tot <= 0:
        return float('nan'), 0
    h *= kf.shape[0] / tot
    obs, _, _ = np.histogram2d(th, ps, bins=[thb, psb])
    m = h >= 5.0
    if m.sum() < 10:
        return float('nan'), int(m.sum())
    from scipy.stats import chi2 as chi2dist
    chi2 = float((((obs[m] - h[m]) ** 2) / h[m]).sum())
    return chi2, chi2dist.sf(chi2, int(m.sum()) - 1)


def chi2_1d(counts, expected, n):
    m = expected >= 5.0
    if m.sum() < 5:
        return float('nan'), float('nan')
    from scipy.stats import chi2 as chi2dist
    c = float((((counts[m] - expected[m]) ** 2) / expected[m]).sum())
    return c, chi2dist.sf(c, int(m.sum()) - 1)


# =============================================================================
# Special cases (closed form) - the validation doctrine in action
# =============================================================================
def make_grid(half, n):
    return np.linspace(-half, half, n)


def case_constant():
    """S1: I = I0 everywhere. sigma = 2*pi*I0 exactly (upper hemisphere!),
    outcomes uniform on the hemisphere (z uniform in [0,1]). The special-case
    result 2*pi (not 4*pi) is what locks the one-sided/hemisphere convention."""
    rows = []
    k, I0 = 1.0, 3.0
    tab = Table2D(make_grid(2.2, 401), make_grid(2.2, 401),
                  np.full((401, 401), I0), 'const')
    ki = (0.0, 0.0, 1.0)
    sig = sigma_table(tab, ki, k)
    rows.append(('S1 constant table: sigma = 4*pi*I0', 4 * math.pi * I0, sig, 'rel', 1e-6))

    cdf = AngularCDF(tab, k)
    rows.append(('S1 sampling-grid sigma consistency (angular CDF vs quad)',
                 sig, cdf.sigma_grid, 'rel', 5e-3))

    rng = np.random.default_rng(11)
    kf = cdf.sample_kf(rng, 200_000)
    rows.append(('S1 angular-CDF outcomes: max ||kf|-1| (elasticity)', 0.0,
                 float(np.abs(np.linalg.norm(kf, axis=1) - 1.0).max()), 'near', 1e-12))
    z = np.clip(kf[:, 2], -1.0, 1.0)
    cnt, _ = np.histogram(z, bins=20, range=(-1, 1))
    c2, p = chi2_1d(cnt, np.full(20, z.size / 20), z.size)
    rows.append(('S1 angular-CDF outcomes uniform on sphere (p)', 0.01, p, 'p', None))

    kf2 = sample_baseline(tab, ki, k, np.random.default_rng(12), 100_000)
    z2 = np.clip(kf2[:, 2], -1.0, 1.0)
    cnt2, _ = np.histogram(z2, bins=20, range=(-1, 1))
    c2b, pb = chi2_1d(cnt2, np.full(20, z2.size / 20), z2.size)
    rows.append(('S1 baseline outcomes uniform on sphere (p)', 0.01, pb, 'p', None))
    return rows


def case_zero_energy():
    """S2: k -> 0. The disc shrinks to the origin, |kf.z| -> 1, so
    sigma -> 4*pi*I(0) for ANY table and ANY beam direction."""
    rows = []
    n = 401
    g = make_grid(0.8, n)
    X, Y = np.meshgrid(g, g)
    vals = 1.0 + np.exp(-((X - 0.25) ** 2 + (Y + 0.3) ** 2) / (2 * 0.12 ** 2))
    tab = Table2D(g, g, vals, 'blob')
    e = 1e-6
    k = k_of_e(e)
    ki = np.array([0.6, -0.5, 1.1]); ki /= np.linalg.norm(ki)
    sig = sigma_table(tab, ki, k)
    rows.append(('S2 E->0 limit (anisotropic table, tilted beam): sigma = 4*pi*I(0)',
                 4 * math.pi * tab.q0(), sig, 'rel', 1e-4))
    return rows


def case_ezero_expansion():
    """S8: the E->0 Taylor expansion as a SHARP anchor. The disc is centred
    at -k*a_perp (not at the origin), so expanding I about the origin and
    angle-integrating (eq.~alphaint) gives
      sigma = 4*pi*I0 - 4*pi*k*a_perp.grad I0
              + k^2*[2*pi*a_perp^T H0 a_perp + (2*pi/3)*lap I0] + O(k^3).
    For a LINEAR ramp the Taylor series terminates at first order, so the
    formula is EXACT for any k (up to quadrature): the strongest closed-form
    check of the first-order tilt coefficient in the model. For a wide
    Gaussian the k^2 term is measurable and pins the curvature coefficient.
    """
    rows = []
    # ---- L1: linear ramp, exact at any k -------------------------------
    b = 0.5
    n = 401
    g = make_grid(0.8, n)
    X, Y = np.meshgrid(g, g)
    vals = 1.0 + b * X          # I >= 0.6 on the table, gradient (b, 0)
    tab = Table2D(g, g, vals, 'ramp')
    for k, beta, tag in [(0.3, 0.0, 'on-axis'),
                         (0.3, 0.03, 'beta=30mrad'),
                         (0.3, -0.03, 'beta=-30mrad (sign)'),
                         (0.3, 0.2, 'beta=200mrad')]:
        a_x = math.sin(beta)
        exact = 4 * math.pi * (1.0 - b * k * a_x)
        sig = sigma_table(tab, (a_x, 0.0, math.cos(beta)), k)
        rows.append((f'S8 ramp: sigma exact linear formula ({tag})',
                     exact, sig, 'rel', 1e-6))
    # ---- L2: wide Gaussian, k^2 curvature coefficient ------------------
    s_q = 1.0
    g = make_grid(2.5, 401)
    X, Y = np.meshgrid(g, g)
    vals = np.exp(-(X ** 2 + Y ** 2) / (2 * s_q ** 2))
    tab = Table2D(g, g, vals, 'widegauss')
    k = k_of_e(0.1)             # the model's lowest energy node
    sig = sigma_table(tab, (0, 0, 1), k)
    approx = 4 * math.pi * (1.0 - k ** 2 / (3 * s_q ** 2))
    rows.append(('S8 wide Gaussian E->0: sigma vs 4pi(1 - k^2/(3 s^2)) '
                 '(O(k^3) residual)', approx, sig, 'rel', 1e-3))
    return rows


def case_single_pixel():
    """S4: one bright pixel. The disc (radius k) is large compared to the
    pixel, so sigma -> I0*dA/(k^2*|kf.z0|) and outcomes cluster at
    kf0 = ki + (Qx0,Qy0,kfz0)/k."""
    rows = []
    k, I0, h = 1.0, 5.0, 0.02
    n = 201
    g = make_grid(2.02, n)
    vals = np.zeros((n, n))
    ix = np.searchsorted(g, 0.3) - 1
    iy = np.searchsorted(g, -0.15) - 1
    vals[iy, ix] = I0
    tab = Table2D(g, g, vals, 'pixel')
    ki = (0.0, 0.0, 1.0)
    kfz0 = math.sqrt(1.0 - (0.3 ** 2 + 0.15 ** 2))
    sig = sigma_table(tab, ki, k, nalpha=1200, nphi=2400)
    rows.append(('S4 single pixel: sigma vs 2*I0*dA/(k^2*|kf.z|) (2 branches)',
                 2 * I0 * h * h / (k * k * kfz0), sig, 'rel', 0.05))
    cdf = AngularCDF(tab, k, nth=384, npsi=768)
    kf = cdf.sample_kf(np.random.default_rng(21), 20_000)
    kf0 = np.array([0.3, -0.15, kfz0])
    # two elastic branches preimage the pixel: kf0 and its z-mirror
    d = np.minimum(np.linalg.norm(kf - kf0, axis=1),
                   np.linalg.norm(kf - kf0 * np.array([1, 1, -1]), axis=1))
    frac = float((d < 3 * h / k).mean())
    rows.append(('S4 pixel: outcomes cluster at the 2 predicted kf (frac<3h/k)',
                 0.98, frac, 'p', None))
    return rows


def case_gaussian():
    """S5: Gaussian table. No closed form for sigma, but three INDEPENDENT
    evaluations must agree (continuum quadrature, table quadrature, sampling
    grid), the fast-path sampler must reproduce the exact plane density
    (shape only), and tilted beams must work through the dilated proposal."""
    rows = []
    s, k = 0.15, 1.0
    n = 401
    g = make_grid(2.0, n)
    X, Y = np.meshgrid(g, g)
    vals = np.exp(-(X ** 2 + Y ** 2) / (2 * s ** 2))
    tab = Table2D(g, g, vals, 'gauss')
    ki = (0.0, 0.0, 1.0)
    func = lambda x, y: np.exp(-(x ** 2 + y ** 2) / (2 * s ** 2))

    sig_c = sigma_cont_generic(func, ki, k)
    sig_t = sigma_table(tab, ki, k)
    rows.append(('S5 Gaussian: continuum sigma vs table sigma (401^2 grid)',
                 sig_c, sig_t, 'rel', 1e-4))
    g2 = make_grid(2.0, 801)
    X2, Y2 = np.meshgrid(g2, g2)
    tab2 = Table2D(g2, g2, np.exp(-(X2 ** 2 + Y2 ** 2) / (2 * s ** 2)), 'gauss fine')
    rows.append(('S5 Gaussian: table sigma 401^2 vs 801^2 (discretisation)',
                 sigma_table(tab2, ki, k), sig_t, 'rel', 1e-4))

    cdf = AngularCDF(tab, k)
    rows.append(('S5 Gaussian: sampling-grid sigma vs table quadrature',
                 sig_t, cdf.sigma_grid, 'rel', 5e-3))

    kf = cdf.sample_kf(np.random.default_rng(31), 200_000)
    c2, p = chi2_angular(tab, ki, k, kf, nbins=24)
    rows.append(('S5 Gaussian: angular-CDF sampler shape chi2 p (200k)', 1e-3, p, 'p', None))
    rows.append(('S5 Gaussian: sampler |kf|-1 max', 0.0,
                 float(np.abs(np.linalg.norm(kf, axis=1) - 1).max()), 'near', 1e-12))

    # tilted beam, tier 2 (dilated proposal, pointwise mask)
    beta = 0.05
    ki_t = (math.sin(beta), 0.0, math.cos(beta))
    tab_dil = tab.dilated(k * math.sin(beta) * 1.001)
    cdf_dil = AngularCDF(tab_dil, k)
    kt, acc = sample_tilted(tab, tab_dil, cdf_dil, ki_t, k,
                            np.random.default_rng(32), 150_000)
    c2, p = chi2_angular(tab, ki_t, k, kt, nbins=24)
    rows.append(('S5 tilted beam (beta=50mrad): shape chi2 p (150k)', 1e-3, p, 'p', None))
    rows.append(('S5 tilted beam: acceptance (report)', acc, acc, 'near', 0.0))
    # sigma of the tilted beam via ROTATION INVARIANCE: rotating the table so
    # the tilted beam becomes z-hat must reproduce sigma(ki_t) (independent
    # evaluation path; residual = table-rotation interpolation error).
    from scipy.ndimage import rotate as ndrotate
    sig_tilt = sigma_table(tab, ki_t, k)
    best = None
    for sgn in (+1, -1):
        vals_rot = ndrotate(vals, sgn * math.degrees(beta), axes=(1, 0),
                            reshape=False, order=1, mode='constant', cval=0.0)
        tab_rot = Table2D(g, g, vals_rot, 'gauss rotated')
        s = sigma_table(tab_rot, (0.0, 0.0, 1.0), k)
        d = abs(s / sig_tilt - 1.0)
        if best is None or d < best:
            best = d
    rows.append(('S5 tilted beam: sigma vs rotated-table on-axis sigma',
                 0.0, best, 'near', 5e-3))

    # tier 1 (small tilt): exact on-axis CDF + ratio mask with scanned envelope
    beta1 = 0.002
    ki_t1 = (math.sin(beta1), 0.0, math.cos(beta1))
    cdf_ex = AngularCDF(tab, k)
    M = 0.0
    th = (np.arange(256) + 0.5) * (0.5 * math.pi / 256)
    ps = (np.arange(512) + 0.5) * (2 * math.pi / 512)
    T, P = np.meshgrid(th, ps, indexing='ij')
    kfx, kfy = np.sin(T) * np.cos(P), np.sin(T) * np.sin(P)
    I_p = tab.eval(k * kfx, k * kfy)
    I_t = tab.eval(k * (kfx - ki_t1[0]), k * (kfy - ki_t1[1]))
    M = 1.01 * float(np.max(I_t[I_p > 0] / I_p[I_p > 0]))
    kf1 = cdf_ex.sample_kf(np.random.default_rng(33), 120_000)
    acc1 = tab.eval(k * (kf1[:, 0] - ki_t1[0]), k * (kf1[:, 1] - ki_t1[1])) / \
        np.maximum(tab.eval(k * kf1[:, 0], k * kf1[:, 1]), 1e-300) / M
    keep = np.random.default_rng(34).random(kf1.shape[0]) < acc1
    kf1 = kf1[keep]
    c2, p = chi2_angular(tab, ki_t1, k, kf1, nbins=24)
    rows.append(('S5 tier-1 tilt (2mrad, exact CDF + envelope): shape chi2 p',
                 1e-3, p, 'p', None))
    return rows


def case_cylinder():
    """S6: radially symmetric continuum reference (infinite cylinder parallel
    to the beam: its form factor at qz=0 is the Airy disc 2J1(x)/x). No table
    needed for the reference; validates continuum vs table vs BOTH samplers."""
    rows = []
    R, k = 8.0, 0.6
    func = lambda x, y: (np.sinc(np.sqrt(x ** 2 + y ** 2) * R / math.pi)) ** 2
    ki = (0.0, 0.0, 1.0)
    n = 501
    g = make_grid(1.3, n)
    X, Y = np.meshgrid(g, g)
    tab = Table2D(g, g, func(X, Y), 'airy')
    sig_c = sigma_cont_generic(func, ki, k, nalpha=800, nphi=1440)
    sig_t = sigma_table(tab, ki, k)
    rows.append(('S6 Airy (cylinder||beam): continuum sigma vs table sigma (501^2)',
                 sig_c, sig_t, 'rel', 1e-4))

    kf = sample_baseline(tab, ki, k, np.random.default_rng(41), 150_000)
    c2, p = chi2_angular(tab, ki, k, kf, nbins=24)
    rows.append(('S6 Airy: baseline sampler shape chi2 p (150k)', 1e-3, p, 'p', None))
    psi = np.arctan2(k * kf[:, 1], k * kf[:, 0])
    cnt, _ = np.histogram(psi, bins=24, range=(-math.pi, math.pi))
    c2, p = chi2_1d(cnt, np.full(24, psi.size / 24), psi.size)
    rows.append(('S6 Airy: psi uniformity p', 1e-3, p, 'p', None))

    cdf = AngularCDF(tab, k)
    kf2 = cdf.sample_kf(np.random.default_rng(42), 150_000)
    c2, p = chi2_angular(tab, ki, k, kf2, nbins=24)
    rows.append(('S6 Airy: angular-CDF sampler shape chi2 p (150k)', 1e-3, p, 'p', None))
    return rows


def case_conventions():
    """S7: orientation conventions. A blob at POSITIVE Qx must send outcomes
    to positive Qx (scattering vector, not momentum transfer), and the plane
    distribution must be independent of the sign of the beam z-component
    (the model depends on the beam only through the disc centre -ka_x,-ka_y)."""
    rows = []
    k = 1.0
    n = 241
    g = make_grid(1.2, n)
    X, Y = np.meshgrid(g, g)
    vals = np.exp(-((X - 0.4) ** 2 + Y ** 2) / (2 * 0.08 ** 2))
    tab = Table2D(g, g, vals, 'xblob')
    cdf = AngularCDF(tab, k)
    kf = cdf.sample_kf(np.random.default_rng(51), 100_000)
    mean_qx = float((k * kf[:, 0]).mean())
    rows.append(('S7 x-blob: modal outcome Qx in [0.34,0.46]', 1.0,
                 float(0.34 < np.median(k * kf[:, 0]) < 0.46), 'near', 0.0))
    ki_dn = (0.0, 0.0, -1.0)
    sig_up = sigma_table(tab, (0, 0, 1), k)
    sig_dn = sigma_table(tab, ki_dn, k)
    rows.append(('S7 beam +z vs -z: same sigma (plane coords beam-sign independent)',
                 sig_up, sig_dn, 'rel', 1e-12))
    return rows


CASES = [('S1 constant table', case_constant),
         ('S2 E->0 limit', case_zero_energy),
         ('S4 single pixel', case_single_pixel),
         ('S5 Gaussian', case_gaussian),
         ('S6 Airy/cylinder', case_cylinder),
         ('S7 conventions', case_conventions),
         ('S8 E->0 expansion', case_ezero_expansion)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--latex', action='store_true', help='emit LaTeX result rows')
    ap.add_argument('--only', default=None, help='run a single case (substring)')
    args = ap.parse_args()

    allrows = []
    np.seterr(all='ignore')
    for name, fn in CASES:
        if args.only and args.only not in name:
            continue
        print(f'--- {name} ...', flush=True)
        try:
            allrows += fn()
        except Exception as exc:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            allrows.append((f'{name}: CRASHED', 0.0, -1.0, 'near', 1.0))

    print('\n' + '=' * 100)
    print(f"{'case':<70}{'expected':>14}{'got':>14}{'rel.err/p':>11}  verdict")
    print('-' * 100)
    nfail = 0
    latex_rows = []
    for name, exp, got, kind, tol in allrows:
        if kind == 'rel':
            err = abs(got / exp - 1.0) if exp != 0 else float('inf')
            ok = err < tol
            metric = f'{err:.2e}'
        elif kind == 'p':
            ok = got >= exp and got == got
            metric = f'{got:.4f}' if got == got else 'nan'
        else:  # near
            err = abs(got - exp)
            ok = err <= tol
            metric = f'{err:.2e}'
        nfail += 0 if ok else 1
        print(f'{name:<70}{exp:>14.6g}{got:>14.6g}{metric:>11}  '
              + ('PASS' if ok else 'FAIL'))
        latex_rows.append((name, exp, got, metric, ok))
    print('=' * 100)
    print('ALL PASS' if nfail == 0 else f'{nfail} FAILURE(S)')

    if args.latex:
        print('\n% ---- LaTeX rows ----')
        for name, exp, got, metric, ok in latex_rows:
            print(f"{name} & {exp:.4g} & {got:.4g} & {metric} & "
                  f"{'pass' if ok else 'FAIL'} \\\\")

    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
