"""Wigner and Husimi analysis for the semiclassical coupled-rotor paper.

Units hbar = I = 1.  Sector M = 0: relative Hamiltonian h = k^2 + lam cos(phi),
eigenvector c_k, rotor-1 Schmidt weights p_k = c_k^2.

Produces
  * relative Wigner function W(phi, l), l in Z/2, for one libration and one
    rotation level at lam = 200 (data for fig5_wigner.pdf), with the checks
    marginal == p_k, half-integer marginal == 0, normalisation, negativity;
  * Husimi (Wehrl) entropies S_W[p](sigma) of the quantum weights and
    S_W[w](sigma) of the continuous classical density, for every level at
    lam = 200 and 800, with the fringe-model prediction D(sigma; E);
  * the deep-rotation (Bessel) table: S_W - S_q against (1 - ln 2) - D(sigma).
Results go to wigner_husimi.npz and are printed.
"""
import numpy as np
from scipy.special import jv, ellipk
from scipy.integrate import quad
from rotorlib import pendulum_parity, classical_cells, shannon

C0 = 1.0 - np.log(2.0)


# ---------------------------------------------------------------- fringe model
def f_contrast(A):
    """<(1 + A cos x) ln(1 + A cos x)>_x  for 0 <= A <= 1."""
    s = np.sqrt(np.clip(1.0 - A * A, 0.0, None))
    return np.log((1.0 + s) / 2.0) + 1.0 - s


def D_rot(sig):
    """Deep-rotation smoothing loss D(sigma) = (2/pi) int_0^{pi/2} f(e^{-2 s^2 b^2}) db."""
    return (2 / np.pi) * quad(lambda b: f_contrast(np.exp(-2 * sig ** 2 * b ** 2)),
                              0, np.pi / 2, limit=200)[0]


def c1_series(nmax=4000):
    """sigma * D(sigma) for sigma >~ 1: (1/2sqrt(pi)) sum_n C(2n,n)/(4^n 2n(2n-1) sqrt n)."""
    tot, t = 0.0, 1.0
    for n in range(1, nmax):
        t *= (2 * n) * (2 * n - 1) / (n * n * 4.0)        # C(2n,n)/4^n recursively
        tot += t / (2 * n * (2 * n - 1) * np.sqrt(n))
    return tot / (2 * np.sqrt(np.pi))


# ---------------------------------------------------------------- smoothing
def smoothed_h(q, sigma, K=None):
    """Differential entropy of sum_j q_j N(x; j, sigma^2), q on consecutive integers."""
    if K is None:
        K = max(8, int(np.ceil(10 / sigma)))
    R = int(np.ceil(8 * sigma)) + 1
    off = np.arange(-R, R + 1)
    H = 0.0
    for r in range(K):
        x = off + r / K
        g = np.exp(-x ** 2 / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma)
        h = np.convolve(q, g[::-1], mode='full')
        h = h[h > 1e-300]
        H += -np.sum(h * np.log(h))
    return H / K


def orbit_samples(E, lam, ns=200001):
    """Momenta p, relative angles phi in [0,pi], and weights of a uniform-in-time
    sampling of the classical level set (both tori above the separatrix)."""
    if E > lam:                                    # rotation, phi = pi - 2s
        k2 = 2 * lam / (E + lam)
        s = (np.arange(ns) + 0.5) * np.pi / ns
        wt = 1 / np.sqrt(1 - k2 * np.sin(s) ** 2)
        p = np.sqrt(E + lam) * np.sqrt(1 - k2 * np.sin(s) ** 2)
        p = np.concatenate([p, -p]); wt = np.concatenate([wt, wt])
    else:                                          # libration, sin(x/2) = kappa sin t
        kap2 = (E + lam) / (2 * lam)
        t = (np.arange(ns) + 0.5) * 2 * np.pi / ns
        wt = 1 / np.sqrt(1 - kap2 * np.sin(t) ** 2)
        p = np.sqrt(E + lam) * np.cos(t)
    wt = wt / wt.sum()
    phi = np.arccos(np.clip((E - p * p) / lam, -1, 1))
    return p, phi, wt


def band_B(E, lam):
    if E > lam:
        return np.sqrt(1 / (E + lam)) * ellipk(2 * lam / (E + lam))
    return np.sqrt(1 / (2 * lam)) * ellipk((E + lam) / (2 * lam))


def classical_diff_entropy(E, lam):
    """S_cl = -int w ln w dp evaluated as an orbit average of -ln w."""
    p, phi, wt = orbit_samples(E, lam)
    B = band_B(E, lam)
    w = 1 / (2 * B * np.sqrt(np.clip(lam ** 2 - (E - p * p) ** 2, 1e-300, None)))
    return float(-np.sum(wt * np.log(w)))


def classical_smoothed_h(E, lam, sigma, dx=0.01):
    """Wehrl entropy of the classical density: h[w * g_sigma]."""
    p, phi, wt = orbit_samples(E, lam)
    L = np.abs(p).max() + 10 * sigma + 1
    edges = np.arange(-L, L + dx, dx)
    hist, _ = np.histogram(p, bins=edges, weights=wt)
    xk = np.arange(-int(8 * sigma / dx) - 1, int(8 * sigma / dx) + 2) * dx
    g = np.exp(-xk ** 2 / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma)
    Q = np.convolve(hist, g, mode='same')
    Q = Q[Q > 1e-300]
    return float(-np.sum(Q * np.log(Q)) * dx)


def D_general(E, lam, sigma):
    """Orbit average of f(exp(-sigma^2 d(2 phi)^2 / 2)), d = distance to 2 pi Z."""
    p, phi, wt = orbit_samples(E, lam, ns=40001)
    d = np.minimum(2 * phi, 2 * np.pi - 2 * phi)
    return float(np.sum(wt * f_contrast(np.exp(-sigma ** 2 * d ** 2 / 2))))


# ---------------------------------------------------------------- Wigner
def signed_vector(lam, N, E_target, parity):
    """Real eigenvector c_k (k = -N..N) of the M = 0 sector with given parity."""
    d = N + 1 if parity > 0 else N
    mm = np.arange(0, N + 1) if parity > 0 else np.arange(1, N + 1)
    H = np.diag(mm.astype(float) ** 2)
    off = np.full(d - 1, lam / 2)
    if parity > 0:
        off[0] = lam / np.sqrt(2)
    H += np.diag(off, 1) + np.diag(off, -1)
    w, v = np.linalg.eigh(H)
    j = int(np.argmin(np.abs(w - E_target)))
    a = v[:, j]
    c = np.zeros(2 * N + 1)
    if parity > 0:
        c[N] = a[0]
        c[N + mm[1:]] = a[1:] / np.sqrt(2)
        c[N - mm[1:]] = a[1:] / np.sqrt(2)
    else:
        c[N + mm] = a / np.sqrt(2)
        c[N - mm] = -a / np.sqrt(2)
    return c, w[j]


def wigner_rel(c, N, nphi=720, kmax=None):
    """W(phi, l) = (1/2pi) sum_{k+k'=2l} c_k c_k' cos((k-k')phi), l in Z/2."""
    ks = np.arange(-N, N + 1)
    nz = np.nonzero(np.abs(c) > 1e-12)[0]
    phi = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    ells = np.arange(-2 * N, 2 * N + 1) / 2.0
    W = np.zeros((len(ells), nphi))
    for a in nz:
        for b in nz:
            li = ks[a] + ks[b] + 2 * N
            W[li] += c[a] * c[b] * np.cos((ks[a] - ks[b]) * phi) / (2 * np.pi)
    return phi, ells, W


def main():
    out = {}
    # ------------------------------------------------------------ Wigner maps
    lam = 200.0
    E, Dl, P, N = pendulum_parity(lam)
    print('=' * 78)
    print('Relative Wigner function, lam = 200')
    for tag, target in (('lib', 0.5), ('rot', 1.5)):
        k0 = int(np.argmin(np.abs(E / lam - target)))
        c, e = signed_vector(lam, N, E[k0], P[k0])
        assert abs(e - E[k0]) < 1e-8
        phi, ells, W = wigner_rel(c, N)
        dphi = phi[1] - phi[0]
        marg = W.sum(1) * dphi
        integer = np.abs(ells - np.round(ells)) < 1e-9
        p = c ** 2
        norm = W.sum() * dphi
        neg = np.abs(W).sum() * dphi - 1
        # position marginal vs |u(phi)|^2
        u = (c[None, :] * np.exp(1j * np.outer(phi, np.arange(-N, N + 1)))).sum(1) / np.sqrt(2 * np.pi)
        posm = W.sum(0)
        print('  %s: E/lam = %.4f  parity %+d' % (tag, e / lam, P[k0]))
        print('     normalisation sum W dphi           = %.12f' % norm)
        print('     max|marginal - p_k| (integer l)     = %.2e' % np.abs(marg[integer] - p).max())
        print('     max|marginal| (half-integer l)      = %.2e' % np.abs(marg[~integer]).max())
        print('     max|sum_l W - |u(phi)|^2|           = %.2e' % np.abs(posm - np.abs(u) ** 2).max())
        print('     negativity volume sum|W| dphi - 1   = %.3f' % neg)
        print('     entropy of W_1 = p_k/2pi (cell h)   = %.6f  ;  S_vN = %.6f'
              % (shannon(p), shannon(p)))
        out['W_' + tag] = W; out['phi_' + tag] = phi; out['ells_' + tag] = ells
        out['E_' + tag] = e; out['neg_' + tag] = neg

    # -------------------------------------------- deep rotation, Bessel weights
    print('=' * 78)
    print('Deep rotation (Bessel weights): Wehrl entropy against the fringe model')
    c1 = c1_series()
    print('  c1 = lim sigma D(sigma) = %.6f   (leading term 1/(8 sqrt pi) = %.6f)'
          % (c1, 1 / (8 * np.sqrt(np.pi))))
    sigs = (0.5, 1.0, 2.0, 4.0)
    rows = []
    print('%8s %5s %9s %9s %9s %9s %9s %9s' % ('gamma', 'sig', 'S_q', 'S_cl', 'S_W[p]',
          'S_W[w]', 'S_W-S_q', 'pred'))
    for gam in (30., 100., 1000., 10000.):
        jm = int(gam + 40 * gam ** (1 / 3) + 50)
        j = np.arange(-jm, jm + 1)
        p = jv(j, gam) ** 2
        Sq, Scl = shannon(p), np.log(np.pi * gam / 2)
        for s in sigs:
            SWp = smoothed_h(p, s)
            # classical arcsine smoothed: orbit average with p = gam cos(beta)
            beta = (np.arange(400001) + 0.5) * np.pi / 400001
            x = gam * np.cos(beta)
            dx = 0.01
            L = gam + 10 * s + 1
            edges = np.arange(-L, L + dx, dx)
            hist, _ = np.histogram(x, bins=edges)
            hist = hist / hist.sum()
            xk = np.arange(-int(8 * s / dx) - 1, int(8 * s / dx) + 2) * dx
            g = np.exp(-xk ** 2 / (2 * s ** 2)) / (np.sqrt(2 * np.pi) * s)
            Q = np.convolve(hist, g, mode='same'); Q = Q[Q > 1e-300]
            SWw = float(-np.sum(Q * np.log(Q)) * dx)
            pred = C0 - D_rot(s)
            rows.append((gam, s, Sq, Scl, SWp, SWw, SWp - Sq, pred))
            print('%8.0f %5.2f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f'
                  % (gam, s, Sq, Scl, SWp, SWw, SWp - Sq, pred))
    out['bessel_rows'] = np.array(rows)
    sgrid = np.linspace(0.5, 8, 60)
    out['sgrid'] = sgrid
    out['Dgrid'] = np.array([D_rot(s) for s in sgrid])
    # Bessel S_W - S_q on the sigma grid for the figure
    for gam in (100., 1000., 10000.):
        jm = int(gam + 40 * gam ** (1 / 3) + 50)
        p = jv(np.arange(-jm, jm + 1), gam) ** 2
        Sq = shannon(p)
        out['bes_gap_%d' % gam] = np.array([smoothed_h(p, s) - Sq for s in sgrid[::4]])
    out['sgrid_pts'] = sgrid[::4]

    # ----------------------------------------- generic energies, lam = 200, 800
    print('=' * 78)
    print('Generic energies: Wehrl entropies (window averages)')
    for lam in (200.0, 800.0):
        E, Dl, P, N = pendulum_parity(lam)
        sel = [i for i in range(len(E)) if -0.9 < E[i] / lam < 3.0]
        res = []
        for i in sel:
            p = Dl[i]
            Sq = shannon(p)
            Scl = classical_diff_entropy(E[i], lam)
            Sbin = shannon(classical_cells(E[i], lam, N))
            row = [E[i] / lam, Sq, Scl, Sbin]
            for s in (0.5, 1.0, 2.0):
                row += [smoothed_h(p, s), classical_smoothed_h(E[i], lam, s),
                        D_general(E[i], lam, s)]
            res.append(row)
        res = np.array(res)
        out['gen_%d' % lam] = res
        for lo, hi, name in ((0.2, 0.7, 'libration'), (1.5, 2.5, 'rotation')):
            m = (res[:, 0] >= lo) & (res[:, 0] <= hi)
            r = res[m].mean(0)
            print('  lam=%4.0f %-9s (%3d levels): S_cl-S_q=%.3f  S_bin-S_q=%.3f'
                  % (lam, name, m.sum(), r[2] - r[1], r[3] - r[1]))
            for a, s in enumerate((0.5, 1.0, 2.0)):
                SWp, SWw, Dp = r[4 + 3 * a], r[5 + 3 * a], r[6 + 3 * a]
                print('        sigma=%.1f  S_W[p]-S_q=%.3f  (1-ln2)-D=%.3f  '
                      'S_W[w]-S_W[p]=%.3f  S_W[w]-S_cl=%.3f'
                      % (s, SWp - r[1], C0 - Dp, SWw - SWp, SWw - r[2]))
    np.savez('wigner_husimi.npz', **out)


if __name__ == '__main__':
    main()
