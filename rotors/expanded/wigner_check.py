"""Does a Wigner (or Husimi) distribution work in place of the classical one?

Relative motion in sector M=0: u(phi) = sum_k c_k e^{ik phi}/sqrt(2 pi), the
rotor-1 Schmidt weights are p_k = |c_k|^2.  Wigner function on the cylinder
(phi in [0,2pi), l in Z/2):
    W(phi, l) = (1/2pi) sum_{k+k'=2l} c_k c_k' e^{i(k-k') phi}.
"""
import numpy as np
from rotorlib import pendulum_parity, classical_cells, shannon

def wigner_rel(c, N, nphi=512):
    ks = np.arange(-N, N + 1)
    phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    ells = np.arange(-2 * N, 2 * N + 1) / 2.0
    W = np.zeros((len(ells), nphi))
    for a, k in enumerate(ks):
        if abs(c[a]) < 1e-14:
            continue
        for b, kp in enumerate(ks):
            if abs(c[b]) < 1e-14:
                continue
            li = int(round(k + kp)) + 2 * N        # index of l = (k+k')/2
            W[li] += c[a] * c[b] * np.cos((k - kp) * phi) / (2 * np.pi)
    return phi, ells, W

def smoothed_entropy(p, ks, sigma, grid):
    """Differential entropy of sum_k p_k g_sigma(x - k)  (Husimi momentum marginal)."""
    h = (p[:, None] * np.exp(-(grid[None, :] - ks[:, None]) ** 2 / (2 * sigma ** 2))).sum(0)
    h /= np.sqrt(2 * np.pi) * sigma
    dx = grid[1] - grid[0]
    h = h / (h.sum() * dx)
    hp = h[h > 1e-300]
    return -np.sum(hp * np.log(hp)) * dx

lam = 200.
E, D, P, N = pendulum_parity(lam)
ks = np.arange(-N, N + 1)
for target in (0.5, 1.5):
    k0 = int(np.argmin(np.abs(E / lam - target)))
    p = D[k0]
    sgn = P[k0]
    # rebuild real amplitudes c_k (sign structure irrelevant for p; needed for W)
    # -> re-diagonalise in the parity sector to get the signed vector
    from rotorlib import sector_H
    H, st, idx = sector_H(0, N, lam)
    w, v = np.linalg.eigh(H)
    # pick the eigenvector of matching energy and parity
    cands = [j for j in range(len(w)) if abs(w[j] - E[k0]) < 1e-6]
    c = None
    for j in cands:
        vec = v[:, j]
        sym = vec + sgn * vec[::-1]
        if np.linalg.norm(sym) > 1e-6:
            c = sym / np.linalg.norm(sym); break
    phi, ells, W = wigner_rel(c, N)
    dphi = phi[1] - phi[0]
    marg = W.sum(1) * dphi
    integer = np.abs(ells - np.round(ells)) < 1e-9
    neg = np.sum(np.abs(W)) * dphi - 1
    Sq = shannon(p)
    w_cl = classical_cells(E[k0], lam, N)
    Scl = shannon(w_cl)
    print('=' * 72)
    print('lam=%.0f  E/lam=%.2f   (%s)' % (lam, E[k0] / lam, 'libration' if E[k0] < lam else 'rotation'))
    print('  phi-marginal of W on integer l vs p_k: max diff = %.2e ; on half-integer l: max |.| = %.2e'
          % (np.abs(marg[integer] - p).max(), np.abs(marg[~integer]).max()))
    print('  Wigner negativity volume  sum|W| dphi - 1 = %.3f   (min W = %.3e)' % (neg, W.min()))
    print('  S_q (von Neumann)                        = %.4f' % Sq)
    print('  entropy of reduced Wigner W_1 = p_k/2pi  = %.4f   (identical by construction)' % Sq)
    print('  S_cl (classical, binned)                 = %.4f' % Scl)
    grid = np.linspace(-N - 6, N + 6, 40 * (2 * N + 13))
    print('  Husimi-smoothed momentum marginals (width sigma in units of hbar):')
    print('     %6s %12s %12s %10s' % ('sigma', 'S[p*g]', 'S[w_cl*g]', 'gap'))
    for s in (0.3, 0.5, 0.7071, 1.0, 1.5, 2.5):
        a, b = smoothed_entropy(p, ks, s, grid), smoothed_entropy(w_cl, ks, s, grid)
        print('     %6.3f %12.4f %12.4f %10.4f' % (s, a, b, b - a))
