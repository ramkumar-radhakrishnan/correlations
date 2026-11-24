"""Numerical checks for the DKS-style coupled-rotor paper (hbar = I = 1).

1. Exact diagonalisation of the sector Hamiltonian h_kk' = k^2 d_kk' + (lam/2)(d_k',k+1 + d_k',k-1)
   for a highly excited branch state, against the WKB/Bessel weights J_j(gamma)^2.
2. Branch entropy S_E(gamma) = -sum J_j^2 ln J_j^2: small-gamma (perturbative) and
   large-gamma (classical) limits; coefficient of the gamma^(-1/3) approach.
3. Classical reduced distribution from a time average over a classical trajectory,
   against the closed form and the Boltzmann entropy formula.
4. Perturbative reduced density matrix (notes, Eq. 33 corrected) against exact numerics.
"""
import numpy as np
from scipy.special import jv
from scipy.integrate import solve_ivp, quad

def shannon(p):
    p = np.asarray(p, float); p = p[p > 1e-300]
    return float(-np.sum(p * np.log(p)))

def S_bessel(g):
    jm = int(g + 40 * g ** (1 / 3) + 50)
    return shannon(jv(np.arange(-jm, jm + 1), g) ** 2)

def branch_state(n, lam, half=None):
    """Exact eigenvector of the sector Hamiltonian nearest to the Wannier-Stark level k0^2,
    in the window k in [k0-half, k0+half] (k0 = n/2; the mirror branch is far away)."""
    k0 = n // 2
    if half is None:
        half = min(k0 - 1, int(lam / n + 40 * (lam / n) ** (1 / 3) + 60))
    ks = np.arange(k0 - half, k0 + half + 1)
    H = np.diag(ks.astype(float) ** 2) + np.diag(np.full(len(ks) - 1, lam / 2), 1) \
        + np.diag(np.full(len(ks) - 1, lam / 2), -1)
    w, v = np.linalg.eigh(H)
    i = int(np.argmax(np.abs(v[half, :]) ** 2 * (np.abs(w - k0 ** 2) < lam)))
    # pick the level whose weight distribution is centred on k0
    cent = (v ** 2 * ks[:, None]).sum(0)
    i = int(np.argmin(np.abs(cent - k0) + 1e3 * (np.abs(w - k0 ** 2) > 2 * lam)))
    return ks - k0, v[:, i] ** 2, w[i]

def main():
    print('=' * 78)
    print('1. Exact diagonalisation vs Bessel weights (branch state, n = 2 k0)')
    for n, g in ((800, 2.0), (800, 5.0), (800, 20.0), (1600, 20.0)):
        lam = g * n
        j, p, e = branch_state(n, lam)
        pb = jv(j, g) ** 2
        print('  n=%5d gamma=%5.1f lam/E_rel=%.3f  max|p-J^2|=%.2e  S_exact=%.5f  S_Bessel=%.5f'
              % (n, g, lam / (n / 2) ** 2, np.abs(p - pb).max(), shannon(p), S_bessel(g)))

    print('=' * 78)
    print('2. Branch entropy limits')
    for g in (0.02, 0.05, 0.1, 0.2):
        F = g * g / 4
        Spt = -(1 - 2 * F) * np.log(1 - 2 * F) - 2 * F * np.log(F)
        print('  gamma=%.2f  S_Bessel=%.6e  S_pert(F=g^2/4)=%.6e  rel.err=%.2e'
              % (g, S_bessel(g), Spt, abs(S_bessel(g) - Spt) / S_bessel(g)))
    gs = np.array([1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6])
    res = []
    for g in gs:
        S = S_bessel(g)
        res.append(S - np.log(np.pi * g / np.e))
        print('  gamma=%8.0f  S_E=%.6f  ln(pi g/2)=%.6f  ln(pi g/e)=%.6f  S_E-ln(pi g/e)=%.5f  x g^(1/3)=%.4f'
              % (g, S, np.log(np.pi * g / 2), np.log(np.pi * g / np.e), res[-1], res[-1] * g ** (1 / 3)))
    res = np.array(res)
    A = np.vstack([gs ** (-1 / 3), gs ** (-2 / 3)]).T
    coef = np.linalg.lstsq(A, res, rcond=None)[0]
    print('  fit S_E - ln(pi g/e) = c1 g^(-1/3) + c2 g^(-2/3):  c1=%.4f  c2=%.4f' % tuple(coef))

    print('=' * 78)
    print('3. Classical trajectory: time average of p1 vs closed form (I = 1)')
    lam, J, E_rel = 50.0, 0.0, 2500.0          # p0 = 50, a = lam I / (2 p0) = 0.5 ... make a larger:
    lam, E_rel = 400.0, 2500.0                  # p0 = 50, a = 4, lam/E_rel = 0.16
    p0 = np.sqrt(E_rel)
    def rhs(t, y):
        th1, th2, p1, p2 = y
        s = np.sin(th1 - th2)
        return [p1, p2, lam * s, -lam * s]
    # initial condition on the torus with p>0: phi = pi/2 (cos phi = 0), p = sqrt(E_rel)
    y0 = [np.pi / 2, 0.0, J / 2 + p0, J / 2 - p0]
    T = 2000.0
    sol = solve_ivp(rhs, (0, T), y0, max_step=0.002, rtol=1e-10, atol=1e-10, dense_output=False,
                    t_eval=np.linspace(0, T, 2_000_001))
    p1 = sol.y[2]
    Echeck = (sol.y[2] ** 2 + sol.y[3] ** 2) / 2 + lam * np.cos(sol.y[0] - sol.y[1])
    print('  energy drift %.2e' % (np.abs(Echeck - Echeck[0]).max()))
    # exact density w(p1) = 1/(2B sqrt(lam^2 - (E_rel - p^2)^2)) on the p>0 band, normalised to 1 there
    pl, ph = np.sqrt(E_rel - lam), np.sqrt(E_rel + lam)
    f = lambda p: 1 / np.sqrt(max(lam ** 2 - (E_rel - p * p) ** 2, 1e-300))
    Bnum = quad(f, pl, ph, limit=200)[0]
    edges = np.arange(np.floor(pl) - 0.5, np.ceil(ph) + 1.5, 1.0)     # cells of width hbar centred on integers
    hist, _ = np.histogram(p1 - J / 2, bins=edges)
    wt = hist / hist.sum()
    wex = np.array([quad(f, max(a, pl), min(b, ph), points=[x for x in (pl, ph) if a < x < b], limit=200)[0]
                    if min(b, ph) > max(a, pl) else 0.0 for a, b in zip(edges[:-1], edges[1:])]) / Bnum
    print('  max |trajectory cell weights - closed form| = %.2e' % np.abs(wt - wex).max())
    # Boltzmann entropy with Delta = h: exact (differential) vs linearised arcsine ln(pi^2 a / Delta)
    Sdiff = quad(lambda p: -f(p) / Bnum * np.log(f(p) / Bnum), pl, ph, points=[], limit=400)[0]
    a = lam / (2 * p0)
    print('  a = %.3f;  S_CE exact (Delta=h) = %.5f ; linearised ln(pi a/2) = %.5f ; binned = %.5f'
          % (a, Sdiff, np.log(np.pi * a / 2), shannon(wex)))

    print('=' * 78)
    print('4. Perturbative reduced density matrix (corrected notes) vs exact')
    for (m1, m2, lam) in ((4, 4, 0.02), (9, 2, 0.02), (9, 2, 0.2)):
        n, M = m1 - m2, m1 + m2
        al = lam / 2
        ep, em = al ** 2 / (n + 1) ** 2, al ** 2 / (n - 1) ** 2
        N = 60
        st = [(a, M - a) for a in range(-N, N + 1)]
        idx = {s: i for i, s in enumerate(st)}
        H = np.zeros((len(st), len(st)))
        for i, (a, b) in enumerate(st):
            H[i, i] = (a * a + b * b) / 2
            for d in (1, -1):
                jj = idx.get((a + d, b - d))
                if jj is not None:
                    H[i, jj] = lam / 2
        w, v = np.linalg.eigh(H)
        k = int(np.argmax(v[idx[(m1, m2)], :] ** 2))
        p = {s[0]: v[i, k] ** 2 for i, s in enumerate(st)}
        S_ex = shannon(list(p.values()))
        S_pt = ep * np.log(1 / ep) + em * np.log(1 / em) + ep + em
        S_bad = ep * np.log(1 / ep) + em * np.log(1 / em)
        print('  (m1,m2)=(%d,%d) lam=%.2f: p0 ex=%.8f pt=%.8f | p+ ex=%.3e pt=%.3e | S ex=%.6e pt=%.6e no-norm=%.6e'
              % (m1, m2, lam, p[m1], 1 - ep - em, p[m1 + 1], ep, S_ex, S_pt, S_bad))


if __name__ == '__main__':
    main()
