"""Shared numerics for the coupled-rotor paper (hbar = I = 1 unless stated)."""
import numpy as np
from scipy.integrate import quad


def sector_H(M, N, lam, hbar=1.0, I=1.0):
    """H in the sector m1+m2 = M, basis (m1, M-m1) with |m_i| <= N."""
    st = [(m1, M - m1) for m1 in range(-N, N + 1) if abs(M - m1) <= N]
    idx = {s: i for i, s in enumerate(st)}
    H = np.zeros((len(st), len(st)))
    for i, (a, b) in enumerate(st):
        H[i, i] = hbar ** 2 / (2 * I) * (a * a + b * b)
        for d in (+1, -1):
            j = idx.get((a + d, b - d))
            if j is not None:
                H[i, j] += lam / 2
    return H, st, idx


def exact_state(m1, m2, N, lam):
    H, st, idx = sector_H(m1 + m2, N, lam)
    w, v = np.linalg.eigh(H)
    i0 = idx[(m1, m2)]
    k = int(np.argmax(np.abs(v[i0, :]) ** 2))
    return v[:, k] * np.sign(v[i0, k]), w[k], st, idx


def rho1_diag(psi, st, N):
    p = np.zeros(2 * N + 1)
    for i, (a, b) in enumerate(st):
        p[a + N] += psi[i] ** 2
    return p


def shannon(p):
    p = np.asarray(p, float)
    p = p[p > 1e-300]
    return float(-np.sum(p * np.log(p)))


def pendulum_parity(lam, N=None):
    """M = 0 sector = quantum pendulum H = m^2 + (lam/2)(shift).  Returns the
    energies and m-distributions of the exact parity eigenstates, sorted."""
    if N is None:
        N = int(8 * np.sqrt(lam)) + 150
    out = []
    for par in (+1, -1):
        d = N + 1 if par > 0 else N
        mm = np.arange(0, N + 1) if par > 0 else np.arange(1, N + 1)
        H = np.diag(mm.astype(float) ** 2)
        off = np.full(d - 1, lam / 2)
        if par > 0:
            off[0] = lam / np.sqrt(2)     # <0|H|(1+-1)/sqrt2>
        H += np.diag(off, 1) + np.diag(off, -1)
        w, v = np.linalg.eigh(H)
        for k in range(len(w)):
            c = v[:, k]
            p = np.zeros(2 * N + 1)
            if par > 0:
                p[N] = c[0] ** 2
                p[N + mm[1:]] = c[1:] ** 2 / 2
                p[N - mm[1:]] = c[1:] ** 2 / 2
            else:
                p[N + mm] = c ** 2 / 2
                p[N - mm] = c ** 2 / 2
            out.append((w[k], p, par))
    out.sort(key=lambda t: t[0])
    return (np.array([o[0] for o in out]), [o[1] for o in out],
            np.array([o[2] for o in out]), N)


def classical_cells(E, lam, N):
    """Microcanonical reduced momentum density f(p) ∝ 1/sqrt(lam^2-(E-p^2)^2),
    integrated over unit cells centred on the integers -N..N."""
    pmin, pmax = np.sqrt(max(E - lam, 0.0)), np.sqrt(E + lam)
    f = lambda p: 1.0 / np.sqrt(max(lam ** 2 - (E - p * p) ** 2, 1e-300))
    w = np.zeros(2 * N + 1)
    for m in range(0, N + 1):
        a, b = max(m - 0.5, pmin), min(m + 0.5, pmax)
        if b <= a:
            continue
        pts = [x for x in (pmin, pmax) if a < x < b]
        val, _ = quad(f, a, b, points=pts or None, limit=400)
        if m == 0:
            w[N] = val
        else:
            w[N + m] = w[N - m] = val / 2
    return w / w.sum()
