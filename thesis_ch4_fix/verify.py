"""Numerical verification of the corrected Chapter-4 wave functions (hbar = m = 1)."""
import numpy as np
from numpy.polynomial.hermite import hermval
from math import factorial, sqrt, pi
from scipy.linalg import eigh

def hfun(m, y):
    c = np.zeros(m + 1); c[m] = 1
    return hermval(y, c) * np.exp(-y * y / 2) / sqrt(2 ** m * factorial(m) * sqrt(pi))

def ymat(D):
    a = np.diag(np.sqrt(np.arange(1, D)), 1)
    return (a + a.T) / np.sqrt(2)

def first_order(n, Hp):
    c = {}
    for m in range(Hp.shape[0] - 8):
        if m != n and abs(Hp[m, n]) > 1e-14:
            c[m] = Hp[m, n] / (n - m)
    return c

def state(n, c, y, mom=False):
    """normalized first-order state on grid y (position: y = sqrt(w) x; momentum: q = k/sqrt(w))."""
    N = 1 / np.sqrt(1 + sum(v * v for v in c.values()))
    if not mom:
        return N * (hfun(n, y) + sum(v * hfun(m, y) for m, v in c.items()))
    return N * ((-1j) ** n * hfun(n, y) + sum(v * (-1j) ** m * hfun(m, y) for m, v in c.items()))

D = 220
Y = ymat(D); Y4 = Y @ Y @ Y @ Y; Y6 = Y4 @ Y @ Y
y = np.linspace(-14, 14, 4001); dy = y[1] - y[0]
k = np.linspace(-8, 8, 801)

def ft(psi, x, kk):
    return (psi[None, :] * np.exp(-1j * np.outer(kk, x))).sum(1) * (x[1] - x[0]) / np.sqrt(2 * np.pi)

print('=' * 72)
print('Structure checks (quartic, b = 0.1; symmetric well, lambda = 1.5)')
for label, Hp, w in (('quartic', 0.1 * Y4, 1.0),
                     ('sym. well', (np.sqrt(10) * 1.5 / 14) * Y4 + (3 * 1.5 ** 2 / 28) * Y6, 1.5 * np.sqrt(2 / 45))):
    for n in range(4):
        c = first_order(n, Hp)
        x = y / np.sqrt(w)
        psi = state(n, c, y) * w ** 0.25            # psi(x), normalized in x
        norm = (psi ** 2).sum() * (x[1] - x[0])
        par = np.abs(psi[::-1] - (-1) ** n * psi).max()
        phi_num = ft(psi, x, k)
        phi_an = state(n, c, k / np.sqrt(w), mom=True) * w ** (-0.25)
        rot = phi_num / (-1j) ** n                   # should be real, parity (-1)^n
        print('  %-9s n=%d: norm=%.12f  |psi(-x)-(-1)^n psi|=%.1e  |FT-analytic|=%.1e  '
              'max|Im[phi/(-i)^n]|=%.1e  |phi(-k)-(-1)^n phi(k)|=%.1e'
              % (label, n, norm, par, np.abs(phi_num - phi_an).max(), np.abs(rot.imag).max(),
                 np.abs(phi_num[::-1] - (-1) ** n * phi_num).max()))

print('=' * 72)
print('Quartic: first-order energy E_n = n + 1/2 + (3b/4)(2n^2+2n+1) and fidelity |<exact|PT>|^2')
for bb in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 1.3, 1.5):
    H = np.diag(np.arange(D) + 0.5) + bb * Y4
    w_, v_ = eigh(H[:200, :200])
    row = []
    for n in range(4):
        c = first_order(n, bb * Y4)
        vec = np.zeros(200); vec[n] = 1
        for m, val in c.items():
            vec[m] = val
        vec /= np.linalg.norm(vec)
        row.append((abs(v_[:, n] @ vec) ** 2, w_[n], n + 0.5 + 0.75 * bb * (2 * n * n + 2 * n + 1)))
    print('  b=%4.2f  ' % bb + '  '.join('n=%d: F=%.4f E=%.4f E1st=%.4f' % (n, *r) for n, r in enumerate(row)))

print('=' * 72)
print('Symmetric well, truncated model (x^6) and full potential: fidelity of first-order state')
def V_full(x, lam):
    z = lam * x
    out = np.empty_like(z)
    small = np.abs(z) < 1e-3
    out[small] = 1 / 3 + z[small] ** 2 / 45
    zz = z[~small]
    out[~small] = (1 - zz / np.tan(zz)) / zz ** 2
    return out
for lam in (0.5, 1.0, 1.5, 2.0):
    w = lam * np.sqrt(2 / 45)
    Hp = (np.sqrt(10) * lam / 14) * Y4 + (3 * lam ** 2 / 28) * Y6
    H = (np.diag(np.arange(D) + 0.5) + Hp)[:200, :200]
    w_, v_ = eigh(H)
    # full potential on a grid inside the walls |x| < pi/lam (finite differences)
    L = np.pi / lam
    xg = np.linspace(-L, L, 3001)[1:-1]; h = xg[1] - xg[0]
    Tm = (np.diag(np.full(len(xg), 1.0)) - 0.5 * np.diag(np.ones(len(xg) - 1), 1)
          - 0.5 * np.diag(np.ones(len(xg) - 1), -1)) / h ** 2
    Hg = Tm + np.diag(V_full(xg, lam))
    wf, vf = eigh(Hg)
    row = []
    for n in range(4):
        c = first_order(n, Hp)
        vec = np.zeros(200); vec[n] = 1
        for m, val in c.items():
            vec[m] = val
        vec /= np.linalg.norm(vec)
        Ft = abs(v_[:, n] @ vec) ** 2
        psi = state(n, c, np.sqrt(w) * xg) * w ** 0.25
        ex = vf[:, n] / np.sqrt(h)
        Ff = (abs((psi * ex).sum() * h) ** 2) / ((psi ** 2).sum() * h)
        row.append((Ft, Ff))
    print('  lambda=%.1f  ' % lam + '  '.join('n=%d: F_trunc=%.3f F_full=%.3f' % (n, *r) for n, r in enumerate(row)))
