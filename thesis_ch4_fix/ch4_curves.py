"""Corrected Figs. 4.7-4.10: exact S_x, S_k, S_x+S_k, F_x, F_k, F_x F_k for n = 0..3
(quartic vs b, symmetric well vs lambda), plus first-order PT (quartic, dashed).  hbar = m = 1."""
import json
import numpy as np
from numpy.polynomial.hermite import hermval
from math import factorial, sqrt, pi
from scipy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 9, 'legend.fontsize': 7, 'axes.linewidth': 0.6})
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#9b59b6']


def hfun(m, y):
    c = np.zeros(m + 1); c[m] = 1
    return hermval(y, c) * np.exp(-y * y / 2) / sqrt(2 ** m * factorial(m) * sqrt(pi))


def aho_coeffs(n, g):
    c = {n: 1.0}
    c[n + 4] = -g * np.sqrt((n + 1) * (n + 2) * (n + 3) * (n + 4)) / 16
    c[n + 2] = -g * (2 * n + 3) * np.sqrt((n + 1) * (n + 2)) / 4
    if n >= 2:
        c[n - 2] = g * (2 * n - 1) * np.sqrt(n * (n - 1)) / 4
    if n >= 4:
        c[n - 4] = g * np.sqrt(n * (n - 1) * (n - 2) * (n - 3)) / 16
    N = 1 / np.sqrt(sum(v * v for v in c.values()))
    return {m: N * v for m, v in c.items()}


def shannon_1d(rho, d):
    r = rho[rho > 1e-300]
    return float(-np.sum(r * np.log(r)) * d)


def sinc_dvr(V, L, N):
    x = np.linspace(-L, L, N); d = x[1] - x[0]
    i = np.arange(N); dif = i[:, None] - i[None, :]
    with np.errstate(divide='ignore'):
        T = np.where(dif == 0, np.pi ** 2 / 3, 2.0 / np.where(dif == 0, 1, dif) ** 2) * (-1.0) ** dif / (2 * d * d)
    e, v = eigh(T + np.diag(V(x)))
    return x, d, e, v


def box_dvr(V, a, b, N):
    i = np.arange(1, N); x = a + (b - a) * i / N
    I, J = np.meshgrid(i, i, indexing='ij')
    with np.errstate(divide='ignore', invalid='ignore'):
        off = 1 / np.sin(np.pi * (I - J) / (2 * N)) ** 2 - 1 / np.sin(np.pi * (I + J) / (2 * N)) ** 2
        diag = (2 * N ** 2 + 1) / 3 - 1 / np.sin(np.pi * i / N) ** 2
    T = np.where(I == J, 0.0, off); T[np.arange(N - 1), np.arange(N - 1)] = diag
    T = T * (-1.0) ** (I - J) * (np.pi ** 2 / (b - a) ** 2) / 4
    e, v = eigh(T + np.diag(V(x)))
    return x, (b - a) / N, e, v


def measures_psi(x, d, psi, k, E=None, Vx=None):
    """[S_x, S_k, F_x, F_k]; F_x = 4<k^2> = 8(E - <V>) when E, Vx are given (exact on the DVR grid)."""
    dk = k[1] - k[0]
    phi = (psi[None, :] * np.exp(-1j * np.outer(k, x))).sum(1) * d / np.sqrt(2 * np.pi)
    rx, rk = psi ** 2, np.abs(phi) ** 2
    sx2 = float(np.sum(x ** 2 * rx) * d); sk2 = float(np.sum(k ** 2 * rk) * dk)
    if E is not None:
        sk2 = float(2 * (E - np.sum(Vx * rx) * d))
    return np.array([shannon_1d(rx, d), shannon_1d(rk, dk), 4 * sk2, 4 * sx2])


def V_well(lam, V0=1.0):
    def V(x):
        z = lam * x; out = np.empty_like(z); s = np.abs(z) < 1e-4
        out[s] = 1 / 3 + z[s] ** 2 / 45
        out[~s] = (1 - z[~s] / np.tan(z[~s])) / z[~s] ** 2
        return V0 * out
    return V


k = np.linspace(-30, 30, 4001)
bs = np.concatenate([[0.0], np.geomspace(0.005, 1.5, 28)])
Q = np.zeros((len(bs), 4, 4)); QP = np.full((len(bs), 4, 4), np.nan)
for i, b in enumerate(bs):
    Vq = lambda xx: xx ** 2 / 2 + b * xx ** 4
    x, d, e, v = sinc_dvr(Vq, 7.0, 351)
    for n in range(4):
        psi = v[:, n] / np.sqrt(d)
        Q[i, n] = measures_psi(x, d, psi, k, E=e[n], Vx=Vq(x))
        if b <= 0.1:
            c = aho_coeffs(n, b)
            QP[i, n] = measures_psi(x, d, sum(cv * hfun(m, x) for m, cv in c.items()), k)
lams = np.linspace(0.6, 2.0, 22)
Wl = np.zeros((len(lams), 4, 4))
for i, lam in enumerate(lams):
    a = np.pi / lam
    x, d, e, v = box_dvr(V_well(lam), -a, a, 360)
    kk = np.linspace(-40, 40, 4001)
    for n in range(4):
        Wl[i, n] = measures_psi(x, d, v[:, n] / np.sqrt(d), kk, E=e[n], Vx=V_well(lam)(x))

json.dump(dict(b=bs.tolist(), quartic=Q.tolist(), quartic_pt=np.nan_to_num(QP, nan=-1).tolist(),
               lam=lams.tolist(), well=Wl.tolist()), open('ch4_curves.json', 'w'))

fig, axs = plt.subplots(4, 3, figsize=(7.2, 8.6), gridspec_kw={'hspace': 0.55, 'wspace': 0.38})
panels = [
    (0, bs, Q, QP, '$b$', 'quartic'), (2, lams, Wl, None, r'$\lambda$', 'symmetric well ($V_0=1$)')]
for r0, par, D, P, lab, name in panels:
    # Shannon row
    for j, (yl, f) in enumerate([(r'$S_x$', lambda A: A[..., 0]), (r'$S_k$', lambda A: A[..., 1]),
                                 (r'$S_x+S_k$', lambda A: A[..., 0] + A[..., 1])]):
        ax = axs[r0, j]
        for n in range(4):
            ax.plot(par, f(D[:, n]), color=COL[n], lw=1.3, label='$n=%d$' % n)
            if P is not None:
                ax.plot(par, f(P[:, n]), color=COL[n], lw=0.9, ls='--')
        if j == 2:
            ax.axhline(1 + np.log(np.pi), color='#8a8984', lw=0.7, ls=':')
        ax.set_xlabel(lab); ax.set_ylabel(yl); ax.set_title(name, fontsize=7.5)
        if name == 'quartic':
            ax.axvspan(1.0, 1.5, color='#eda100', alpha=0.12, lw=0)
    for j, (yl, f, lg) in enumerate([(r'$F_x$', lambda A: A[..., 2], False), (r'$F_k$', lambda A: A[..., 3], False),
                                     (r'$F_xF_k$', lambda A: A[..., 2] * A[..., 3], True)]):
        ax = axs[r0 + 1, j]
        for n in range(4):
            ax.plot(par, f(D[:, n]), color=COL[n], lw=1.3, label='$n=%d$' % n)
            if P is not None:
                ax.plot(par, f(P[:, n]), color=COL[n], lw=0.9, ls='--')
        if lg:
            ax.set_yscale('log')
        ax.set_xlabel(lab); ax.set_ylabel(yl); ax.set_title(name, fontsize=7.5)
        if name == 'quartic':
            ax.axvspan(1.0, 1.5, color='#eda100', alpha=0.12, lw=0)
axs[0, 0].legend(frameon=False, loc='upper right', ncol=2)
fig.savefig('figE_ch4_exact.pdf', bbox_inches='tight'); plt.close(fig)
print('quartic b=0..1.5, n=0 Sx+Sk range %.4f-%.4f' % ((Q[:, 0, 0] + Q[:, 0, 1]).min(), (Q[:, 0, 0] + Q[:, 0, 1]).max()))
for n in range(4):
    s = Wl[:, n, 0] + Wl[:, n, 1]; f = Wl[:, n, 2] * Wl[:, n, 3]
    print('well n=%d: Sx+Sk %.4f..%.4f   FxFk %.3f..%.3f' % (n, s.min(), s.max(), f.min(), f.max()))
