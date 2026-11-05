import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import jv
from rotorlib import sector_H, shannon

mpl.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
                     'mathtext.fontset': 'dejavuserif', 'font.size': 9,
                     'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'legend.fontsize': 7.5, 'xtick.direction': 'in',
                     'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                     'figure.dpi': 200, 'savefig.bbox': 'tight'})

def S_bessel(g):
    jmax = int(g + 40 * max(g, 1) ** (1 / 3) + 40)
    j = np.arange(-jmax, jmax + 1)
    return shannon(jv(j, g) ** 2)

def S_cl_binned(g):
    jmax = int(g) + 2
    e = np.clip(np.arange(-jmax - 0.5, jmax + 1.5), -g, g)
    return shannon((np.arcsin(e[1:] / g) - np.arcsin(e[:-1] / g)) / np.pi)

def branch_entropy(n, lam):
    """Exact branch entropy in sector M=0 for the Stark state centred on m1=n/2."""
    m1 = n // 2
    N = m1 + int(lam / n) + 60
    H, st, idx = sector_H(0, N, lam)
    w, v = np.linalg.eigh(H)
    k = int(np.argmin(np.abs(w - m1 ** 2)))
    p = np.zeros(2 * N + 1)
    for i, (a, b) in enumerate(st):
        p[a + N] += v[i, k] ** 2
    pos, neg = p[N + 1:], p[:N][::-1]  # the two exchange branches
    pb = pos if pos.sum() >= neg.sum() else neg
    return shannon(pb / pb.sum())

g = np.logspace(-2, 3, 160)
SB = np.array([S_bessel(x) for x in g])
pts = []
for n in (100, 200, 400):
    for gam in (0.1, 0.3, 1.0, 3.0, 10.0, 20.0):
        if gam / n < 0.06:
            pts.append((n, gam, branch_entropy(n, gam * n)))
pts = np.array(pts)

gd = np.logspace(1, 6, 21)
D = np.array([S_cl_binned(x) - S_bessel(x) for x in gd])

fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9)); fig.subplots_adjust(wspace=0.32)
ax[0].semilogx(g, SB, 'k-', lw=1.3, label=r'$-\sum_j J_j^2\ln J_j^2$')
gs = g[g < 1.5]
ax[0].semilogx(gs, gs ** 2 * (np.log(2 / gs) + 0.5), 'C0--', lw=1.0,
               label=r'perturbative $\gamma^2[\ln(2/\gamma)+\frac{1}{2}]$')
gl = g[g > 2]
ax[0].semilogx(gl, np.log(np.pi * gl / 2), 'C1:', lw=1.2, label=r'classical $\ln(\pi\gamma/2)$')
ax[0].semilogx(gl, np.log(np.pi * gl / 2) - (1 - np.log(2)), 'C3-.', lw=1.0,
               label=r'$\ln(\pi\gamma/2)-(1-\ln 2)$')
mk = {100: 'o', 200: 's', 400: '^'}
for n in (100, 200, 400):
    sel = pts[:, 0] == n
    ax[0].plot(pts[sel, 1], pts[sel, 2], mk[n], color='C2', ms=4, mfc='none',
               label=r'exact, $n=%d$' % n)
ax[0].set_xlabel(r'$\gamma=\lambda I/\hbar^2 n$'); ax[0].set_ylabel(r'$S_{\rm branch}$')
ax[0].set_title('(a) deep rotation: universal crossover', fontsize=9)
ax[0].legend(frameon=False, loc='upper left', fontsize=6.5); ax[0].set_ylim(-0.1, 7.5)
ax[1].semilogx(gd, D, 'ko-', ms=3, lw=1, label=r'$S_{\rm cl}-S_q$ (Bessel)')
ax[1].axhline(1 - np.log(2), color='C3', ls=':', lw=1)
ax[1].text(15, 1 - np.log(2) + 0.006, r'$1-\ln 2$', fontsize=8, color='C3')
c = (1 - np.log(2) - D[-1]) * gd[-1] ** (1 / 3)
ax[1].semilogx(gd, 1 - np.log(2) - c * gd ** (-1 / 3), 'C0--', lw=0.9,
               label=r'$1-\ln 2-c\,\gamma^{-1/3}$')
ax[1].set_xlabel(r'$\gamma$'); ax[1].set_ylabel(r'$S_{\rm cl}-S_q$')
ax[1].set_title(r'(b) approach to $1-\ln 2$', fontsize=9)
ax[1].legend(frameon=False, loc='lower right'); ax[1].set_ylim(0, 0.34)
fig.savefig('fig2_rotation.pdf'); plt.close(fig)
print('collapse points (n, gamma, S_exact, S_Bessel):')
for n, gam, s in pts:
    print('   %4d %6.2f %.5f %.5f' % (n, gam, s, S_bessel(gam)))
print('c =', c)
np.savez('fig2_data.npz', g=g, SB=SB, pts=pts, gd=gd, D=D)
