import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.linalg import eigh
plt.rcParams.update({'font.size': 9, 'axes.labelsize': 10, 'legend.fontsize': 8, 'axes.linewidth': 0.6})
BLUE, ORANGE, AQUA, YEL, INK, GREY = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#2b2b29', '#8a8984'
R = json.load(open('results.json'))

# Fig 1: corrected Fig. 3.2 (W and Q for n = 0, 1, 5)
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.9), gridspec_kw={'hspace': 0.45, 'wspace': 0.42})
for j, n in enumerate((0, 1, 5)):
    d = np.load('ps_n%d.npz' % n); x = d['x']
    W, Q = d['W'], d['Q']
    ax = axs[0, j]
    im = ax.pcolormesh(x, x, W.T, cmap='RdBu_r', norm=TwoSlopeNorm(0, -1 / np.pi, 1 / np.pi), shading='nearest', rasterized=True)
    ax.set_title(r'$W_{%d}(x,p)$, range [%.3f, %.3f]' % (n, W.min(), W.max()), fontsize=8)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=7)
    ax = axs[1, j]
    im = ax.pcolormesh(x, x, Q.T, cmap='Blues', shading='nearest', rasterized=True)
    th = np.linspace(0, 2 * np.pi, 300)
    if n > 0:
        ax.plot(np.sqrt(2 * n) * np.cos(th), np.sqrt(2 * n) * np.sin(th), color=ORANGE, lw=0.8, ls='--')
    ax.set_title(r'$H_{%d}(x,p)$, max %.4f' % (n, Q.max()), fontsize=8)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=7)
    for a_ in axs[:, j]:
        a_.set_xlim(-5, 5); a_.set_ylim(-5, 5); a_.set_aspect('equal')
        a_.set_xlabel('$x$'); a_.set_ylabel('$p$')
fig.savefig('figA_wigner_husimi.pdf', bbox_inches='tight', dpi=300); plt.close(fig)

# Fig 2: corrected Fig. 3.4 (position marginals)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.3), gridspec_kw={'wspace': 0.32})
for j, n in enumerate((0, 1, 5)):
    d = np.load('ps_n%d.npz' % n); x = d['x']
    ax = axs[j]
    ax.plot(x, d['Wx'], color=BLUE, lw=1.4, label=r'$\rho^{W}_{x}=|\psi|^{2}$')
    ax.plot(x, d['Qx'], color=ORANGE, lw=1.2, ls='--', label=r'$\rho^{H}_{x}$')
    ax.set_xlim(-5.5, 5.5); ax.set_xlabel('$x$'); ax.set_title('$n=%d$' % n, fontsize=9)
    if j == 0:
        ax.legend(frameon=False)
fig.savefig('figB_marginals.pdf', bbox_inches='tight'); plt.close(fig)

# Fig 3: Chapter 3 measures vs lambda
c3 = R['ch3']
lam = [r['lam'] for r in c3['0']]
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.4), gridspec_kw={'hspace': 0.5, 'wspace': 0.45})
spec = [(('Sx', 'Sp'), r'$S_x$ (solid), $S_p$ (dashed)'), (('SxH', 'SpH'), r'$S^{H}_x$, $S^{H}_p$'),
        (('SH',), r'Wehrl $S_H$'), (('Fx', 'Fp'), r'$F_x$, $F_p$'), (('R2H', 'R4H'), r'$R^{H}_2$, $R^{H}_4$'),
        (('KL', 'DCS'), r'$S_{KL}(\rho^W_x\Vert\rho^H_x)$, $D_{CS}(W,H)$')]
for ax, (keys, title) in zip(axs.ravel(), spec):
    for n, col in (('0', BLUE), ('1', ORANGE)):
        for kk, ls in zip(keys, ('-', '--')):
            ax.plot(lam, [r[kk] for r in c3[n]], ls, color=col, marker='o', ms=3, lw=1.1)
    ax.set_title(title, fontsize=8); ax.set_xlabel(r'$\lambda$')
axs[0, 0].plot([], [], color=BLUE, label='$n=0$'); axs[0, 0].plot([], [], color=ORANGE, label='$n=1$')
axs[0, 0].legend(frameon=False, fontsize=7)
fig.savefig('figC_ch3_measures.pdf', bbox_inches='tight'); plt.close(fig)

# Fig 4: validity of first-order PT (fidelity) and S_x+S_k for many states
import importlib.util, sys
from numpy.polynomial.hermite import hermval
def ymat(D):
    a = np.diag(np.sqrt(np.arange(1, D)), 1); return (a + a.T) / np.sqrt(2)
D = 200; Y = ymat(D + 20); Y4 = (Y @ Y @ Y @ Y)[:D, :D]
bs = np.linspace(0.001, 1.5, 60)
fid = {n: [] for n in range(4)}
for bb in bs:
    w_, v_ = eigh(np.diag(np.arange(D) + 0.5) + bb * Y4)
    for n in range(4):
        vec = np.zeros(D); vec[n] = 1
        for m in range(D):
            if m != n and abs(Y4[m, n]) > 1e-14:
                vec[m] = bb * Y4[m, n] / (n - m)
        vec /= np.linalg.norm(vec)
        fid[n].append(abs(v_[:, n] @ vec) ** 2)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw={'wspace': 0.3})
ax = axs[0]
for n, col in zip(range(4), (BLUE, ORANGE, AQUA, YEL)):
    ax.plot(bs, fid[n], color=col, lw=1.4, label='$n=%d$' % n)
ax.axvspan(1.0, 1.5, color=GREY, alpha=0.15, lw=0)
ax.text(1.25, 0.95, 'thesis\nrange', ha='center', va='top', fontsize=7)
ax.set_xlabel('$b$'); ax.set_ylabel(r'$|\langle\psi_{\rm exact}|\psi_{\rm 1st}\rangle|^2$')
ax.set_title('(a) quartic: overlap of first-order state with exact', fontsize=8.5)
ax.legend(frameon=False, loc='lower left'); ax.set_ylim(0, 1.02)
ax = axs[1]
nn = np.arange(0, 101, 2)
ax.plot(nn, R['many']['quartic_b1'], 'o-', color=BLUE, ms=2.5, lw=1, label=r'quartic, $b=1$')
ax.plot(nn, R['many']['well_lam1'], 's-', color=ORANGE, ms=2.5, lw=1, label=r'symmetric well, $\lambda=1$')
ax.axhline(1 + np.log(np.pi), color=GREY, lw=0.7, ls=':')
ax.text(98, 1 + np.log(np.pi) - 0.17, r'$1+\ln\pi$', ha='right', fontsize=7)
ax.set_xlabel('$n$'); ax.set_ylabel(r'$S_x+S_k$'); ax.legend(frameon=False, loc='center right')
ax.set_title('(b) exact diagonalization, even $n$ up to 100', fontsize=8.5)
fig.savefig('figD_validity_many.pdf', bbox_inches='tight'); plt.close(fig)
print('figures written')
