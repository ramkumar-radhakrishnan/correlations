"""Figures for the DKS-style coupled-rotor paper (hbar = I = 1)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import jv
from checks import S_bessel, branch_state, shannon

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 10, 'legend.fontsize': 8,
                     'axes.linewidth': 0.6, 'lines.linewidth': 1.4})
BLUE, ORANGE, AQUA, GREY, INK = '#2a78d6', '#eb6834', '#1baf7a', '#8a8984', '#2b2b29'
C0 = 1 - np.log(2)

# ---------------------------------------------------------------- Fig. 1: weights
g, n = 30.0, 20000
j, p, e = branch_state(n, g * n)
jj = np.arange(-45, 46)
edges = np.clip(np.arange(-45.5, 46.5), -g, g)
wcl = (np.arcsin(edges[1:] / g) - np.arcsin(edges[:-1] / g)) / np.pi
fig, ax = plt.subplots(figsize=(6.2, 2.7))
ax.bar(jj, wcl, width=0.9, color='#d9d8d3', label=r'classical cell weights (arcsine law)')
ax.plot(jj, jv(jj, g) ** 2, '-', color=BLUE, lw=1.1, label=r'WKB: $J_{j}^{2}(\gamma)$')
sel = np.abs(j) <= 45
ax.plot(j[sel], p[sel], 'o', ms=3.2, mfc='none', mec=ORANGE, mew=0.9,
        label=r'exact diagonalization, $n=20000$')
ax.set_xlabel(r'$j=m_{1}^{0}-m_{1}$')
ax.set_ylabel(r'eigenvalue of $\hat\rho_{1}$')
ax.set_xlim(-42, 42); ax.set_ylim(0, 0.105)
ax.legend(loc='upper center', frameon=False, ncol=1)
ax.set_title(r'$\gamma=\lambda/\hbar|\omega_{1}-\omega_{2}|=30$', fontsize=9)
fig.savefig('fig1_weights.pdf', bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- Fig. 2: entropies
gg = np.logspace(-2, 3, 160)
SB = np.array([S_bessel(x) for x in gg])
F = gg ** 2 / 4
Spt = np.where(F < 0.2, -(1 - 2 * F) * np.log(np.clip(1 - 2 * F, 1e-300, None)) - 2 * F * np.log(F), np.nan)
gpts = np.array([0.5, 1, 2, 5, 10, 20])
Sex = []
for x in gpts:
    jx, px, _ = branch_state(1600, x * 1600)
    Sex.append(shannon(px))
G = np.logspace(1, 6, 26)
gap = np.array([np.log(np.pi * x / 2) - S_bessel(x) for x in G])

fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={'wspace': 0.32})
ax = axs[0]
ax.plot(gg, SB, color=INK, lw=1.5, label=r'quantum, $-\sum_j J_j^2\ln J_j^2$')
ax.plot(gpts, Sex, 'o', ms=4, mfc='none', mec=ORANGE, mew=1.0, label=r'exact diag., $n=1600$')
ax.plot(gg, Spt, '--', color=BLUE, lw=1.2, label=r'$-(1-2F)\ln(1-2F)-2F\ln F$')
m = gg > 0.6
ax.plot(gg[m], np.log(np.pi * gg[m] / 2), ':', color=ORANGE, lw=1.5, label=r'classical, $\Delta=h$')
ax.plot(gg[m], np.log(np.pi * gg[m] / np.e), '-.', color=AQUA, lw=1.3, label=r'classical, $\Delta=eh/2$')
ax.set_xscale('log'); ax.set_xlim(1e-2, 1e3); ax.set_ylim(-0.3, 7.5)
ax.set_xlabel(r'$\gamma$'); ax.set_ylabel('entanglement entropy (one branch)')
ax.legend(loc='upper left', frameon=False, fontsize=7.2)
ax.set_title('(a)', fontsize=9)
ax = axs[1]
ax.axhline(C0, color=GREY, lw=0.8, ls=(0, (4, 2)))
ax.text(1.3e6, C0 + 0.005, r'$1-\ln2$', ha='right', va='bottom', fontsize=8)
ax.plot(G, gap, 'o-', color=INK, ms=3, lw=1.0, label=r'$S_{CE}(\Delta=h)-S_{E}$')
ax.plot(G, C0 - 0.965 * G ** (-1 / 3), '--', color=BLUE, lw=1.1,
        label=r'$1-\ln2-0.965\,\gamma^{-1/3}$')
ax.set_xscale('log'); ax.set_xlim(8, 1.5e6); ax.set_ylim(0.0, 0.34)
ax.set_xlabel(r'$\gamma$'); ax.set_ylabel('entropy difference')
ax.legend(loc='lower right', frameon=False)
ax.set_title('(b)', fontsize=9)
fig.savefig('fig2_entropy.pdf', bbox_inches='tight'); plt.close(fig)
print('done', list(np.round(Sex, 5)))
