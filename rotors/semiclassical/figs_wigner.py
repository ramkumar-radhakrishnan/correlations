"""Figures 5 (Wigner) and 6 (Husimi) of the semiclassical paper; reads wigner_husimi.npz."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 10, 'legend.fontsize': 8,
                     'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'lines.linewidth': 1.4})
BLUE, ORANGE, AQUA, GREY, INK = '#2a78d6', '#eb6834', '#1baf7a', '#8a8984', '#2b2b29'
C0 = 1 - np.log(2)
d = np.load('wigner_husimi.npz')
lam = 200.0

# ------------------------------------------------------------------ Figure 5
fig = plt.figure(figsize=(7.0, 2.55))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.42)
for a, (tag, title) in enumerate((('lib', '(a) libration'), ('rot', '(b) rotation'))):
    ax = fig.add_subplot(gs[0, a])
    W, phi, ells, E = d['W_' + tag], d['phi_' + tag], d['ells_' + tag], float(d['E_' + tag])
    sel = np.abs(ells) <= 26
    vm = np.quantile(np.abs(W[sel]), 0.995)
    dphi = phi[1] - phi[0]
    ax.pcolormesh(phi, ells[sel], W[sel], cmap='RdBu_r', shading='nearest',
                  norm=TwoSlopeNorm(0, -vm, vm), rasterized=True)
    ph = np.linspace(0, 2 * np.pi, 600)
    pp = E - lam * np.cos(ph)
    ok = pp >= 0
    for s in (1, -1):
        ax.plot(ph[ok], s * np.sqrt(pp[ok]), color=INK, lw=0.7, ls=(0, (3, 2)))
    ax.set_xlim(0, 2 * np.pi); ax.set_ylim(-26, 26)
    ax.set_xticks([0, np.pi, 2 * np.pi]); ax.set_xticklabels(['0', r'$\pi$', r'$2\pi$'])
    ax.set_xlabel(r'$\varphi$'); ax.set_ylabel(r'$\ell=p/\hbar$')
    ax.set_title(r'%s, $E_{\rm rel}/\lambda=%.2f$' % (title, E / lam), fontsize=9)

ax = fig.add_subplot(gs[0, 2])
W, phi, ells, E = d['W_rot'], d['phi_rot'], d['ells_rot'], float(d['E_rot'])
for l, col, lab in ((16.0, BLUE, r'$\ell=16$ (bright fringe)'),
                    (14.0, ORANGE, r'$\ell=14$ (dark fringe)')):
    i = int(np.argmin(np.abs(ells - l)))
    ax.plot(phi, 1e2 * W[i], color=col, lw=1.0, label=lab)
    pl = np.arccos((E - l * l) / lam)
    for x in (pl, 2 * np.pi - pl):
        ax.axvline(x, color=col, lw=0.6, ls=(0, (1, 1.5)))
ax.axhline(0, color=GREY, lw=0.5)
ax.set_xlim(0, 2 * np.pi)
ax.set_xticks([0, np.pi, 2 * np.pi]); ax.set_xticklabels(['0', r'$\pi$', r'$2\pi$'])
ax.set_xlabel(r'$\varphi$'); ax.set_ylabel(r'$10^{2}\,W(\varphi,\ell)$')
ax.set_title('(c) cuts of (b)', fontsize=9)
ax.legend(loc='upper center', frameon=False, handlelength=1.4, borderaxespad=0.1, ncol=1)
ax.set_ylim(-3.2, 8.6)
fig.savefig('fig5_wigner.pdf', bbox_inches='tight', dpi=300)
plt.close(fig)

# ------------------------------------------------------------------ Figure 6
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={'wspace': 0.32})
ax = axs[0]
sg, Dg = d['sgrid'], d['Dgrid']
ax.axhline(C0, color=GREY, lw=0.8, ls=(0, (4, 2)))
ax.text(7.9, C0 + 0.006, r'$1-\ln 2$', ha='right', va='bottom', color=INK, fontsize=8)
ax.plot(sg, C0 - Dg, color=INK, lw=1.2, label=r'$(1-\ln2)-D(\sigma)$')
sp = d['sgrid_pts']
ax.plot(sp, d['bes_gap_1000'], 'o', ms=4.5, mfc=BLUE, mec='white', mew=0.6,
        label=r'Bessel, $\gamma=10^{3}$')
ax.plot(sp, d['bes_gap_100'], 's', ms=4.0, mfc='none', mec=ORANGE, mew=1.0,
        label=r'Bessel, $\gamma=10^{2}$')
ax.plot(sp, d['bes_gap_10000'], '^', ms=4.5, mfc=AQUA, mec=INK, mew=0.5,
        label=r'Bessel, $\gamma=10^{4}$')
ax.set_xlabel(r'coherent-state width $\sigma$ (units of $\hbar$)')
ax.set_ylabel(r'$S_{W}-S_{q}$')
ax.set_xlim(0.4, 8.1); ax.set_ylim(0.12, 0.36)
ax.set_title('(a) deep rotation', fontsize=9)
ax.legend(loc='lower right', frameon=False)

ax = axs[1]
r = d['gen_800']
x = r[:, 0]
ax.axhline(0, color=GREY, lw=0.5)
ax.axvline(1, color=GREY, lw=0.6, ls=(0, (4, 2)))
ax.plot(x, r[:, 3] - r[:, 1], '-', color=GREY, lw=1.0, label=r'unsmoothed: $S^{\rm bin}_{\rm cl}-S_{q}$')
for a, (s, col, mk) in enumerate(((0.5, ORANGE, 's'), (1.0, AQUA, '^'), (2.0, BLUE, 'o'))):
    SWp, SWw = r[:, 4 + 3 * a], r[:, 5 + 3 * a]
    ax.plot(x, SWw - SWp, '-', color=col, lw=1.0, marker=mk, ms=3.2, mfc=col, mec='white',
            mew=0.4, label=r'$S_{W}[w]-S_{W}[p]$, $\sigma=%g$' % s)
ax.set_xlabel(r'$E_{\rm rel}/\lambda$'); ax.set_ylabel('entropy difference')
ax.set_title(r'(b) generic energies, $\lambda I/\hbar^{2}=800$', fontsize=9)
ax.set_xlim(-0.9, 3.0); ax.set_ylim(-0.12, 0.8)
ax.legend(loc='upper left', frameon=False, fontsize=7.2, ncol=1, handlelength=1.8)
fig.savefig('fig6_husimi.pdf', bbox_inches='tight')
plt.close(fig)
print('figures written')
