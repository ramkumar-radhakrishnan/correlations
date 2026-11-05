import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from rotorlib import pendulum_parity, classical_cells, shannon

mpl.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
                     'mathtext.fontset': 'dejavuserif', 'font.size': 9,
                     'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'legend.fontsize': 7.5, 'xtick.direction': 'in',
                     'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                     'figure.dpi': 200, 'savefig.bbox': 'tight'})

# ---- fig 4: fringes + deficit vs lambda at generic energies
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9)); fig.subplots_adjust(wspace=0.32)
lam = 800.
E, D, P, N = pendulum_parity(lam)
k = int(np.argmin(np.abs(E / lam - 1.5)))
mm = np.arange(-N, N + 1)
pc = classical_cells(E[k], lam, N)
w = np.abs(mm) < 60
ax[0].plot(mm[w], D[k][w], 'C3-', lw=0.9, label=r'quantum $|c_m|^2$')
ax[0].plot(mm[w], pc[w], 'k-', lw=2.2, alpha=0.35, label=r'classical $w_m$')
ax[0].set_xlabel(r'$m_1$'); ax[0].set_ylabel(r'$p_{m_1}$')
ax[0].set_title(r'(a) $\lambda I/\hbar^2=800$, $E/\lambda=%.2f$' % (E[k] / lam), fontsize=9)
ax[0].legend(frameon=False, loc='upper center')
d = np.load('fig4_data.npz')
for key, c, lab in (('lib', 'C3', r'libration, $E/\lambda\in[0.2,0.7]$'),
                    ('rot', 'C0', r'rotation, $E/\lambda\in[1.5,2.5]$')):
    m = d[key]
    ax[1].errorbar(d['lamsD'], -m[:, 0], yerr=m[:, 1], fmt='o-', color=c, ms=3.5,
                   lw=1, label=lab)
ax[1].axhline(1 - np.log(2), color='k', ls=':', lw=0.9)
ax[1].text(260, 1 - np.log(2) + 0.008, r'$1-\ln 2$', fontsize=7.5)
ax[1].set_xscale('log'); ax[1].set_ylim(0, 0.35)
ax[1].set_xlabel(r'$\lambda I/\hbar^2$'); ax[1].set_ylabel(r'$\langle S_{\rm cl}-S_q\rangle$')
ax[1].set_title('(b) deficit at generic energies', fontsize=9)
ax[1].legend(frameon=False, loc='lower right')
fig.savefig('fig4_fringes.pdf'); plt.close(fig)

# ---- fig 5: deep libration
d3 = np.load('fig3_data.npz')
lamsB, S0 = d3['lamsB'], d3['S0']
c0 = 0.5 * np.log(np.pi * np.e / np.sqrt(2))
ll = np.logspace(0.9, 4.1, 200)
fig, ax = plt.subplots(1, 1, figsize=(3.7, 2.8))
ax.semilogx(lamsB, S0 - 0.25 * np.log(lamsB), 'ko', ms=4, label='exact ground state')
ax.semilogx(ll, c0 + 0 * ll, 'C0--', lw=0.9, label=r'harmonic, $\frac{1}{2}\ln(\pi e/\sqrt{2})$')
ax.semilogx(ll, c0 - 1 / (8 * np.sqrt(2 * ll)), 'C3-', lw=1.0, label='with anharmonic correction')
ax.set_xlabel(r'$\lambda I/\hbar^2$'); ax.set_ylabel(r'$S_0-\frac{1}{4}\ln(\lambda I/\hbar^2)$')
ax.legend(frameon=False, loc='lower right')
fig.savefig('fig5_libration.pdf'); plt.close(fig)
print('done; fringe-state S_q=%.4f S_cl=%.4f L1=%.3f' %
      (shannon(D[k]), shannon(pc), 0.5 * np.abs(D[k] - pc).sum()))
