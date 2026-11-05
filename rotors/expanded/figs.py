import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from rotorlib import *

mpl.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
                     'mathtext.fontset': 'dejavuserif', 'font.size': 9,
                     'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'legend.fontsize': 7.5, 'xtick.direction': 'in',
                     'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                     'figure.dpi': 200, 'savefig.bbox': 'tight'})
COL = {50: '#1b6ca8', 200: '#c0392b', 800: '#2e7d32'}

# ------------------------------------------------------------------ fig 1
lams = np.logspace(-2.3, -0.3, 40)
Sx = []
for l in lams:
    psi, E, st, idx = exact_state(4, 4, 30, l)
    Sx.append(shannon(rho1_diag(psi, st, 30)))
Sx = np.array(Sx)
a = lams / 2
ep = em = a ** 2
p0 = 1 - ep - em
Sw = -p0 * np.log(p0) - ep * np.log(ep) - em * np.log(em)
Swo = -ep * np.log(ep) - em * np.log(em)
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9)); fig.subplots_adjust(wspace=0.30)
ax[0].loglog(lams, Sx, 'k-', lw=2.6, alpha=0.32, label='exact diagonalisation')
ax[0].loglog(lams, Sw, 'C3-', label=r'with the $\mathcal{N}^2$ term')
ax[0].loglog(lams, Swo, 'C0--', label=r'without the $\mathcal{N}^2$ term')
ax[0].set_xlabel(r'$\lambda I/\hbar^2$'); ax[0].set_ylabel(r'$S$')
ax[0].set_title(r'(a) weak coupling, $n=0$', fontsize=9)
ax[0].legend(frameon=False, loc='upper left')
ax[1].semilogx(lams, 100 * (Sw - Sx) / Sx, 'C3-', label=r'with $\mathcal{N}^2$')
ax[1].semilogx(lams, 100 * (Swo - Sx) / Sx, 'C0--', label=r'without $\mathcal{N}^2$')
ax[1].axhline(0, color='0.6', lw=0.6)
ax[1].set_xlabel(r'$\lambda I/\hbar^2$'); ax[1].set_ylabel('relative error in $S$ [%]')
ax[1].set_title('(b) effect of the normalisation term', fontsize=9)
ax[1].legend(frameon=False, loc='center left'); ax[1].set_ylim(-16, 4)
fig.savefig('fig1_weak.pdf'); plt.close(fig)
np.savez('fig1_data.npz', lams=lams, Sx=Sx, Sw=Sw, Swo=Swo)
print('fig1')

# ------------------------------------------------------------------ fig 2
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9)); fig.subplots_adjust(wspace=0.28)
peaks = {}
for lam in (50., 200., 800.):
    E, D, P, N = pendulum_parity(lam)
    sel = np.where((E / lam > -1.0) & (E / lam < 3.0))[0]
    allS = np.array([shannon(D[k]) for k in sel])
    peaks[int(lam)] = E[sel][np.argmax(allS)] / lam
    sel = sel[:: max(1, len(sel) // 55)]
    x = E[sel] / lam
    Sq = np.array([shannon(D[k]) for k in sel])
    Sc = np.array([shannon(classical_cells(E[k], lam, N)) for k in sel])
    ax[0].plot(x, Sq, '-', color=COL[int(lam)], lw=1.2, label=r'$\lambda I/\hbar^2=%d$' % lam)
    ax[0].plot(x, Sc, ':', color=COL[int(lam)], lw=1.4)
    ax[1].plot(x, Sq - Sc, '-', color=COL[int(lam)], lw=1.1)
for a_ in ax:
    a_.axvline(1, color='0.5', ls='--', lw=0.8)
ax[0].text(1.04, 1.9, 'separatrix', fontsize=7, rotation=90, color='0.35')
ax[0].set_xlabel(r'$E/\lambda$'); ax[0].set_ylabel(r'$S$')
ax[0].set_title('(a) quantum (solid) and classical (dotted)', fontsize=9)
ax[0].legend(frameon=False, loc='lower right')
ax[1].axhline(0, color='0.6', lw=0.6)
ax[1].axhline(-(1 - np.log(2)), color='k', ls=':', lw=0.8)
ax[1].text(2.1, -(1 - np.log(2)) + 0.03, r'$-(1-\ln 2)$', fontsize=7)
ax[1].set_xlabel(r'$E/\lambda$'); ax[1].set_ylabel(r'$S_q-S_{\rm cl}$')
ax[1].set_title('(b) absolute deviation', fontsize=9); ax[1].set_ylim(-0.6, 0.5)
fig.savefig('fig3_separatrix.pdf'); plt.close(fig)
print('fig2, entropy maxima at E/lam =', peaks)

# ------------------------------------------------------------------ fig 3
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9)); fig.subplots_adjust(wspace=0.32)
lam = 800.
E, D, P, N = pendulum_parity(lam)
k = int(np.argmin(np.abs(E / lam - 1.5)))
mm = np.arange(-N, N + 1)
pc = classical_cells(E[k], lam, N)
L1 = 0.5 * np.sum(np.abs(D[k] - pc))
w = np.abs(mm) < 60
ax[0].plot(mm[w], D[k][w], 'C3-', lw=0.9, label=r'quantum $|c_m|^2$')
ax[0].plot(mm[w], pc[w], 'k-', lw=2.2, alpha=0.35, label=r'classical $w_m$')
ax[0].set_xlabel(r'$m$'); ax[0].set_ylabel(r'$p_m$')
ax[0].set_title(r'(a) $\lambda I/\hbar^2=800$, $E/\lambda=%.2f$' % (E[k] / lam), fontsize=9)
ax[0].legend(frameon=False, loc='upper right')
lamsB = np.array([10., 40., 160., 640., 2560., 10240.])
S0 = []
for l in lamsB:
    Eb, Db, Pb, Nb = pendulum_parity(l)
    S0.append(shannon(Db[0]))
S0 = np.array(S0)
c0 = 0.5 * np.log(np.pi * np.e / np.sqrt(2))
ll = np.logspace(0.9, 4.1, 200)
ax[1].semilogx(lamsB, S0 - 0.25 * np.log(lamsB), 'ko', ms=4, label='exact ground state')
ax[1].semilogx(ll, c0 + 0 * ll, 'C0--', lw=0.9, label=r'harmonic: $\frac{1}{2}\ln(\pi e/\sqrt{2})$')
ax[1].semilogx(ll, c0 - 1 / (8 * np.sqrt(2 * ll)), 'C3-', lw=1.0,
               label=r'with anharmonic correction')
ax[1].set_xlabel(r'$\lambda I/\hbar^2$'); ax[1].set_ylabel(r'$S_0-\frac{1}{4}\ln(\lambda I/\hbar^2)$')
ax[1].set_title('(b) deep libration', fontsize=9)
ax[1].legend(frameon=False, loc='lower right')
plt.close(fig)   # panels re-drawn in figs_rest.py
print('fig3, L1 =', L1, ' S_q =', shannon(D[k]), ' S_cl =', shannon(pc))
np.savez('fig3_data.npz', lamsB=lamsB, S0=S0)

# ------------------------------------------------------------------ fig 4
lamsD = np.array([200., 800., 3200., 12800., 51200.])
res = {'lib': [], 'rot': []}
for lam in lamsD:
    E, D, P, N = pendulum_parity(lam)
    for key, lo, hi in (('lib', 0.2, 0.7), ('rot', 1.5, 2.5)):
        ks = np.where((E / lam > lo) & (E / lam < hi))[0]
        ks = ks[:: max(1, len(ks) // 20)]
        d = [shannon(D[k]) - shannon(classical_cells(E[k], lam, N)) for k in ks]
        res[key].append((np.mean(d), np.std(d) / np.sqrt(len(d))))
    print('  lam', lam, res['lib'][-1], res['rot'][-1])
fig, ax = plt.subplots(1, 1, figsize=(3.7, 2.8))
for key, c, lab in (('lib', 'C3', r'libration, $E/\lambda\in[0.2,0.7]$'),
                    ('rot', 'C0', r'rotation, $E/\lambda\in[1.5,2.5]$')):
    m = np.array(res[key])
    ax.errorbar(lamsD, -m[:, 0], yerr=m[:, 1], fmt='o-', color=c, ms=3.5, lw=1, label=lab)
ax.axhline(1 - np.log(2), color='k', ls=':', lw=0.9)
ax.text(260, 1 - np.log(2) + 0.008, r'$1-\ln 2$', fontsize=7.5)
ax.set_xscale('log'); ax.set_ylim(0, 0.35)
ax.set_xlabel(r'$\lambda I/\hbar^2$'); ax.set_ylabel(r'$\langle S_{\rm cl}-S_q\rangle$')
ax.legend(frameon=False, loc='lower right')
plt.close(fig)   # re-drawn in figs_rest.py
np.savez('fig4_data.npz', lamsD=lamsD, lib=np.array(res['lib']), rot=np.array(res['rot']))
print('fig4 done')
