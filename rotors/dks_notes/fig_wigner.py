"""Fig. 3 of the notes: closed-form Wigner function of the deep-rotation branch state."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.special import jv

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 10, 'legend.fontsize': 8,
                     'axes.linewidth': 0.6, 'lines.linewidth': 1.3})
BLUE, ORANGE, INK = '#2a78d6', '#eb6834', '#2b2b29'
g = 20.0
phi = np.linspace(0, 2 * np.pi, 800)
n2 = np.arange(-64, 65)                 # n2 = 2(k0 - l), i.e. l - k0 = -n2/2
W = jv(n2[:, None], 2 * g * np.cos(phi)[None, :]) / (2 * np.pi)
dl = -n2 / 2.0                          # l - k0
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={'wspace': 0.32, 'width_ratios': [1.1, 1]})
ax = axs[0]
vm = np.quantile(np.abs(W), 0.995)
ax.pcolormesh(phi, dl, W, cmap='RdBu_r', shading='nearest', norm=TwoSlopeNorm(0, -vm, vm), rasterized=True)
ax.plot(phi, -g * np.cos(phi), color=INK, lw=0.8, ls=(0, (3, 2)), label=r'classical: $\ell-k_0=-\gamma\cos\varphi$')
ax.plot(phi, g * np.cos(phi), color=INK, lw=0.6, ls=(0, (1, 1.5)), label=r'ghost: $\varphi\to\varphi+\pi$')
ax.set_xlim(0, 2 * np.pi); ax.set_ylim(-30, 30)
ax.set_xticks([0, np.pi, 2 * np.pi]); ax.set_xticklabels(['0', r'$\pi$', r'$2\pi$'])
ax.set_xlabel(r'$\varphi$'); ax.set_ylabel(r'$\ell-k_{0}$')
ax.set_title(r'(a) $W(\varphi,\ell)=J_{2(k_0-\ell)}(2\gamma\cos\varphi)/2\pi$, $\gamma=20$', fontsize=8.5)
ax.legend(loc='upper center', frameon=True, framealpha=0.85, fontsize=7, borderaxespad=0.2)
ax = axs[1]
G = 60.0
j = np.arange(-50, 51)
beta = np.arccos(j / G)
Phi = 2 * np.sqrt(G * G - j * j) - 2 * j * beta - np.pi / 2
ax.plot(j, 1e3 * jv(2 * j, 2 * G) / (2 * np.pi), 'o', ms=2.6, mfc='none', mec=ORANGE, mew=0.8,
        label=r'exact $W(0,k_0-j)$')
jj = np.linspace(-59.5, 59.5, 2000)
bb = np.arccos(jj / G)
PP = 2 * np.sqrt(G * G - jj * jj) - 2 * jj * bb - np.pi / 2
ax.plot(j, 1e3 * np.sqrt(1 / (np.pi * G * np.sin(beta))) * np.cos(Phi + np.pi / 4) / (2 * np.pi), '.',
        color=BLUE, ms=3.0, label=r'Debye: $\frac{1}{2\pi}(\pi\gamma\sin\beta)^{-1/2}\cos(\Phi_j+\frac{\pi}{4})$')
ax.axhline(0, color='#8a8984', lw=0.5)
ax.set_xlabel(r'$j=k_{0}-\ell$'); ax.set_ylabel(r'$10^{3}\,W(\varphi=0,\ell)$')
ax.set_title(r'(b) chord term at $\varphi=0$, $\gamma=60$', fontsize=8.5)
ax.set_xlim(-52, 52); ax.set_ylim(-13, 19)
ax.legend(loc='upper center', frameon=True, framealpha=0.9, fontsize=7, borderaxespad=0.2)
fig.savefig('fig3_wigner.pdf', bbox_inches='tight', dpi=300)
print('ok')
