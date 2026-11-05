import numpy as np
from scipy.special import jv
from rotorlib import exact_state, rho1_diag, shannon

print('(1) Free rotation: exact Schmidt weights vs J_j(gamma)^2, gamma = 2 alpha/n = lam I/(hbar^2 n)')
for n, lam in ((200, 400.), (400, 2000.)):
    m1, m2 = n // 2, -n // 2                    # M = 0 sector, E/lam = n^2/(4 lam)
    psi, E, st, idx = exact_state(m1, m2, n // 2 + 60, lam)
    p = rho1_diag(psi, st, n // 2 + 60)
    N = n // 2 + 60
    g = lam / n
    j = np.arange(-25, 26)
    pj = p[N + m1 + j]
    print('   n=%d lam=%.0f gamma=%.2f  E/lam=%.1f  max|p_j - J_j^2| = %.2e   S_exact=%.5f  S_Bessel=%.5f'
          % (n, lam, g, n * n / (4 * lam), np.abs(pj - jv(j, g) ** 2).max(), shannon(p),
             shannon(jv(np.arange(-int(g) - 60, int(g) + 61), g) ** 2)))

print()
print('(2) Bessel entropy vs classical arcsine entropy ln(pi gamma/2):  deficit -> 1 - ln 2 = %.5f ?'
      % (1 - np.log(2)))
for g in (10., 30., 100., 300., 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6):
    jmax = int(g + 40 * g ** (1 / 3) + 50)
    j = np.arange(-jmax, jmax + 1)
    pj = jv(j, g) ** 2
    Sq = shannon(pj)
    # classical arcsine density binned onto integer cells
    edges = np.clip(np.arange(-jmax - 0.5, jmax + 1.5), -g, g)
    w = (np.arcsin(edges[1:] / g) - np.arcsin(edges[:-1] / g)) / np.pi
    Sc = shannon(w)
    print('   gamma=%9.0f  sum p=%.12f  S_q=%.6f  S_cl(binned)=%.6f  ln(pi g/2)=%.6f  '
          'S_cl-S_q=%.5f' % (g, pj.sum(), Sq, Sc, np.log(np.pi * g / 2), Sc - Sq))
