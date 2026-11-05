import numpy as np
from scipy.special import jv
from rotorlib import sector_H, shannon

for n, lam in ((200, 400.), (400, 2000.), (1000, 20000.)):
    N = n // 2 + 80
    H, st, idx = sector_H(0, N, lam)
    w, v = np.linalg.eigh(H)
    m1 = n // 2
    E0 = m1 ** 2                       # hbar = I = 1, M = 0: E0 = m^2
    k = int(np.argmin(np.abs(w - E0)))  # Stark state centred on j = 0
    p = np.zeros(2 * N + 1)
    for i, (a, b) in enumerate(st):
        p[a + N] += v[i, k] ** 2
    g = lam / n
    j = np.arange(-int(g) - 30, int(g) + 31)
    pj = p[N + m1 + j]
    print('n=%5d lam=%7.0f gamma=%5.1f  E/lam=%6.1f  max|p_j-J_j^2|=%.2e  '
          'S_exact=%.5f  S_Bessel=%.5f'
          % (n, lam, g, n * n / (4 * lam), np.abs(pj - jv(j, g) ** 2).max(),
             shannon(p), shannon(jv(j, g) ** 2)))
