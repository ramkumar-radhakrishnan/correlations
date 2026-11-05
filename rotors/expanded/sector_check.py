import numpy as np
from scipy.special import mathieu_a, mathieu_b
from rotorlib import sector_H, pendulum_parity

lam = 7.3; q = 2 * lam
for M in (0, 2, 1, 3):
    H, st, idx = sector_H(M, 80, lam)
    w = np.sort(np.linalg.eigvalsh(H))[:8] - M * M / 4.0     # subtract hbar^2 M^2/(4I)
    if M % 2 == 0:
        ref = sorted([mathieu_a(2 * r, q) for r in range(6)] + [mathieu_b(2 * r + 2, q) for r in range(6)])
    else:
        ref = sorted([mathieu_a(2 * r + 1, q) for r in range(6)] + [mathieu_b(2 * r + 1, q) for r in range(6)])
    print('M=%d : E - M^2/4 =' % M, np.round(w, 6))
    print('      a/4 Mathieu =', np.round(np.array(ref[:8]) / 4, 6))

print('\nlibrating states below the separatrix vs (4/pi) sqrt(q), q = 2 lam I/hbar^2')
for lam in (50., 200., 800., 3200.):
    E, D, P, N = pendulum_parity(lam)
    print('   lam=%6.0f  N_lib = %4d   (4/pi)sqrt(2 lam) = %.1f' % (lam, (E < lam).sum(),
          4 / np.pi * np.sqrt(2 * lam)))
