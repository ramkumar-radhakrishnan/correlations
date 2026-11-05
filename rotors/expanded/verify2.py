import numpy as np
from rotorlib import *

print('=' * 76)
print('(b) Large-n asymptotics of S (branch state; doublet splitting ~ lam^|n|')
print('    is below machine precision here, so eigh returns the branch state)')
print('   %4s %6s %12s %12s %14s %14s' % ('n', 'lam', 'S exact', 'S_pert full',
                                         'ln(2n/lam)+1/2', 'ln(n/lam)+1/2'))
for n in (12, 20, 30):
    m1, m2 = n, 0
    for lam in (0.05, 0.2):
        psi, E, st, idx = exact_state(m1, m2, 80, lam)
        Sx = shannon(rho1_diag(psi, st, 80))
        a = lam / 2
        ep, em = a ** 2 / (n + 1) ** 2, a ** 2 / (n - 1) ** 2
        p0 = 1 - ep - em
        Sp = -p0 * np.log(p0) - ep * np.log(ep) - em * np.log(em)
        print('   %4d %6.2f %12.4e %12.4e %14.4e %14.4e' % (
            n, lam, Sx, Sp, lam ** 2 / n ** 2 * (np.log(2 * n / lam) + 0.5),
            lam ** 2 / n ** 2 * (np.log(n / lam) + 0.5)))

print('=' * 76)
print('(e) Fringe deficit: <S_q - S_cl> averaged over an energy window,')
print('    compared with -(1 - ln 2) = %.4f' % -(1 - np.log(2)))
for lam in (200., 800., 3200., 12800.):
    E, D, P, N = pendulum_parity(lam)
    for lo, hi, lab in ((0.2, 0.7, 'libration'), (1.5, 2.5, 'rotation ')):
        ks = np.where((E / lam > lo) & (E / lam < hi))[0]
        ks = ks[:: max(1, len(ks) // 25)]
        d = [shannon(D[k]) - shannon(classical_cells(E[k], lam, N)) for k in ks]
        print('   lam=%6.0f  %s  E/lam in [%.1f,%.1f]: mean diff = %+.4f  (sd %.4f, %d states)'
              % (lam, lab, lo, hi, np.mean(d), np.std(d), len(d)))
