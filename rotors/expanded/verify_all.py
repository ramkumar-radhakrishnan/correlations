import numpy as np
from scipy.special import mathieu_a, mathieu_b
from rotorlib import *

LN2 = np.log(2)
print('=' * 76)
print('(a) Second-order energy  E2 = lam^2 I / (2 hbar^2 (n^2-1))')
for (m1, m2) in ((9, 2), (12, 2), (4, 4)):
    n = m1 - m2
    for lam in (0.02, 0.1):
        psi, E, st, idx = exact_state(m1, m2, 40, lam)
        E0 = 0.5 * (m1 ** 2 + m2 ** 2)
        print('   n=%2d lam=%.2f : (E-E0) exact = %+.6e   formula = %+.6e'
              % (n, lam, E - E0, lam ** 2 / (2 * (n * n - 1))))

print('=' * 76)
print('(b) Large-n entropy: exact branch entropy vs asymptotic forms')
print('   %4s %6s %12s %12s %12s %12s' % ('n', 'lam', 'S_exact-ln2', 'S_pert(full)',
                                         'corrected', 'draft (38)'))
for n in (12, 20, 30):
    m1, m2 = n // 2 + 5, -(n - n // 2) + 5
    m1, m2 = m1, m1 - n
    for lam in (0.05, 0.2):
        psi, E, st, idx = exact_state(m1, m2, 80, lam)
        Sx = shannon(rho1_diag(psi, st, 80)) - LN2
        a = lam / 2
        ep, em = a ** 2 / (n + 1) ** 2, a ** 2 / (n - 1) ** 2
        p0 = 1 - ep - em
        Sp = -p0 * np.log(p0) - ep * np.log(ep) - em * np.log(em)
        corr = lam ** 2 / n ** 2 * (np.log(2 * n / lam) + 0.5)
        draft = lam ** 2 / n ** 2 * (np.log(n / lam) + 0.5)
        print('   %4d %6.2f %12.4e %12.4e %12.4e %12.4e' % (n, lam, Sx, Sp, corr, draft))

print('=' * 76)
print('(c) M=0 sector vs Mathieu characteristic values, q = 4 alpha = 2 lam')
lam = 7.3
E, D, P, N = pendulum_parity(lam, N=60)
q = 2 * lam
mat = sorted([mathieu_a(2 * r, q) for r in range(6)] +
             [mathieu_b(2 * r + 2, q) for r in range(6)])
print('   E (exact diag)  :', np.round(E[:8], 6))
print('   a/4 (Mathieu)   :', np.round(np.array(mat[:8]) / 4, 6))

print('=' * 76)
print('(d) Deep libration: S0 - [1/4 ln lam + 1/2 ln(pi e/sqrt2)] vs -1/(8 sqrt(2 lam))')
c0 = 0.5 * np.log(np.pi * np.e / np.sqrt(2))
for lam in (10., 40., 160., 640., 2560., 10240.):
    E, D, P, N = pendulum_parity(lam)
    S0 = shannon(D[0])
    print('   lam=%7.0f  residual=%+.6f  prediction=%+.6f  ratio=%.4f'
          % (lam, S0 - 0.25 * np.log(lam) - c0, -1 / (8 * np.sqrt(2 * lam)),
             (S0 - 0.25 * np.log(lam) - c0) / (-1 / (8 * np.sqrt(2 * lam)))))
