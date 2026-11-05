import numpy as np
from rotorlib import *
for lam in (51200.,):
    E, D, P, N = pendulum_parity(lam)
    for lo, hi, lab in ((0.2, 0.7, 'libration'), (1.5, 2.5, 'rotation ')):
        ks = np.where((E / lam > lo) & (E / lam < hi))[0]
        ks = ks[:: max(1, len(ks) // 20)]
        d = [shannon(D[k]) - shannon(classical_cells(E[k], lam, N)) for k in ks]
        print('lam=%6.0f  %s: mean diff = %+.4f  (sd %.4f, %d states)'
              % (lam, lab, np.mean(d), np.std(d), len(d)))
