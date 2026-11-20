# Refinement of order B/C of c7_order.py: the kappa integral extended to K = 300 and the remaining tail added
# analytically (Hu(kappa) -> -pi f(z) sgn(kappa) for m = 0 at large |kappa|, Hu -> 0 for m = 1).
import numpy as np
from scipy.integrate import quad
from scipy.special import sici
exec(open(__file__.replace('c7_refine.py', 'c7_order.py')).read().split('for m in (0, 1):')[0])
K = 300.0
for m in (0, 1):
    def Hu(kap, a=1.0, U=4000.0):
        return quad(lambda u: h(m, kap, u), -U, a * U, points=[p for p in (-3, -1, 0, 1, 3) if -U < p < a * U], limit=1000, epsabs=1e-13)[0]
    for a in (1.0, 3.0):
        d = np.pi * 0.5 * (Hu(1e-9, a) + Hu(-1e-9, a))
        g = lambda kap, part: (np.real if part == 0 else np.imag)((np.exp(1j * kap) * Hu(kap, a) - np.exp(-1j * kap) * Hu(-kap, a)) / kap)
        brk = list(np.arange(0, K + 1e-9, 2.0))
        pv = sum(quad(g, brk[i], brk[i + 1], args=(0,), limit=200)[0] + 1j * quad(g, brk[i], brk[i + 1], args=(1,), limit=200)[0] for i in range(len(brk) - 1))
        if m == 0:   # tail: (e^{ik}(-pi fz) - e^{-ik}(+pi fz))/k = -2 pi fz cos(k)/k ; int_K^inf cos k/k = -Ci(K)
            pv += -2 * np.pi * fz * (-sici(K)[1])
        v = d + 1j * pv
        print('m=%d  order B with u in (-U, %gU), kappa to %g + analytic tail: %.6f%+.6fi' % (m, a, K, v.real, v.imag), flush=True)
