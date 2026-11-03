"""Numerical check of PV int_{p1}^{p2} dp A/(pA - k B) = log|(p2 A - k B)/(p1 A - k B)|,
and of the i*pi from the +-i eps prescriptions.   python3 check.py"""
import numpy as np
from scipy.integrate import quad

k = 1.0
rows = []
for A, B, p1, p2 in [(0.7, 1.9, 0.0, 10.0), (0.7, 1.9, 1.0, 10.0), (1.3, 0.4, 0.0, 1.0), (0.5, 0.2, 1.0, 10.0)]:
    p0 = k * B / A
    exact = np.log(abs((p2 * A - k * B) / (p1 * A - k * B)))
    if p1 < p0 < p2:
        pv = quad(lambda p: 1.0, p1, p2, weight='cauchy', wvar=p0)[0]
        eps = 1e-3     # +i eps:  int A dp/(pA - kB + i eps) = log(u2 + i eps) - log(u1 + i eps), u = pA - kB
        z = np.log(complex(p2 * A - k * B, eps)) - np.log(complex(p1 * A - k * B, eps))
        re, im = z.real, z.imag
    else:
        pv = quad(lambda p: A / (p * A - k * B), p1, p2)[0]
        re, im = pv, 0.0
    rows.append((A, B, p1, p2, p0, pv, exact, re, im))
    print('A=%.1f B=%.1f (%g,%g) p0=%.3f inside=%s  PV=%.10f  log|.|=%.10f  +i eps: Re=%.6f Im/pi=%.6f'
          % (A, B, p1, p2, p0, p1 < p0 < p2, pv, exact, re, im / np.pi))
