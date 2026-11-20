# Eq. (F.12) (eq:app_IVfin): order of the xi integral and of the y - y' Fourier integral.
# Toy with the structure of the finite term of Row II.IV: h(r) = [cos(2 alpha) + ln(r^2)] e^{-r^2/2}  (the traceless
# tensor term and a collinear logarithm), phase e^{-i k.r/xi}, weight (1-xi)/xi.
#   A: Fourier transform in r at fixed xi first (integrand O(xi) at xi -> 0), then xi.
#   B: xi integral at fixed r first, G(a) = int_0^1 dxi (1-xi)/xi e^{-ia/xi} = (1+ia) E1(ia) - e^{-ia}  (log on k.r = 0), then r.
#   C: B with the r integral cut at |r| < R (absolute convergence check) and with an anisotropic Abel factor.
import numpy as np
from scipy.integrate import quad, dblquad
from scipy.special import exp1, jv
import mpmath as mp
def htil(q):   # 2D FT int d^2r e^{-i q.r} h(r), q along the first axis, |q| = q   (closed forms)
    # cos(2 alpha) part: int dal cos(2al) e^{-i q r cos al} = -2 pi J_2(qr), int r J_2(qr) e^{-r^2/2} = 2(1-e^{-q^2/2})/q^2 - e^{-q^2/2}
    # ln r^2 part: 2 pi d/dnu [2^nu Gamma(1+nu) 1F1(1+nu;1;-q^2/2)] at nu = 0
    q = mp.mpf(q)
    a = -2 * mp.pi * (2 * (1 - mp.exp(-q * q / 2)) / (q * q) - mp.exp(-q * q / 2))
    b = 2 * mp.pi * mp.diff(lambda nu: 2 ** nu * mp.gamma(1 + nu) * mp.hyp1f1(1 + nu, 1, -q * q / 2), 0)
    return float(a + b)
# order A
IA = quad(lambda xi: (1 - xi) / xi * htil(1 / xi) if xi > 1e-6 else 0.0, 0, 1, limit=500, epsabs=1e-12, points=[0.01, 0.1])[0]
# order B
def G(a):
    if abs(a) < 1e-300: a = 1e-300
    return (1 + 1j * a) * exp1(1j * a) - np.exp(-1j * a)
def IB(R=40.0, aniso=0.0, eps=0.0):
    def f(al, part):
        def g(r):
            v = r * np.exp(-r * r / 2) * (np.cos(2 * al) + np.log(r * r)) * G(r * np.cos(al)) * np.exp(-eps * r * (1 + aniso * np.cos(al - 0.4)))
            return v.real if part == 0 else v.imag
        return quad(g, 0, R, limit=400, epsabs=1e-12)[0]
    return sum((1 if p == 0 else 1j) * quad(f, 0, 2 * np.pi, args=(p,), points=[np.pi / 2, 3 * np.pi / 2], limit=400, epsabs=1e-11)[0] for p in (0, 1))
print('order A (Fourier first, then xi):   %.10f' % IA)
vB = IB(); print('order B (xi first, then r):          %.10f %+.2ei' % (vB.real, vB.imag))
vC = IB(aniso=0.7, eps=1e-3); print('order B, anisotropic Abel factor:    %.10f %+.2ei' % (vC.real, vC.imag))
print('A - B = %.2e' % abs(IA - vB))
