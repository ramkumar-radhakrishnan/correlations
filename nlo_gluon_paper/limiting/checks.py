r"""Checks for Limiting_terms.tex.   python3 checks.py -> checks.json"""
import json
import os

import mpmath as mp
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# 1. limits that are taken analytically
a = sp.symbols('a', positive=True)
out['log(1-a)log(1/a) at a=1e-4,1e-8'] = [float((sp.log(1 - a) * sp.log(1 / a)).subs(a, v)) for v in (1e-4, 1e-8)]
out['limit log(1-a)log(1/a)'] = str(sp.limit(sp.log(1 - a) * sp.log(1 / a), a, 0, '+'))
p2, k2 = sp.symbols('p2 k2', positive=True)
out['R^i log: log(p^2/(k^2+a p^2)) - log(p^2/k^2) at a->0'] = str(sp.limit(sp.log(p2 / (k2 + a * p2)) - sp.log(p2 / k2), a, 0))

# 2. Group I Row II part 1 + c.c. in units of LO: g^4 Nc/(32 pi^4) / (g^2/(2 pi^2)) x 2 (c.c.) x (-1) (colour order)
g, Nc = sp.symbols('g N_c', positive=True)
al = g**2 / (4 * sp.pi)
ratio = (g**4 * Nc / (32 * sp.pi**4)) / (g**2 / (2 * sp.pi**2)) * 2
out['G1 II(1)+cc / LO  = alpha_s Nc/(2 pi) ?'] = str(sp.simplify(ratio / (al * Nc / (2 * sp.pi))))
K = Nc * (sp.Rational(67, 18) - sp.pi**2 / 6)
out['Nc(67/9 - pi^2/3) = 2 K_CMW'] = str(sp.simplify(Nc * (sp.Rational(67, 9) - sp.pi**2 / 3) - 2 * K))

# 3. Group III: p+ -> zeta = k+/(p+ + k+), the plus prescription and the measures (test function F)
kp = 1.3
F = lambda p: np.cos(0.7 * p / kp) * np.exp(-0.2 * p / kp) + 0.3 / (1 + p / kp)   # noqa: E731
G = lambda z: F(kp * (1 - z) / z)                                                     # noqa: E731
from scipy.integrate import quad
lhs = quad(lambda p: kp / (p * (p + kp)) * (F(p) - F(0.0)), 0, np.inf, limit=500)[0]
rhs = quad(lambda z: (G(z) - G(1.0)) / (1 - z), 0, 1, limit=500)[0]
out['[k+/(p+(p++k+))]_+ dp+  =  [1/(1-zeta)]_+ dzeta'] = [lhs, rhs, abs(lhs - rhs)]
lhs = quad(lambda p: F(p) / (p + kp), 0, np.inf, limit=500)[0]
rhs = quad(lambda z: G(z) / z, 0, 1, limit=500)[0]
out['dp+/(p++k+) = dzeta/zeta'] = [lhs, rhs, abs(lhs - rhs)]
F2 = lambda p: np.cos(0.7 * p / kp) * np.exp(-0.2 * p / kp)       # noqa: E731  (must decay: no 1/p+ tail)
lhs = quad(lambda p: F2(p) / kp, 0, np.inf, limit=500)[0]
rhs = quad(lambda z: F2(kp * (1 - z) / z) / z**2, 0, 1, limit=500)[0]
out['dp+/k+ = dzeta/zeta^2'] = [lhs, rhs, abs(lhs - rhs)]

# 4. phases: Group III and Group I Row V as one phase with a weighted centre
rng = np.random.default_rng(3)
k, yp, w, z, x, wp = (rng.normal(size=2) for _ in range(6))
errs = []
for zeta in (0.13, 0.5, 0.91):
    xi = (1 - zeta) / zeta
    l = np.exp(-1j * k @ (yp - w)) * np.exp(-1j * xi * k @ (yp - z))
    r = np.exp(-1j * (k / zeta) @ (yp - zeta * w - (1 - zeta) * z))
    errs.append(abs(l - r))
for xi in (0.13, 0.5, 0.91):
    l = np.exp(-1j * k @ (wp - z)) * np.exp(-1j * xi * k @ (z - x))
    r = np.exp(-1j * k @ (wp - (1 - xi) * z - xi * x))
    errs.append(abs(l - r))
out['phase identities (G3, G1 V)'] = max(errs)

# 5. Group I Row V first part: closed forms with Ein(s) = int_0^1 (1 - e^{-s t})/t dt
def ein(s):
    return mp.quad(lambda t: (1 - mp.e**(-s * t)) / t, [0, 1])
errs = []
for al_ in (0.3, -1.7, 4.2):
    s = 1j * al_
    num1 = mp.quad(lambda t: mp.e**(-1j * al_ * t), [0, 1])
    errs.append(abs(num1 - (1 - mp.e**(-s)) / s))
    num2 = mp.quad(lambda t: (mp.e**(-1j * al_ * t) - 1) / t, [0, 1])
    errs.append(abs(num2 + ein(s)))
    series = mp.nsum(lambda n: (-1)**(n + 1) * s**n / (n * mp.factorial(n)), [1, mp.inf])
    errs.append(abs(series - ein(s)))
    errs.append(abs(ein(s) - (mp.euler + mp.log(s) + mp.e1(s))))
out['G1 V(1) closed forms (Ein)'] = float(max(errs))

# 6. Group III Row III part 2, last line: int_0^inf dp B/(A (p+k)(pA + kB)) = B log(B/A)/(k A (B - A))
errs = []
for A, B in ((0.7, 1.9), (2.0, 0.3), (1.0, 1.0000001)):
    num = quad(lambda p: B / (A * (p + kp) * (p * A + kp * B)), 0, np.inf)[0]
    cf = B * np.log(B / A) / (kp * A * (B - A))
    errs.append(abs(num - cf) / abs(cf))
out['G3 III(2) phase-free line closed form (rel)'] = max(errs)

json.dump(out, open(os.path.join(HERE, 'checks.json'), 'w'), indent=1, default=str)
for kk, v in out.items():
    print(kk, ':', v)

# 7. q+ lines -> xi = 1 - q+/k+: phases, denominators, measures, plus prescription
X = sp.symbols('X1 X2 W1 W2 Z1 Z2 Y1 Y2 K1 K2', real=True)
xv, wv, zv, yv, kv = (sp.Matrix(X[i:i + 2]) for i in (0, 2, 4, 6, 8))
eta, xi_, kp_ = sp.symbols('eta xi kplus', positive=True)
ph_q2a = -(kv.dot(wv - xv)) - kv.dot(xv - zv) * eta                     # G1 V 2a and V(1) line 2 (w' -> w here)
ph_p = -(kv.dot(wv - zv)) - kv.dot(zv - xv) * xi_
out['q-line phase = p-line phase at xi = 1 - eta'] = str(sp.simplify(sp.expand(ph_q2a.subs(eta, 1 - xi_) - ph_p)))
A2, B2 = sp.symbols('A2 B2', positive=True)                               # (y-x)^2, (y-z)^2
q = kp_ * (1 - xi_)
out['denominator (k+-q+)(y-x)^2 + q+(y-z)^2 = k+ D_xi'] = str(sp.simplify((kp_ - q) * A2 + q * B2 - kp_ * (xi_ * A2 + (1 - xi_) * B2)))
Ff = lambda t: np.cos(1.1 * t) + t**2                                   # noqa: E731
lhs = quad(lambda qq: (Ff(1 - qq) - Ff(1.0)) / qq, 0, 1)[0]              # (1/k+) int dq+ [1/q+]_+ f, k+ = 1
rhs = quad(lambda t: (Ff(t) - Ff(1.0)) / (1 - t), 0, 1)[0]
out['(1/k+) int_0^k+ dq+ [1/q+]_+ = int_0^1 dxi [1/(1-xi)]_+'] = abs(lhs - rhs)
errs = []
for al_ in (0.3, -1.7, 4.2):
    num = mp.quad(lambda t: (mp.e**(-1j * al_ * t) - mp.e**(-1j * al_)) / (1 - t), [0, 1])
    errs.append(abs(num + mp.e**(-1j * al_) * ein(-1j * al_)))
out['int_0^1 [1/(1-xi)]_+ e^{-i al xi} = -e^{-i al} Ein(-i al)'] = float(max(errs))

# 8. colour: 2a colour factor vanishes at z = x (so 2a + 2b has the x->z behaviour of 2b alone)
import sys
sys.path.insert(0, os.path.join(HERE, '..', 'three_rho'))
import numerics as nm  # noqa: E402
from scipy.linalg import expm  # noqa: E402
TS, Fs = nm.structure_constants(3)
th = rng.normal(size=8)
Uad = nm.adjoint(expm(1j * sum(th[i] * TS[i] for i in range(8))), TS)
col = -np.einsum('adc,fbe,db,ce->af', Fs, Fs, Uad, Uad) + 3 * Uad
out['-f^{adc} f^{fbe} U^{db}(x) U^{ce}(x) + N_c U^{af}(x) = 0'] = float(np.abs(col).max())

json.dump(out, open(os.path.join(HERE, 'checks.json'), 'w'), indent=1, default=str)
for kk in list(out)[-6:]:
    print(kk, ':', out[kk])
