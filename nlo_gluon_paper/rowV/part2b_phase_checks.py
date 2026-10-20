"""Checks for Part2b_UV_phase.pdf: UV term of part 2b with the phase exp(-ik.(z-x)p+/k+) kept.
 1. J(q) = (mu^2 e^gamma/4pi)^(eps/2) int d^dx e^{iq.(x-z)} c_d^2 [(x-z)^2]^(1-d) = pi[2/eps - log(q^2/mu^2)] + O(eps)  (sympy)
 2. 2D: int d^2r [e^{iq.r}/r^2 - (Y-r).Y/((Y-r)^2 r^2)] = -pi log(q^2 Y^2 e^{2 gamma}/4)  (numerical, Bessel J0)
 3. A(lam) = int_lam^{1-lam} (xi/(1-xi) + (1-xi)/xi) log xi,  B(lam) = int_lam^{1-lam} xi(1-xi) log xi: closed forms (numerical)
 4. lam -> 0 limits: -2A -> log^2 lam - 4 + pi^2/3,  -2B -> 5/18
 5. finite term with the phase in the counterterm: r^2 <2b - CT_phase> -> 0 at small |x-z| (2 dims)
 6. large |x-z|: int d^2r e^{iq.r}/r^2 over R1<r<R2 converges (oscillation) for q != 0
-> part2b_phase_checks.json"""
import json
import os

import numpy as np
import sympy as sp
from mpmath import mp, polylog
from scipy.integrate import quad
from scipy.special import j0

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'part2b_checks.py')).read()
exec(src[:src.index('# 1. short distance')])          # sq, brackets, coeffs, phihat, V, S, sphere
out = {}

# 1. dim reg with the phase: int d^dr e^{iq.r} (r^2)^(-a) = pi^(d/2) 2^(d-2a) Gamma(d/2-a)/Gamma(a) (q^2)^(a-d/2), a = d-1
eps, mu, q = sp.symbols('epsilon mu q', positive=True)
d = 2 - eps
a = d - 1
cd2 = (sp.gamma(d / 2) * sp.pi**(1 - d / 2))**2
J = (mu**2 * sp.exp(sp.EulerGamma) / (4 * sp.pi))**(eps / 2) * cd2 * sp.pi**(d / 2) * 2**(d - 2 * a) * sp.gamma(d / 2 - a) / sp.gamma(a) * (q**2)**(a - d / 2)
ser = sp.series(J, eps, 0, 1).removeO()
out['J_minus_pi(2/eps-log(q2/mu2))'] = str(sp.simplify(sp.expand_log(ser - sp.pi * (2 / eps - sp.log(q**2 / mu**2)), force=True)))

# 2. the 2D difference between the two counterterms (angular averages: J0(qr) and theta(Y-r)/r^2)
rows = []
for qv, Yv in ((1.3, 0.7), (0.2, 2.0), (3.0, 0.4)):
    f = lambda r: (j0(qv * r) - (1.0 if r < Yv else 0.0)) / r
    val = 2 * np.pi * (quad(f, 0, Yv, limit=400)[0] + quad(lambda r: j0(qv * r) / r, Yv, 200 / qv, limit=4000)[0]
                       + quad(lambda r: j0(qv * r) / r, 200 / qv, 4000 / qv, limit=20000)[0])
    exp_ = -np.pi * np.log(qv**2 * Yv**2 * np.exp(2 * np.euler_gamma) / 4)
    rows.append(dict(q=qv, Y=Yv, numeric=val, closed=exp_))
out['difference_of_counterterms_2D'] = rows

# 3./4. A, B closed forms
mp.dps = 30
def A_closed(l):
    return float(polylog(2, l) - polylog(2, 1 - l)) + 0.5 * np.log(1 - l)**2 - 0.5 * np.log(l)**2 - 2 * (1 - l) * np.log(1 - l) + 2 * l * np.log(l) + 2 * (1 - 2 * l)
def B_closed(l):
    F = lambda x: x**2 / 2 * np.log(x) - x**2 / 4 - x**3 / 3 * np.log(x) + x**3 / 9
    return F(1 - l) - F(l)
ab = []
for l in (0.3, 0.1, 1e-2, 1e-4):
    An = quad(lambda x: (x / (1 - x) + (1 - x) / x) * np.log(x), l, 1 - l, limit=400)[0]
    Bn = quad(lambda x: x * (1 - x) * np.log(x), l, 1 - l, limit=400)[0]
    ab.append(dict(lam=l, A_rel=abs(An - A_closed(l)) / abs(An), B_rel=abs(Bn - B_closed(l)) / abs(Bn),
                   m2A_minus_limit=-2 * A_closed(l) - (np.log(l)**2 - 4 + np.pi**2 / 3),
                   m2B_minus_limit=-2 * B_closed(l) - 5 / 18))
out['A_B'] = ab

# 5. finite term with the phase kept in the counterterm (pure 1/(x-z)^2 kernel), d_perp = 2, 3
kv = np.array([1.3, 0.4]); Yv = np.array([0.8, -0.5]); t = np.array([1.0, 0.6])
G = lambda X: np.cos(X @ np.array([0.7, -0.3])) + 0.2 * (X @ X)
dirs, wts = sphere(2, 256)
fin = []
for dperp in (2.0, 3.0):
    for xi in (0.02, 0.37, 0.95):
        for r in (1e-2, 1e-3):
            o = sum(wq * np.exp(1j * xi * kv @ (u * r)) * (V(xi, u * r, Yv, np.zeros(2), dperp) @ t) * (G(Yv) - G(u * r)) for u, wq in zip(dirs, wts))
            c = sum(wq * np.exp(1j * xi * kv @ (u * r)) * ((dperp / 2) * xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi) * (Yv @ t) / sq(Yv)
                    / (r * r) * (G(Yv) - G(np.zeros(2))) for u, wq in zip(dirs, wts))
            fin.append(dict(dperp=dperp, xi=xi, r=r, r2_orig=float(abs(r * r * o)), r2_rest=float(abs(r * r * (o - c)))))
out['finite_term_small_r'] = fin

# 6. convergence at large r of int d^2r e^{iq.r}/r^2 = 2 pi int dr J0(qr)/r
tails = []
for xi in (0.02, 0.37):
    qv = xi * np.linalg.norm(kv)
    tails.append(dict(xi=xi, tail_R_to_10R=[2 * np.pi * quad(lambda r: j0(qv * r) / r, R, 10 * R, limit=20000)[0] for R in (10.0, 100.0, 1000.0)]))
out['large_r_tail'] = tails

json.dump(out, open(os.path.join(HERE, 'part2b_phase_checks.json'), 'w'), indent=1)
for k, v in out.items():
    print(k, ':', v)
