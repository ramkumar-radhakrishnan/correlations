"""Checks for Part2b_UV_finite.pdf: part 2b split into a UV term (dim reg, MS-bar) and a finite term, d_perp kept general.

 1. short distance x -> z with d_perp != d: r^2 <V> -> [(d_perp/d) xi(1-xi) + xi/(1-xi) + (1-xi)/xi] (y-z)/(y-z)^2, d = 2, 3
 2. exact p+ integral with d_perp/d (sympy)
 3. d = 2 - eps: (d_perp/d)[2/eps + 2 gamma + L] = (d_perp/2)[2/eps + 2 gamma + L + 1] + O(eps) (sympy)
 4. d_perp = 2 - eps reproduces -11/6 + 4 lam - lam^2 + 2 lam^3/3 (sympy)
 5. the finite term with d_perp general: r^2 <rest> -> 0 at small |x-z| (2 dims, phase kept, smooth test colour factor)
 6. prefactor (sympy)
Functions for the bracket are taken from part2b_checks.py. -> part2b_dperp_checks.json"""
import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'part2b_checks.py')).read()
exec(src[:src.index('# 1. short distance')])          # sq, brackets, coeffs, phihat, V, S, sphere
out = {}

# 1. short distance with d_perp != d
rng = np.random.default_rng(5)
worst = 0.0
for n in (2, 3):
    Y = rng.normal(size=n)
    z = np.zeros(n)
    dirs, wts = sphere(n, 64)
    for dperp in (1.0, 2.0, 3.0, 5.0):
        for xi in (0.1, 0.37, 0.8):
            r = 1e-4
            avg = sum(wq * V(xi, u * r, Y, z, dperp) for u, wq in zip(dirs, wts))
            got = r * r * avg * sq(Y) / Y
            exp = (dperp / n) * xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi
            worst = max(worst, float(np.max(np.abs(got - exp)) / exp))
out['short_distance_dperp_ne_d_max_rel'] = worst

# 2. exact p+ integral, lam = Lambda/k+
xi, lam, dp, d, eps, gam, L = sp.symbols('xi lambda d_perp d epsilon gamma_E L', positive=True)
I = sp.integrate((dp / d) * xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi, (xi, lam, 1 - lam))
target = (dp / d) * (sp.Rational(1, 6) - lam**2 + sp.Rational(2, 3) * lam**3) + 2 * sp.log((1 - lam) / lam) - 2 + 4 * lam
diffs = [abs(complex((I - target).subs({lam: lv, dp: dv, d: dd}).evalf())) for lv in (0.01, 0.1, 0.3, 0.45)
         for dv in (2, 2.7) for dd in (1.9, 2.0)]
out['p_integral_minus_closed_form_max'] = max(diffs)

# 3. d = 2 - eps
expr = (dp / (2 - eps)) * (2 / eps + 2 * gam + L) - (dp / 2) * (2 / eps + 2 * gam + L + 1)
out['dperp_over_d_expansion_O(eps^0)'] = str(sp.simplify(sp.series(expr, eps, 0, 1).removeO()))

# 4. d_perp = 2 - eps (CDR)
E = 2 / eps + 2 * gam + L
full = (2 * sp.log((1 - lam) / lam) - 2 + 4 * lam) * E + ((2 - eps) / 2) * (sp.Rational(1, 6) - lam**2 + sp.Rational(2, 3) * lam**3) * (E + 1)
cdr = (2 * sp.log((1 - lam) / lam) - sp.Rational(11, 6) + 4 * lam - lam**2 + sp.Rational(2, 3) * lam**3) * E
out['CDR_minus_previous_result_O(eps^0)'] = str(sp.simplify(sp.series(full - cdr, eps, 0, 1).removeO()))

# 5. finite term with d_perp general (2 dims): r^2 <2b - CT> -> 0, CT coefficient (d_perp/2) xi(1-xi) + xi/(1-xi) + (1-xi)/xi
kv = np.array([1.3, 0.4])
Yv = np.array([0.8, -0.5])
G = lambda X: np.cos(X @ np.array([0.7, -0.3])) + 0.2 * (X @ X)
dirs, wts = sphere(2, 256)
rows = []
for dperp in (2.0, 3.0):
    for xi in (0.02, 0.37, 0.95):
        for r in (1e-2, 1e-3):
            o = sum(wq * np.exp(1j * xi * kv @ (u * r)) * (V(xi, u * r, Yv, np.zeros(2), dperp) @ np.array([1.0, 0.6])) * (G(Yv) - G(u * r))
                    for u, wq in zip(dirs, wts))
            c = sum(wq * ((dperp / 2) * xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi) * (Yv @ np.array([1.0, 0.6])) / sq(Yv)
                    * S(u * r, Yv, np.zeros(2)) * (G(Yv) - G(np.zeros(2))) for u, wq in zip(dirs, wts))
            rows.append(dict(dperp=dperp, xi=xi, r=r, r2_orig=float(abs(r * r * o)), r2_rest=float(abs(r * r * (o - c)))))
out['rest_small_r'] = rows

# 6. prefactor: (g^4/8 pi^4)(1/k+) x (1/2pi) x pi x N_c = g^4 N_c/(16 pi^4)(1/k+)
g, Nc = sp.symbols('g N_c', positive=True)
out['prefactor_minus_g4Nc_16pi4'] = str(sp.simplify(g**4 / (8 * sp.pi**4) / (2 * sp.pi) * sp.pi * Nc - g**4 * Nc / (16 * sp.pi**4)))

json.dump(out, open(os.path.join(HERE, 'part2b_dperp_checks.json'), 'w'), indent=1)
for k, v in out.items():
    print(k, ':', v)
