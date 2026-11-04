r"""Checks for G2IV_full_derivation.tex.   python3 checks.py -> checks.json"""
import json
import numpy as np
import sympy as sp
from scipy.integrate import quad

out = {}
eps, gE, L, a, xi, Nc, al, b0 = sp.symbols('epsilon gamma_E L a xi N_c alpha_s beta_0', positive=True)
g = sp.symbols('g', positive=True)
d = 2 - eps
# 1. bracket
B = -d / eps - d * gE - d / 2 * L - 1
out['bracket = -2/eps - (d/2)(2 gamma_E + L)'] = sp.simplify(B - (-2 / eps - d / 2 * (2 * gE + L))) == 0
# 2. prefactor: g^2/(4 pi^2) N_c = (alpha_s/2pi) 2 N_c
out['g^2 Nc/(4pi^2) = (alpha_s/2pi) 2Nc'] = sp.simplify((g**2 * Nc / (4 * sp.pi**2)).subs(g, sp.sqrt(4 * sp.pi * al)) - al / (2 * sp.pi) * 2 * Nc) == 0
# 3. (xi/(1-xi)_+ = xi [1/(1-xi)]_+ acting on xi G, subtraction at xi = 1, no log(1-a))
# where line 2's log(k+/Lambda) + log(1-a) come from: the cut xi <= k+/(k+ + Lambda) on xi/(1-xi)
kp, Lam, vee = 1.0, 1e-5, 50.0
av, xmax = kp / vee, kp / (kp + Lam)
G = lambda x: np.cos(2.3 * x) * x**0.5 + 0.4                                  # noqa: E731  (smooth test function)
cut = quad(lambda x: (x * G(x) - G(1.0)) / (1 - x), av, xmax, limit=400)[0] + G(1.0) * np.log((1 - av) / (1 - xmax))  # exact split of the cut integral
yours = quad(lambda x: (x * G(x) - G(1.0)) / (1 - x), av, 1, limit=400)[0] + G(1.0) * (np.log(kp / Lam) + np.log(1 - av))
out['cut integral vs your plus + G(1)[log(k+/Lambda)+log(1-a)] (Lambda/k+=1e-5; diff = O(Lambda/k+))'] = [cut, yours, abs(cut - yours)]
# 4. your plus + log(1-a) G(1) = standard plus on [a,1]: identity by definition; ∫_a^1 (1-xi)/xi
out['int_a^1 (1-xi)/xi'] = str(sp.simplify(sp.integrate((1 - xi) / xi, (xi, a, 1))))
# 5. Fourier transform of log r^2 (k != 0): -4 pi/q^2, from Laplacian(log r^2) = 4 pi delta;
#    numerical: 2 pi int r J0(qr) log(r^2) e^{-delta r^2} dr -> -4 pi/q^2
from scipy.special import j0
q = 1.7
vals = []
for dl in (1e-2, 3e-3, 1e-3):
    f = lambda r: 2 * np.pi * r * j0(q * r) * np.log(r * r) * np.exp(-dl * r * r)   # noqa: E731
    vals.append(quad(f, 0, 12 / np.sqrt(dl), limit=4000)[0])
out['FT log r^2 at q=1.7 (delta -> 0) vs -4pi/q^2'] = [vals, -4 * np.pi / q**2]
# 6. H4 normalisation: (alpha_s/2pi) 2Nc (4pi/k^2) = g^2 Nc/(pi k^2)
k = sp.symbols('k', positive=True)
out['(alpha_s/2pi) 2Nc 4pi/k^2 = g^2 Nc/(pi k^2)'] = sp.simplify((al / (2 * sp.pi) * 2 * Nc * 4 * sp.pi / k**2).subs(al, g**2 / (4 * sp.pi)) - g**2 * Nc / (sp.pi * k**2)) == 0
# 7. beta0/2 delta part + coupling counterterm -> running coupling; mu_R independence
muR, r, c0 = sp.symbols('mu_R r c_0', positive=True)
rc = al * b0 / (4 * sp.pi) * (2 / eps + sp.log(r**2 * muR**2 / c0**2)) - al * b0 / (4 * sp.pi) * 2 / eps
out['RC after counterterm'] = str(sp.simplify(rc))
out['d/dlog muR^2 [LO(alpha(muR)) + RC] = 0 at O(alpha^2)'] = sp.simplify(-b0 * al**2 / (4 * sp.pi) + al * sp.diff(rc.subs(muR, sp.exp(sp.Symbol('t') / 2)), sp.Symbol('t'))) == 0
json.dump(out, open('checks.json', 'w'), indent=1, default=str)
for kk, v in out.items():
    print(kk, ':', v)
