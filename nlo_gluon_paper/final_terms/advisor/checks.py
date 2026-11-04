r"""Checks for Advisor_reply.tex.   python3 checks.py -> checks.json"""
import json
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import j0, j1

out = {}
al, b0, X, mu, eps, gE, muF, muR = sp.symbols('alpha_s beta_0 X mu epsilon gamma_E mu_F mu_R', positive=True)
g = sp.sqrt(4 * sp.pi * al)
c0 = 2 * sp.exp(-gE)
muhat2 = mu**2 * sp.exp(2 * gE) / 4
# 1. bare RC + MS-bar coupling counterterm = (alpha_s^2 beta0/4pi) ln(X^2 muhat^2) per unit LO/alpha_s
rc = g**2 * b0 / (16 * sp.pi**2) * (2 / eps + 2 * gE + sp.log(X**2 * mu**2 / 4)) - al * b0 / (4 * sp.pi) * 2 / eps
adv = al * b0 / (4 * sp.pi) * sp.log(X**2 * muhat2)
out['RC (after counterterm) - advisor alpha_s beta0/(4pi) ln(X^2 muhat^2)'] = str(sp.simplify(sp.expand_log(rc - adv, force=True)))
out['muhat^2 = mu^2/c0^2'] = str(sp.simplify(muhat2 - mu**2 / c0**2))
# 2. beta0 bookkeeping: delta part of P_gg in DGLAP at mu_F + RC at mu_R = (alpha_s beta0/4pi) ln(muR^2/muF^2)
d_part = -al / (2 * sp.pi) * b0 / 2 * sp.log(X**2 * muF**2 / c0**2)
rcR = al * b0 / (4 * sp.pi) * sp.log(X**2 * muR**2 / c0**2)
out['delta-part(mu_F) + RC(mu_R)'] = str(sp.simplify(sp.expand_log(d_part + rcR, force=True)))
# 3. w, w' integrals of H4: |int d^2w e^{ik.w} (w-z)^i/(w-z)^2|^2 = 4 pi^2/k^2 ; radial form: 2 pi i k^i/k^2
k = 1.3
for dl in (1e-3,):
    rad = quad(lambda r: 2 * np.pi * j1(k * r) * np.exp(-dl * r * r), 0, 40 / np.sqrt(dl), limit=4000)[0]
out['int d^2r e^{ik.r} r^i/r^2 = 2 pi i k^i/k^2 : radial 2pi int J1(kr)dr vs 2pi/k'] = [rad, 2 * np.pi / k]
# 4. y' integral of the subtraction: int d^2r e^{-iq.r} log r^2 = -4pi/q^2
q = 1.3
vals = [quad(lambda r: 2 * np.pi * r * j0(q * r) * np.log(r * r) * np.exp(-dl * r * r), 0, 12 / np.sqrt(dl), limit=4000)[0] for dl in (3e-3, 1e-3)]
out['int d^2r e^{-iq.r} log r^2 (delta->0) vs -4pi/q^2'] = [vals, -4 * np.pi / q**2]
# 5. normalisation: H4 prefactor (g^4 Nc/8pi^5)(4pi^2/k^2) = (g^2 Nc/(pi k^2)) (g^2/2pi^2)
gg, Nc, kk = sp.symbols('g N_c k', positive=True)
out['H4 = g^2 Nc/(pi k^2) LO^y'] = str(sp.simplify(gg**4 * Nc / (8 * sp.pi**5) * 4 * sp.pi**2 / kk**2 - gg**2 * Nc / (sp.pi * kk**2) * gg**2 / (2 * sp.pi**2)))
# 6. subtraction term = -[log(1/a) - (1-a)] H4 ; plus -(1-a) H4 line = -log(1/a) H4
a, xi = sp.symbols('a xi', positive=True)
sub = -sp.integrate((1 - xi) / xi, (xi, a, 1))
out['subtraction + (1-a)-line, in units of H4'] = str(sp.simplify(sub - (1 - a)))
json.dump(out, open('checks.json', 'w'), indent=1, default=str)
for kq, v in out.items():
    print(kq, ':', v)
