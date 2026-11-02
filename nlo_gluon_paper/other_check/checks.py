r"""Checks for Other_terms_analysis.tex.  python3 checks.py  ->  checks.json

1. mu table: coefficient of log mu^2 (= coefficient of 2/eps) of every source, in units of alpha_s/2pi.
2. prefactors: G1 II part 1, G1 V 2b UV, G2 IV lines in units of LO; G2 IV line 2 and H4.
3. plus prescriptions: your [xi + xi^2/(1-xi)_+] = standard xi/(1-xi)_+ on [a,1] minus log(1-a) delta(1-xi).
4. Fourier identity int d^2r e^{-iq.r} log r^2 = -4 pi/q^2 (q != 0), and int_a^1 (1-xi)/xi = log(1/a) - (1-a).
5. fragmentation convolution: int_a^1 dzeta/zeta^3 D(zeta) int_{a/zeta}^1 dxi/xi^3 P(xi) F(z=zeta xi)
   = int_a^1 dz/z^3 F(z) int_z^1 dzeta/zeta P(z/zeta) D(zeta), with P the [0,1] plus distribution cut at xi >= a/zeta.
"""
import json
import numpy as np
import sympy as sp
from scipy import integrate, special

out = {}
g, Nc, eps, a, L = sp.symbols('g N_c epsilon a Lambda_over_kplus', positive=True)
pi = sp.pi
alpha = g**2 / (4 * pi)
LOpref = g**2 / (2 * pi**2)

# ---- 2. prefactors (ratio to the LO prefactor 1/(2pi)^3 g^2/(2pi^2) 1/k+), with the colour sign
r_G1II = sp.simplify(-(g**4 * Nc / (32 * pi**4)) / LOpref / (alpha * Nc / (4 * pi)))      # colour bracket = -LO colour
r_2bUV = sp.simplify((g**4 * Nc / (16 * pi**4)) / LOpref / (alpha * Nc / (2 * pi)))       # colour bracket = +LO colour
r_G2IV = sp.simplify((g**2 / (4 * pi**2)) / (alpha / pi))
# line 2 after the y' integral: -(g^2/4pi^2) N_c (g^2/2pi^2) (4 pi/k^2) [log(vee/k+) - (1-a)];  H4: g^4 N_c/(8 pi^5) (2pi)^2/k^2
kk = sp.symbols('k2', positive=True)
line2 = (g**2 / (4 * pi**2)) * Nc * (g**2 / (2 * pi**2)) * 4 * pi / kk
H4 = g**4 * Nc / (8 * pi**5) * (2 * pi)**2 / kk
out['prefactors'] = {
    'G1 II part 1 / [(alpha_s N_c/4pi) LO], incl. colour sign -1': str(r_G1II),
    'G1 V 2b UV / [(alpha_s N_c/2pi) LO], colour sign +1': str(r_2bUV),
    'G2 IV lines 1,4,5: g^2/(4pi^2) / (alpha_s/pi)': str(r_G2IV),
    'G2 IV line 2 coefficient / H4 coefficient': str(sp.simplify(line2 / H4)),
}
assert r_G1II == -1 and r_2bUV == 1 and r_G2IV == 1 and sp.simplify(line2 / H4) == 1

# ---- 1. mu table (units alpha_s/2pi x LO; P = int dxi/xi^2 P^(r)_{gg,+} LO_xi)
b0 = sp.Rational(11, 3) * Nc
l1 = sp.log(1 - a)
lam = sp.log(1 - L)
P = sp.Symbol('P_conv')
# G1 II part 1 + cc = (alpha_s N_c/4pi)[11/3+l1] LO[E_X+E_X'];  each E has one log mu^2  -> 2 x (1/2) in alpha_s/2pi units
src = {
    'G1 II part 1 + c.c.': Nc * (sp.Rational(11, 3) + l1),
    'G1 V 2b UV + c.c.': 2 * Nc * (2 * lam - sp.Rational(11, 6)),
    'G2 IV lines 1+4': -P,
    'G2 IV line 2': sp.Integer(0),
    'coupling counterterm (mu_R)': -b0 / 2,
    'fragmentation counterterm (mu_F)': P + b0 / 2,
}
tot = sp.simplify(sum(src.values()))
out['mu table (coefficient of log mu^2 = coefficient of 2/eps, units alpha_s/2pi x LO)'] = {k: str(sp.expand(v)) for k, v in src.items()}
out['mu table: sum'] = str(sp.expand(tot))
assert sp.simplify(tot - (Nc * l1 + 4 * Nc * lam)) == 0
out['mu table: beta0 part of G1 II + G1 V 2b UV'] = str(sp.simplify(src['G1 II part 1 + c.c.'].subs(a, 0) + src['G1 V 2b UV + c.c.'].subs(L, 0)))

# ---- 3. plus prescriptions, with a test function G(xi) (= LO_xi)
def G(x):
    return np.cos(3.1 / x) * np.exp(-0.4 / x) / (1 + x * x) + 0.3 * x


res = []
for av in (0.05, 0.3, 0.7):
    user = integrate.quad(lambda x: G(x) / x, av, 1, limit=400)[0] + integrate.quad(lambda x: (G(x) - G(1)) / (1 - x), av, 1, limit=400)[0]
    # standard: int_a^1 dxi/xi^2 xi/(1-xi)_+ G = int_a^1 [G/xi - G(1)]/(1-xi) + G(1) log(1-a)
    std = integrate.quad(lambda x: (G(x) / x - G(1)) / (1 - x), av, 1, limit=400)[0] + G(1) * np.log(1 - av)
    res.append(abs(user - (std - np.log(1 - av) * G(1))))
out['user [xi+xi^2/(1-xi)_+] - (standard xi/(1-xi)_+ - log(1-a) G(1)), a=0.05,0.3,0.7'] = res
assert max(res) < 1e-9

# ---- 4. Fourier identity (Gaussian regulator), and the xi integral
q = 1.7
vals = []
for d in (1e-2, 1e-3, 1e-4):
    f = lambda r: 2 * np.pi * r * special.j0(q * r) * np.exp(-d * r * r) * np.log(r * r)
    v = 0.0
    edges = np.linspace(0, 12 / np.sqrt(d), 4000)
    for lo, hi in zip(edges[:-1], edges[1:]):
        v += integrate.quad(f, lo, hi, limit=100)[0]
    vals.append(v)
out['int d^2r e^{-iq.r} e^{-d r^2} log r^2, d=1e-2,1e-3,1e-4 (target -4pi/q^2)'] = vals + [-4 * np.pi / q**2]
assert abs(vals[-1] + 4 * np.pi / q**2) < 2e-2
xi = sp.symbols('xi', positive=True)
assert sp.simplify(sp.integrate((1 - xi) / xi, (xi, a, 1)) - (-sp.log(a) - (1 - a))) == 0
out['int_a^1 (1-xi)/xi'] = 'log(1/a) - (1-a)'

# ---- 5. convolution with the [0,1] plus distribution restricted to xi >= a/zeta
def Dz(z):
    return z**0.7 * (1 - z)**2.3


def F(z):
    return np.exp(-1.3 / z) * (1 + np.sin(4 * z))


def Preg(x):            # 2N_c[(1-x)/x + x(1-x)] with N_c = 3, plus the plus part handled below
    return 6 * ((1 - x) / x + x * (1 - x))


def conv_gluon(zeta, av):
    """int_{a/zeta}^1 dxi/xi^3 P(xi) F(zeta xi), P = Preg + 6 xi/(1-xi)_+ ([0,1] plus, cut at a/zeta)"""
    lo = av / zeta
    h = lambda x: F(zeta * x) / x**3                     # dxi/xi^3 x test function
    reg = integrate.quad(lambda x: Preg(x) * h(x), lo, 1, limit=400)[0]
    plus = integrate.quad(lambda x: (6 * x * h(x) - 6 * h(1)) / (1 - x), lo, 1, limit=400)[0] + 6 * h(1) * np.log(1 - lo)
    return reg + plus


def PD(z):
    """[P x D](z) = int_z^1 dzeta/zeta P(z/zeta) D(zeta), same plus distribution, in xi = z/zeta: int_z^1 dxi/xi P(xi) D(z/xi)"""
    hh = lambda x: Dz(z / x) / x
    reg = integrate.quad(lambda x: Preg(x) * hh(x), z, 1, limit=400)[0]
    plus = integrate.quad(lambda x: (6 * x * hh(x) - 6 * hh(1)) / (1 - x), z, 1, limit=400)[0] + 6 * hh(1) * np.log(1 - z)
    return reg + plus


cv = []
for av in (0.1, 0.4):
    lhs = integrate.quad(lambda zt: Dz(zt) / zt**3 * conv_gluon(zt, av), av, 1, limit=200)[0]
    rhs = integrate.quad(lambda z: F(z) / z**3 * PD(z), av, 1, limit=200)[0]
    cv.append([lhs, rhs, abs(lhs - rhs) / abs(rhs)])
out['convolution: hadron level = int dz/z^3 F(z) [P x D](z), a=0.1,0.4: [lhs, rhs, rel]'] = cv
assert max(c[2] for c in cv) < 1e-6
# the same without the log(1-a/zeta) term (i.e. without your log(1-k+/vee) line) does not factorize:
def conv_gluon_nolog(zeta, av):
    lo = av / zeta
    h = lambda x: F(zeta * x) / x**3
    reg = integrate.quad(lambda x: Preg(x) * h(x), lo, 1, limit=400)[0]
    return reg + integrate.quad(lambda x: (6 * x * h(x) - 6 * h(1)) / (1 - x), lo, 1, limit=400)[0]


bad = []
for av in (0.1, 0.4):
    lhs = integrate.quad(lambda zt: Dz(zt) / zt**3 * conv_gluon_nolog(zt, av), av, 1, limit=200)[0]
    rhs = integrate.quad(lambda z: F(z) / z**3 * PD(z), av, 1, limit=200)[0]
    bad.append([lhs, rhs, abs(lhs - rhs) / abs(rhs)])
out['convolution without the log(1-a) line: [lhs, rhs, rel] (must differ)'] = bad
assert min(b[2] for b in bad) > 1e-3
json.dump(out, open('checks.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
