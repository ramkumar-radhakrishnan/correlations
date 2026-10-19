"""Checks for RowIV_FF_DGLAP.tex: Row IV as a convolution with LO, the plus prescription with exact limits,
the MS-bar counterterm, DGLAP, and the delta(1-xi) poles of Group I (Row II and Row V 2b, with c.c.).
Results in frag_checks.json."""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad, dblquad

HERE = os.path.dirname(os.path.abspath(__file__))
R = {}

xi, d, Eb, E, g, Nc, kp, al, eps, ell = sp.symbols('xi d_perp Ebar E g N_c kp alpha_s epsilon ell')
XX, XD = sp.symbols('XX XD')
Phi = xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi

# 1. the user's xi form vs. RowIV_dperp.pdf Eq. (1.5) (E form), both on int_{xi1}^{xi0}
user = (g**4 * d / (4 * sp.pi**3)) * Nc / (2 * sp.pi) / xi**2 * Phi * (E / 2) * XX \
    + (g**4 / (4 * sp.pi**3)) * Nc / (2 * sp.pi) / xi**2 * (XD * d * xi * (1 - xi) + XX * ((1 - xi) / xi + xi / (1 - xi)))
ours = (g**4 * Nc / (16 * sp.pi**4)) / xi**2 * (XX * (d * Phi * E + 2 * (xi / (1 - xi) + (1 - xi) / xi)) + 2 * d * xi * (1 - xi) * XD)
R['user_xi_form_minus_ours'] = str(sp.simplify(user - ours))

# 2. conversion to LO(k+/xi, k/xi): (g^4 Nc/16 pi^4)(1/k+)(1/xi^2) * (2 pi^2 k+/(g^2 xi)) = (alpha_s/2pi) Nc / xi^3
conv = (g**4 * Nc / (16 * sp.pi**4)) / kp / xi**2 * (2 * sp.pi**2 * kp / (g**2 * xi))
R['LO_conversion_minus_alpha_Nc_over_2pi_xi3'] = str(sp.simplify(conv.subs(g, sp.sqrt(4 * sp.pi * al)) - al / (2 * sp.pi) * Nc / xi**3))

# 3. convolution form {P_r (d Ebar-2)/2 XX + 2 Nc xi(1-xi)(XX - d XD)}/xi^3 vs. compact form (before the endpoint split)
Pr = 2 * Nc * Phi
conv_form = (Pr * (d * Eb - 2) / 2 * XX + 2 * Nc * xi * (1 - xi) * (XX - d * XD)) / xi**3
compact = Nc / xi * ((d * Eb - 2) * XX * (xi**-2 / (1 - xi) + (1 - 2 * xi) / xi**3) + d * (1 - xi) / xi * (Eb * XX - 2 * XD))
R['convolution_form_minus_compact'] = str(sp.simplify(conv_form - compact))

# 4. CDR: d = 2 - eps, 1/epsbar = -2/eps, Ebar = 1/epsbar - ell  =>  (d Ebar - 2)/2 = 1/epsbar - ell + O(eps)
w = ((2 - eps) * (-2 / eps - ell) - 2) / 2
R['CDR_weight_minus_(1/epsbar-ell)'] = str(sp.series(w - (-2 / eps - ell), eps, 0, 1).removeO())
# d_perp = 2: (d Ebar - 2)/2 = 1/epsbar - ell - 1
R['dperp2_weight_minus_(1/epsbar-ell-1)'] = str(sp.simplify(((2) * (-2 / eps - ell) - 2) / 2 - (-2 / eps - ell - 1)))

# 5. plus prescription with exact limits:
#    int_{x0}^{x1} Phi g = int_{x0}^1 Phi_+ g + g(1) log(1/(1-x1)) - int_{x1}^1 [Phi g - g(1)/(1-x)]
#    with int_{x0}^1 Phi_+ g = int_{x0}^1 [(x g - g(1))/(1-x) + ((1-x)/x + x(1-x)) g] + g(1) log(1-x0)
Phin = lambda x: x * (1 - x) + x / (1 - x) + (1 - x) / x
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(12):
    x0 = rng.uniform(0.01, 0.4)
    x1 = rng.uniform(0.6, 0.98)
    a = rng.uniform(0.2, 3.0)
    for part in (np.cos, np.sin):
        gf = lambda x: part(a / x) / x**2
        lhs = quad(lambda x: Phin(x) * gf(x), x0, x1, limit=400)[0]
        plus = quad(lambda x: (x * gf(x) - gf(1)) / (1 - x) + ((1 - x) / x + x * (1 - x)) * gf(x), x0, 1, limit=400)[0] + gf(1) * np.log(1 - x0)
        rem = quad(lambda x: Phin(x) * gf(x) - gf(1) / (1 - x), x1, 1, limit=400)[0]
        rhs = plus + gf(1) * np.log(1 / (1 - x1)) - rem
        worst = max(worst, abs(lhs - rhs) / max(abs(lhs), 1e-12))
R['plus_prescription_exact_limits_max_rel'] = worst

# 5b. R_Lambda/delta -> -2 Nc [g(1) + g'(1)] as delta -> 0 (Nc = 1 here)
a = 1.3
gf = lambda x: np.cos(a / x) / x**2
gp1 = a * np.sin(a) - 2 * np.cos(a)          # d/dx [cos(a/x)/x^2] at x = 1
lim = []
for dl in (1e-2, 1e-3, 1e-4):
    x1 = 1 / (1 + dl)
    rem = quad(lambda x: 2 * (Phin(x) * gf(x) - gf(1) / (1 - x)), x1, 1, limit=200)[0]
    lim.append(rem / dl)
R['R_Lambda_over_delta'] = lim
R['R_Lambda_over_delta_expected'] = -2 * (gf(1) + gp1)

# 6. momentum sum rule: int_0^1 xi P_gg = 0 fixes the delta(1-xi) coefficient to 11 Nc/6 (n_f = 0)
z = sp.symbols('z', positive=True)
reg = 2 * Nc * ((1 - z) / z + z * (1 - z))
plus_part = 2 * Nc * sp.integrate((z * z - 1) / (1 - z), (z, 0, 1))      # z * [1/(1-z)]_+ against f(z) = z
mom = sp.simplify(plus_part + sp.integrate(z * reg, (z, 0, 1)))
R['int_xi_Preal_plus'] = str(mom)
R['delta_coefficient_from_sum_rule'] = str(sp.simplify(-mom))

# 7. Jacobian: int_{a}^1 dxi/xi^3 [int_xi^1 dz/z P(z) D(xi/z)] F(k/xi) = int_a^1 dzeta/zeta^3 D(zeta) int_{a/zeta}^1 dz/z^3 P(z) F(k/(zeta z))
Pt = lambda zz: (1 - zz) / zz + zz * (1 - zz)
Dt = lambda zz: zz**2 * (1 - zz)**3 + 0.3 * zz
Ft = lambda q: np.exp(-0.2 * q) * q / (1 + q**2)
kk, a0 = 1.7, 0.15
lhs = quad(lambda x: Ft(kk / x) / x**3 * quad(lambda zz: Pt(zz) * Dt(x / zz) / zz, x, 1)[0], a0, 1)[0]
rhs = dblquad(lambda zz, ze: Dt(ze) / ze**3 * Pt(zz) / zz**3 * Ft(kk / (ze * zz)), a0, 1, lambda ze: a0 / ze, lambda ze: 1)[0]
R['jacobian_dz_over_z_with_dxi_over_xi3_rel'] = abs(lhs - rhs) / abs(lhs)

# 8. mu dependence: d/dlog mu^2 [(alpha_s/2pi)(1/epsbar)] with d alpha_s/dlog mu^2 = -(eps/2) alpha_s, 1/epsbar = -2/eps
R['dlogmu2_of_alpha_over_epsbar_over_alpha'] = str(sp.simplify((-eps / 2) * (-2 / eps)))

# 9. delta(1-xi) poles, units (alpha_s/2pi)(1/epsbar) Nc LO(k), exact in delta = Lambda/k+ and xi0 = k+/vee
dl, x0 = sp.symbols('delta xi0', positive=True)
rowIV = 2 * sp.log((1 + dl) / dl)                                 # 2 Nc log((k+ + Lambda)/Lambda), [0,1] plus cut at xi0
B_II = sp.Rational(11, 3) + (2 * sp.log(dl) - sp.log(1 - dl) - sp.log(1 + dl)) + sp.log(1 - x0) - 8 * dl + 2 * dl**2 - sp.Rational(4, 3) * dl**3
B_V = 2 * sp.log(1 / dl) + 2 * sp.log(1 - dl) - sp.Rational(11, 6) + 4 * dl - dl**2 + sp.Rational(2, 3) * dl**3
II_cc = -B_II                       # Row II: -(alpha_s Nc/4pi) Ebar B_II LO per ordering; x2 for c.c.
V_cc = -2 * B_V                     # Row V 2b UV: -(alpha_s Nc/2pi) Ebar B_V LO per ordering; x2 for c.c.
CT = sp.Rational(11, 6)             # charge renormalization, (alpha_s/2pi)(beta0/2)(1/epsbar) LO
tot = rowIV + II_cc + V_cc + CT
target = sp.Rational(11, 6) + 3 * sp.log((1 + dl) / (1 - dl)) - sp.log(1 - x0)
diff = sp.expand_log(sp.expand(tot - target), force=True)
R['delta_tally_minus_(11/6+3log((1+d)/(1-d))-log(1-xi0))'] = str(sp.simplify(diff))
R['tally_without_CT_at_delta0_xi00'] = str(sp.limit(sp.limit(sp.expand_log(rowIV + II_cc + V_cc, force=True), dl, 0), x0, 0))
R['leftover_series'] = str(sp.series(3 * sp.log((1 + dl) / (1 - dl)), dl, 0, 4).removeO())
# per-ordering beta0 and soft-log content
R['II_plus_V_per_ordering_11/3_part'] = str(sp.Rational(11, 3) + 2 * (-sp.Rational(11, 6)))
# prefactor conversions to (alpha_s/2pi) LO
cII = (g**4 * Nc / (32 * sp.pi**4)) / kp * (2 * sp.pi**2 * kp / g**2)
cV = (g**4 * Nc / (16 * sp.pi**4)) / kp * (2 * sp.pi**2 * kp / g**2)
R['RowII_prefactor_over_(alpha_s/2pi)'] = str(sp.simplify(cII.subs(g, sp.sqrt(4 * sp.pi * al)) / (al / (2 * sp.pi))))
R['RowV2b_prefactor_over_(alpha_s/2pi)'] = str(sp.simplify(cV.subs(g, sp.sqrt(4 * sp.pi * al)) / (al / (2 * sp.pi))))

# 10. colour signs: Row II has [U(x')-U(w')][U(w)-U(x)] = -O_LO; Row V 2b has [U(x')-U(w')][U(y)-U(z)] = +O_LO(y->x, z->y, w'->y')
rng = np.random.default_rng(3)
def rot(n=8):
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    return q
Ux, Uxp, Uw, Uwp, Uy, Uz = (rot() for _ in range(6))
r1, r2 = rng.normal(size=8), rng.normal(size=8)
OLO = ((Uwp - Uxp) @ r2) @ ((Uw - Ux) @ r1)
OII = ((Uxp - Uwp) @ r2) @ ((Uw - Ux) @ r1)
OLOV = ((Uwp - Uxp) @ r2) @ ((Uz - Uy) @ r1)
OV = ((Uxp - Uwp) @ r2) @ ((Uy - Uz) @ r1)
R['colour_RowII_over_LO'] = OII / OLO
R['colour_RowV2b_over_LO'] = OV / OLOV

json.dump(R, open(os.path.join(HERE, 'frag_checks.json'), 'w'), indent=1)
for k, v in R.items():
    print(k, ':', v)
