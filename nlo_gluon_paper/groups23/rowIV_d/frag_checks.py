"""Checks for RowIV_FF_DGLAP.tex: Row IV in fragmentation form with the phase exp(-i k.(y'-y)/xi) kept explicit,
the plus prescription with exact limits, the MS-bar counterterm, DGLAP, and the delta(1-xi) poles of Group I
(Row II and Row V 2b, with c.c.). Results in frag_checks.json.

Throughout, a = k.(y'-y), and the phase at momentum fraction xi is exp(-i a/xi)."""
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
Pr = 2 * Nc * Phi


def cquad(f, a, b, **kw):
    """complex quad"""
    return quad(lambda x: f(x).real, a, b, limit=400, **kw)[0] + 1j * quad(lambda x: f(x).imag, a, b, limit=400, **kw)[0]


# 1. the user's xi form vs. RowIV_dperp.pdf Eq. (1.5) (E form), both on int_{xi1}^{xi0}, same phase and 1/xi^2
user = (g**4 * d / (4 * sp.pi**3)) * Nc / (2 * sp.pi) / xi**2 * Phi * (E / 2) * XX \
    + (g**4 / (4 * sp.pi**3)) * Nc / (2 * sp.pi) / xi**2 * (XD * d * xi * (1 - xi) + XX * ((1 - xi) / xi + xi / (1 - xi)))
ours = (g**4 * Nc / (16 * sp.pi**4)) / xi**2 * (XX * (d * Phi * E + 2 * (xi / (1 - xi) + (1 - xi) / xi)) + 2 * d * xi * (1 - xi) * XD)
R['user_xi_form_minus_ours'] = str(sp.simplify(user - ours))

# 2. prefactor: g^4/(16 pi^4) = (alpha_s/2pi) * g^2/(2 pi^2)   (nothing depends on xi)
R['prefactor_g4_16pi4_minus_alpha_2pi_g2_2pi2'] = str(sp.simplify(
    (g**4 / (16 * sp.pi**4) - al / (2 * sp.pi) * g**2 / (2 * sp.pi**2)).subs(al, g**2 / (4 * sp.pi))))

# 3. braces: Nc{XX[d Phi Ebar - 2(xi/(1-xi)+(1-xi)/xi)] - 2 d xi(1-xi) XD}
#           = P_r (d Ebar-2)/2 XX + 2 Nc xi(1-xi)(XX - d XD)       (same 1/xi^2 and same phase on both sides)
lhs = Nc * (XX * (d * Phi * Eb - 2 * (xi / (1 - xi) + (1 - xi) / xi)) - 2 * d * xi * (1 - xi) * XD)
rhs = Pr * (d * Eb - 2) / 2 * XX + 2 * Nc * xi * (1 - xi) * (XX - d * XD)
R['braces_explicit_form_minus_RowIV_dperp'] = str(sp.simplify(lhs - rhs))
# and against the compact form of RowIV_compact.pdf (which carries xi^-2 and the phase inside its brackets)
compact = Nc * ((d * Eb - 2) * XX * (1 / (1 - xi) + (1 - 2 * xi) / xi) + d * (1 - xi) * xi * (Eb * XX - 2 * XD))
R['braces_explicit_form_minus_compact'] = str(sp.simplify(rhs - compact))

# 4. CDR: d = 2 - eps, 1/epsbar = -2/eps, Ebar = 1/epsbar - ell  =>  (d Ebar - 2)/2 = 1/epsbar - ell + O(eps)
w = ((2 - eps) * (-2 / eps - ell) - 2) / 2
R['CDR_weight_minus_(1/epsbar-ell)'] = str(sp.series(w - (-2 / eps - ell), eps, 0, 1).removeO())
R['dperp2_weight_minus_(1/epsbar-ell-1)'] = str(sp.simplify(((2) * (-2 / eps - ell) - 2) / 2 - (-2 / eps - ell - 1)))

# 5. plus prescription with exact limits, phase explicit (Nc = 1):
#    int_{x0}^{x1} P_r xi^-2 e^{-ia/xi}
#      = 2 int_{x0}^1 [ (xi^-1 e^{-ia/xi} - e^{-ia})/(1-xi) + ((1-xi)/xi^3 + (1-xi)/xi) e^{-ia/xi} ] + 2 log(1-x0) e^{-ia}
#        + 2 log(1/(1-x1)) e^{-ia} - int_{x1}^1 [ P_r xi^-2 e^{-ia/xi} - 2 e^{-ia}/(1-xi) ]
Prn = lambda x: 2 * (x * (1 - x) + x / (1 - x) + (1 - x) / x)
ph = lambda a, x: np.exp(-1j * a / x)
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(24):
    x0 = rng.uniform(0.01, 0.4)
    x1 = rng.uniform(0.6, 0.98)
    a = rng.uniform(-3.0, 3.0)
    lhs = cquad(lambda x: Prn(x) * ph(a, x) / x**2, x0, x1)
    plus = 2 * cquad(lambda x: (ph(a, x) / x - ph(a, 1)) / (1 - x) + ((1 - x) / x**3 + (1 - x) / x) * ph(a, x), x0, 1) \
        + 2 * np.log(1 - x0) * ph(a, 1)
    rem = cquad(lambda x: Prn(x) * ph(a, x) / x**2 - 2 * ph(a, 1) / (1 - x), x1, 1)
    rhs = plus + 2 * np.log(1 / (1 - x1)) * ph(a, 1) - rem
    worst = max(worst, abs(lhs - rhs) / abs(lhs))
R['plus_prescription_explicit_phase_max_rel'] = worst

# 5b. R_Lambda/delta -> 2 Nc (1 - i a) e^{-ia} as delta -> 0 (Nc = 1)
a = 1.3
lim = []
for dl in (1e-2, 1e-3, 1e-4):
    x1 = 1 / (1 + dl)
    rem = cquad(lambda x: Prn(x) * ph(a, x) / x**2 - 2 * ph(a, 1) / (1 - x), x1, 1)
    lim.append(rem / dl)
R['R_Lambda_over_delta'] = [[v.real, v.imag] for v in lim]
exp_ = 2 * (1 - 1j * a) * np.exp(-1j * a)
R['R_Lambda_over_delta_expected'] = [exp_.real, exp_.imag]
R['R_Lambda_a'] = a

# 5c. the full Row IV in plus form (Eq. 3.6) vs. the direct xi integral (Eq. 2.2), random d, Ebar, XX, XD, a (Nc = 1)
worst = 0.0
for _ in range(12):
    x0 = rng.uniform(0.01, 0.4)
    dl = rng.uniform(0.02, 0.8)
    x1 = 1 / (1 + dl)
    a = rng.uniform(-3.0, 3.0)
    dd, Ebv, xx, xd = rng.uniform(1.5, 3.5), rng.normal(), rng.normal(), rng.normal()
    direct = cquad(lambda x: ph(a, x) / x**2 * (Prn(x) * (dd * Ebv - 2) / 2 * xx + 2 * x * (1 - x) * (xx - dd * xd)), x0, x1)
    plus = 2 * cquad(lambda x: (ph(a, x) / x - ph(a, 1)) / (1 - x) + ((1 - x) / x**3 + (1 - x) / x) * ph(a, x), x0, 1) \
        + 2 * np.log(1 - x0) * ph(a, 1)
    rem = cquad(lambda x: Prn(x) * ph(a, x) / x**2 - 2 * ph(a, 1) / (1 - x), x1, 1)
    form = (dd * Ebv - 2) / 2 * xx * (plus + 2 * np.log((1 + dl) / dl) * ph(a, 1) - rem) \
        + 2 * cquad(lambda x: (1 - x) / x * ph(a, x), x0, x1) * (xx - dd * xd)
    worst = max(worst, abs(direct - form) / abs(direct))
R['full_RowIV_plus_form_vs_direct_max_rel'] = worst

# 6. momentum sum rule: int_0^1 xi P_gg = 0 fixes the delta(1-xi) coefficient to 11 Nc/6 (n_f = 0)
z = sp.symbols('z', positive=True)
reg = 2 * Nc * ((1 - z) / z + z * (1 - z))
plus_part = 2 * Nc * sp.integrate((z * z - 1) / (1 - z), (z, 0, 1))      # z * [1/(1-z)]_+ against f(z) = z
mom = sp.simplify(plus_part + sp.integrate(z * reg, (z, 0, 1)))
R['int_xi_Preal_plus'] = str(mom)
R['delta_coefficient_from_sum_rule'] = str(sp.simplify(-mom))

# 7. composition with the phase explicit (the 1/k+ of the parton spectrum at k+/zeta gives zeta/k+):
#    int_{x0}^1 dzeta/zeta^3 D(zeta) * zeta * int_{x0/zeta}^1 dz/z^2 P(z) e^{-ia/(zeta z)}
#      = int_{x0}^1 dxi/xi^2 e^{-ia/xi} int_xi^1 dz/z P(z) D(xi/z)
Pt = lambda zz: (1 - zz) / zz + zz * (1 - zz)
Dt = lambda zz: zz**2 * (1 - zz)**3 + 0.3 * zz
worst = 0.0
for a, a0 in ((1.7, 0.15), (-0.6, 0.05), (3.1, 0.3)):
    out = []
    for part in (np.cos, np.sin):
        lhs = dblquad(lambda zz, ze: Dt(ze) / ze**2 * Pt(zz) / zz**2 * part(a / (ze * zz)), a0, 1, lambda ze: a0 / ze, lambda ze: 1)[0]
        rhs = quad(lambda x: part(a / x) / x**2 * quad(lambda zz: Pt(zz) * Dt(x / zz) / zz, x, 1)[0], a0, 1)[0]
        out.append((lhs, rhs))
    l = out[0][0] - 1j * out[1][0]
    r = out[0][1] - 1j * out[1][1]
    worst = max(worst, abs(l - r) / abs(r))
R['composition_explicit_phase_max_rel'] = worst

# 8. mu dependence: d/dlog mu^2 [(alpha_s/2pi)(1/epsbar)] with d alpha_s/dlog mu^2 = -(eps/2) alpha_s, 1/epsbar = -2/eps
R['dlogmu2_of_alpha_over_epsbar_over_alpha'] = str(sp.simplify((-eps / 2) * (-2 / eps)))

# 9. delta(1-xi) poles (all with the phase exp(-i k.(y'-y)), no 1/xi), units
#    (alpha_s/2pi)(1/epsbar) Nc (1/(2pi)^3) (g^2/(2 pi^2 k+)) int e^{-ik.(y'-y)} X.X' O; exact in delta and xi0
dl, x0 = sp.symbols('delta xi0', positive=True)
rowIV = 2 * sp.log((1 + dl) / dl)                                 # 2 Nc log((k+ + Lambda)/Lambda), [0,1] plus cut at xi0
B_II = sp.Rational(11, 3) + (2 * sp.log(dl) - sp.log(1 - dl) - sp.log(1 + dl)) + sp.log(1 - x0) - 8 * dl + 2 * dl**2 - sp.Rational(4, 3) * dl**3
B_V = 2 * sp.log(1 / dl) + 2 * sp.log(1 - dl) - sp.Rational(11, 6) + 4 * dl - dl**2 + sp.Rational(2, 3) * dl**3
II_cc = -B_II                       # Row II: -(alpha_s/2pi)(Nc/2) Ebar B_II per ordering; x2 for c.c.
V_cc = -2 * B_V                     # Row V 2b UV: -(alpha_s/2pi) Nc Ebar B_V per ordering; x2 for c.c.
CT = sp.Rational(11, 6)             # charge renormalization, (alpha_s/2pi)(beta0/2)(1/epsbar)
tot = rowIV + II_cc + V_cc + CT
target = sp.Rational(11, 6) + 3 * sp.log((1 + dl) / (1 - dl)) - sp.log(1 - x0)
diff = sp.expand_log(sp.expand(tot - target), force=True)
R['delta_tally_minus_(11/6+3log((1+d)/(1-d))-log(1-xi0))'] = str(sp.simplify(diff))
R['tally_without_CT_at_delta0_xi00'] = str(sp.limit(sp.limit(sp.expand_log(rowIV + II_cc + V_cc, force=True), dl, 0), x0, 0))
R['leftover_series'] = str(sp.series(3 * sp.log((1 + dl) / (1 - dl)), dl, 0, 4).removeO())
R['II_plus_V_per_ordering_11/3_part'] = str(sp.Rational(11, 3) + 2 * (-sp.Rational(11, 6)))
cII = (g**4 * Nc / (32 * sp.pi**4)) / kp
cV = (g**4 * Nc / (16 * sp.pi**4)) / kp
unit = al / (2 * sp.pi) * g**2 / (2 * sp.pi**2 * kp)
R['RowII_prefactor_over_unit'] = str(sp.simplify((cII / unit).subs(al, g**2 / (4 * sp.pi))))
R['RowV2b_prefactor_over_unit'] = str(sp.simplify((cV / unit).subs(al, g**2 / (4 * sp.pi))))

# 10. colour signs: Row II has [U(x')-U(w')][U(w)-U(x)] = -O; Row V 2b has [U(x')-U(w')][U(y)-U(z)] = +O (renamed)
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
R['colour_RowII_over_O'] = OII / OLO
R['colour_RowV2b_over_O'] = OV / OLOV

json.dump(R, open(os.path.join(HERE, 'frag_checks.json'), 'w'), indent=1)
for k, v in R.items():
    print(k, ':', v)
