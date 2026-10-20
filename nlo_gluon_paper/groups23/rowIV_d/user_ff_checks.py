"""Checks for RowIV_FF_userform.tex: fragmentation function and DGLAP from the user's form of Row IV
(subtraction exp(-ik.(y'-y))/(xi(1-xi)), delta(1-xi) coefficient log(vee/Lambda)+log(1-k+/vee)).
a = k.(y'-y); the phase at momentum fraction xi is exp(-i a/xi). Results in user_ff_checks.json."""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
R = {}


def cq(f, a, b):
    return quad(lambda x: f(x).real, a, b, limit=500)[0] + 1j * quad(lambda x: f(x).imag, a, b, limit=500)[0]


Phi = lambda x: x * (1 - x) + x / (1 - x) + (1 - x) / x
rest = lambda x: x * (1 - x) + (1 - x) / x
ph = lambda a, x: np.exp(-1j * a / x)

# 1. the user's form vs. the exact xi integral of Row IV (coefficient of X.X' times the weight)
#    upper limit xi1: exact; upper limit 1: off by Rhat = -i (Lambda/k+) a e^{-ia} + O(Lambda^2)
rows = []
for kp, vee, a in ((1.0, 20., 1.3), (2.0, 5., -0.7), (0.7, 3., 2.2)):
    x0 = kp / vee
    U = lambda x: (ph(a, x) - ph(a, 1)) / (x * (1 - x)) + rest(x) * ph(a, x) / x**2
    for Lam in (0.3 * kp, 1e-2 * kp, 1e-4 * kp):
        x1 = kp / (kp + Lam)
        direct = cq(lambda x: Phi(x) * ph(a, x) / x**2, x0, x1)
        dterm = (np.log(vee / Lam) + np.log(1 - kp / vee)) * ph(a, 1)
        exact = cq(U, x0, x1) + dterm
        up1 = cq(U, x0, 1) + dterm
        Rhat = cq(U, x1, 1)
        rows.append(dict(kp=kp, vee=vee, a=a, Lam_over_kp=Lam / kp, exact_rel=abs(exact - direct) / abs(direct),
                         upper1_minus_direct=abs(up1 - direct),
                         Rhat_over_delta=[(Rhat / (Lam / kp)).real, (Rhat / (Lam / kp)).imag],
                         expected=[(-1j * a * ph(a, 1)).real, (-1j * a * ph(a, 1)).imag]))
R['user_form_vs_exact'] = rows
R['user_form_exact_max_rel'] = max(r['exact_rel'] for r in rows)

# 2. dxi/(xi(1-xi)) integral: int_{k+/vee}^{k+/(k++Lambda)} = log(vee/Lambda) + log(1-k+/vee), exactly
kp, vee, Lam = sp.symbols('k_p vee Lambda', positive=True)
xs = sp.symbols('xi', positive=True)
F = sp.log(xs / (1 - xs))
val = F.subs(xs, kp / (kp + Lam)) - F.subs(xs, kp / vee)
R['int_dxi_over_xi(1-xi)_minus_user_bracket'] = str(sp.simplify(sp.expand_log(val - (sp.log(vee / Lam) + sp.log(1 - kp / vee)), force=True)))

# 3. conversion to the standard plus distribution (exact):
#    U_+ = int_{x0}^1 {[e(x)-e(1)]/(x(1-x)) + rest e/x^2}
#    S_+ = int_{x0}^1 {[e(x)/x - e(1)]/(1-x) + rest e/x^2} + log(1-x0) e(1)      (standard [0,1] plus, cut at x0)
#    U_+ = S_+ - [log(vee/k+) + log(1-k+/vee)] e(1)
worst = 0.0
rng = np.random.default_rng(11)
for _ in range(20):
    x0 = rng.uniform(0.01, 0.6)
    a = rng.uniform(-3, 3)
    Up = cq(lambda x: (ph(a, x) - ph(a, 1)) / (x * (1 - x)) + rest(x) * ph(a, x) / x**2, x0, 1)
    Sp = cq(lambda x: (ph(a, x) / x - ph(a, 1)) / (1 - x) + rest(x) * ph(a, x) / x**2, x0, 1) + np.log(1 - x0) * ph(a, 1)
    rhs = Sp - (np.log(1 / x0) + np.log(1 - x0)) * ph(a, 1)
    worst = max(worst, abs(Up - rhs) / abs(Up))
R['user_plus_to_standard_plus_max_rel'] = worst

# 4. weight W = -d/eps - d gamma - (d/2) log(D^2 mu^2/4) - 1 with d = 2 - eps  ->  -2/eps - 2 gamma - log(D^2 mu^2/4) + O(eps)
eps, gam, lg = sp.symbols('epsilon gamma_E lg')
W = -(2 - eps) / eps - (2 - eps) * gam - (2 - eps) / 2 * lg - 1
R['W_CDR_minus_(-2/eps-2gamma-log)'] = str(sp.series(W - (-2 / eps - 2 * gam - lg), eps, 0, 1).removeO())
dd = sp.symbols('d_perp')
Wd = -dd / eps - dd * gam - dd / 2 * lg - 1
R['W_dperp2_minus_(-2/eps-2gamma-log-1)'] = str(sp.simplify(Wd.subs(dd, 2) - (-2 / eps - 2 * gam - lg - 1)))

# 5. prefactor: g^2 Nc/(4 pi^2) = (alpha_s/2pi) 2 Nc
g, Nc, al = sp.symbols('g N_c alpha_s', positive=True)
R['g2Nc_4pi2_minus_alpha_2pi_2Nc'] = str(sp.simplify((g**2 * Nc / (4 * sp.pi**2) - al / (2 * sp.pi) * 2 * Nc).subs(al, g**2 / (4 * sp.pi))))

# 6. B = A with the standard plus distribution and cut supports, phase explicit:
#    A = int_{x0}^1 dxi/xi^2 e^{-ia/xi} [ int_xi^1 dz/z [1/(1-z)]_+ D(xi/z) ]
#    B = int_{x0}^1 dzeta/zeta^2 D(zeta) [ int_{x0/zeta}^1 dz/z^2 [1/(1-z)]_+ e^{-ia/(zeta z)} ]
Dt = lambda t: t**2 * (1 - t)**3 + 0.3 * t
worst = 0.0
for a, x0 in ((1.7, 0.15), (-0.6, 0.05), (2.4, 0.3)):
    def inner_A(x):
        return quad(lambda z: (Dt(x / z) / z - Dt(x)) / (1 - z), x, 1, limit=200)[0] + Dt(x) * np.log(1 - x)
    A = cq(lambda x: ph(a, x) / x**2 * inner_A(x), x0, 1)

    def inner_B(ze):
        c = x0 / ze
        h = lambda z: ph(a, ze * z) / z**2
        return cq(lambda z: (h(z) - h(1)) / (1 - z), c, 1) + h(1) * np.log(1 - c)
    B = cq(lambda ze: Dt(ze) / ze**2 * inner_B(ze), x0, 1)
    worst = max(worst, abs(A - B) / abs(A))
R['B_equals_A_with_plus_and_cuts_max_rel'] = worst

# 7. delta(1-xi) tally at Lambda -> 0, in the user's logs (units (alpha_s/2pi)(-2/eps) Nc e^{-ik.(y'-y)}):
dY, Lvk, x0s = sp.symbols('deltaY logvk xi0')           # deltaY = log(vee/Lambda), logvk = log(vee/k+)
rowIV_std = 2 * (dY - Lvk)                             # 2 log(k+/Lambda), after conversion to the standard plus
rowIV_user = 2 * (dY + sp.log(1 - x0s))                # your form, before conversion
II_cc = -(sp.Rational(11, 3) - 2 * (dY - Lvk) + sp.log(1 - x0s))
V_cc = -2 * (2 * dY - 2 * Lvk - sp.Rational(11, 6))
CT = sp.Rational(11, 6)
R['tally_std'] = str(sp.simplify(rowIV_std + II_cc + V_cc + CT))
R['conversion_shift'] = str(sp.simplify(rowIV_user - rowIV_std))

# 8. mu dependence of (alpha_s/2pi)(-2/eps) with d alpha_s/dlog mu^2 = -(eps/2) alpha_s
R['dlogmu2_of_alpha_times_minus2_over_eps_over_alpha'] = str(sp.simplify((-eps / 2) * (-2 / eps)))

# 9. momentum sum rule
z = sp.symbols('z', positive=True)
mom = 2 * Nc * sp.integrate((z * z - 1) / (1 - z), (z, 0, 1)) + sp.integrate(z * 2 * Nc * ((1 - z) / z + z * (1 - z)), (z, 0, 1))
R['delta_coefficient_from_sum_rule'] = str(sp.simplify(-mom))

json.dump(R, open(os.path.join(HERE, 'user_ff_checks.json'), 'w'), indent=1, default=float)
for k, v in R.items():
    if k != 'user_form_vs_exact':
        print(k, ':', v)
for r in R['user_form_vs_exact']:
    print(r)
