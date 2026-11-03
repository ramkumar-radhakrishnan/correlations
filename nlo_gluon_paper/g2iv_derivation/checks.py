r"""Checks of the derivation of Group II Row IV from C_1 (user's steps).  python3 checks.py"""
import itertools
import numpy as np
import sympy as sp
from scipy import integrate, special

p, k, d, xi = sp.symbols('p k d xi', positive=True)
P = p + k
out = {}

# 1. sum_{ij} T_{kmij} T_{k'm'ij} x pk/P^4 x (P/k)^4 (k/P)^2   vs the bracket of step 3 (with 1/k+)
A, B = P / p, P / k
cKM_KM = d - 2 * A - 2 * B            # coefficient of delta_km delta_k'm'
cMM_KK = A**2 + B**2                  # delta_mm' delta_kk'
cMK_KM = 2 * A * B                    # delta_mk' delta_m'k
fac = p * k / P**4 * (P / k)**4 * (k / P)**2
user = {'KM_KM': d * p / P**2 / k - 2 / k**2, 'MM_KK': (p / k**2 + 1 / p) / k, 'MK_KM': 2 / k**2}
mine = {'KM_KM': fac * cKM_KM, 'MM_KK': fac * cMM_KK, 'MK_KM': fac * cMK_KM}
out['step3 bracket'] = {key: sp.simplify(mine[key] - user[key]) == 0 for key in user}
# explicit index check of the contraction in d=2 and d=3
for dim in (2, 3):
    Av, Bv = 1.7, 2.9
    I = np.eye(dim)
    T = (np.einsum('km,ij->kmij', I, I) - Av * np.einsum('jm,ik->kmij', I, I) - Bv * np.einsum('im,kj->kmij', I, I))
    lhs = np.einsum('kmij,lnij->kmln', T, T)            # k m k' m'
    rhs = ((dim - 2 * Av - 2 * Bv) * np.einsum('km,ln->kmln', I, I) + (Av**2 + Bv**2) * np.einsum('mn,kl->kmln', I, I)
           + 2 * Av * Bv * np.einsum('ml,nk->kmln', I, I))
    out['contraction d=%d' % dim] = float(abs(lhs - rhs).max())

# 2. xi = k/(p+k): brackets and measure
sub = {p: k * (1 - xi) / xi}
b1 = sp.simplify((p / P**2 + p / k**2 + 1 / p).subs(sub) * k - (xi * (1 - xi) + (1 - xi) / xi + xi / (1 - xi)))
b2 = sp.simplify((p / P**2).subs(sub) * k - xi * (1 - xi))
jac = sp.simplify(sp.diff(k * (1 - xi) / xi, xi) + k / xi**2)
out['step 7: brackets, dp = -k dxi/xi^2'] = [b1 == 0, b2 == 0, jac == 0]

# 3. z integral: traceless part -pi (r^m r^m'/r^2 - delta/2), trace -pi[2/eps + ell] (2d check of the r-dependence)
#    cos(2phi) component: int d^2k e^{ik.r} cos2phi_k / k^2 = -2 pi int_0^inf J2(kr) dk/k = -pi
val = -2 * np.pi * integrate.quad(lambda t: special.jv(2, t) / t, 0, np.inf, limit=2000)[0]
out['step 5: traceless coefficient (expect -pi)'] = [val, -np.pi]

# 4. plus prescription with the exact upper limit k/(k+Lambda)
K, V, L = 1.0, 9.0, 0.02
a, x1 = K / V, K / (K + L)
G = lambda x: np.cos(2.3 / x) * np.exp(-0.3 / x) + 0.4 * x
lhs = integrate.quad(lambda x: G(x) / (x * (1 - x)), a, x1, limit=400)[0]                    # your step 6 term xi/(1-xi)
user_plus = integrate.quad(lambda x: (G(x) - G(1)) / (x * (1 - x)), a, x1, limit=400)[0]      # your prescription on [a, x1]
blue = G(1) * (np.log(K / L) + np.log(V / K - 1))
ext = integrate.quad(lambda x: (G(x) - G(1)) / (x * (1 - x)), x1, 1, limit=400)[0]            # extension x1 -> 1
out['step 9: exact on [k/vee, k/(k+L)]'] = abs(lhs - (user_plus + blue))
out['step 9: upper limit k/(k+L) -> 1 changes by (O(L/k))'] = ext
print(out)
