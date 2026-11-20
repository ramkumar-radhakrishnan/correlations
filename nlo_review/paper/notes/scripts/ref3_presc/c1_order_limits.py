# Ref.3 / error item "leftover pole ~ ell_vee/2 removed by vee -> oo before eps -> 0".
# Origin of ell_vee = ln(1-k+/vee) in the pole factor (11/3 - 2L + ell_vee) of the one-rho part of A^(3)  (app_wavefunctions,
# "Sum of the one-rho terms"): the xi-integral of (2-xi)/(xi(1-xi)) = 2/xi + 1/(1-xi) over k+/vee < xi < k+/(k++Lambda),
# where xi k+ ... the parent of the observed gluon carries k+/xi.  The 2/xi part gives 2 ln(vee/k+) (BFKL, cancels in the
# pole coefficient against the ln vee of the other diagrams); the 1/(1-xi) part, which is NOT enhanced at xi -> 0, gives
# L + ell_vee, and ell_vee = -int_0^{k+/vee} dxi/(1-xi): the contribution that is missing because parents faster than vee
# (valence modes) are excluded.  It is a power correction, O(k+/vee).
# With a smooth cut-off f(parent/vee) (f(0)=1, f -> 0 at oo) instead of theta(vee - parent):
#   2/xi part   -> 2[ln(vee/k+) + c_f] + O(a),          c_f = int du/u [f(u) - theta(1-u)]      (universal shift of ln vee)
#   1/(1-xi)    -> L + ell_f,   ell_f = -int_0^1 dxi [1 - f(a/xi)]/(1-xi) = -kappa_f a + o(a),  a = k+/vee,
#                  kappa_f = int_0^oo dv [1-f(v)]/v^2  (sharp: 1, exp(-v^2): sqrt(pi), 1/(1+v^2): pi/2; f with f'(0)!=0: O(a L))
# i.e. the coefficient of the leftover pole is cut-off-shape dependent: it is not a physical quantity of the effective theory.
import numpy as np
from scipy.integrate import quad

fs = {"sharp": lambda u: 1.0 if u < 1 else 0.0, "exp(-u^2)": lambda u: np.exp(-u * u),
      "1/(1+u^2)": lambda u: 1.0 / (1 + u * u)}

def ell_f(a, f):
    # -int_0^1 dxi [1 - f(a/xi)]/(1-xi); integrand is regular at xi->1 since f(a)->1 for a->0 (up to O(a) handled numerically)
    g = lambda xi: (1.0 - f(a / xi)) / (1.0 - xi) if xi < 1 else 0.0
    pts = [a, 10 * a] if a < 0.1 else None
    return -quad(g, 0, 1, points=pts, limit=400)[0]

def kappa(f):
    return quad(lambda v: (1 - f(v)) / v ** 2, 0, 1, limit=400)[0] + quad(lambda v: (1 - f(v)) / v ** 2, 1, np.inf, limit=400)[0]

print("ell_vee for the sharp cut-off: ln(1-a) vs numerical:")
for a in (1e-1, 1e-2, 1e-3):
    print(f"   a=k+/vee={a:.0e}:  ln(1-a)={np.log(1-a): .6e}   -int_0^a dxi/(1-xi) = {ell_f(a, fs['sharp']): .6e}")
print("\nsmooth cut-offs: ell_f(a)/a  (-> -kappa_f) ;  kappa_f")
for name, f in fs.items():
    try:
        kf = kappa(f)
    except Exception:
        kf = float('nan')
    vals = [ell_f(a, f) / a for a in (1e-2, 1e-3, 1e-4)]
    print(f"   {name:10s}: ell_f/a at a=1e-2,1e-3,1e-4: {vals[0]: .5f} {vals[1]: .5f} {vals[2]: .5f}   kappa_f = {kf: .5f}")
print("   (for f with f'(0) != 0, e.g. exp(-u), also the observed gluon itself gets the weight f(k+/vee) = 1 - O(k+/vee),")
print("    which multiplies L as well: ell_f ~ -a [ln(k+/Lambda) + O(1)], again O(k+/vee) times logarithms; not shown.)")

# The leftover pole (single-eps bookkeeping):  S(eps, a) = -(ell/2)(2/eps).  Order of limits:
print("\norder of limits for S = -(ell_vee/2)(2/eps) ~ (a/2)(2/eps):")
for eps in (1e-1, 1e-2, 1e-3):
    row = [-(np.log(1 - a) / 2) * (2 / eps) for a in (1e-2, 1e-4, 1e-6)]
    print(f"   eps={eps:.0e}: a=1e-2,1e-4,1e-6 -> {row[0]: .3e} {row[1]: .3e} {row[2]: .3e}")
print("   vee->oo at fixed eps: 0;  eps->0 at fixed vee: diverges.  At leading power in a = k+/vee the term is absent.")
