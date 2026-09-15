"""Does this row contribute to DGLAP?  Compare its longitudinal structure with P_gg."""
import sympy as sp, numpy as np
xi = sp.symbols('xi')                      # no positivity assumption: we want ALL poles

print("A. LONGITUDINAL STRUCTURE OF THIS ROW   (xi = p+/k+, in units of 1/k+)")
for nm, t in [("1/(1+xi)", 1/(1+xi)), ("-1/xi", -1/xi), ("-1", sp.Integer(-1))]:
    p = sp.solve(sp.denom(sp.together(t)), xi)
    inside = [q for q in p if q.is_real and 0 <= q <= 1]
    print("   %-10s poles at xi = %-10s inside [0,1]: %s" % (nm, str(p), str(inside)))
print("   -> the only singularity in the physical range is xi -> 0, the SOFT endpoint.")
print("      No xi -> 1 pole: no collinear endpoint.")
print()

print("B. WHY -- the vertex's collinear denominator is switched off by the routing")
kp, pp = sp.symbols('kplus pplus', positive=True)
print("   C's third term carries k+/(k+ - p+), a genuine 1/(1-xi) collinear pole.")
print(f"   But C is used with first argument K+ = p+ + k+, so it becomes")
print(f"      K+/(K+ - p+) = {sp.simplify((pp+kp)/kp)}  -- finite for every p+ > 0.")
print()

print("C. COMPARE WITH THE SPLITTING FUNCTION")
Cuv = xi*(1-xi) + xi/(1-xi) + (1-xi)/xi              # = P_gg/(2 N_c)
p = sp.solve(sp.denom(sp.together(Cuv)), xi)
print(f"   C_UV(xi) = P_gg/(2N_c) = xi(1-xi) + xi/(1-xi) + (1-xi)/xi")
print(f"      poles at xi = {p}   <- BOTH endpoints, 0 AND 1")
print(f"   this row      = -1/xi          poles at xi = [0]   <- soft endpoint only")
print("   the xi(1-xi) and xi/(1-xi) pieces of P_gg are simply absent here.")
print()

print("D. THE EARLIER (UV-projected) ROW, FOR CONTRAST")
print("   there the light-cone vertex supplied a SECOND 1/(x-z)^2, giving a true 1/r^2,")
print("   a true transverse UV log, and its coefficient WAS C_UV(xi) = P_gg/(2N_c).")
f = sp.lambdify(xi, Cuv, 'numpy')
from scipy.integrate import quad
for lam in (1e-6, 1e-9):
    I = quad(f, lam, 1-lam, limit=400)[0]
    lk = np.log(1/lam)
    print(f"   lam={lam:.0e}:  int C_UV dxi = {I:16.9f}    2 l_k - 11/6 = {2*lk-11/6:16.9f}")
print("   and (2 l_k - 11/6) x 2N_c = 4 N_c l_k - 11 N_c/3 = 4 N_c l_k - b  at N_f = 0:")
print("   the constant IS the beta-function coefficient -> coupling renormalisation.")
print("   That is the DGLAP / running-coupling structure.  This row has none of it.")
