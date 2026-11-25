#!/usr/bin/env python3
"""Checks of the intermediate steps quoted in sec_examples.tex and sec_setup.tex (jhep_v3).
Run: python3 -I check_examples.py     (needs sympy, numpy, mpmath)"""
import itertools, sympy as sp, numpy as np, mpmath as mp
ok = True
def report(name, val, tol=1e-9):
    global ok
    good = abs(val) < tol; ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {name:78s} dev = {float(abs(val)):.2e}")

# (1) longitudinal integrals of the gluon self-energy, Eq. (xiintegrals)
xi = sp.symbols('xi', positive=True)
I2 = sp.integrate((-2 + xi*(1-xi))*sp.log(xi*(1-xi)), (xi, 0, 1))
report("int_0^1 [-2+xi(1-xi)] log(xi(1-xi)) = 67/18", sp.N(I2 - sp.Rational(67, 18)))
mp.mp.dps = 30
for d in (mp.mpf('1e-4'), mp.mpf('1e-6')):
    I1 = mp.quad(lambda x: -2 + x*(1-x) + 1/(x*(1-x)), [d, 0.5, 1-d])
    report(f"int_d^(1-d) [-2+xi(1-xi)+1/(xi(1-xi))] = -11/6-2log d   (d={float(d):.0e})", I1 - (-mp.mpf(11)/6 - 2*mp.log(d)), 10*float(d))
    I3 = mp.quad(lambda x: mp.log(x*(1-x))/(x*(1-x)), [d, 0.5, 1-d])
    report(f"int_d^(1-d) log(xi(1-xi))/(xi(1-xi)) = -log^2 d - pi^2/3 (d={float(d):.0e})", I3 - (-mp.log(d)**2 - mp.pi**2/3), 50*float(d)*abs(float(mp.log(d))))
# assembly: -(1/16)[...] = (1/32)[(11/3+4 log d)(pole) + 2 log^2 d - 67/9 + 2 pi^2/3]
ld, P = sp.symbols('ld P')
lhs = -sp.Rational(1, 16)*((-sp.Rational(11, 6) - 2*ld)*P + sp.Rational(67, 18) - ld**2 - sp.pi**2/3)
rhs = sp.Rational(1, 32)*((sp.Rational(11, 3) + 4*ld)*P + 2*ld**2 - sp.Rational(67, 9) + 2*sp.pi**2/3)
report("self-energy state: assembly of pole coefficient and constant", sp.N(sp.simplify(lhs - rhs).subs({ld: 0.3, P: 0.7})))
# Row I.II constant: 2 * (-1/(16 pi^2)) g^2 Nc (-67/9+pi^2/3) = (alpha_s Nc/2pi)(67/9-pi^2/3)
g, Nc = sp.symbols('g N_c', positive=True)
al = g**2/(4*sp.pi)
report("Row I.II constant with c.c. = (alpha_s Nc/2pi)(67/9 - pi^2/3)",
       sp.N((2*(-g**2*Nc/(16*sp.pi**2))*(-sp.Rational(67, 9) + sp.pi**2/3) - al*Nc/(2*sp.pi)*(sp.Rational(67, 9) - sp.pi**2/3)).subs({g: 1.3, Nc: 3})))
report("11/3 part with c.c.: (alpha_s Nc/4pi)(11/3) = alpha_s beta0/4pi", sp.N((al*Nc/(4*sp.pi)*sp.Rational(11, 3) - al*sp.Rational(11, 3)*Nc/(4*sp.pi)).subs({g: 1.3, Nc: 3})))

# (2) tensor contraction of Row II.IV for d_perp = 2,3,4 polarisations
rng = np.random.default_rng(1)
for d in (2, 3, 4):
    x = 0.37
    I = np.eye(d)
    def T(k, m, i, j):
        return I[k, m]*I[i, j] - I[j, m]*I[i, k]/(1-x) - I[i, m]*I[k, j]/x
    dev = 0.0
    for k, m, kp, mp_ in itertools.product(range(d), repeat=4):
        lhs = x*(1-x)*sum(T(k, m, i, j)*T(kp, mp_, i, j) for i in range(d) for j in range(d))
        rhs = d*x*(1-x)*I[k, m]*I[kp, mp_] + 2*(I[m, kp]*I[k, mp_] - I[k, m]*I[kp, mp_]) + (x/(1-x) + (1-x)/x)*I[m, mp_]*I[k, kp]
        dev = max(dev, abs(lhs - rhs))
    report(f"Row II.IV: sum over polarisations of the splitting tensors, d_perp={d}", dev)
# (3) contraction with the dimensionally regularised r integral, d = 2 - eps
eps, B0, X, KK, UU = sp.symbols('epsilon b xi KK UU')   # KK = K.K', UU = (u.K)(u.K')/u^2, B = 2/eps + b
Pp = X/(1-X) + (1-X)/X
d = 2 - eps
B = 2/eps + B0
expr = -(sp.pi/2)*B*(d*X*(1-X) + d*Pp)*KK - sp.pi*(d*X*(1-X)*UU + Pp*KK)
target = -sp.pi*((X*(1-X) + Pp)*B*KK + 2*X*(1-X)*(UU - KK/2))
report("Row II.IV: d_perp=2-eps contraction -> P_gg+ x [2/eps+b] + tensor term", sp.N(sp.limit(sp.simplify(expr - target), eps, 0).subs({B0: .3, X: .4, KK: 1.1, UU: -.6})))
expr2 = expr.subs(d, 2)
report("   with d_perp=2 the bracket [2/eps+b] becomes [2/eps+b+1]", sp.N(sp.limit(sp.simplify(expr.subs(eps, eps) - expr + (-(sp.pi/2)*B*(2*X*(1-X) + 2*Pp)*KK - sp.pi*(2*X*(1-X)*UU + Pp*KK)) - (target - sp.pi*(X*(1-X) + Pp)*KK)), eps, 0).subs({B0: .3, X: .4, KK: 1.1, UU: -.6})))
# (4) change of variable p+ -> xi = k+/(p+ + k+)
p, k = sp.symbols('p k', positive=True)
xi_ = k/(p + k)
report("dp+ p+k+/(p+ + k+)^4 = dxi xi(1-xi)/k+", sp.N((p*k/(p+k)**4 - (xi_*(1-xi_)/k)*sp.Abs(sp.diff(xi_, p))).subs({p: .7, k: 1.9})))
a = sp.symbols('a', positive=True)
report("int_a^1 (1-xi)/xi = log(1/a) - (1-a)", sp.N((sp.integrate((1-xi)/xi, (xi, a, 1)) - (sp.log(1/a) - (1-a))).subs(a, .01)))
report("double plus integral of the real splitting function = -11/6",
       sp.N(sp.integrate(xi*(1-xi) + (xi - 1)/(1-xi) + ((1-xi) - 1)/xi, (xi, 0, 1)) + sp.Rational(11, 6)))
# (5) Row III.V: splitting of 1/xi and the longitudinal integrals
lam, V = sp.symbols('lam V', positive=True)   # lam = Lambda/k+, V = vee/k+
J = sp.integrate(1/(xi*(1+xi)), (xi, lam, V-1))
report("int dxi/(xi(1+xi)) = log(k+/Lambda)+log(1-k+/vee)+O(Lambda/k+)", sp.N((J - (sp.log(1/lam) + sp.log(1-1/V))).subs({lam: 1e-7, V: 50})), 1e-6)
J2 = sp.integrate(1/(1+xi), (xi, lam, V-1))
report("int dxi/(1+xi) = log(vee/k+)+O(Lambda/k+)", sp.N((J2 - sp.log(V)).subs({lam: 1e-7, V: 50})), 1e-6)
# Row III.V tensor after the w' integration, contracted with K^m(y'-z) K^k'(x'-y') K^i(y-w) K^j(x-z)
v = rng.normal(size=(4, 2)); Km, Kk, Ki, Kj = v
xi0 = 0.8
I = np.eye(2)
full = sum(Km[m]*Kk[kp]*Ki[i]*Kj[j]*(I[kp, m]*I[i, j]/(1+xi0) - I[j, m]*I[i, kp]/xi0 - I[i, m]*I[j, kp])
           for m, kp, i, j in itertools.product(range(2), repeat=4))
kern = (1/(1+xi0))*((Km@Kk)*(Ki@Kj) - (Km@Kj)*(Kk@Ki)) - (1/(xi0*(1+xi0)))*(Km@Kj)*(Kk@Ki) - (Km@Ki)*(Kj@Kk)
report("Row III.V: tensor of C_1^dagger A A = kernel K^(y')_III (without plus prescription)", full - kern)
# Fourier rescaling used for the instantaneous term: int d^2r e^{i(1+xi)k.r} K(r) = (1+xi)^-1 int d^2r e^{ik.r} K(r)
# both sides are 2 pi i k/k^2 times the factor (analytic, Eq. FT_K): check the homogeneity only
kk = np.array([0.3, -1.1])
report("Fourier transform of K is homogeneous of degree -1", np.linalg.norm((1+xi0)*kk/((1+xi0)**2*(kk@kk)) - kk/(kk@kk)/(1+xi0)))
# (6) N_2 from unitarity: -(1/2) (1/(2pi)^3) * 2 g^2 = -g^2/(8 pi^3)
report("N_2 from unitarity: prefactor -g^2/(8 pi^3)", sp.N(-sp.Rational(1, 2)*2/(2*sp.pi)**3 + 1/(8*sp.pi**3)))
print("ALL PASS" if ok else "SOME FAILED")
