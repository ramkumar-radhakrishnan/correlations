"""MSbar for  G^{mm'} = pi^{1-e}(r^2)^{-e} mu^{2e} N(e) [ delta/d_perp (1/e - 1) - r^m r^m'/r^2 ]."""
import sympy as sp
e=sp.symbols('epsilon'); mu,r=sp.symbols('mu r',positive=True); GE=sp.EulerGamma
dperp=2-2*e
print("1.  the delta coefficient is EXACTLY 1/(2 eps) -- no expansion needed:")
print("    (1/d_perp)(1/eps - 1) =", sp.simplify((1/dperp)*(1/e-1)))
print("    (if you set d_perp = 2 first you get (1/2)(1/eps - 1) -- wrong by -1/2, the d_perp trap)")

N=sp.gamma(1-e)**2/sp.gamma(2-2*e)*sp.gamma(1+e)
print("\n2.  N(eps) =", sp.series(N,e,0,2))
pref=sp.pi**(1-e)*(r**2)**(-e)*mu**(2*e)*N
ser=sp.series(sp.expand(pref/sp.pi),e,0,2).removeO()
print("    prefactor/pi = 1 + eps*[ ... ] with [...] =", sp.simplify(sp.expand(ser-1)/e))

B=sp.simplify(sp.expand(sp.series(pref/(2*e),e,0,1).removeO()*2/sp.pi))
print("\n3.  bracket (delta part, x 2/pi) =", sp.logcombine(sp.expand(B),force=True))
Bexp=1/e+2-GE+sp.log(mu**2/(sp.pi*r**2))
print("    equals  1/eps + 2 - gammaE + log(mu^2/(pi r^2)) ? ", sp.simplify(sp.expand(B-Bexp)))

MS=1/e-GE+sp.log(4*sp.pi)
print("\n4.  MSbar regrouping:")
print("    B - [1/eps_MSbar + 2 + log(mu^2/(4 pi^2 r^2))] =",
      sp.simplify(sp.expand(Bexp-(MS+2+sp.log(mu**2/(4*sp.pi**2*r**2))))))
print("\n5.  by rescaling mu^2 -> mu^2 e^gE/(4 pi)  (correct direction for a mu^{+2eps}):")
resc=Bexp.subs(mu,sp.sqrt(mu**2*sp.exp(GE)/(4*sp.pi)))
print("    ->", sp.simplify(sp.logcombine(sp.expand(resc),force=True)))
print("    equals 1/eps + 2 + log(mu^2/(4 pi^2 r^2)) ? ",
      sp.simplify(sp.expand(resc-(1/e+2+sp.log(mu**2/(4*sp.pi**2*r**2))))))
print("\n6.  for comparison, the mu^{-2eps} version (dimensionally consistent):")
pref2=sp.pi**(1-e)*(r**2)**(-e)*mu**(-2*e)*N
B2=sp.simplify(sp.expand(sp.series(pref2/(2*e),e,0,1).removeO()*2/sp.pi))
resc2=B2.subs(mu,sp.sqrt(4*sp.pi*mu**2*sp.exp(-GE)))
print("    after mu^2 -> 4 pi mu^2 e^-gE :", sp.simplify(sp.logcombine(sp.expand(resc2),force=True)))
print("    equals 1/eps + 2 - log(4 pi^2 mu^2 r^2) ? ",
      sp.simplify(sp.expand(resc2-(1/e+2-sp.log(4*sp.pi**2*mu**2*r**2)))))
