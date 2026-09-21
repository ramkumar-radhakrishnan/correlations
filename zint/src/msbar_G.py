"""MSbar for  G^{mm'} = (pi/2) delta [1/eps + 2 - gammaE - log(pi mu^2 r^2)] - pi r^m r^m'/r^2 ."""
import sympy as sp
eps,gE,mu,r,pi=sp.symbols('epsilon gamma_E mu r pi',positive=True)
GE=sp.EulerGamma

# 0.  where the expression comes from: N(eps) (pi mu^2 r^2)^{-eps} [ 1/(2 eps) delta - rhat rhat ]
N=sp.gamma(1+eps)*sp.gamma(1-eps)**2/sp.gamma(2-2*eps)
print("N(eps) = Gamma(1+e)Gamma(1-e)^2/Gamma(2-2e) =", sp.series(N,eps,0,2))
pref=sp.series(N*(sp.pi*mu**2*r**2)**(-eps)/(2*eps),eps,0,1).removeO()
print("N (pi mu^2 r^2)^-eps /(2 eps) =", sp.simplify(sp.expand(pref)))
print("   i.e. (1/2)[1/eps + 2 - gammaE - log(pi mu^2 r^2)]   <-- your bracket, confirmed")

B = 1/eps + 2 - GE - sp.log(sp.pi*mu**2*r**2)          # your bracket
MS = 1/eps - GE + sp.log(4*sp.pi)                       # the MSbar pole combination
print()
print("STEP 1  (scheme-independent regrouping):")
print("   bracket - [1/eps_MSbar] - 2 + log(4 pi^2 mu^2 r^2) =",
      sp.simplify(B - MS - 2 + sp.log(4*sp.pi**2*mu**2*r**2)))
print("   => bracket = 1/eps_MSbar + 2 - log(4 pi^2 mu^2 r^2)")
print()
print("STEP 2  (the rescaling that makes mu the MSbar scale):")
wrong = B.subs(mu, sp.sqrt(mu**2*sp.exp(GE)/(4*sp.pi)))
right = B.subs(mu, sp.sqrt(4*sp.pi*mu**2*sp.exp(-GE)))
print("   mu^2 -> mu^2 e^gE /(4 pi)  gives :", sp.simplify(sp.expand(sp.logcombine(sp.expand(wrong),force=True))))
print("   mu^2 -> 4 pi mu^2 e^-gE    gives :", sp.simplify(sp.expand(sp.logcombine(sp.expand(right),force=True))))
print()
print("   check the first one equals  1/eps + 2 - log( mu^2 r^2 e^{2 gE} / 4 ) :",
      sp.simplify(sp.expand(wrong - (1/eps + 2 - sp.log(mu**2*r**2*sp.exp(2*GE)/4)))))
print("   check the second equals     1/eps + 2 - log( 4 pi^2 mu^2 r^2 )      :",
      sp.simplify(sp.expand(right - (1/eps + 2 - sp.log(4*sp.pi**2*mu**2*r**2)))))
