"""Group II Row IV with d_perp kept symbolic everywhere (never set to 2)."""
import sympy as sp
e,P,K,S=sp.symbols('epsilon pplus kplus S'); XX,Xr,r2=sp.symbols('XX Xr rsq')
gE=sp.EulerGamma; mu,r=sp.symbols('mu r',positive=True); pi=sp.pi
dp=2-2*e                                  # d_perp, kept exact
N=sp.gamma(1+e)*sp.gamma(1-e)**2/sp.gamma(2-2*e)
Om=sp.series(sp.expand(pi**(-e)*(r**2)**(-e)*mu**(-2*e)*N),e,0,2).removeO()   # the mu^{-2eps} choice
print("="*74); print("1.  G^{mm'} WITH d_perp KEPT"); print("="*74)
print("   Omega(eps) = pi^-e (r^2)^-e mu^-2e N(e) = 1 + eps*[%s]"
      %sp.simplify(sp.expand(Om-1)/e))
coef=sp.series(sp.expand(Om*(1/e-1)),e,0,1).removeO()
print("   Omega * (1/eps - 1)  =  %s"%sp.simplify(coef))
M=1/e+1-gE-sp.log(sp.pi*mu**2*r**2)
print("   i.e.  M = 1/eps + 1 - gammaE - log(pi mu^2 r^2)   ->  check: %s"%sp.simplify(coef-M))
MS=1/e-gE+sp.log(4*sp.pi)
print("   in MSbar :  M = 1/eps_MSbar + 1 - log(4 pi^2 mu^2 r^2)  ->  check: %s"
      %sp.simplify(M-(MS+1-sp.log(4*sp.pi**2*mu**2*r**2))))
print()
print("   G^{mm'} = pi [ delta^{mm'} M / d_perp  -  r^m r^m' / r^2 ]")
print("   sanity, expanding 1/d_perp = (1+eps)/2 :  pi delta/2 * (1+eps)(1/eps+1-...) = pi delta/2 (1/eps+2-...)")
Mc=sp.Symbol('Mc')                      # M as a symbol for the algebra below
print("   trace:  pi[ M - 1 ] = pi[1/eps - log(4 pi^2 mu^2 r^2)]  ->  check: %s"
      %sp.simplify((M-1)-(MS-sp.log(4*sp.pi**2*mu**2*r**2))))

print(); print("="*74); print("2.  THE CONTRACTION, d_perp NEVER SET TO 2"); print("="*74)
# contraction tables:  structure -> (against delta^{mm'} , against r^m r^m')
tabD={'a':XX,'b':dp*XX,'c':XX,'d':XX}
tabR={'a':Xr,'b':r2*XX,'c':Xr,'d':Xr}
def con(t): return dp*P/S**2*t['a'] + 2/K*(t['c']-t['d']) + (P/K**2+1/P)*t['b']
Gd=pi*Mc/dp*con(tabD)          # delta piece
Gr=-pi/r2*con(tabR)            # r^m r^m' piece
tot=sp.simplify(sp.expand(Gd+Gr))
target=pi*( P/S**2*(Mc*XX - dp*Xr/r2) + (P/K**2+1/P)*(Mc-1)*XX )
print("   delta piece  = %s"%sp.simplify(sp.expand(Gd)))
print("     -> every d_perp CANCELS: pi M [ p+/S^2 + p+/k+^2 + 1/p+ ] (X.X')")
print("   r r   piece  = %s"%sp.simplify(sp.expand(Gr)))
print("   total - pi[ p+/S^2 (M (X.X') - d_perp (X.r)(X'.r)/r^2) + (p+/k+^2+1/p+)(M-1)(X.X') ] = %s"
      %sp.simplify(sp.expand(tot-target)))
print()
print("   NOTE: no eps*L -> 1 manoeuvre was needed anywhere.  Keeping d_perp does the")
print("   bookkeeping automatically -- that is exactly why it is the safer convention.")
print("   The only surviving d_perp multiplies a FINITE term, so you may leave it symbolic.")
