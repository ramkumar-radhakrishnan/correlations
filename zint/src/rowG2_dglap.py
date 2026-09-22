"""Group II Row IV: the p+ bracket IS P_gg, and the p+ integral IS a DGLAP convolution."""
import sympy as sp, mpmath as mp
P,K,S,z=sp.symbols('pplus kplus S zeta',positive=True); Nc=sp.Symbol('N_c',positive=True)
print("="*76); print("1.  THE BRACKET IS THE GLUON SPLITTING FUNCTION"); print("="*76)
print("   the C_1 vertex splits a parent of  S = p+ + k+  into the measured gluon (k+)")
print("   and an unmeasured one (p+).  The measured momentum fraction is")
print("        zeta = k+/(k+ + p+) ,   1 - zeta = p+/(k+ + p+)")
sub={P:K*(1-z)/z}                                   # p+ = k+ (1-z)/z
def inzeta(expr): return sp.simplify(sp.expand(sp.simplify(expr.subs(P,K*(1-z)/z))))
for name,expr in (("p+/(k+ + p+)^2", P/(K+P)**2), ("p+/k+^2", P/K**2), ("1/p+", 1/P)):
    print("   %-16s = %s"%(name, sp.simplify(K*inzeta(expr))), " / k+")
tot=inzeta(P/(K+P)**2 + P/K**2 + 1/P)
Cuv=z*(1-z)+(1-z)/z+z/(1-z)
print()
print("   sum x k+  = %s"%sp.simplify(K*tot))
print("   zeta zbar + zbar/zeta + zeta/zbar = %s"%sp.simplify(Cuv))
print("   difference = %s      <-- ZERO"%sp.simplify(sp.expand(K*tot-Cuv)))
Pgg=2*Nc*Cuv
print()
print("   and  P_gg(zeta) = 2 N_c [ zeta/(1-zeta) + (1-zeta)/zeta + zeta(1-zeta) ] , so")
print("        p+/S^2 + p+/k+^2 + 1/p+  =  P_gg(zeta) / (2 N_c k+)")

print(); print("="*76); print("2.  THE MEASURE AND THE HARD SCALE"); print("="*76)
print("   p+ = k+ (1-zeta)/zeta  =>  dp+ = - k+ dzeta / zeta^2 , so")
print("        int_Lambda^{V-k+} dp+  =  k+ int_{k+/V}^{k+/(k++Lambda)} dzeta / zeta^2")
print("   -> the FRAGMENTATION measure  dzeta/zeta^2 .")
print("   the phase is  exp[-i k.(y'-y) (k+ + p+)/k+] = exp[-i (k/zeta).(y'-y)]")
print("   -> the hard factor sits at  k_perp / zeta , not k_perp .")
print("   the Wilson lines sit at x, x', y, y' -- the PARENT positions, zeta-independent.")

print(); print("="*76); print("3.  CONSISTENCY WITH THE sigma = 1/zeta MASTERS"); print("="*76)
s=sp.Symbol('sigma',positive=True)
lhs=sp.simplify(sp.expand(Cuv.subs(z,1/s)))
rhs=(s-1)/s**2 + (s-1) + 1/(s-1)
print("   C_UV(1/sigma) - [ (sigma-1)/sigma^2 + (sigma-1) + 1/(sigma-1) ] = %s"
      %sp.simplify(sp.expand(lhs-rhs)))
print("   i.e.  J1 + J2  =  int dzeta/zeta^2 e^{-i kappa/zeta} P_gg(zeta)/(2 N_c) .")

print(); print("="*76); print("4.  THE PLUS PRESCRIPTION IS THE DGLAP ONE"); print("="*76)
print("   p+ -> 0  <=>  zeta -> 1 , and the pole there is  zeta/(1-zeta) .")
print("   its companion in the integrand is  h(zeta) = e^{-i kappa/zeta} / zeta , h(1) = e^{-i kappa} : FINITE.")
print("   so   1/(1-zeta) -> [1/(1-zeta)]_+ + delta(1-zeta) log(k+/Lambda)")
mp.mp.dps=20
def lhsI(kap,lam,z0): return mp.quad(lambda t: mp.e**(-1j*kap/t)/(t*(1-t)),[z0,1-lam])
def rhsI(kap,lam,z0):
    h=lambda t: mp.e**(-1j*kap/t)/t
    return h(1)*mp.log((1-z0)/lam)+mp.quad(lambda t:(h(t)-h(1))/(1-t),[z0,1-1e-12])
print("   %-22s %26s %26s"%("(kappa,lam,zeta0)","direct","plus form"))
for kap,lam,z0 in ((0.7,1e-6,0.02),(2.0,1e-7,0.05)):
    print("   k=%.1f lam=%.0e z0=%.2f %26s %26s"%(kap,lam,z0,mp.nstr(lhsI(kap,lam,z0),10),mp.nstr(rhsI(kap,lam,z0),10)))

print(); print("="*76); print("5.  THE DGLAP NORMALISATION"); print("="*76)
g,pi=sp.symbols('g pi',positive=True); al=sp.Symbol('alpha_s',positive=True)
print("   NLO pole coefficient  :  g^4/(8 pi^4 k+)  x  P_gg/2")
print("   LO  (2/(2pi)^3)|A^(1)|^2 :  g^2/(pi^2 k+)")
ratio=sp.simplify((g**4/(8*pi**4))/(g**2/pi**2))
print("   ratio = %s  =  g^2/(8 pi^2)  =  alpha_s/(2 pi)   [alpha_s = g^2/4pi]"%ratio)
print("   so   d^3N_NLO = (alpha_s/2pi)(1/eps) int dzeta/zeta^2 [P_gg(zeta)/2] d^3N_LO(k/zeta) + ...")
print("   i.e. exactly the collinear-factorisation form; the residual factor 1/2 on P_gg")
print("   depends on the LO normalisation convention and on whether the hermitian-conjugate")
print("   row is counted separately -- fix it against your own LO expression.")
