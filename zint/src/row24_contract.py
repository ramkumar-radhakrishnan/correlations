"""Row 24: the T^{mm'} contraction with the bracket, step by step."""
import numpy as np, sympy as sp, itertools
print("="*80); print("1.  THE THREE CONTRACTIONS, DONE INDEX BY INDEX (numerically, d=2)"); print("="*80)
rng=np.random.default_rng(101); d=2
T=rng.normal(size=(d,d)); T=(T+T.T)/2          # T^{mm'} is symmetric
X=rng.normal(size=d); Xp=rng.normal(size=d); dd=np.eye(d)
c=lambda A: np.einsum('mp,k,q->',T,X,Xp)*0+sum(
    A[m,k,mp,kp]*T[m,mp]*X[k]*Xp[kp] for m,k,mp,kp in itertools.product(range(d),repeat=4))
B1=np.zeros((d,d,d,d)); B2=np.zeros((d,d,d,d)); B3=np.zeros((d,d,d,d))
for m,k,mp,kp in itertools.product(range(d),repeat=4):
    B1[m,k,mp,kp]=dd[k,m]*dd[kp,mp]; B2[m,k,mp,kp]=dd[m,mp]*dd[k,kp]; B3[m,k,mp,kp]=dd[m,kp]*dd[k,mp]
TXX=X@T@Xp; trT=np.trace(T); TXpX=Xp@T@X
print("   delta_km delta_k'm'  ->  T^{mm'} X^m X'^m'  = T(X,X') :  %14.9f  vs  %14.9f"%(c(B1),TXX))
print("   delta_mm' delta_kk'  ->  (tr T)(X.X')                 :  %14.9f  vs  %14.9f"%(c(B2),trT*(X@Xp)))
print("   delta_mk' delta_km'  ->  T^{mm'} X'^m X^m' = T(X',X)  :  %14.9f  vs  %14.9f"%(c(B3),TXpX))
print("   and T symmetric => T(X',X) = T(X,X') :  diff %.2e"%abs(TXX-TXpX))
print("   ==> the b1 and b3 terms carry the SAME structure and simply ADD.")

print(); print("="*80); print("2.  b1 + b3 : THE 2/k+ CANCELS (the momentum sum rule)"); print("="*80)
P,K=sp.symbols('pplus kplus',positive=True); S=P+K; dp=sp.Symbol('d_perp')
c1=dp-2*S/P-2*S/K; c2=S**2*(1/P**2+1/K**2); c3=2*S**2/(P*K)
b1,b2,b3=[sp.simplify(P/S**2*q) for q in (c1,c2,c3)]
print("   b1 = (p+/S^2) c1 = %s"%sp.simplify(sp.expand(b1)))
print("        = d_perp p+/S^2 - 2/S - 2p+/(S k+)   and  2/S + 2p+/(S k+) = (2/k+)(zeta + zetabar) = 2/k+")
print("   b3 = (p+/S^2) c3 = %s"%b3)
print("   b1 + b3          = %s"%sp.simplify(b1+b3))
print("   check  b1+b3 - d_perp p+/S^2 = %s"%sp.simplify(b1+b3-dp*P/S**2))

print(); print("="*80); print("3.  THE TWO T-CONTRACTIONS IN d = 2-2eps"); print("="*80)
eps=sp.Symbol('epsilon'); L,Xc,Rc,cc=sp.symbols('Lcal Xcal Rcal c')
Ce=sp.pi*(1+eps*cc)                      # C_eps ;  C_eps/eps -> pi(Lcal+2)  <=>  1/eps + c = Lcal+2
sub=lambda e: sp.simplify(sp.expand(sp.series(e,eps,0,1).removeO()).subs(cc,L+2-1/eps))
TXXe=Ce*(Xc/(2*eps)-Rc)
trTe=Ce*((2-2*eps)/(2*eps)-1)
print("   (i)  X_m X'_m' T^{mm'} = C_eps [ Xcal/(2eps) - Rcal ]        (X, X' carry NO eps)")
print("          ->  %s"%sp.simplify(sub(TXXe)))
print("   (ii) tr T = delta_mm' T^{mm'} = C_eps [ d_perp/(2eps) - 1 ] ,  d_perp = 2-2eps")
print("          = C_eps/eps (1-eps) - C_eps = C_eps/eps - C_eps - C_eps")
print("          ->  %s        <-- the +2 is EATEN"%sp.simplify(sub(trTe)))

print(); print("="*80); print("4.  ASSEMBLE"); print("="*80)
first=sp.simplify(sub((dp*P/S**2).subs(dp,2-2*eps)*TXXe))
second=sp.simplify(sub(b2*trTe*Xc))
print("   (b1+b3) T(X,X') = (2-2eps) p+/S^2 x C_eps[Xcal/2eps - Rcal]")
print("        = (p+/S^2)[ (C_eps/eps - C_eps) Xcal - 2 C_eps Rcal ]  ->  (pi p+/S^2)[ (Lcal+1) Xcal - 2 Rcal ]")
print("        symbolic: %s"%sp.factor(first))
print("   b2 (tr T) Xcal  ->  pi (p+/k+^2 + 1/p+) Lcal Xcal")
print("        symbolic: %s"%sp.factor(second))
Q=sp.simplify(first+second)

print(); print("="*80); print("5.  TO zeta VARIABLES"); print("="*80)
zt=sp.Symbol('zeta',positive=True); zb=1-zt
print("   S = k+/zeta ,  p+ = (zetabar/zeta) k+  =>")
print("     p+/S^2   = zeta zetabar / k+")
print("     p+/k+^2  = zetabar/(zeta k+)")
print("     1/p+     = zeta/(zetabar k+)")
Cuv=zt/(1-zt)+(1-zt)/zt+zt*(1-zt)
target=sp.pi/K*( Cuv*L*Xc + zt*zb*(Xc-2*Rc) )
print("   Q - (pi/k+){ C_UV(zeta) Lcal Xcal + zeta zetabar [Xcal - 2 Rcal] }  =  %s"
      %sp.simplify(sp.expand(Q - target.subs(zt,K/S))))
print("   (the '+1' left over in (Lcal+1) is exactly what supplies the zeta zetabar Xcal)")
print()
print("   C_UV(zeta) = zeta/zetabar + zetabar/zeta + zeta zetabar = P_gg(zeta)/(2 N_c) :")
for z in (0.2,0.5,0.8):
    print("     zeta=%.1f :  C_UV = %.6f   =  zeta zetabar + 1/(zeta zetabar) - 2 = %.6f"
          %(z,z/(1-z)+(1-z)/z+z*(1-z), z*(1-z)+1/(z*(1-z))-2))

print(); print("="*80); print("6.  THE p+ -> zeta MEASURE"); print("="*80)
V,Lam=sp.symbols('vee Lambda',positive=True)
print("   p+ = k+(1-zeta)/zeta  =>  dp+/dzeta = %s"%sp.simplify(sp.diff(K*(1-zt)/zt,zt)))
print("   |dp+| = k+ dzeta/zeta^2 ;  p+=Lambda -> zeta = k+/(k+ + Lambda) ;  p+=V-k+ -> zeta = k+/V")
print("   int_Lambda^{V-k+} dp+  =  k+ int_{k+/V}^{k+/(k+ + Lambda)} dzeta/zeta^2")
print()
print("   overall constant:  (1/(2pi)^3)(g^4/4pi^4)(N_c/k+) x (1/2pi) x (pi/k+) = %s"
      %sp.simplify(1/(2*sp.pi)**3*sp.Symbol('g')**4/(4*sp.pi**4)*sp.Symbol('N_c')/K*(1/(2*sp.pi))*sp.pi/K))
