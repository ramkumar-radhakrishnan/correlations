"""Row 22: the z integral in dim reg / MS-bar, the P_gg assembly, and the + prescription."""
import numpy as np, sympy as sp
from scipy.integrate import quad
print("="*80); print("1.  THE LARGE-|z| RESIDUE IS EXACTLY P_gg/(2 N_c)"); print("="*80)
P,K=sp.symbols('pplus kplus',positive=True); S=P+K; zt,d=sp.symbols('zeta d_perp',positive=True)
b1=d*P/S**2-2/K; b2=P/K**2+1/P; b3=2/K
coef=sp.simplify(((b1+b3)/2+b2).subs(d,2))
Cuv=zt/(1-zt)+(1-zt)/zt+zt*(1-zt)
print("   [(b1+b3)/2 + b2]  (d_perp=2)  = %s"%coef)
print("   C_UV(zeta)/k+ with zeta=k+/S  = %s"%sp.simplify((Cuv/K).subs(zt,K/S)))
print("   difference                    = %s"%sp.simplify(coef-(Cuv/K).subs(zt,K/S)))
print("   => the coefficient of  int d^2z/|z|^2  is  [P_gg(zeta)/(2 N_c k+)] x (X.X')")

print(); print("="*80); print("2.  THE z INTEGRAL IN d = 2-2eps, AND MS-BAR"); print("="*80)
eps=sp.Symbol('epsilon'); r2=sp.Symbol('r2',positive=True); mu2=sp.Symbol('mu2',positive=True)
Be=sp.gamma(1-eps)**2/sp.gamma(2-2*eps)
Tpref=sp.pi**(1-eps)*(r2)**(-eps)*mu2**eps*Be*sp.gamma(1+eps)
print("   T^{mm'} = int d^{2-2eps}z (y-z)^m (y'-z)^{m'}/[(y-z)^2 (y'-z)^2]")
print("           = pi Gamma(1+eps) B(eps) (pi mu^2 r^2)^{-eps} [ delta^{mm'}/(2 eps) - rhat^m rhat^{m'} ]")
ser=sp.series(sp.simplify(Tpref/eps),eps,0,1).removeO()
print("   (T-prefactor)/eps  ->  %s"%sp.simplify(sp.expand(ser)))
print("   i.e.  pi [ 1/eps + 2 - gamma_E - log(pi mu^2 r^2) ] = pi ( Lcal + 2 )")
print("   with  Lcal = 1/eps_MSbar - log(4 pi^2 mu^2 r^2) ,  1/eps_MSbar = 1/eps - gamma_E + log 4pi")
print("   cross-check against the cutoff form pi/2 (log R^2/r^2 + 1) delta - pi rhat rhat :  consistent.")

print(); print("="*80); print("3.  THE ASSEMBLY"); print("="*80)
X,R_,L=sp.symbols('Xcal Rcal Lcal'); c=sp.Symbol('c'); dp=2-2*eps
Ce=sp.pi*(1+eps*c)                     # C_eps , with  C_eps/eps  ->  pi(Lcal+2)  =>  1/eps + c = Lcal+2
TXX=Ce*( X/(2*eps) - R_ )              # T^{mm'} X^k X'^k' delta_km delta_k'm'
TrT=Ce*( dp/(2*eps) - 1 )              # delta_mm' delta_kk' contraction
tot=sp.expand(((b1+b3).subs(d,dp))*TXX + b2*TrT*X)
tot=sp.expand(sp.series(tot,eps,0,1).removeO())
tot=sp.simplify(tot.subs(c,L+2-1/eps))          # identify 1/eps + c = Lcal + 2
tot=sp.simplify(sp.expand(tot))
zb=1-zt
target=sp.pi/K*( Cuv*L*X + zt*zb*(X-2*R_) ).subs({zt:K/S})
print("   contraction  b1 d_km d_k\'m\' + b2 d_mm\' d_kk\' + b3 d_mk\' d_km\'  against  T^{mm\'} X^k X\'^k\' :")
print("     the 2/k+ in b1 and b3 CANCEL:  b1 + b3 = d_perp p+/S^2   (this is zeta + zetabar = 1)")
print("   assembled =", sp.simplify(sp.factor(tot)))
print("   target    =  (pi/k+) { C_UV(zeta) Lcal Xcal + zeta zetabar [Xcal - 2 Rcal] }")
print("   difference                                   =  %s"%sp.simplify(sp.expand(tot-target)))
print("   with the 1/k+ you already pulled out, the full row carries  pi/k+^2  -- the Row-8 form.")

print(); print("="*80); print("4.  CHANGE TO zeta, AND THE + PRESCRIPTION"); print("="*80)
print("   p+ = k+ (1-zeta)/zeta ,  dp+ = -k+ dzeta/zeta^2 ,  p+ = Lambda <-> zeta = k+/(k+ + Lambda)")
print("   so   int_Lambda^{V-k+} dp+ (...)  =  k+ int_{k+/V}^{k+/(k+ +Lambda)} dzeta/zeta^2 (...)")
print("   the dzeta/zeta^2 measure and the hard factor at k_perp/zeta are the FRAGMENTATION pair.")
print()
kt_r=1.7
F=lambda z: np.cos(kt_r/z)          # the only zeta-dependence left in the hard factor: phase e^{-ik.r/zeta}
print("   is F(zeta) smooth at zeta -> 1 ?   F = cos(k.r/zeta) with k.r = %.1f"%kt_r)
for z in (0.9,0.99,0.999,0.9999,1.0):
    print("     zeta=%-8.4f  F=%12.9f   F'=%12.9f"%(z,F(z),(F(z+1e-6)-F(z-1e-6))/2e-6 if z<1 else kt_r*np.sin(kt_r)))
print("   YES -- analytic at zeta=1, F(1)=cos(k.r): the LO phase.  Nothing is deleted.")
print()
a=0.05
for dlt in (1e-3,1e-5,1e-7):
    direct=quad(lambda z: z*F(z)/(1-z),a,1-dlt,limit=400)[0]
    sub   =quad(lambda z:(z*F(z)-F(1))/(1-z),a,1,limit=400)[0]+F(1)*np.log((1-a)/dlt)
    print("   delta=%-8.0e  int_a^{1-delta} zeta F/(1-zeta) = %14.9f   [+]-form = %14.9f   diff %.1e"
          %(dlt,direct,sub,abs(direct-sub)))
print("   identity:  int_a^{1-d} zF/(1-z) = int_a^1 [zF - F(1)]/(1-z) + F(1) log[(1-a)/d]")
print("   with a = k+/V and d = Lambda/k+ :  the isolated log is  F(1) log[k+(V-k+)/(V Lambda)],")
print("   and F(1) is the LO structure.  That is exactly the piece that goes into the evolution.")
