"""From the z integral to P_gg: every step."""
import sympy as sp, numpy as np
eps=sp.symbols('epsilon'); z_=sp.symbols('zeta',positive=True); zb=1-z_
X,R_,L=sp.symbols('X R Lcal'); S,Kp,Pp=sp.symbols('S kplus pplus',positive=True)
mu,r2=sp.symbols('mu r2',positive=True)

print("="*76); print("STEP 6.  THE TRANSVERSE LOG IN MS-BAR"); print("="*76)
C=sp.pi**(1-eps)*sp.gamma(1+eps)*sp.gamma(1-eps)**2/sp.gamma(2-2*eps)
Sms=(sp.exp(sp.EulerGamma)/(4*sp.pi))**eps
expr=sp.expand(sp.series(sp.simplify(C*Sms*(r2*mu**2)**(-eps)/eps),eps,0,1).removeO())
print("   attach the MS-bar factor (e^gamma/4pi)^eps that rides on g^2 mu^{2eps}:")
print("      C(eps)(e^gamma/4pi)^eps (r^2 mu^2)^{-eps}/eps =",sp.simplify(expr))
print("      = pi [ 1/eps + 2 - log(4 pi^2 mu^2 r^2) ] :",
      sp.simplify(expr - sp.pi*(1/eps+2-sp.log(4*sp.pi**2*mu**2*r2)))==0)
print()
print("   define   Lcal = 1/eps - log(4 pi^2 mu^2 r^2)  ,  r = y' - y .")
print("   cutoff dictionary:  (pi/2)(log R^2/r^2 + 1) = (pi/2)[Lcal + 2]  =>  R = e^{1/2}/(2 pi mu)")
Rc=sp.Symbol('R',positive=True)
print("      solve:",[sp.simplify(s) for s in
      sp.solve(sp.Eq(sp.log(Rc**2/r2)+1, 2-sp.log(4*sp.pi**2*mu**2*r2)),Rc) if s.is_positive is not False])
print()
print("   WHY IT IS A COLLINEAR LOG:  r = y'-y is conjugate to the measured momentum")
print("   (the surviving phase is e^{-ik.r}), so Lcal is the logarithm between the")
print("   factorisation scale mu and the hard transverse scale ~ 1/r ~ k_perp.")

print(); print("="*76); print("STEP 7.  CONTRACT INTO THE ROW, IN THE SPLITTING VARIABLE"); print("="*76)
print("   zeta = k+/(p+ + k+) , zetabar = p+/(p+ + k+) :  the C vertex splits ONE gluon")
print("   of S = p+ + k+ into the measured k+ and the unobserved p+.  Then")
print("      a = S/p+ = 1/zetabar ,   b = S/k+ = 1/zeta .")
a,b=1/zb,1/z_; dp=sp.Symbol('d_perp')
c1=dp-2*a-2*b; c2=a**2+b**2; c3=2*a*b
print("   the two C brackets contract to")
print("      c1 = d_perp - 2/zetabar - 2/zeta   (delta_km delta_k'm')")
print("      c2 = 1/zeta^2 + 1/zetabar^2        (delta_mm' delta_kk')")
print("      c3 = 2/(zeta zetabar)              (delta_km' delta_k'm)")
print()
print("   T^{mm'} IS SYMMETRIC (step 4), so c1 and c3 contract identically and only c1+c3 acts:")
print("      c1 + c3 =",sp.simplify(c1+c3)," = d_perp exactly, because")
print("      2/(zeta zetabar) - 2/zeta - 2/zetabar = 2[1 - zeta - zetabar]/(zeta zetabar) =",
      sp.simplify(2*(1-z_-zb)/(z_*zb)))
print("   -> the 2/k+ term of the bracket dies.  It is the momentum sum rule.")
print()
print("   with d_perp = 2-2eps :  d_perp/(2eps) = 1/eps - 1 , (tr T)/Tc = d/(2eps) - 1 = 1/eps - 2,")
print("   and  Tc/eps -> pi(Lcal+2) , Tc -> pi :")
brace=((L+1)*X - 2*R_) + c2*L*X
net=z_*zb/Kp**2
full=sp.simplify(sp.expand(net*brace))
Cuv=z_/zb+zb/z_+z_*zb
target=(Cuv*L*X + z_*zb*(X-2*R_))/Kp**2
print("      result  minus  (pi/k+^2){ C_UV(zeta) Lcal X + zeta zetabar [X - 2R] } :",
      sp.simplify(sp.expand(sp.together(full-target))))

print(); print("="*76); print("STEP 8.  THAT COEFFICIENT IS P_gg"); print("="*76)
print("   C_UV(zeta) = zeta/zetabar + zetabar/zeta + zeta zetabar = zeta zetabar + 1/(zeta zetabar) - 2 :",
      sp.simplify(Cuv-(z_*zb+1/(z_*zb)-2))==0)
print("   P_gg(zeta) = 2 Nc [zeta/zetabar + zetabar/zeta + zeta zetabar] = 2 Nc C_UV(zeta)")
print("   the colour reduction of the row supplies exactly one Nc, so  Nc C_UV = P_gg/2 .")
print()
print("   dp+ = k+ dzeta/zeta^2 (fixed k+, S = k+/zeta), and dzeta/zeta^2 = (1/zeta) dS/S :")
zz=sp.Symbol('z',positive=True)
prim=2*sp.log(zz)-sp.log(1-zz)-zz+1/zz-1/(2*zz**2)
print("      primitive of C_UV/zeta^2 :",prim,"   check dF/dz - C_UV/z^2 =",
      sp.simplify(sp.diff(prim,zz)-(Cuv/z_**2).subs(z_,zz)))
print("      zeta -> 1 (p+ -> 0):  -log(1-zeta) -> log(k+/Lambda)   the RAPIDITY log")
print("      zeta -> 0 (p+ -> V):  1/zeta^3      -> V^2/(2k+^2)     the kinematic end")
print("      restricted to p+ <= k+ (zeta >= 1/2):  int = log(2 k+/Lambda), clean")
f=sp.lambdify(zz,Cuv.subs(z_,zz),'numpy')
from scipy.integrate import quad
for kp,lam in [(1.0,1e-6),(2.0,1e-5)]:
    b_=kp/(kp+lam)
    print("        k+=%.1f Lam=%.0e :  quad %12.6f   log(2k+/Lam) = %12.6f"
          %(kp,lam,quad(lambda t:f(t)/t**2,0.5,b_,limit=400)[0],np.log(2*kp/lam)))
