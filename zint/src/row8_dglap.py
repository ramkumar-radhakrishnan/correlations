"""Row 8: the DGLAP content.  z integrated FIRST, then the bracket rewritten in the
splitting variable zeta = k+/(p+ + k+)."""
import sympy as sp, numpy as np

z_, eps, dp = sp.symbols('zeta epsilon d_perp')
zb = 1 - z_
S, Kp, Pp = sp.symbols('S kplus pplus', positive=True)
X, R_, L = sp.symbols('X R Lcal')          # Xcal, Rcal, Lcal = 1/eps - log(4 pi^2 mu^2 r^2)

print("="*76); print("1.  THE SPLITTING VARIABLE"); print("="*76)
print("   A^dagger creates ONE gluon of light-cone momentum  S = p+ + k+  from rho.")
print("   C splits it into the measured gluon (k+) and the unobserved one (p+).")
print("   So the DGLAP variable is the measured daughter's fraction of the PARENT:")
print("       zeta = k+/(p+ + k+),   zetabar = p+/(p+ + k+),   zeta + zetabar = 1")
print("   NOT xi = p+/k+.  In these variables  a = S/p+ = 1/zetabar,  b = S/k+ = 1/zeta.")

print(); print("="*76); print("2.  THE BRACKET IN zeta -- and why the 2/k+ term dies"); print("="*76)
a, b = 1/zb, 1/z_
c1 = dp - 2*a - 2*b          # delta_km delta_k'm'
c2 = a**2 + b**2             # delta_mm' delta_kk'
c3 = 2*a*b                   # delta_km' delta_k'm
print("   c1 = d_perp - 2/zetabar - 2/zeta ,  c2 = 1/zeta^2 + 1/zetabar^2 ,  c3 = 2/(zeta zetabar)")
print("   after int d^2z the tensor T^{mm'} is SYMMETRIC, so c1 and c3 contract identically:")
print("       c1 + c3 =", sp.simplify(c1+c3), "   <- exactly d_perp, because")
print("       2/(zeta zetabar) - 2/zeta - 2/zetabar = 2[1 - zeta - zetabar]/(zeta zetabar) =",
      sp.simplify(2*(1-z_-zb)/(z_*zb)))
print("   the cancellation IS the momentum sum rule zeta + zetabar = 1.")

print(); print("="*76); print("3.  z INTEGRATED FIRST, IN d = 2 - 2 eps"); print("="*76)
print("   T^{mm'} = C(eps)(r^2 mu^2)^{-eps} [ delta^{mm'}/(2 eps) - rhat^m rhat^{m'} ]")
print("   delta_km delta_k'm' , delta_km' delta_k'm  ->  T^{kk'} X^k X'^k'  = Tc [X/(2eps) - R]")
print("   delta_mm' delta_kk'                        ->  (tr T) (X.X')      = Tc (d/(2eps) - 1) X")
print("   with d_perp = 2 - 2 eps :  d_perp/(2 eps) = 1/eps - 1 ,  d/(2eps) - 1 = 1/eps - 2")
tot_eps = (c1+c3)*(X/(2*eps) - R_) + c2*(dp/(2*eps) - 1)*X
tot_eps = tot_eps.subs(dp, 2-2*eps)
print("   longitudinal x transverse, before expanding:")
print("     ", sp.simplify(tot_eps))

print(); print("="*76); print("4.  EXPANDING:  Tc/eps -> pi(Lcal + 2),  Tc -> pi"); print("="*76)
# substitute 1/eps -> Lcal + 2 in the pole pieces, d_perp -> 2 elsewhere
expr = ( 2*((L+2)/2*X - R_)                       # (c1+c3) = d_perp ; d_perp/(2eps)X = (1/eps -1)X
         - X                                       # the -1 from d_perp/(2 eps) = 1/eps - 1
         + c2*((L+2) - 2)*X )                      # c2 (1/eps - 2) X
expr = sp.simplify(sp.expand(expr))
print("   { } =", sp.factor(sp.simplify(expr)))
pref = zb/z_                                       # the net k+ factor p+/(k+ S^2) = (zb/z)/S^2
full = sp.simplify(sp.expand(pref*expr))
print("   x (zetabar/zeta)  [the net k+ factor is (zetabar/zeta)/S^2] :")
print("     ", sp.simplify(full))

print(); print("="*76); print("5.  THE SPLITTING FUNCTION"); print("="*76)
Cuv = z_/zb + zb/z_ + z_*zb
print("   C_UV(zeta) = zeta/zetabar + zetabar/zeta + zeta zetabar =",
      sp.simplify(sp.expand(Cuv)), " = P_gg/(2 Nc)")
print("   also      = zeta zetabar + 1/(zeta zetabar) - 2 :",
      sp.simplify(Cuv - (z_*zb + 1/(z_*zb) - 2)) == 0)
coefL = sp.simplify(sp.expand(full.coeff(L)))
print()
print("   coefficient of Lcal   :", sp.factor(coefL))
print("   C_UV(zeta)/zeta^2     :", sp.factor(sp.simplify(Cuv/z_**2)))
print("   DIFFERENCE            :", sp.simplify(sp.expand(coefL - Cuv/z_**2*X)))
rem = sp.simplify(sp.expand(full - coefL*L))
print("   non-log remainder     :", sp.factor(rem))
print("   compare zeta zetabar (X - 2R)/zeta^2 :",
      sp.simplify(sp.expand(rem - z_*zb*(X-2*R_)/z_**2)))

print(); print("="*76); print("6.  S^2 = k+^2/zeta^2  ABSORBS THE 1/zeta^2"); print("="*76)
print("   net k+ factor = (zetabar/zeta)/S^2 = (zetabar/zeta)(zeta^2/k+^2) = zeta zetabar/k+^2")
print("   => the whole z-integrated longitudinal x transverse structure is")
print()
print("        (pi/k+^2) { C_UV(zeta) Lcal X  +  zeta zetabar [ X - 2 R ] }")
print()
chk = sp.simplify(sp.expand( (z_**2)*full/1 - (Cuv*L*X + z_*zb*(X-2*R_)) ))
print("   check against the direct expansion  (should be 0) :", chk)

print(); print("="*76); print("7.  CROSS-CHECK AGAINST THE p+ FORM ALREADY ESTABLISHED"); print("="*76)
old = (L*X/Pp + Pp/Kp**2*L*X + Pp/(Kp+Pp)**2*((L+1)*X - 2*R_))/Kp
new = (Cuv*L*X + z_*zb*(X-2*R_))/Kp**2
sub = {z_: Kp/(Kp+Pp)}
print("   old (p+ variables, x 1/k+) minus new (zeta variables) :",
      sp.simplify(sp.expand(sp.together(old - new.subs(sub).subs(zb, Pp/(Kp+Pp))))))

print(); print("="*76); print("8.  THE zeta INTEGRAL AND ITS TWO ENDS"); print("="*76)
print("   dp+ = k+ d zeta / zeta^2   (at fixed k+, since S = k+/zeta)")
zz = sp.Symbol('z', positive=True)
integ = sp.simplify((Cuv/z_**2).subs(z_, zz))
prim = sp.simplify(sp.integrate(sp.apart(integ, zz), zz))
print("   integrand C_UV/zeta^2 =", sp.apart(integ, zz))
print("   primitive             =", prim)
zp = Kp/(Kp+sp.Symbol('Lambda',positive=True)); zm = Kp/sp.Symbol('V',positive=True)
Lam, V = sp.Symbol('Lambda',positive=True), sp.Symbol('V',positive=True)
val = sp.simplify(prim.subs(zz,zp) - prim.subs(zz,zm))
print("   evaluated on  zeta in [k+/V , k+/(k+ + Lambda)] :")
print("     ", sp.simplify(sp.expand(val)))
print()
print("   zeta -> 1 end (p+ -> 0):  -log(1-zeta) -> log(k+/Lambda) = the RAPIDITY log")
print("   zeta -> 0 end (p+ -> V):  1/zeta^3 -> V^2/(2 k+^2), the known kinematic growth")
print()
print("   and the classic sum rule, C_UV integrated over the full range:")
print("      int_0^1 [ zeta zetabar + 1/(zeta zetabar) - 2 ] d zeta = 1/6 - 2 + (rapidity logs)")
print("      1/6 - 2 =", sp.Rational(1,6)-2, " ->  2 l_k - 11/6 ;  x 2Nc  ->  4 Nc l_k - 11Nc/3 = -b0")
