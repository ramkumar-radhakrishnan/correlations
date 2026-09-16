"""Row 8 DGLAP: clean extraction of the Lcal coefficient and the zeta endpoints."""
import sympy as sp, numpy as np
z_ = sp.symbols('zeta', positive=True); zb = 1-z_
X, R_, L = sp.symbols('X R Lcal')
Kp, Pp = sp.symbols('kplus pplus', positive=True)

# --- the z-integrated structure, built term by term (no shortcuts) -------------
c1 = 2 - 2/zb - 2/z_          # d_perp -> 2 in the non-pole parts
c13 = sp.simplify((2 - 2/zb - 2/z_) + 2/(z_*zb))
c2  = 1/z_**2 + 1/zb**2
print("c1 + c3 =", sp.simplify(c13), "   (= d_perp; the 2/k+ term is gone)")

# (c1+c3)=d_perp : d_perp/(2eps) X - d_perp R  ->  (Lcal+2-1) X - 2R
# c2            : (d_perp/(2eps) - 1) X       ->  c2 (Lcal+2-2) X
brace = ((L+1)*X - 2*R_) + c2*L*X
net   = z_*zb/Kp**2                       # p+/(k+ S^2) with S = k+/zeta
full  = sp.simplify(sp.expand(net*brace))
Cuv   = z_/zb + zb/z_ + z_*zb
target= (Cuv*L*X + z_*zb*(X - 2*R_))/Kp**2
print("full - (C_UV Lcal X + zeta zetabar [X-2R])/k+^2 =",
      sp.simplify(sp.expand(sp.together(full - target))))
print("coefficient of Lcal in full, x k+^2 :",
      sp.factor(sp.simplify(sp.expand(full*Kp**2).coeff(L))), " = C_UV(zeta) X :",
      sp.simplify(sp.expand(full*Kp**2).coeff(L) - Cuv*X) == 0)
print("P_gg(zeta)/(2 Nc) = C_UV :", sp.factor(sp.expand(Cuv)))

# --- the zeta integral --------------------------------------------------------
print()
print("PRIMITIVE of C_UV(zeta)/zeta^2 :")
prim = 2*sp.log(z_) - sp.log(1-z_) - z_ + 1/z_ - 1/(2*z_**2)
print("   F(zeta) =", prim, "     dF/dzeta - C_UV/zeta^2 =",
      sp.simplify(sp.diff(prim,z_) - Cuv/z_**2))
Lam, V = sp.symbols('Lambda V', positive=True)
zp, zm = Kp/(Kp+Lam), Kp/V
up = sp.simplify(sp.series(prim.subs(z_, zp), Lam, 0, 1).removeO())
lo = sp.simplify(prim.subs(z_, zm))
print()
print("   upper end zeta+ = k+/(k+ + Lambda), small Lambda :", sp.simplify(up))
print("   lower end zeta- = k+/V                           :", sp.simplify(lo))
print()
print("   => int = log(k+/Lambda) - 1/2 + 2 log(V/k+) - V/k+ + V^2/(2 k+^2)")
half = sp.simplify(prim.subs(z_, sp.Rational(1,2)))
print()
print("   if instead p+ <= k+ (i.e. zeta >= 1/2, the measured gluon the harder one):")
print("      F(1/2) =", sp.simplify(half))
print("      int_{1/2}^{zeta+} = log(k+/Lambda) + log 2 = log(2 k+/Lambda)   -- finite, no V")

# --- numerical confirmation ---------------------------------------------------
print()
print("NUMERICAL CHECK of the two endpoint formulae")
f = sp.lambdify(z_, Cuv/z_**2, 'numpy')
from scipy.integrate import quad
for kp,lam,Vv in [(1.0,1e-6,50.0),(2.0,1e-5,30.0),(0.5,1e-7,80.0)]:
    a,b = kp/Vv, kp/(kp+lam)
    num = quad(f,a,b,limit=400)[0]
    ana = np.log(kp/lam) - 0.5 + 2*np.log(Vv/kp) - Vv/kp + Vv**2/(2*kp**2)
    num2 = quad(f,0.5,b,limit=400)[0]; ana2 = np.log(2*kp/lam)
    print("   k+=%-4.1f Lam=%.0e V=%-5.1f | full: quad %14.6f  formula %14.6f  (%.1e)"
          % (kp,lam,Vv,num,ana,abs(num-ana)))
    print("   %24s | p+<k+: quad %14.6f  formula %14.6f  (%.1e)" % ("",num2,ana2,abs(num2-ana2)))
