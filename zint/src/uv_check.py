"""Independent re-derivation of the UV (r -> 0) coefficient, three ways."""
import numpy as np, sympy as sp

# ---------------------------------------------------------------- 1. symbolic, general d
print("1. ANGULAR AVERAGE IN d TRANSVERSE DIMENSIONS")
print("   <r^i r^j> = delta^ij r^2 / d   ->   < (P.r)(s.r)/r^4 > = (s.P)/(d r^2)")
print("   c_1 = d_perp xi xib   (the instantaneous delta^im term of Theta_1 carries d_perp)")
d, dp, xi = sp.symbols('d d_perp xi', positive=True)
xb = 1-xi
C_UV = dp*xi*xb/d + xi/xb + xb/xi
print("   C_UV = c_1/d + c_2 =", sp.simplify(C_UV))
print("   consistent dim reg, d = d_perp :", sp.simplify(C_UV.subs(d, dp)))
print("   Pgg/(2Nc)                      :", sp.simplify(xi*xb + xi/xb + xb/xi))
print("   equal?", sp.simplify(C_UV.subs(d, dp) - (xi*xb+xi/xb+xb/xi)) == 0)
print("   hybrid (vertex in d_perp, 2D angular average), d = 2 :", sp.simplify(C_UV.subs(d,2)))
print()

# ------------------------------------------- 2. numeric angular average of the RAW bracket
print("2. DIRECT ANGULAR AVERAGE OF THE BRACKET OF EQ. (1.1), built from its tensors")
s = np.array([0.7,-0.4]); P = np.array([1.3,0.5]); s2 = s@s; P2 = P@P; sP = s@P
def bracket(r, xi, dperp=2.0):
    """P^i r^m / (P^2 r^2) contracted with the six structures of (1.1), verbatim."""
    xb = 1-xi
    r2 = np.sum(r*r,-1); u = s-r; u2 = np.sum(u*u,-1)           # u = y-z = s-r
    Pr = r@P; sr = r@s; Pu = u@P; su = u@s
    dd = Pr/(P2*r2)                                             # delta^im contracted
    G1 = dd*dperp/(2*r2)*(s2-u2)*xi*xb
    G2 = -(xi/xb)*( (s@P)*(-r2)/r2 - (s@P)*(np.sum(u*r,-1))/(2*u2) )/(P2*r2)
    G3 = xb*dd*( np.sum(u*(-r),-1)/r2 + su/(2*s2) - (s2-u2)/(2*r2) )
    G4 = -(xb/xi)*( Pu*(-r2)/r2 + sr*Pu/(2*s2) )/(P2*r2)
    G5 = xi*dd*( -sr/r2 - su/(2*u2) - (s2-u2)/(2*r2) )
    G6 = -( sr*(-Pr)/r2 - Pu*sr/(2*u2) + np.sum(u*r,-1)*(-Pr)/r2 + (s@P)*np.sum(u*r,-1)/(2*s2) )/(P2*r2)
    return G1+G2+G3+G4+G5+G6
th = 2*np.pi*np.arange(20000)/20000
for xi_ in (0.3, 0.65, 0.5):
    xb_ = 1-xi_
    for R in (1e-3, 1e-4):
        rv = np.stack([R*np.cos(th), R*np.sin(th)], -1)
        avg = bracket(rv, xi_).mean()*R**2/(sP/P2)
        print(f"   xi={xi_:<5} R={R:.0e}: <bracket> r^2 P^2/(s.P) = {avg:.9f}"
              f"   C_UV = {xi_*xb_+xi_/xb_+xb_/xi_:.9f}   ratio {avg/(xi_*xb_+xi_/xb_+xb_/xi_):.7f}")
print()

# ------------------------------------------------- 3. the S1 term on its own, both factors
print("3. THE S_1 TERM ON ITS OWN")
print("   S_1 = (P.r)(s.r)/r^4 ,  total coefficient c_1 = d_perp xi xib = 2 xi xib  (d_perp=2)")
for R in (1e-3,):
    rv = np.stack([R*np.cos(th), R*np.sin(th)], -1)
    S1 = ((rv@P)*(rv@s)/np.sum(rv*rv,-1)**2).mean()
    print(f"   <S_1> r^2 / (s.P) = {S1*R**2/sP:.9f}    predicted 1/d = 1/2 = 0.5")
    print(f"   so  c_1 <S_1> = 2 xi xib * (s.P)/(2 r^2) = xi xib (s.P)/r^2   -> weight xi(1-xi) in Pgg")
