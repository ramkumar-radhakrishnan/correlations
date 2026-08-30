"""The p+ (xi) integration of the C1 x B2 row, in the style of Lublinsky-Mulian
   1610.03453 section 4.2 / appendix H.2.

Dictionary to the paper (Sec. 4.2):   z_paper <-> x ,  z'_paper <-> z ,  y <-> y
   Y  = y - z_paper  <-> y - x = s          Y' = y - z'_paper <-> y - z = s - r
   Z  = z_paper - z'_paper <-> x - z = -r   so Z^2 = r^2
   =>  xibar*D = xibar (s-r)^2 + xi s^2 = (1-xi)(Y')^2 + xi Y^2   <-- Eq. (4.16) denominator
"""
import sympy as sp

r1,r2,s1,s2,P1,P2,xi,dperp,kap = sp.symbols('r1 r2 s1 s2 P1 P2 xi d_perp kappa', real=True)
xb = 1-xi
r=sp.Matrix([r1,r2]); s=sp.Matrix([s1,s2]); P=sp.Matrix([P1,P2])
dot=lambda a,b:(a.T*b)[0,0]
r2_=dot(r,r); s2_=dot(s,s); P2_=dot(P,P); sr=dot(s,r); sP=dot(s,P); Pr=dot(P,r)
u=s-r; u2=dot(u,u)
d=lambda i,j: sp.Integer(1) if i==j else sp.Integer(0)

# ---- the six vertex structures V_i, i.e. the user's (1.4)-(1.10) with the
#      xi-weight stripped off.  The bracket of the source is
#         sum_i w_i(xi) V_i ,   w = (xi*xib, xi/xib, xib, xib/xi, xi, 1)
V = {}
V[1] = dperp/2*(2*Pr*sr/(P2_*r2_**2) - Pr/(P2_*r2_))
V[2] = sP/(P2_*r2_) + sP*(sr-r2_)/(2*P2_*r2_*u2)
V[3] = Pr/(P2_*r2_)*(2 - 2*sr/r2_ - sr/(2*s2_))
V[4] = (sP-Pr)/(P2_*r2_)*(1 - kap*sr/s2_)
V[5] = -xi*0 - Pr/P2_*(2*sr/r2_**2 + (s2_-sr)/(2*r2_*u2) - 1/(2*r2_))
V[6] = (2*Pr*sr/(P2_*r2_**2) - Pr/(P2_*r2_) + sP/(2*P2_*s2_) - sP*sr/(2*P2_*r2_*s2_)
        + sr*(sP-Pr)/(2*P2_*r2_*u2))
w = {1: xi*xb, 2: xi/xb, 3: xb, 4: xb/xi, 5: xi, 6: sp.Integer(1)}
bracket = sum(w[i]*V[i] for i in V)

# ---- cross-check against the contraction of Eq. (1.1) itself -------------------
def G(i,m):
    G1 = d(i,m)*dperp/(2*r2_)*(s2_-u2)*xi*xb
    G2 = -(xi/xb)*( s[i]*(-r[m])/r2_ - s[i]*u[m]/(2*u2) )
    G3 = xb*d(i,m)*( dot(u,-r)/r2_ + dot(s,u)/(2*s2_) - (s2_-u2)/(2*r2_) )
    G4 = -(xb/xi)*( u[i]*(-r[m])/r2_ + s[m]*u[i]/(2*s2_) )
    G5 = xi*d(i,m)*( dot(s,-r)/r2_ - dot(s,u)/(2*u2) - (s2_-u2)/(2*r2_) )
    G6 = -( s[m]*(-r[i])/r2_ - s[m]*u[i]/(2*u2) + u[m]*(-r[i])/r2_ + s[i]*u[m]/(2*s2_) )
    return [G1,G2,G3,G4,G5,G6]
src = sum(P[i]*r[m]*sum(G(i,m))/(P2_*r2_) for i in range(2) for m in range(2))
print("bracket == contraction of (1.1)  [kappa=1/2] :",
      sp.simplify(bracket.subs(kap,sp.Rational(1,2)) - src) == 0)
print("bracket == contraction of (1.1)  [kappa=1]   :",
      sp.simplify(bracket.subs(kap,1) - src) == 0)

# ---- small-r (UV) coefficient of each V_i, from the angular average ------------
th = sp.symbols('theta'); R = sp.symbols('R', positive=True)
ang = {r1: R*sp.cos(th), r2: R*sp.sin(th)}
print("\nUV coefficient of each V_i  (coefficient of (s.P)/(P^2 r^2) as r->0):")
for i in range(1,7):
    e = sp.simplify(V[i].subs(ang))
    avg = sp.integrate(e, (th, 0, 2*sp.pi))/(2*sp.pi)
    lead = sp.simplify(sp.limit(sp.expand(sp.simplify(avg))*R**2, R, 0))
    print(f"   V{i}: {sp.simplify(lead/(sP/P2_))}")
