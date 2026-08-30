"""Symbolic reduction of the C1 x B2 (two-rho) bracket to the r-space structure T(r).

Everything is done with explicit 2-component vectors r, s, P and the light-cone
fraction xi.  Nothing is assumed: T(r) is built by contracting the six structures
of Eq. (1.1) with the prefactor (x'-w')^i (z-x)^m / [(x'-w')^2 (z-x)^2], and then
matched against
   (a) the user's hand-reduced Eqs. (1.4)-(1.10),
   (b) the 9-structure basis  S1..S9  used in the note.
"""
import sympy as sp

r1, r2, s1, s2, P1, P2, xi, dperp = sp.symbols('r1 r2 s1 s2 P1 P2 xi d_perp', real=True)
xb = 1 - xi
r = sp.Matrix([r1, r2]); s = sp.Matrix([s1, s2]); P = sp.Matrix([P1, P2])
dot = lambda a, b: (a.T*b)[0, 0]
r2_ = dot(r, r); s2_ = dot(s, s); P2_ = dot(P, P)
sr = dot(s, r); sP = dot(s, P); Pr = dot(P, r)
u  = s - r                      # y - z = s - r
u2 = dot(u, u)                  # (y-z)^2 = (r-s)^2

d = lambda i, j: sp.Integer(1) if i == j else sp.Integer(0)

def G(i, m):
    """The six structures of Eq. (1.1), as the tensor multiplying P^i r^m/(P^2 r^2)."""
    G1 = d(i, m)*dperp/(2*r2_)*(s2_ - u2)*xi*xb
    G2 = -(xi/xb)*( s[i]*(-r[m])/r2_ - s[i]*u[m]/(2*u2) )
    G3 = xb*d(i, m)*( dot(u, -r)/r2_ + dot(s, u)/(2*s2_) - (s2_ - u2)/(2*r2_) )
    G4 = -(xb/xi)*( u[i]*(-r[m])/r2_ + s[m]*u[i]/(2*s2_) )
    G5 = xi*d(i, m)*( dot(s, -r)/r2_ - dot(s, u)/(2*u2) - (s2_ - u2)/(2*r2_) )
    G6 = -( s[m]*(-r[i])/r2_ - s[m]*u[i]/(2*u2) + u[m]*(-r[i])/r2_ + s[i]*u[m]/(2*s2_) )
    return [G1, G2, G3, G4, G5, G6]

# T_n = P^i r^m G_n^{im} / r^2   (the universal 1/P^2 is stripped off)
Tn = [sp.together(sp.simplify(sum(P[i]*r[m]*G(i, m)[n] for i in range(2) for m in range(2))/r2_))
      for n in range(6)]
T_total = sp.simplify(sum(Tn))

# ---- (a) the user's Eqs. (1.4)-(1.10), P^2 stripped off ------------------------
U1 = dperp*xi*xb/2*(2*Pr*sr/r2_**2 - Pr/r2_)
U2 = (xi/xb)*(sP/r2_ + sP*(sr - r2_)/(2*r2_*u2))
U3 = xb*Pr/r2_*(2 - 2*sr/r2_ - sr/(2*s2_))
U4_user = (xb/xi)*(sP - Pr)/r2_*(1 - sr/s2_)          # as printed in the user's PDF
U4_fix  = (xb/xi)*(sP - Pr)/r2_*(1 - sr/(2*s2_))      # with the 1/2 restored
U5 = -xi*Pr*(2*sr/r2_**2 + (s2_ - sr)/(2*r2_*u2) - 1/(2*r2_))
U6 = (2*Pr*sr/r2_**2 - Pr/r2_ + sP/(2*s2_) - sP*sr/(2*r2_*s2_)
      + sr*(sP - Pr)/(2*r2_*u2))

print("Theta_1 == (1.4):", sp.simplify(Tn[0] - U1) == 0)
print("Theta_2 == (1.5):", sp.simplify(Tn[1] - U2) == 0)
print("Theta_3 == (1.6)/(1.7):", sp.simplify(Tn[2] - U3) == 0)
print("Theta_4 == (1.8) as printed:", sp.simplify(Tn[3] - U4_user) == 0)
print("Theta_4 == (1.8) with s^2 -> 2 s^2:", sp.simplify(Tn[3] - U4_fix) == 0)
print("Theta_5 == (1.9):", sp.simplify(Tn[4] - U5) == 0)
print("Theta_6 == (1.10):", sp.simplify(Tn[5] - U6) == 0)

# ---- (b) the 9-structure basis -----------------------------------------------
S = {1: Pr*sr/r2_**2, 2: sP/r2_, 3: Pr/r2_, 4: Pr*sr/(s2_*r2_), 5: sP*sr/(s2_*r2_),
     6: sP/s2_, 7: sP*dot(r, r - s)/(r2_*u2), 8: Pr*dot(s, r - s)/(r2_*u2),
     9: sr*dot(P, r - s)/(r2_*u2)}
c = {1: 2*xi*xb, 2: xi/xb + xb/xi, 3: xi**2 - sp.Rational(5, 2)*xi + 2 - 1/xi,
     4: xb**2/(2*xi), 5: -1/(2*xi), 6: sp.Rational(1, 2), 7: -xi/(2*xb),
     8: xi/2, 9: sp.Rational(-1, 2)}
T_basis = sum(c[k]*S[k] for k in S)
print("T (d_perp=2) == sum_k c_k S_k :",
      sp.simplify((T_total.subs(dperp, 2) - T_basis)) == 0)

# ---- UV split:  T = C_UV (s.P)/r^2 + T_fin ------------------------------------
C_UV = xi*xb + xi/xb + xb/xi
T_fin = sp.simplify(T_total.subs(dperp, 2) - C_UV*sP/r2_)
Tfin_basis = (2*xi*xb*(S[1] - S[2]/2) + c[3]*S[3] + c[4]*S[4] + c[5]*S[5]
              + c[6]*S[6] + c[7]*S[7] + c[8]*S[8] + c[9]*S[9])
print("T_fin basis form:", sp.simplify(T_fin - Tfin_basis) == 0)
print("C_UV =", sp.simplify(C_UV), " = Pgg/(2Nc)")
print("c3 factorised:", sp.factor(sp.together(c[3])))
