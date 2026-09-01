"""Explicit reduction of the Mathematica ConditionalExpression for

      I(a,lam) = int_lam^{1-lam} dxi [ xi xib + xi/xib + xib/xi ] log( a/(1-xi) )

down to   (2 l_k - 11/6) log a + l_k^2/2 + pi^2/6 - 67/36.

Every regrouping below is checked as an exact polynomial identity in
(A, l, w, lam, Li2a, Li2b), where A = log a, l = log(1/lam), w = log(1-lam),
Li2a = Li2(lam), Li2b = Li2(1-lam).  Nothing is expanded or approximated.
"""
import sympy as sp

A, l, w, lam, Li2a, Li2b = sp.symbols('A l w lambda Li2a Li2b')
u = A - w                                   # Log[(a k)/(k - Lambda)]
v = A + l                                   # Log[(a k)/Lambda]
P = 12 - 12*lam + 3*lam**2 - 2*lam**3       # bracket in the -6 Log[a/(1-lam)] term
Q =  1 + 12*lam - 3*lam**2 + 2*lam**3       # bracket in the +6 Log[a/lam]  term

# --- the Mathematica result, verbatim, times 36 -------------------------------
T0 = -67 + 138*lam - 12*lam**2 + 8*lam**3
T1 = -18*u**2 + 18*v**2
T2 = -6*l                                   # 6 Log[Lambda/k]
T3 = -6*u*(P - 6*l)                         # -6 Log[a/(1-lam)] (P + 6 Log[Lambda/k])
T4 = -6*w                                   # -6 Log[1 - Lambda/k]
T5 =  6*v*(Q + 6*w)                         # +6 Log[a/lam] (Q + 6 Log[1-Lambda/k])
T6 = -36*Li2a + 36*Li2b
M36 = T0 + T1 + T2 + T3 + T4 + T5 + T6

# --- step 1 : difference of squares -------------------------------------------
step1 = 36*A*(l + w) + 18*(l**2 - w**2)
print("step 1  T1 = 36 A (l+w) + 18(l^2 - w^2)                   :", sp.expand(T1 - step1) == 0)

# --- step 2 : the two mixed products; the l*w cross terms cancel --------------
step2 = 6*A*(Q - P) + 36*A*(l + w) + 6*w*P + 6*l*Q
print("step 2  T3+T5 = 6A(Q-P) + 36A(l+w) + 6wP + 6lQ            :", sp.expand(T3 + T5 - step2) == 0)

# --- step 3 : absorb the lone -6l and -6w -------------------------------------
print("step 3a T2 + 6lQ = 6l(Q-1) = 6l(12lam-3lam^2+2lam^3)      :",
      sp.expand(T2 + 6*l*Q - 6*l*(12*lam - 3*lam**2 + 2*lam**3)) == 0)
print("step 3b T4 + 6wP = 6w(P-1) = 6w(11-12lam+3lam^2-2lam^3)   :",
      sp.expand(T4 + 6*w*P - 6*w*(11 - 12*lam + 3*lam**2 - 2*lam**3)) == 0)

# --- step 4 : the coefficient of log a ----------------------------------------
print("step 4a P + Q = 13                                        :", sp.expand(P + Q - 13) == 0)
print("step 4b Q - P = -11 + 24lam - 6lam^2 + 4lam^3             :",
      sp.expand(Q - P - (-11 + 24*lam - 6*lam**2 + 4*lam**3)) == 0)
Acoef = 2*l + 2*w - sp.Rational(11,6) + 4*lam - lam**2 + sp.Rational(2,3)*lam**3
print("step 4c coefficient of A is  2l + 2w - 11/6 + 4lam - lam^2 + 2lam^3/3 :",
      sp.expand(sp.together(72*(l+w) + 6*(Q-P) - 36*Acoef)) == 0)

# --- step 5 : the remainder ----------------------------------------------------
Brem = sp.Rational(1,36)*(18*l**2 - 67 + 36*(Li2b - Li2a)
        + 66*w - 18*w**2 + 138*lam - 12*lam**2 + 8*lam**3
        + 6*(l - w)*(12*lam - 3*lam**2 + 2*lam**3))
print("step 5  M/36 = Acoef * A + Brem   (exact identity)        :",
      sp.expand(M36/36 - (Acoef*A + Brem)) == 0)

# --- step 6 : the limit --------------------------------------------------------
lim = Acoef.subs({w:0, lam:0})*A + Brem.subs({w:0, lam:0, Li2a:0, Li2b:sp.pi**2/6})
print("step 6  lambda -> 0 :", sp.simplify(lim))
print("        target      :", sp.simplify((2*l - sp.Rational(11,6))*A + l**2/2 + sp.pi**2/6 - sp.Rational(67,36)))
print("        equal?      :",
      sp.simplify(lim - ((2*l - sp.Rational(11,6))*A + l**2/2 + sp.pi**2/6 - sp.Rational(67,36))) == 0)
