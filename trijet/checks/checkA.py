import sympy as sp
z,x = sp.symbols('zeta xi', positive=True)
zb = 1-z
eta = zb - x                       # eta = 1 - zeta - xi

# exact prefactor * [Phi.Phi*]_{delta delta} for the three contributions
A_aq = 8*z**2      *(eta**2+zb**2)*(z**2+zb**2)
A_q  = 8*eta**2    *(z**2+(1-eta)**2)*(eta**2+(1-eta)**2)
A_x  = -8*eta*z    *(1-x-2*z*eta)*(x*(1-x)+2*z*eta)
soft = 16*z**2*zb**2*(z**2+zb**2)

for n,A in [("A_qbar",A_aq),("A_q",A_q),("A_qqbar(interf.)",A_x)]:
    at0 = sp.factor(sp.simplify(A.subs(x,0)))
    print(f"{n:18s} at xi=0 : {at0}    equals soft? {sp.simplify(at0-soft)==0}")

print("\nexact subtracted combinations  [A(xi) - A(0)]/xi   (must be regular at xi=0):")
for n,A in [("A_qbar",A_aq),("A_q",A_q),("A_qqbar",A_x)]:
    d = sp.simplify(sp.cancel(sp.expand(A - A.subs(x,0))/x))
    print(f"  {n:8s}: {sp.factor(sp.simplify(d))}")
    print(f"           value at xi=0 : {sp.factor(sp.simplify(d.subs(x,0)))}")

print("\nclosed form for the antiquark channel:")
print("  A_qbar(xi) =", sp.factor(sp.expand(A_aq)))
print("  8 zeta^2 (zeta^2+zetabar^2)(2 zetabar^2 - 2 zetabar xi + xi^2) - A_qbar =",
      sp.simplify(8*z**2*(z**2+zb**2)*(2*zb**2-2*zb*x+x**2) - A_aq))
