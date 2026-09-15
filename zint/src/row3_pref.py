import sympy as sp
g, pi = sp.symbols('g pi', positive=True)
kp, pp = sp.symbols('kplus pplus', positive=True); I = sp.I
print("PREFACTOR of row 3")
Adag = -I*g/(sp.sqrt(2)*pi*sp.sqrt(kp))                        # A^dag(k+, -w')
Cpre = -g/(2*pi*sp.sqrt(2*kp))*sp.sqrt(pp*(kp-pp))/kp          # C(k+, w-z; p+, z-x)
# B_2(k+-p+, -z; p+, -x): with K = k+-p+, P = p+ :  sqrt(P K)/(4 pi^2 [...] (P+K)),  P+K = k+
Bpre = I*g**2*sp.sqrt(pp*(kp-pp))/(4*pi**2*kp)
overall = sp.Integer(2)/(2*pi)**3
print("   signs/i : (-i) from A^dag, (-1) from C, (+i) from B  ->  (-i)(-1)(i) = -1  (real)")
print("   numeric : (2/(2pi)^3) x 1/(sqrt2 pi) x 1/(2 sqrt2 pi) x 1/(4 pi^2) =",
      sp.simplify(overall/(sp.sqrt(2)*pi)/(2*sp.sqrt(2)*pi)/(4*pi**2)),
      " vs quoted (1/(2pi)^3)(1/(8 pi^4)) =", sp.simplify(1/(2*pi)**3/(8*pi**4)))
kpow = sp.simplify(1/sp.sqrt(kp) * sp.sqrt(pp*(kp-pp))/(kp*sp.sqrt(kp)) * sp.sqrt(pp*(kp-pp))/kp)
print("   k+ powers:", kpow, "  = p+ (k+ - p+) / k+^3      <-- exactly the quoted form")
tot = sp.simplify(overall*Adag*Cpre*Bpre)
print("   full  :", tot)
print("   quoted: (1/(2pi)^3) (g^4/(8 pi^4)) p+(k+-p+)/k+^3 =",
      sp.simplify(g**4/(2*pi)**3/(8*pi**4)*pp*(kp-pp)/kp**3))
print("   ratio :", sp.simplify(tot/(g**4/(2*pi)**3/(8*pi**4)*pp*(kp-pp)/kp**3)),
      "  (-1 = the overall sign carried by the bracket structure)")
print()
print("B_2's Feynman denominator with K = k+-p+, w -> z, P = p+, z -> x, internal x -> y:")
print("   P(y-x)^2 + K(y-z)^2 = p+ (y-x)^2 + (k+-p+)(y-z)^2     <-- matches")
print("C's delta with z_arg = w-z, x_arg = z-x:")
print("   delta^(2)( w-z + (p+/k+)(z-x) )                        <-- matches, no Jacobian needed")
print()
print("C's index bracket here has FIRST argument k+ and SECOND p+, so it carries")
print("   k+/p+   AND   k+/(k+-p+)   -- BOTH endpoints, unlike the earlier rows.")
print("Consistently, C's Theta(k+ - p+ - Lambda) gives the limits int_Lambda^(k+-Lambda).")
