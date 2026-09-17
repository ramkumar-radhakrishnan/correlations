"""Row 9: check the N1 / N2 split of |B_2|^2 against the quoted vertex."""
import sympy as sp
g,pi = sp.symbols('g pi',positive=True); P,K = sp.symbols('pplus kplus',positive=True)
S = P+K; I = sp.I
over = sp.Integer(4)/(2*pi)**3

B1 = I*g**2*sp.sqrt(P*K)/(4*pi**2*S)          # the one-rho piece of B_2
B2 = -g**2/(8*pi**2*sp.sqrt(P*K))             # the two-rho piece of B_2

print("="*74); print("PREFACTORS"); print("="*74)
n2 = sp.simplify(over*B1*sp.conjugate(B1).subs(sp.conjugate(g),g).subs(sp.conjugate(pi),pi))
n2 = sp.simplify(over*(I*g**2*sp.sqrt(P*K)/(4*pi**2*S))*(-I*g**2*sp.sqrt(P*K)/(4*pi**2*S)))
print("N2 = (1rho)x(1rho):", sp.simplify(n2), " =  (1/(2pi)^3)(g^4/4pi^4) x",
      sp.simplify(n2*(2*pi)**3*4*pi**4/g**4))
n1 = sp.simplify(over*B1*B2)
print("N1 = (1rho)x(2rho):", sp.simplify(n1), " = -(1/(2pi)^3)(i g^4/8pi^4) x",
      sp.simplify(n1/(-(1/(2*pi)**3)*I*g**4/(8*pi**4))))
print()
print("  sqrt(p+k+) x 1/sqrt(p+k+) = 1  -> N1 carries NO p+k+ factor   MATCHES")
print("  the vertex's overall 1/S distributes as")
print("     (1/S)[ d^ij(..)/2(z-w)^2 + (S/k+)(..) + (S/p+)(..) ]")
print("   = d^ij(..)/(2 S (z-w)^2) + (1/k+)(..) + (1/p+)(..)      MATCHES the N1 bracket")
print()
print("="*74); print("WHAT IS MISSING"); print("="*74)
print("  |B_2|^2 has THREE pieces, not two:")
print("     (1rho)x(1rho) -> 2 rho   = your N2")
print("     (1rho)x(2rho) -> 3 rho   = your N1   (two cross terms, the two you wrote)")
print("     (2rho)x(2rho) -> 4 rho   = NOT in your expression")
n4 = sp.simplify(over*B2*B2)
print("  the missing 4-rho term would carry", sp.simplify(n4),
      " = (1/(2pi)^3)(g^4/16 pi^4)/(p+ k+)")
print("  and is a pure product of four WW kernels:")
print("     (y-w)^i(x-z)^j(y'-w')^i(x'-z)^j / [(y-w)^2 (x-z)^2 (y'-w')^2 (x'-z)^2]")
print("  -> flag it; it has a 1/(p+ k+) which is a DOUBLE soft pole.")
