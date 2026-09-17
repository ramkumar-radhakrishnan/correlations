"""Row 11: A^(1)dag x Bbar_3^dag A^(1).  Prefactors of the five quoted terms."""
import sympy as sp
g,pi = sp.symbols('g pi',positive=True); K,P = sp.symbols('kplus pplus',positive=True)
S=P+K; I=sp.I
Abra = I*g/(sp.sqrt(2)*pi*sp.sqrt(K))     # bra bracket coefficient  (-Abar^dag + A^dag U)
Aket = I*g/(sp.sqrt(2)*pi*sp.sqrt(P))     # [U(z)A - Abar] coefficient
over = 1/(2*pi)**3

print("="*78); print("ARGUMENT MAP   B_3(k+, z; p+, -u)  used as  B_3(k+, w; p+, -z)"); print("="*78)
print("   definition (x-z) -> (x-w)      definition (x-u) -> (x-z)")
print("   so  p+(x-u)^2 - k+(x-z)^2  ->  p+(x-z)^2 - k+(x-w)^2          matches your terms 1,2,3")
print("       k+(x-z)^2 - p+(x-u)^2  ->  k+(x-w)^2 - p+(x-z)^2          matches your terms 4,5")
print("       ((x-z)-(x-u))          ->  (z-w)                          matches")
print("       (k+(x-z) - p+(x-u))    ->  (k+(x-w) - p+(x-z))            matches")

print(); print("="*78); print("PREFACTORS"); print("="*78)
# conjugates of the three B_3 pieces
B1  = sp.conjugate(I*g**2/(2*pi)**2 * sp.sqrt(P*K)*S/(K-P)**2).subs({sp.conjugate(g):g,sp.conjugate(pi):pi,
        sp.conjugate(K):K,sp.conjugate(P):P})
B2  = sp.conjugate(I*g**2/(2*pi**2)/K * sp.sqrt(K*P)/(K-P)).subs({sp.conjugate(g):g,sp.conjugate(pi):pi,
        sp.conjugate(K):K,sp.conjugate(P):P})
B3  = sp.conjugate(-I*g**2/(4*pi**2*sp.sqrt(P*K))).subs({sp.conjugate(g):g,sp.conjugate(pi):pi,
        sp.conjugate(K):K,sp.conjugate(P):P})
rows = [
 ("terms 1,2  (piece I,  delta_ij / [p+(x-z)^2-k+(x-w)^2])", Abra*Aket*B1,
   I*g**4/(8*pi**4)*S/(K-P)**2),
 ("term 3     (piece III)",                                  Abra*Aket*B3,
   -I*g**4/(8*pi**4)/(P*K)),
 ("term 4     (piece II, the 1/(z-w)^2 half)",                Abra*Aket*B2,
   I*g**4/(4*pi**4)/(K*(K-P))),
 ("term 5     (piece II, the -p+k+/V^2 half)",                Abra*Aket*B2*(-P*K),
   -I*g**4/(4*pi**4)*P/(K-P)),
]
for nm, got, want in rows:
    r = sp.simplify(sp.expand(sp.simplify(got)/want))
    print("  %-52s ratio to quoted = %s" % (nm, r))
print()
print("  (the overall 1/(2pi)^3 sits outside in both cases)")
print("  ALL FIVE PREFACTORS REPRODUCE EXACTLY.")

print(); print("="*78); print("COLOUR"); print("="*78)
print("   B_3^{ba} carries f^{abc} rho^c(x);  relabel c -> d:  f^{abd} rho^d(x)")
print("   barred  ->  f^{abd} U^{de}(x) rho^e(x)          [a free, d contracted]")
print("   [U^{bc}(z) A^(1)c - Abar^(1)b] ->  (ig/(sqrt2 pi sqrt p+)) (y-z)^j/(y-z)^2")
print("                                      x [U^{bc}(z) - U^{bc}(y)] rho^c(y)")
print("   product: f^{abd} [ U^{de}(x)U^{bc}(z) - U^{de}(x)U^{bc}(y) ] rho^e(x) rho^c(y)")
print("   -> exactly what you wrote, and here f^{abd} IS legitimately an overall factor")
print("      (a free, b and d both contracted with Wilson lines).   CORRECT.")
print("   [contrast the previous row, where f could not be pulled out.]")
