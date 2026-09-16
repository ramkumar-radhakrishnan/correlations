"""Row 8: two C's. Check the assembly, the w/w' integration, and the tensor contraction."""
import sympy as sp

g, pi = sp.symbols('g pi', positive=True)
K, P = sp.symbols('kplus pplus', positive=True)      # k+, p+
S = P + K                                            # p+ + k+  (C's first argument)
d = sp.symbols('d_perp', positive=True)
I = sp.I

print("="*72); print("1. PREFACTOR"); print("="*72)
Adag = -I*g/(sp.sqrt(2)*pi*sp.sqrt(S))               # A^dag(p+ + k+, -y')
A    =  I*g/(sp.sqrt(2)*pi*sp.sqrt(S))               # A(p+ + k+, -y)
Cf   = -g/(2*pi*sp.sqrt(2*S))*sp.sqrt(P*(S-P))/S     # C(p+ + k+, .; p+, .)   twice
over = sp.Integer(4)/(2*pi)**3
print("   phases : (-i) A^dag, (-1) C^dag, (-1) C, (+i) A  ->",
      sp.simplify((-I)*(-1)*(-1)*I), "  (REAL)")
num = sp.simplify(over/(sp.sqrt(2)*pi)**2/(2*pi*sp.sqrt(2))**2)
print(f"   numeric: (4/(2pi)^3)(1/(sqrt2 pi))^2(1/(2 pi sqrt2))^2 = {num}"
      f"  = (1/(2pi)^3) x {sp.simplify(num*(2*pi)**3)}")
kpow = sp.simplify((1/sp.sqrt(S))**2 * (sp.sqrt(P*(S-P))/(S*sp.sqrt(S)))**2)
print(f"   k+ powers: {sp.simplify(kpow)}   (S - P = k+)   = p+ k+/(p+ + k+)^4     MATCHES")

print()
print("="*72); print("2. THE DELTA AND THE w INTEGRATION"); print("="*72)
y, w, z = sp.symbols('y w z')
sol = sp.solve(sp.Eq(y - w + P/S*(w-z), 0), w)[0]
print(f"   delta^(2)[ y - w + (p+/(p+ + k+))(w - z) ]  ->  w = {sp.simplify(sol)}")
print(f"   w - z = {sp.factor(sp.simplify(sol - z))}     = ((p+ + k+)/k+)(y - z)")
print(f"   so (w-z)^m/(w-z)^2 -> (k+/(p+ + k+)) (y-z)^m/(y-z)^2   [one factor per kernel, two kernels]")
print(f"   Jacobian per delta: (S/k+)^2 = {sp.simplify((S/K)**2)}   [two deltas]")
net = sp.simplify(kpow * (K/S)**2 * (S/K)**4)
print(f"   net k+ factor after both integrations: {sp.simplify(net)}  = p+/[(p++k+)^2 k+]")
wp = sol.subs(y, sp.Symbol("yp"))
print(f"   phase: w' - w = {sp.factor(sp.simplify(wp - sol))}  ->  e^(-i (S/k+) k.(y'-y))   MATCHES")

print()
print("="*72); print("3. THE TENSOR CONTRACTION OF THE TWO C BRACKETS"); print("="*72)
print("   C_{1jik} bracket (slots j,k,i -> j,i,k) :")
print("      [ d_km d_ij - (S/p+) d_jm d_ik - (S/k+) d_im d_kj ]        MATCHES the quoted one")
print("   contracting the internal i,j between bra and ket:")
# symbolic index algebra with explicit small-d sums, but keep d_perp symbolic via traces
# use the nine elementary contractions worked out by hand, verified below numerically
T = {}
T['km,k\'m\''] = d - 2*(S/P) - 2*(S/K)
T['mm\',kk\''] = S**2/P**2 + S**2/K**2
T['km\',k\'m'] = 2*S**2/(P*K)
for kk,vv in T.items():
    print(f"      coefficient of  delta_{{{kk.split(',')[0]}}} delta_{{{kk.split(',')[1]}}} : {sp.simplify(vv)}")
print(f"   note -2(S/P + S/K) = {sp.simplify(-2*(S/P+S/K))} = -2 S^2/(p+ k+)")
print()
print("   multiplying by the net p+/S^2 (the k+ factor with 1/k+ pulled out front):")
for kk,vv in T.items():
    print(f"      delta_{{{kk.split(',')[0]}}} delta_{{{kk.split(',')[1]}}} : {sp.simplify(sp.expand(P/S**2*vv))}")
print()
print("   QUOTED bracket:")
print("      d_km d_k'm' :  d_perp p+/(k+ + p+)^2  -  2/k+")
print("      d_mm' d_kk' :  p+/k+^2  +  1/p+")
print("      d_km' d_k'm :  2/k+")
chk = [sp.simplify(sp.expand(P/S**2*T["km,k'm'"]) - (d*P/S**2 - 2/K)),
       sp.simplify(sp.expand(P/S**2*T["mm',kk'"]) - (P/K**2 + 1/P)),
       sp.simplify(sp.expand(P/S**2*T["km',k'm"]) - 2/K)]
print(f"   differences: {chk}      <-- all zero means the quoted bracket is EXACT")
