"""Check the three p+ integrals of the four-kernel row, and the small-Lambda structure.

Let  kap = k.(y'-z),  xi = p+/k+,  lam = Lambda/k+,  Xi = (V - k+)/k+.
The collapsed bracket is (1/k+)[ A/(1+xi) - B/xi - C ]  with
  A = (K1.K2)(K3.K4) = (y'-z).(x'-y') (y-w).(x-z)
  B = (K1.K4)(K2.K3) = (y'-z).(x-z)  (x'-y').(y-w)
  C = (K1.K3)(K2.K4) = (y'-z).(y-w)  (x'-y').(x-z)
and the p+ integral is  (1/k+) int dp+ (...) = int_lam^Xi dxi (...) e^{-i kap xi}.
"""
import mpmath as mp
mp.mp.dps = 30
Ei = mp.ei
g = mp.euler

def num(f, a, b):
    return mp.quad(f, [a, (a+b)/2, b])

kap = mp.mpf('1.7'); lam = mp.mpf('1e-4'); Xi = mp.mpf('37.0')
ph = lambda xi: mp.e**(-1j*kap*xi)

print("branch check:  Ei(-i eps)  vs  gamma + log(eps) - i pi/2   and   + i pi/2")
for eps in ['1e-3','1e-5','1e-7']:
    e = mp.mpf(eps); v = Ei(-1j*e)
    print(f"   eps={eps}:  Ei = {mp.nstr(v,10)}    g+log e - i pi/2 = "
          f"{mp.nstr(g+mp.log(e)-1j*mp.pi/2,10)}")
print("   (so Ei(-i eps) = gamma + log(eps) - i pi/2 + O(eps) for eps>0)")
print()

print("TERM  1/(1+xi)   [user's third line]")
lhs = num(lambda xi: ph(xi)/(1+xi), lam, Xi)
rhs = mp.e**(1j*kap)*(Ei(-1j*kap*(1+Xi)) - Ei(-1j*kap*(1+lam)))
print(f"   int_lam^Xi dxi e^-i kap xi/(1+xi) = {mp.nstr(lhs,12)}")
print(f"   e^{{i kap}}[Ei(-i kap(1+Xi)) - Ei(-i kap(1+lam))] = {mp.nstr(rhs,12)}")
print(f"   rel diff {mp.nstr(abs(lhs-rhs)/abs(lhs),3)}      <-- MATCHES (incl. the e^{{ik.(y'-z)}} prefactor)")
print()

print("TERM  -1/xi   [user's second line]")
lhs = -num(lambda xi: ph(xi)/xi, lam, Xi)
rhs = -(Ei(-1j*kap*Xi) - Ei(-1j*kap*lam))
print(f"   -int_lam^Xi dxi e^-i kap xi/xi = {mp.nstr(lhs,12)}")
print(f"   -[Ei(-i kap Xi) - Ei(-i kap lam)] = {mp.nstr(rhs,12)}")
print(f"   rel diff {mp.nstr(abs(lhs-rhs)/abs(lhs),3)}      <-- MATCHES")
print()

print("TERM  -1   [user's first line]")
lhs = -num(ph, lam, Xi)
rhs = (mp.e**(-1j*kap*Xi) - mp.e**(-1j*kap*lam))/(1j*kap)
print(f"   -int_lam^Xi dxi e^-i kap xi = {mp.nstr(lhs,12)}")
print(f"   [e^-i kap Xi - e^-i kap lam]/(i kap) = {mp.nstr(rhs,12)}")
print(f"   rel diff {mp.nstr(abs(lhs-rhs)/abs(lhs),3)}      <-- MATCHES")
print("   NOTE the denominator is i*kap = i k.(y'-z), with NO 1/k+ : the k+ of dp+ = k+ dxi")
print("   has already been spent producing it.")
print()

print("SMALL-Lambda STRUCTURE  (lam = Lambda/k+ -> 0)")
print("  1/(1+xi) term : Ei(-i kap (1+lam)) -> Ei(-i kap)          finite, NO rapidity log")
for L in ['1e-2','1e-4','1e-6','1e-8']:
    l = mp.mpf(L)
    print(f"     lam={L}:  Ei(-i kap(1+lam)) = {mp.nstr(Ei(-1j*kap*(1+l)),12)}   "
          f"Ei(-i kap) = {mp.nstr(Ei(-1j*kap),12)}")
print("  -1 term       : e^{-i kap lam} -> 1                        finite, NO rapidity log")
print("  -1/xi term    : Ei(-i kap lam) -> gamma + log(kap) + log(lam) - i pi/2")
print("                  ==> the ONLY rapidity divergence, coefficient of log(lam) is +1")
for L in ['1e-2','1e-4','1e-6','1e-8']:
    l = mp.mpf(L)
    ex = Ei(-1j*kap*l); ap = g + mp.log(kap*l) - 1j*mp.pi/2
    print(f"     lam={L}:  exact {mp.nstr(ex,12)}   asymptotic {mp.nstr(ap,12)}   "
          f"diff {mp.nstr(abs(ex-ap),3)}")
print()
print("  so, with l_k = log(k+/Lambda) = -log(lam) :")
print("     -[Ei(-i kap Xi) - Ei(-i kap lam)]  =  -l_k + log|k.(y'-z)| + gamma")
print("                                            - i(pi/2)sgn(k.(y'-z)) - Ei(-i kap Xi)")
print()
print("UPPER CUTOFF  V -> infinity  (Xi -> infinity): all three are oscillatory, NOT power-growing")
for X in [10, 100, 1000, 10000]:
    X = mp.mpf(X)
    print(f"   Xi={float(X):8.0f}:  -1 term {mp.nstr(abs((mp.e**(-1j*kap*X)-1)/(1j*kap)),8):>12}"
          f"   Ei(-i kap Xi) {mp.nstr(abs(Ei(-1j*kap*X)),8):>12}")
print("   -> bounded; and after the transverse integration Riemann-Lebesgue kills them.")
