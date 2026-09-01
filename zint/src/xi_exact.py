"""The UV xi-integral at finite lambda = Lambda/k+, exactly, and its lambda -> 0 limit.

  I(a,lam) = int_lam^{1-lam} dxi C_UV(xi) log(a/(1-xi)),   C_UV = Pgg/(2Nc),
             a = (y-x)^2/m^2 ,  l_k = log(1/lam) = log(k+/Lambda).

  Exact:   I = A(lam) log a + B(lam)
  Limit:   I -> (2 l_k - 11/6) log a + l_k^2/2 + pi^2/6 - 67/36
"""
import mpmath as mp
mp.mp.dps = 40
Li2 = lambda z: mp.polylog(2, z)
C   = lambda x: x*(1-x) + x/(1-x) + (1-x)/x

def A(lam):
    lk = mp.log(1/lam)
    return 2*lk + 2*mp.log(1-lam) - mp.mpf(11)/6 + 4*lam - lam**2 + mp.mpf(2)/3*lam**3

def B(lam):
    l, L1 = mp.log(lam), mp.log(1-lam)
    Ga = lambda t: -(t**2/2)*mp.log(t) + t**2/4 + (t**3/3)*mp.log(t) - t**3/9
    a_ = Ga(1-lam) - Ga(lam)
    b_ = (-L1**2/2 + (1-lam)*L1 - (1-lam)) - (-l**2/2 + lam*l - lam)
    c_ = Li2(1-lam) - Li2(lam) - lam*l + (1-lam)*L1 - 1 + 2*lam
    return a_ + b_ + c_

def exact(a, lam):  return A(lam)*mp.log(a) + B(lam)
def limit(a, lam):
    lk = mp.log(1/lam)
    return (2*lk - mp.mpf(11)/6)*mp.log(a) + lk**2/2 + mp.pi**2/6 - mp.mpf(67)/36

def mathematica(a, lam):
    """The user's ConditionalExpression with k = 1, Lambda = lam."""
    L1, L2, Ll, L3 = mp.log(a/(1-lam)), mp.log(a/lam), mp.log(lam), mp.log(1-lam)
    return (mp.mpf(1)/36)*(-67 + 138*lam - 12*lam**2 + 8*lam**3
            - 18*L1**2 + 18*L2**2 + 6*Ll
            - 6*L1*(12 - 12*lam + 3*lam**2 - 2*lam**3 + 6*Ll) - 6*L3
            + 6*L2*(1 + 12*lam - 3*lam**2 + 2*lam**3 + 6*L3)
            - 36*Li2(lam) + 36*Li2(1-lam))

if __name__ == "__main__":
    a = mp.mpf(100)
    print("a = (y-x)^2/m^2 = 100")
    print(f"{'lambda':>8} {'direct quadrature':>22} {'Mathematica':>22} {'my exact A,B':>22} "
          f"{'lambda->0 limit':>22} {'limit - exact':>14}")
    for e in range(2, 9):
        lam = mp.mpf(10)**(-e)
        q  = mp.quad(lambda x: C(x)*mp.log(a/(1-x)), [lam, 1-lam])
        mm, ex, li = mathematica(a, lam), exact(a, lam), limit(a, lam)
        print(f"{'1e-%d'%e:>8} {mp.nstr(q,16):>22} {mp.nstr(mm,16):>22} {mp.nstr(ex,16):>22} "
              f"{mp.nstr(li,16):>22} {mp.nstr(li-ex,6):>14}")
    print()
    print("the terms dropped in the limit,  exact - limit  =")
    print("   log(a) * [ 2 log(1-lam) + 4 lam - lam^2 + 2 lam^3/3 ]  +  O(lam log^2 lam)")
    for e in (2, 4, 6):
        lam = mp.mpf(10)**(-e)
        print(f"   lam=1e-{e}:  exact-limit = {mp.nstr(exact(a,lam)-limit(a,lam),8)}"
              f"   [ relative to the answer: {mp.nstr((exact(a,lam)-limit(a,lam))/exact(a,lam),3)} ]")
