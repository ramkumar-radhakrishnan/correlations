"""Is the split an EXACT rearrangement for finite Lambda, or only to O(Lambda)?"""
import mpmath as mp
mp.mp.dps = 30
Ei = mp.ei
kap = mp.mpf('1.7')
print("identity:  int_lam^Xi dxi g/xi  ==  int_lam^Xi dxi [g-g(0)]/xi  +  g(0) log(Xi/lam)")
print("with g = e^{-i kap xi}, g(0) = 1.   NOTE the + piece keeps the lower limit lam.")
Xi = mp.mpf('37.0')
for lam in ['1e-2','1e-4','1e-6']:
    l = mp.mpf(lam)
    lhs = Ei(-1j*kap*Xi) - Ei(-1j*kap*l)
    plus = mp.quad(lambda u: (mp.e**(-1j*kap*u)-1)/u, [l, 1/kap, Xi])
    rhs = plus + mp.log(Xi/l)
    plus0 = mp.quad(lambda u: (mp.e**(-1j*kap*u)-1)/u, [0, 1/kap, Xi])
    print(f"   lam={lam}:  lhs {mp.nstr(lhs,12):>30}   rhs(lower=lam) {mp.nstr(rhs,12):>30}"
          f"   diff {mp.nstr(abs(lhs-rhs),3)}")
    print(f"              {'':>30}   rhs(lower=0)   {mp.nstr(plus0+mp.log(Xi/l),12):>30}"
          f"   diff {mp.nstr(abs(lhs-plus0-mp.log(Xi/l)),3)}")
print("-> keeping the lower limit at Lambda makes the split EXACT (diff ~ 1e-30);")
print("   sending it to 0 costs only O(Lambda).  Your int_Lambda^(V-k+) is exactly right.")
