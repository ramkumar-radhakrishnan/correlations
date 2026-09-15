"""Why Ei(-i kappa lambda) is undefined at Lambda -> 0 OR at y'-z -> 0,
and what that does (and does not) mean for the transverse integral.

kappa = k.(y'-z),  lam = Lambda/k+,  Xi = (V-k+)/k+ .
p+ integral of the pole term:  -[ Ei(-i kappa Xi) - Ei(-i kappa lam) ].

Fast route:  Ei(-i u) = -E_1(i u) - i pi sgn(u)   (checked against mpmath to 1e-10);
the branch term cancels in the difference since Xi, lam > 0.
"""
import numpy as np
from scipy.special import exp1
from scipy.integrate import quad
g = 0.5772156649015329

def Ei_mi(u):                       # Ei(-i u), u real
    u = np.asarray(u, float)
    return -exp1(1j*u) - 1j*np.pi*np.sign(u)

print("A. Ei(-i kappa lam) DEPENDS ONLY ON THE PRODUCT kappa*lam")
print("   Ei(-i u) = gamma + log(i u) + O(u):  same value for very different (kappa, lam)")
for kap, lam in [(1.0, 1e-6), (1e-3, 1e-3), (1e-6, 1.0), (2.0, 5e-7)]:
    u = kap*lam
    print(f"   kappa={kap:<8g} lam={lam:<8g} product={u:.2e}:  Ei = {Ei_mi(u):+.8f}"
          f"   gamma+log(i u) = {g+np.log(1j*u):+.8f}")
print("   -> the blow-up is ONE logarithm, log|kappa| + log(lam).")
print("      log(lam)   = -l_k            : the RAPIDITY divergence (real; this is what you subtract)")
print("      log|kappa| = log|k.(y'-z)|   : a log SINGULARITY of the integrand at the point y'=z")
print()

rng = np.random.default_rng(3)
zz, xx, xpp, yy, ww = rng.normal(size=(5,2))
kk = np.array([0.9,0.6]); d, e = xx-zz, yy-ww
lam, Xi = 1e-6, 37.0

def integrand(r):
    """the p+-integrated pole structure as a function of r = y'-z"""
    r2 = np.sum(r*r,-1); C = xpp-zz-r; C2 = np.sum(C*C,-1)
    T2 = (r@d)/(r2*(d@d)) * (C@e)/(C2*(e@e))
    kap = r@kk
    return -(Ei_mi(kap*Xi) - Ei_mi(kap*lam)) * T2

def ring(rho, nth=2048):
    th = np.linspace(0,2*np.pi,nth,endpoint=False)
    u = np.stack([np.cos(th),np.sin(th)],-1)
    return rho*np.mean(integrand(u*rho))*2*np.pi

def cum(a, b, n=80):
    ed = np.geomspace(a, b, n)
    return sum(quad(lambda x: ring(x).real, p, q, limit=60)[0]
               + 1j*quad(lambda x: ring(x).imag, p, q, limit=60)[0]
               for p, q in zip(ed[:-1], ed[1:]))

print("B. SHORT DISTANCE  y' -> z :  the log is INTEGRABLE. No transverse divergence.")
print("   integrand ~ (1/rho) x log rho ;  int rho drho (log rho)/rho = int drho log rho -> 0")
prev = None
for eps in (1e-1, 1e-2, 1e-3, 1e-4):
    v = cum(1e-9, eps)
    r = f"    ratio {abs(v)/abs(prev):.4f}" if prev is not None else ""
    print(f"   |int_{{|y'-z| < {eps:.0e}}} d^2y'| = {abs(v):.4e}{r}")
    prev = v
print("   ~0.01 per decade => F ~ eps^2.  Compare: log(y'-z)^2 is 'undefined' at y'=z too,")
print("   yet int d^2(y'-z) log(y'-z)^2 is perfectly finite.  Same thing here.")
print()

print("C. LARGE DISTANCE  |y'-z| -> inf :  HERE the same log does hurt.")
print("   Ei(-i kappa Xi) -> -i pi sgn(kappa), so the bracket -> -l_k + gamma + log(i k.(y'-z)),")
print("   growing like log|r| while the kernel structure falls only as 1/r^2:")
th = np.linspace(0,2*np.pi,8192,endpoint=False)
u = np.stack([np.cos(th),np.sin(th)],-1)
for rho in (1e1,1e2,1e3,1e4):
    r = u*rho; r2=np.sum(r*r,-1); C=xpp-zz-r; C2=np.sum(C*C,-1)
    T2 = (r@d)/(r2*(d@d)) * (C@e)/(C2*(e@e))
    print(f"   rho={rho:8.0e}:  rho^2 x <T2> = {np.mean(T2)*rho**2:+.8f}"
          f"    predicted -d.e/(2 d^2 e^2) = {-(d@e)/(2*(d@d)*(e@e)):+.8f}")
print()
print("   the truncated integral and its slope:")
prev = None; R0 = 4.0
for R in (1e2, 1e3, 1e4, 1e5):
    v = cum(1e-9, R)
    s = f"    dF/dlogR = {(v-prev).real/np.log(10):+10.4f}" if prev is not None else ""
    print(f"   F(R={R:.0e}) = {v.real:+12.5f} {v.imag:+12.5f}i{s}")
    prev = v
print("   the slope GROWS (it is -l_k + log R), i.e. F ~ -l_k log R + (1/2) log^2 R.")
print("   That second piece is exactly the double logarithm l_k^2.")
