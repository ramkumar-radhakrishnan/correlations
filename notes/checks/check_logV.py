"""Numerical check: does the y'=y point make the integral diverge as V -> infinity?

Toy model of the pasted expression, keeping exactly the structures that matter:
  y = 0  (so C = y'),  A = x - y fixed,  x' fixed OUTSIDE the integration disk
  (the y'=x' singularity is separately protected by [U(y')-U(x')] in the real
  expression, so we keep it out of the domain),
  f(C) = smooth bounded stand-in for the U,rho factors (nonzero at C=0),
  FV(a,V) = int_{1/V}^1 dxi (1-xi)/xi * exp(-i a/xi),  a = k.C   (exact, via Si/Ci)

Checks:
  1. FV(0,V) = log V - 1 + 1/V              -> diverges   (the y'=y slice)
  2. I(V) = int_{|C|<0.5} d2C br(C) f(C) FV(k.C, V), k=(1,0)  -> converges
  3. same with k=0 (phase removed)           -> diverges like log V
  4. angular integral of the bracket at radius r -> 0 as r -> 0
"""
import numpy as np
from scipy.special import sici
from scipy.integrate import quad

A  = np.array([1.0, 0.2])
X  = np.array([-0.6, 0.8])          # x' ; |X| = 1.0 > disk radius 0.5
sh = np.array([0.13, 0.21])         # center of the smooth factor
R  = 0.5                            # disk radius

def br(cx, cy):
    """the square bracket:  A.B/(A^2 B^2) - 2 (A.C)(B.C)/(A^2 B^2 C^2),  B = x'-C"""
    bx, by = X[0]-cx, X[1]-cy
    A2 = A @ A
    B2 = bx*bx + by*by
    C2 = cx*cx + cy*cy
    if C2 == 0.0:
        return 0.0   # bracket is bounded but direction-ambiguous at C=0; measure zero
    t1 = (A[0]*bx + A[1]*by) / (A2*B2)
    t2 = 2.0*(A[0]*cx + A[1]*cy)*(bx*cx + by*cy) / (A2*B2*C2)
    return t1 - t2

def f(cx, cy):
    return np.exp(-((cx-sh[0])**2 + (cy-sh[1])**2))

def FV(a, V):
    """int_{1/V}^1 dxi (1-xi)/xi exp(-i a / xi)   exactly, via t = 1/xi:
       int_1^V dt (1/t - 1/t^2) e^{-iat} = g1 - g2
       g1 = Ci(|a|V)-Ci(|a|) - i sgn(a) (Si(|a|V)-Si(|a|))
       g2 = e^{-ia} - e^{-iaV}/V - i a g1
    """
    if abs(a) < 1e-14:
        return np.log(V) - 1.0 + 1.0/V
    s1, c1 = sici(abs(a))
    sV, cV = sici(abs(a)*V)
    g1 = (cV - c1) - 1j*np.sign(a)*(sV - s1)
    g2 = np.exp(-1j*a) - np.exp(-1j*a*V)/V - 1j*a*g1
    return g1 - g2

# sanity check of the closed form against direct quadrature of the xi integral
for a, V in [(0.3, 100.0), (-1.7, 1e4), (0.01, 1e6)]:
    re = quad(lambda x: (1-x)/x*np.cos(a/x), 1/V, 1, limit=4000)[0]
    im = quad(lambda x: -(1-x)/x*np.sin(a/x), 1/V, 1, limit=4000)[0]
    cf = FV(a, V)
    assert abs(cf - (re+1j*im)) < 1e-6, (a, V, cf, re+1j*im)
print("closed form FV verified against direct xi-quadrature\n")

def H(cx):
    """inner integral over C_y along the chord of the disk (smooth, fast)"""
    ylim = np.sqrt(max(R*R - cx*cx, 0.0))
    if ylim == 0.0:
        return 0.0
    val, _ = quad(lambda cy: br(cx, cy)*f(cx, cy), -ylim, ylim, limit=400)
    return val

# cache H on a fine grid once (it is V-independent) and interpolate
grid = np.linspace(-R, R, 4001)
Hg = np.array([H(c) for c in grid])
from scipy.interpolate import CubicSpline
Hs = CubicSpline(grid, Hg)

def I_of_V(V, k):
    if k == 0.0:
        c = FV(0.0, V)
        tot = quad(Hs, -R, R, limit=400)[0]
        return c*tot
    fre = quad(lambda a: (FV(k*a, V)).real * Hs(a), -R, R,
               limit=20000, points=[0.0], epsabs=1e-9, epsrel=1e-9)[0]
    fim = quad(lambda a: (FV(k*a, V)).imag * Hs(a), -R, R,
               limit=20000, points=[0.0], epsabs=1e-9, epsrel=1e-9)[0]
    return fre + 1j*fim

print("1) pointwise value on the slice y'=y   (= log V - 1 + 1/V):")
for V in [1e2, 1e4, 1e6, 1e8]:
    print(f"   V={V:8.0e}   FV(0,V) = {FV(0.0,V):.6f}")

print("\n2) FULL integral, k=(1,0)  (phase active):")
for V in [1e2, 1e3, 1e4, 1e5, 1e6, 1e8]:
    r = I_of_V(V, 1.0)
    print(f"   V={V:8.0e}   I(V) = {r.real:+.8f} {r.imag:+.8f}i")

print("\n3) FULL integral, k=0  (phase removed):")
for V in [1e2, 1e4, 1e6, 1e8]:
    r = I_of_V(V, 0.0)
    print(f"   V={V:8.0e}   I(V) = {r:+.8f}")

print("\n4) angular integral of the bracket at radius r (-> 0 as r -> 0):")
for r0 in [0.2, 0.02, 0.002, 0.0002]:
    val = quad(lambda t: br(r0*np.cos(t), r0*np.sin(t)), 0, 2*np.pi, limit=400)[0]
    print(f"   r={r0:7.4f}   loop integral = {val:+.6e}")
