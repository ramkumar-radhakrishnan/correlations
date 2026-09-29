"""Row V 2a, p+ -> 0 plus prescription: may the denominator be set to its p+ = 0 value?

f(p+) = Phi(p+) B4 / D(p+),  D(p+) = p+ (y-x)^2 + (k+ - p+)(y-z)^2 ,  Phi = e^{-ik.(w'-z)} e^{-ik.(z-x) p+/k+}.
Exact identity (k+ = 1):
  -Int_0^1 dp [1/p]_+ f  =  -Int_0^1 dp [1/p]_+ Phi B4/D(0)  +  Int_0^1 dp Phi B4 [(y-x)^2-(y-z)^2] / ((y-z)^2 D(p))
"""
import numpy as np
from scipy import integrate

sq = lambda v: float(v @ v)
def B4c(x, y, z, xp, wp):   # bracket of T4 without -(xibar/xi), contracted with P^i r^m/(P^2 r^2)
    Pv = (xp-wp)/sq(xp-wp); K = (z-x)/sq(z-x)
    return Pv@(np.outer(y-z, x-z)/sq(x-z) + np.outer(y-z, y-x)/(2*sq(y-x)))@K
def phase(k, x, z, wp, p): return np.exp(-1j*(k@(wp-z)) - 1j*p*(k@(z-x)))
def den(x, y, z, p): return p*sq(y-x) + (1-p)*sq(y-z)
def cq(f, a, b, pts=None):
    re = integrate.quad(lambda t: f(t).real, a, b, points=pts, limit=1000, epsabs=1e-13, epsrel=1e-12)[0]
    im = integrate.quad(lambda t: f(t).imag, a, b, points=pts, limit=1000, epsabs=1e-13, epsrel=1e-12)[0]
    return re + 1j*im

rng = np.random.default_rng(11)
print("1. EXACT SPLIT AT FIXED POSITIONS (k+ = 1)")
print("config   exact [1/p]_+ term              D -> D(0) only                  D(0) + extra term               extra term")
for trial in range(4):
    x, y, z, xp, wp = rng.normal(size=(5, 2)); k = rng.normal(size=2)
    b4 = B4c(x, y, z, xp, wp)
    f = lambda p: phase(k, x, z, wp, p)*b4/den(x, y, z, p)
    exact = -cq(lambda p: (f(p) - f(0))/p, 0, 1)
    D0 = den(x, y, z, 0.0)
    approx = -cq(lambda p: (phase(k, x, z, wp, p) - phase(k, x, z, wp, 0))*b4/D0/p, 0, 1)
    extra = cq(lambda p: phase(k, x, z, wp, p)*b4*(sq(y-x) - sq(y-z))/(sq(y-z)*den(x, y, z, p)), 0, 1)
    print(f"  {trial}   {exact:.8f}   {approx:.8f}   {approx+extra:.8f}   {extra:.6f}")
print("-> 'D -> D(0) only' misses a finite piece; adding the extra term restores the exact value.")

print("\n2. y -> z (u = y - z): |u|^2 x angular average of the p+-integrated extra term (-> 0 means integrable in d^2y)")
x = np.array([0.3, -0.2]); z = np.array([1.0, 0.4]); xp = np.array([-0.5, 0.9]); wp = np.array([0.7, 1.3]); k = np.array([0.8, -0.3])
for u in (1e-1, 1e-2, 1e-3):
    vals = []
    for ph in np.linspace(0, 2*np.pi, 64, endpoint=False):
        y = z + u*np.array([np.cos(ph), np.sin(ph)]); b4 = B4c(x, y, z, xp, wp)
        vals.append(cq(lambda p: phase(k, x, z, wp, p)*b4*(sq(y-x) - sq(y-z))/(sq(y-z)*den(x, y, z, p)), 0, 1,
                       pts=[u*u, 10*u*u]))
    print(f"   |u| = {u:.0e}:  |u|^2 <extra> = {abs(u*u*np.mean(vals)):.3e}")
