# integral identities used to rewrite UV brackets / logs as z-integrals
import sympy as sp, numpy as np
from scipy import integrate
e, g, L = sp.symbols("epsilon gamma L")   # L = log(r^2 mu^2)
def half_I(d):
    # (1/2) mu^{d-2}... int d^dz r^2/((z^2)(z-r)^2) = (1/2) pi^{d/2} Gamma(d/2-1)^2 Gamma(2-d/2)/Gamma(d-2) (r^2)^{d/2-1}
    return sp.Rational(1, 2) * sp.pi ** (d / 2) * sp.gamma(d / 2 - 1) ** 2 * sp.gamma(2 - d / 2) / sp.gamma(d - 2)
for sgn, nm in ((+1, "d=2+eps"), (-1, "d=2-eps")):
    d = 2 + sgn * e
    expr = half_I(d) * sp.exp(sgn * e / 2 * L)          # (r^2 mu^2)^{(d-2)/2}
    ser = sp.series(expr / sp.pi, e, 0, 1).removeO()
    ser = sp.simplify(ser.subs(sp.EulerGamma, g))
    # MS-bar: mu^2 = mubar^2 e^{gamma}/(4 pi) for d=2+eps (the scale enters as (mu^2)^{+eps/2}); inverse for d=2-eps
    Lb = sp.symbols("Lbar")   # log(r^2 mubar^2/4)
    print(nm, ": (1/pi) int d^dz (1/2) r^2/((w-z)^2 (x-z)^2) =", sp.expand(ser), "  [L=log(r^2 mu^2)]")
    print("      with mu^2 = mubar^2 e^gamma/(4pi):", sp.simplify(sp.expand(ser).subs(L, Lb + g + sp.log(4) - sp.log(4 * sp.pi))), " [Lbar=log(r^2 mubar^2/4)]")
# numerical: cutoff version  int_{|z-w|>a} [1/(w-z)^2 - K(x-z).K(w-z)] = pi log(r^2/a^2)  and  log identity
rng = np.random.default_rng(3); x, w, y = rng.normal(size=(3, 2))
def K(v): return v / (v @ v)
def polar(f, c, a, R=np.inf):
    # integrate f over |z-c|>a in polar coordinates around c
    val, err = integrate.dblquad(lambda r, th: r * f(c + r * np.array([np.cos(th), np.sin(th)])), 0, 2 * np.pi, a, R, epsabs=1e-7, epsrel=1e-7)
    return val
a = 1e-3; r2 = (x - w) @ (x - w)
f1 = lambda z: 1 / ((w - z) @ (w - z)) - K(x - z) @ K(w - z)
# split domain: inner disc up to R0 and outer
v = polar(f1, w, a, 50.0)
print("cutoff: int_{a<|z-w|<50}[1/(w-z)^2-K(x-z).K(w-z)] = %.5f ; pi log(r^2/a^2) = %.5f (tail beyond 50 ~ O(r/50))" % (v, np.pi * np.log(r2 / a ** 2)))
f2 = lambda z: (K(x - z) - K(w - z)) @ K(y - z)
v2 = polar(f2, y, 0.0, 200.0)
print("log identity: int [K(x-z)-K(w-z)].K(y-z) = %.5f ; pi log((w-y)^2/(x-y)^2) = %.5f" % (v2, np.pi * np.log(((w - y) @ (w - y)) / ((x - y) @ (x - y)))))
