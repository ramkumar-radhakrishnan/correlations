"""Direct test: is the FIRST term's z-integral convergent?

Compute the truncated z-integral  F(R) = int_{|z|<R} d^2z f(z)  for each structure and
look at dF/dlogR.  A convergent integral has dF/dlogR -> 0; a log divergence has a
constant nonzero slope.

  S1 = (y'-z).a/(y'-z)^2 * (x-z).b/(x-z)^2        on  1/(p+ + k+)     keeps full phase
  S2 = (y'-z).(x-z)/[(y'-z)^2 (x-z)^2]            on  [1/p+]_+        phase -> phase - 1
  S3 = (y'-z).b/(y'-z)^2 * (x-z).a/(x-z)^2        on  1/k+            keeps full phase
  phase = exp(-i xi k.(y'-z))
"""
import numpy as np
from scipy.integrate import quad

rng = np.random.default_rng(3)
yp, xx, xpp, yy, ww = rng.normal(size=(5, 2))
a, b = xpp - yp, yy - ww
kk = np.array([0.9, 0.6])

def parts(z, xi):
    A = yp - z; B = xx - z
    A2 = np.sum(A*A, -1); B2 = np.sum(B*B, -1)
    ph = np.exp(-1j*xi*(A @ kk))
    S1 = (A @ a)/A2 * (B @ b)/B2
    S2 = np.sum(A*B, -1)/(A2*B2)
    S3 = (A @ b)/A2 * (B @ a)/B2
    return {"S1 x ph": ph*S1, "S2 x ph": ph*S2, "S3 x ph": ph*S3,
            "S2 x (ph-1)": (ph-1)*S2, "S2 bare (= g(0))": S2+0j}

def cumulative(xi, R, ntheta=6001):
    """2pi int_0^R rho drho <f>_theta, for every structure at once."""
    th = np.linspace(0, 2*np.pi, ntheta, endpoint=False)
    u = np.stack([np.cos(th), np.sin(th)], -1)
    keys = list(parts(np.zeros((1, 2)), xi).keys())
    lamosc = np.pi/(xi*np.hypot(*kk))                     # half oscillation length
    edges = [0.0, 0.5, 1.0, 2.0, 4.0]
    edges += list(np.arange(1, 400)*lamosc + 4.0)
    edges = [e for e in edges if e < R] + [R]
    out = {k: 0j for k in keys}
    for k in keys:
        def g(rho, k=k):
            return rho*np.mean(parts(u*rho, xi)[k])*2*np.pi
        tot = 0j
        for p, q in zip(edges[:-1], edges[1:]):
            tot += (quad(lambda r: g(r).real, p, q, limit=200)[0]
                    + 1j*quad(lambda r: g(r).imag, p, q, limit=200)[0])
        out[k] = tot
    return out

for xi in (0.3, 0.05):
    print(f"=== xi = {xi} ===")
    Rs = [1e2, 1e3, 1e4, 1e5]
    vals = [cumulative(xi, R) for R in Rs]
    keys = list(vals[0].keys())
    print(f"   {'structure':<20}" + "".join(f"{'F(1e%d)'%int(np.log10(R)):>26}" for R in Rs))
    for k in keys:
        print(f"   {k:<20}" + "".join(f"{v[k].real:+13.5f}{v[k].imag:+12.5f}i" for v in vals))
    print(f"   {'dF/dlogR (real)':<20}" + " "*26, end="")
    print()
    for k in keys:
        sl = [ (vals[i+1][k]-vals[i][k]).real/np.log(10) for i in range(len(Rs)-1) ]
        print(f"   {k:<20} dF/dlogR = " + "".join(f"{s:+12.5f}" for s in sl))
    print(f"   (-2 pi = {-2*np.pi:+.5f})")
    print()
