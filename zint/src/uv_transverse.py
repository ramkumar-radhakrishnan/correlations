"""Is there a transverse ULTRAVIOLET (short-distance) divergence anywhere?

Integrand of the row:
   (y'-z)^m (x'-y')^k' (y-w)^i (x-z)^j / [ (y'-z)^2 (x'-y')^2 (y-w)^2 (x-z)^2 ]
Every integration variable w, y', z, y, x, x' is a coincidence point of exactly ONE
kernel, so the test is whether each small-disc integral vanishes with the radius.
"""
import numpy as np
from scipy.integrate import quad

rng = np.random.default_rng(3)
yp0, xx, xpp, yy, ww = rng.normal(size=(5, 2))
zz0 = rng.normal(size=2)
kk = np.array([0.9, 0.6])
a, b = xpp - yp0, yy - ww
xi = 0.3

# ---- the z-integrand, five variants (as in the previous note) -------------
def zf(z):
    A = yp0 - z; B = xx - z
    A2 = np.sum(A*A, -1); B2 = np.sum(B*B, -1)
    ph = np.exp(-1j*xi*(A @ kk))
    S1 = (A @ a)/A2 * (B @ b)/B2
    S2 = np.sum(A*B, -1)/(A2*B2)
    S3 = (A @ b)/A2 * (B @ a)/B2
    return {"S1 x ph": ph*S1, "S2 x ph": ph*S2, "S3 x ph": ph*S3,
            "S2 x (ph-1)": (ph-1)*S2, "S2 bare (g(0))": S2+0j}

# ---- the y'-integrand, the three index structures -------------------------
def yf(yq):
    A = yq - zz0; C = xpp - yq
    A2 = np.sum(A*A, -1); C2 = np.sum(C*C, -1)
    ph = np.exp(-1j*xi*(-A @ kk))
    d, e = xx - zz0, yy - ww
    T1 = np.sum(A*C, -1)/(A2*C2)                      # d_k'm d_ij : dipole kernel in y'
    T2 = (A @ d)/A2 * (C @ e)/C2                      # d_jm d_ik'
    T3 = (A @ e)/A2 * (C @ d)/C2                      # d_im d_jk'
    return {"T1 x ph": ph*T1, "T2 x ph": ph*T2, "T3 x ph": ph*T3,
            "T2 x (ph-1)": (ph-1)*T2, "T2 bare (g(0))": T2+0j}

def disc(f, centre, eps, ntheta=4001):
    """int_{|v-centre|<eps} d^2v f(v)  -- should vanish like eps if f ~ 1/rho."""
    th = np.linspace(0, 2*np.pi, ntheta, endpoint=False)
    u = np.stack([np.cos(th), np.sin(th)], -1)
    keys = list(f(centre + u[:1]*1e-3).keys())
    out = {}
    for k in keys:
        def g(rho, k=k):
            return rho*np.mean(f(centre + u*rho)[k])*2*np.pi
        tot = 0j
        edges = np.concatenate(([0.0], np.geomspace(eps*1e-6, eps, 40)))
        for p, q in zip(edges[:-1], edges[1:]):
            tot += (quad(lambda r: g(r).real, p, q, limit=100)[0]
                    + 1j*quad(lambda r: g(r).imag, p, q, limit=100)[0])
        out[k] = tot
    return out

print("A. z-INTEGRAL, small disc about the two singular points")
for nm, c in (("z -> y'", yp0), ("z -> x", xx)):
    print(f"   {nm}")
    res = [disc(zf, c, e) for e in (1e-1, 1e-2, 1e-3)]
    for k in res[0]:
        vals = "".join(f"{abs(r[k]):13.3e}" for r in res)
        rat = abs(res[1][k])/abs(res[0][k]) if abs(res[0][k]) else float('nan')
        rat2 = abs(res[2][k])/abs(res[1][k]) if abs(res[1][k]) else float('nan')
        print(f"      {k:<18} |F(eps)| for eps=1e-1,1e-2,1e-3:{vals}   ratios {rat:.3f} {rat2:.3f}")
print("   (ratio 0.1 per decade  <=>  F ~ eps  <=>  integrand ~ 1/rho : INTEGRABLE)")
print()

print("B. y'-INTEGRAL, small disc about the two singular points")
for nm, c in (("y' -> z", zz0), ("y' -> x'", xpp)):
    print(f"   {nm}")
    res = [disc(yf, c, e) for e in (1e-1, 1e-2, 1e-3)]
    for k in res[0]:
        vals = "".join(f"{abs(r[k]):13.3e}" for r in res)
        rat = abs(res[1][k])/abs(res[0][k]) if abs(res[0][k]) else float('nan')
        rat2 = abs(res[2][k])/abs(res[1][k]) if abs(res[1][k]) else float('nan')
        print(f"      {k:<18} |F(eps)| for eps=1e-1,1e-2,1e-3:{vals}   ratios {rat:.3f} {rat2:.3f}")
print()

print("C. THE SUBTRACTION MAKES SHORT DISTANCE BETTER, NEVER WORSE")
print("   near z -> y' the phase is analytic:  ph - 1 = -i xi k.(y'-z) + O(rho^2),")
print("   so S2 x (ph-1) ~ (1/rho) x rho = O(1) : the 1/rho is REMOVED altogether.")
for e in (1e-2, 1e-3, 1e-4):
    r1 = disc(zf, yp0, e)
    print(f"      eps={e:.0e}:  |S2 x ph| = {abs(r1['S2 x ph']):.3e}    "
          f"|S2 x (ph-1)| = {abs(r1['S2 x (ph-1)']):.3e}   (ratio {abs(r1['S2 x (ph-1)'])/abs(r1['S2 x ph']):.2e})")
print("   -> the + prescription cannot create a transverse UV divergence.")
print()

print("D. DOES THE COLOUR STRUCTURE KILL THE LARGE-|z| LOG?   (U(z) -> 1 as |z| -> inf)")
print("   first bracket  -> delta^{bd'} [ U^{e'c'}(y') - U^{c'e'}(x') ] rho^{e'}(x')")
print("   second bracket -> U^{ad}(y) rho^d(y) [ U^{be}(x) rho^e(x) - rho^b(x) ]")
print("   both are z-INDEPENDENT, not zero.  So the coefficient of the large-|z|")
print("   logarithm survives; the colour structure does not regulate it.")
