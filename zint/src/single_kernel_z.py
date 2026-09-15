"""int d^2z of the BK real kernel alone:  (y'-z).(x-z)/[(y'-z)^2 (x-z)^2].
Is there a UV (short-distance) divergence?"""
import numpy as np
from scipy.integrate import quad
g = 0.5772156649015329
yp = np.array([-0.5, 0.8]); xx = np.array([0.3, -0.2])
d = np.hypot(*(yp-xx))
def K(z):
    A = yp - z; B = xx - z
    return np.sum(A*B,-1)/(np.sum(A*A,-1)*np.sum(B*B,-1))

th = np.linspace(0, 2*np.pi, 8192, endpoint=False)
u = np.stack([np.cos(th), np.sin(th)], -1)

print("A. SHORT DISTANCE.  small disc about each singular point")
print("   near z = y' :  numerator -> (y'-z).(x-y') = O(rho), denominator = rho^2 (x-y')^2")
print("   so the integrand is O(1/rho) and, the leading piece being odd, the disc")
print("   integral should go as eps^2:")
for nm, c in (("z -> y'", yp), ("z -> x", xx)):
    prev = None
    for eps in (1e-1, 1e-2, 1e-3):
        f = lambda r: r*np.mean(K(c + u*r))*2*np.pi
        ed = np.concatenate(([0.0], np.geomspace(eps*1e-6, eps, 40)))
        v = sum(quad(f, p, q, limit=80)[0] for p, q in zip(ed[:-1], ed[1:]))
        rat = f"   ratio {abs(v)/abs(prev):.4f}" if prev is not None else ""
        print(f"   {nm:8s} eps={eps:.0e}:  |F| = {abs(v):.4e}{rat}")
        prev = v
print("   -> ~0.01 per decade, F ~ eps^2.   NO ULTRAVIOLET DIVERGENCE.")
print()

print("B. EXACT CLOSED FORM.  In complex coordinates K = Re 1/[(zbar-ybar')(z-x)], and")
print("   int_{|z|<R} d^2z / [(zbar-abar)(z-b)] = pi log(R^2/|a-b|^2), so")
print("      int_{|z|<R} d^2z K  =  pi log( R^2 / (y'-x)^2 )   exactly, no extra constant.")
for R in (20.0, 100.0, 500.0, 2000.0):
    f = lambda r: r*np.mean(K(u*r))*2*np.pi
    ed = np.concatenate(([1e-9], np.geomspace(1e-8, R, 600)))
    v = sum(quad(f, p, q, limit=120)[0] for p, q in zip(ed[:-1], ed[1:]))
    print(f"   R={R:8.1f}:  quad {v:+.6f}    pi log(R^2/(y'-x)^2) = {np.pi*np.log(R**2/d**2):+.6f}")
print(f"   |y'-x| = {d:.6f}")
print()

print("C. READING IT OFF")
print("   * the ONLY short-distance scale in the answer is |y'-x|, an EXTERNAL separation.")
print("     No regulator appears -> no UV divergence, nothing to renormalise.")
print("   * the R dependence is 2 pi log R : the divergence is INFRARED, coefficient 2 pi,")
print("     matching the slope +6.28319 measured earlier for this structure.")
print("   * the log(1/(y'-x)^2) is the y' -> x singularity.  At exactly y' = x the kernel")
print("     degenerates to 1/(x-z)^2 and the z-integral IS log divergent at z -> x;")
print("     but that is a set of measure zero and log^n is integrable against d^2(y'-x):")
for e in (1e-2, 1e-4, 1e-6):
    print(f"      |y'-x| = {e:.0e} :  pi log(R^2/(y'-x)^2) at R=100 = {np.pi*np.log(100**2/e**2):+10.4f}"
          f"   (grows only logarithmically)")
