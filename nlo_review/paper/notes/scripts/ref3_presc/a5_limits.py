# Ref.3 / Major 5(a).  Coincidence and large-distance behaviour of the on-shell kernel R of Group I Row III
# (R = P^i Q^{ij} K^j(y-z), Q = 2 Brace, see a2/a3), and an independent Lorentzian check near z = w and at large z.
import numpy as np
exec(open(__file__.replace('a5_limits.py', 'a2_rowIII_delta.py')).read().split("rng = np.random.default_rng(7)")[0])
k = 1.0
x, y, z, w = np.array([0.3, 0.1]), np.array([-0.8, 0.5]), np.array([0.9, -0.6]), np.array([0., 0.]); P = np.array([0.4, 1.1])
phis = np.linspace(0, 2 * np.pi, 9)[:-1] + 0.1
def scan(label, f, power, rs=(1e-2, 1e-3, 1e-4)):
    for r in rs:
        vals = [f(r * np.array([np.cos(p), np.sin(p)])) * r ** power for p in phis]
        print(f"{label:28s} r={r:8.0e}  R*r^{power} over 8 directions:", np.round(vals, 5), " mean", round(float(np.mean(vals)), 6))
print("1) across the line (x-z)^2 = (x-w)^2 (B = A(1+d)): R stays finite (double and single poles cancel)")
A = (x - w) @ (x - w); e = np.array([np.cos(2.0), np.sin(2.0)])
for d in (1e-2, 1e-4, 1e-6, -1e-6, -1e-4, -1e-2):
    zz = x - np.sqrt(A * (1 + d)) * e
    print(f"   d={d:+.0e}: R = {R_analytic(k, x, y, zz, w, P): .8f}")
print("2) z -> w: R -> P.K(y-w)/(z-w)^2, isotropic  [P.K(y-w) = %.6f]" % (P @ ((y - w) / ((y - w) @ (y - w)))))
scan("z->w", lambda d: R_analytic(k, x, y, w + d, w, P), 2)
print("3) other coincidences: bounded or 1/r (integrable)")
scan("z->x", lambda d: R_analytic(k, x, y, x + d, w, P), 0)
scan("w->x", lambda d: R_analytic(k, x, y, z, x + d, P), 0)
scan("z->y", lambda d: R_analytic(k, x, y, y + d, w, P), 1)
print("4) large distances: |z|^3 R bounded (no 1/z^2 tail), |w|^2 R bounded")
scan("z->oo  (r = 1/|z|)", lambda d: R_analytic(k, x, y, d / (d @ d), w, P), -3, rs=(1e-1, 1e-2, 1e-3))
scan("w->oo  (r = 1/|w|)", lambda d: R_analytic(k, x, y, z, d / (d @ d), P), -2, rs=(1e-1, 1e-2, 1e-3))
print("5) Lorentzian (direct p+ integral, D -> D +- i eps) vs analytic R at |z-w| = 1e-2 and at |z| = 20:")
for zz in (w + 1e-2 * np.array([np.cos(0.7), np.sin(0.7)]), 20 * np.array([np.cos(0.4), np.sin(0.4)])):
    print("   analytic", R_analytic(k, x, y, zz, w, P), " Lorentz eps=1e-5,1e-6,1e-7:",
          [R_lorentz(k, x, y, zz, w, P, 1e-9, 1e9, ee) for ee in (1e-5, 1e-6, 1e-7)])
print("6) compact form: Brace^{ij} = [cos(2phi) delta^{ij} + sin(2phi) eps^{ij}] / (2 (z-w)^2), phi = angle from x-w to x-z,")
print("   so R = (1/k+) P^i [cos2phi d^{ij} + sin2phi eps^{ij}] K^j(y-z) / (z-w)^2  -- manifestly regular on (x-z)^2=(x-w)^2")
def R_compact(k, x, y, z, w, P):
    a, b, c = x - w, x - z, z - w
    phi = np.arctan2(a[0] * b[1] - a[1] * b[0], a @ b)
    Kv = (y - z) / ((y - z) @ (y - z))
    M = np.array([[np.cos(2 * phi), np.sin(2 * phi)], [-np.sin(2 * phi), np.cos(2 * phi)]])
    return P @ M @ Kv / (c @ c) / k
rng = np.random.default_rng(11); worst = 0
for _ in range(20000):
    xx, yy, zz, ww, PP = (rng.normal(size=2) for _ in range(5))
    ra, rc = R_analytic(k, xx, yy, zz, ww, PP), R_compact(k, xx, yy, zz, ww, PP)
    worst = max(worst, abs(ra - rc) / (abs(rc) + 1e-300))
print("   max relative difference compact vs residue sum over 20000 random points:", worst)
for d in (1e-2, 1e-6, 1e-10, -1e-10, -1e-6, -1e-2):
    zz = x - np.sqrt(A * (1 + d)) * e
    print(f"   across the line, d={d:+.0e}: R_compact = {R_compact(k, x, y, zz, w, P): .10f}")
