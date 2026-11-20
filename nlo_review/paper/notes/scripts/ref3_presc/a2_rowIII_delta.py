# Ref.3 / Major 5(a).  Group I Row III (LO^dagger x Bbar_3^dagger A): the extra term of an i*eta*0 prescription.
# Row III before the p+ integral (Group_I_checked.tex, Row III, first display), in units of
#   (1/(2pi)^3) (i g^4 f^{abd}/8pi^4) (dp+/2pi)  x  colour  x  sum_n T_n(p+),
# with A=(x-w)^2, B=(x-z)^2, D = p+ B - k+ A, Y=y-z, a=x-w, b=x-z, c=z-w, P=K(x'-w'):
#   T1 =  (P.Y/Y^2) (p+k)/(k-p)^2 / D                                       p in (L,k-L) u (k+L,V)
#   T3 = -(P.a)(b.Y)/(A B Y^2) (kA+pB)/(p k D)                              p in (L,k-L) u (k+L,V)
#   T4 = (2/k) Ten[V'=k a - p b] /(Y^2 c^2 (k-p) D)                          p in (k+L,V)
#   T5 = 2 p Ten[V'] /(Y^2 |k a - p b|^2 (k-p) D)                            p in (L,k-L)
#   Ten[V'] = P^i Y^j c^k V'^m [d_ik d_jm - (k/p) d_jk d_im - k/(k-p) d_ij d_km]
# 1/(D + i eta 0) = P(1/D) - i eta pi delta(D): the extra term is (-i eta pi) * R_theta with
#   R_theta = sum_n [D T_n]_{p=p0} / B * theta(p0 in range_n),   p0 = k A / B.
# Analytic form (this note):  R = P.Rvec,
#   Rvec = Y (A+B)/(k Y^2 (B-A)^2) - (2/k) a (b.Y)/(A B Y^2)
#        + (2/k) [ c (Y.V) - (B/A)(Y.c) V - B/(B-A) Y (c.V) ] / (Y^2 c^2 (B-A)),   V = a - (A/B) b,
# the same expression on both sides of A=B (T4 for A>B, T5 for A<B).
import numpy as np
from scipy.integrate import quad

def terms_times_D(p, k, x, y, z, w, P):
    """D*T_n(p) for n=1,3,4,5 (zero outside the range of the branch is handled by caller)."""
    a, b, c, Y = x - w, x - z, z - w, y - z
    A, B, C, Y2 = a @ a, b @ b, c @ c, Y @ Y
    t1 = (P @ Y) / Y2 * (p + k) / (k - p) ** 2
    t3 = -(P @ a) * (b @ Y) / (A * B * Y2) * (k * A + p * B) / (p * k)
    Vp = k * a - p * b
    ten = (P @ c) * (Y @ Vp) - (k / p) * (Y @ c) * (P @ Vp) - k / (k - p) * (P @ Y) * (c @ Vp)
    t4 = (2.0 / k) * ten / (Y2 * C * (k - p))
    t5 = 2.0 * p * ten / (Y2 * (Vp @ Vp) * (k - p))
    return t1, t3, t4, t5

def R_analytic(k, x, y, z, w, P):
    a, b, c, Y = x - w, x - z, z - w, y - z
    A, B, C, Y2 = a @ a, b @ b, c @ c, Y @ Y
    V = a - (A / B) * b
    R = Y * (A + B) / (k * Y2 * (B - A) ** 2) - (2.0 / k) * a * (b @ Y) / (A * B * Y2) \
        + (2.0 / k) * (c * (Y @ V) - (B / A) * (Y @ c) * V - B / (B - A) * Y * (c @ V)) / (Y2 * C * (B - A))
    return P @ R

def R_lorentz(k, x, y, z, w, P, lam, vee, eps):
    """(1/pi) int dp sum_n [D T_n] eps/(D^2+eps^2) over the ranges: -> R_theta for eps->0.
    The interval is split at p0 +- (1,10,100,1000) Lorentz widths and on a log grid, so that quad resolves the peak."""
    A, B = (x - w) @ (x - w), (x - z) @ (x - z)
    p0 = k * A / B
    wid = eps / B
    def f(p, which):
        D = p * B - k * A
        t1, t3, t4, t5 = terms_times_D(p, k, x, y, z, w, P)
        val = {'13': t1 + t3, '4': t4, '5': t5}[which]
        return val * eps / (D * D + eps * eps) / np.pi
    tot = 0.0
    for lo, hi, which in [(lam, k - lam, '13'), (k + lam, vee, '13'), (k + lam, vee, '4'), (lam, k - lam, '5')]:
        br = set(np.geomspace(lo, hi, 60))
        for m in (0, 1, 10, 100, 1000):
            for sgn in (-1, 1):
                q = p0 + sgn * m * wid
                if lo < q < hi:
                    br.add(q)
        br = sorted(br | {lo, hi})
        for u, v in zip(br[:-1], br[1:]):
            val, err = quad(f, u, v, args=(which,), limit=200, epsabs=1e-14, epsrel=1e-12)
            tot += val
    return tot

rng = np.random.default_rng(7)
k = 1.0
print("check of the residue algebra: R_theta(eps) from the Lorentzian vs analytic R")
for trial in range(8):
    x, y, z, w = (rng.normal(size=2) for _ in range(4))
    P = rng.normal(size=2)
    A, B = (x - w) @ (x - w), (x - z) @ (x - z)
    Ra = R_analytic(k, x, y, z, w, P)
    rows = []
    for eps in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
        rows.append(R_lorentz(k, x, y, z, w, P, 1e-7, 1e7, eps))
    print(f"  A/B={A/B:7.4f}  analytic R={Ra: .10f}  Lorentz eps=1e-3..1e-7: " + "  ".join(f"{r: .10f}" for r in rows))

# R does not vanish: magnitude relative to its separate pieces on a sample
x, y, z, w = np.array([0.3, 0.1]), np.array([-0.8, 0.5]), np.array([0.9, -0.6]), np.array([0.0, 0.0])
P = np.array([0.4, 1.1])
print("sample point: R =", R_analytic(k, x, y, z, w, P))
