# (C7): instantaneous Group III terms.  Toy with the structure of the subtraction bracket:
#   h^m(r) = K^m(r) [f(z + r) - f(z)],  r = y' - z,  f(y) = exp(-|y - c|^2)   (h ~ -f(z) K^m(r) at large r: degree -1)
#   S = int_0^inf dt int d^2r e^{i(1+t) k.r} h^m(r)          (k = (1,0); y' Fourier integral at fixed p+ first)
# A: y' (2D Fourier) at fixed t first, then t.     B: t first -> pi delta(r_par) + i P/r_par, then the perpendicular
# component u symmetric, then r_par (principal value).  C: as B with an asymmetric cut u in (-U, a U).
import numpy as np
from scipy.integrate import quad
z = np.array([0.0, 0.0]); c = np.array([0.7, -0.4]); fz = np.exp(-(z - c) @ (z - c))
def f(y1, y2): return np.exp(-((y1 - c[0]) ** 2 + (y2 - c[1]) ** 2))
def h(m, r1, r2):
    r2s = r1 * r1 + r2 * r2
    return (r1 if m == 0 else r2) / r2s * (f(z[0] + r1, z[1] + r2) - fz)
# A: H(q) = int d^2r e^{i q r1} K^m(r) f(z+r) - fz * 2 pi i q^m/q^2  (q along k)
rr = np.linspace(0, 9, 9001)[1:]; th = np.linspace(0, 2 * np.pi, 721)[:-1]
R, TH = np.meshgrid(rr, th, indexing='ij')
FV = f(z[0] + R * np.cos(TH), z[1] + R * np.sin(TH))
def H(m, q):
    e = np.cos(TH) if m == 0 else np.sin(TH)
    integrand = np.exp(1j * q * R * np.cos(TH)) * e * FV           # K^m r dr = e^m dr
    val = np.trapezoid(np.trapezoid(integrand, rr, axis=0), dx=th[1] - th[0]) + (integrand[:, 0].sum() * 0)
    val = np.sum(np.trapezoid(integrand, rr, axis=0)) * (th[1] - th[0])
    return val - fz * 2j * np.pi * (1.0 if m == 0 else 0.0) / q
for m in (0, 1):
    ts = np.concatenate([np.linspace(0, 20, 801), np.linspace(20.05, 60, 800)])
    Hs = np.array([H(m, 1 + tt) for tt in ts])
    SA = np.trapezoid(Hs, ts)
    # tail beyond t = 60: fit H ~ a/q^2 + b/q^3
    q = 1 + ts[-200:]; Mfit = np.vstack([1 / q ** 2, 1 / q ** 3]).T
    coef = np.linalg.lstsq(Mfit, Hs[-200:], rcond=None)[0]
    SA += coef[0] / (1 + ts[-1]) + coef[1] / (2 * (1 + ts[-1]) ** 2)
    # B and C
    def Hu(kap, a=1.0, U=2000.0):
        g = lambda u, part: (np.real if part == 0 else np.imag)(h(m, kap, u))
        pts = [-3, -1, 0, 1, 3]
        val = quad(lambda u: h(m, kap, u), -U, a * U, points=[p for p in pts if -U < p < a * U], limit=800, epsabs=1e-12)[0]
        if m == 1:   # analytic tail of -fz * u/(kap^2+u^2) beyond the cut is O(kap^2/U^2): drop
            pass
        return val
    def SB(a=1.0, U=2000.0):
        d = np.pi * Hu(0.0 + 1e-12, a, U) if m == 1 else np.pi * 0.5 * (Hu(1e-9, a, U) + Hu(-1e-9, a, U))
        pv = quad(lambda kap: (np.exp(1j * kap) * Hu(kap, a, U) - np.exp(-1j * kap) * Hu(-kap, a, U)).real / kap, 0, 60, limit=400)[0] \
           + 1j * quad(lambda kap: (np.exp(1j * kap) * Hu(kap, a, U) - np.exp(-1j * kap) * Hu(-kap, a, U)).imag / kap, 0, 60, limit=400)[0]
        return d + 1j * pv
    vB = SB(); vC = SB(a=3.0)
    print('m=%d (%s)  A: %.6f%+.6fi   B (symmetric u): %.6f%+.6fi   C (u in (-U,3U)): %.6f%+.6fi' % (
        m, 'parallel to k' if m == 0 else 'perpendicular', SA.real, SA.imag, vB.real, vB.imag, vC.real, vC.imag), flush=True)
