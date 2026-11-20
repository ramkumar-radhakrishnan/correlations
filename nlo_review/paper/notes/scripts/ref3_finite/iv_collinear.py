# Row I.V, Eq. (F.5)/(sIV2): local structure of the kernel K^V_2b = P^i K^m(z-x) F^{im} + (11/6) e^{-is} P.K(y-x) K(y-z).(y-x)/(x-z)^2
# at the collinear point z -> x (the colour factor of K^V_2b does not vanish there; that of K^V_2a does).
# Prints r^2 K(phi) on a circle |z-x| = r: its angular mean (-> 0 by the 11/6 term) and its anisotropic part.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, lib
from scipy.integrate import quad
def Nim(xi, x, y, z):
    yx, yz, xz = y - x, y - z, x - z
    YX, YZ, XZ = yx @ yx, yz @ yz, xz @ xz
    B1 = np.outer(yx, xz) / XZ - np.outer(yx, yz) / (2 * YZ)
    B2 = np.outer(yz, xz) / XZ + np.outer(yz, yx).T * 0 + np.outer(yz, yx) * 0 + np.outer(yz, yx) * 0
    B2 = np.outer(yz, xz) / XZ + np.outer(yz, yx) / (2 * YX)           # B2^{im} = (y-z)^i (x-z)^m/XZ + (y-x)^m (y-z)^i/(2YX)
    T1 = (yz @ xz) / XZ + (yx @ yz) / (2 * YX) - (YX - YZ) / (2 * XZ)
    T2 = (yx @ xz) / XZ - (yx @ yz) / (2 * YZ) - (YX - YZ) / (2 * XZ)
    Q = (np.outer(xz, yx) / XZ - np.outer(yz, yx) / (2 * YZ)        # Q^{im}: (y-x)^m (x-z)^i/XZ - (y-x)^m (y-z)^i/(2YZ)
         + np.outer(xz, yz) / XZ + np.outer(yx, yz) / (2 * YX))     #        + (y-z)^m (x-z)^i/XZ + (y-x)^i (y-z)^m/(2YX)
    d = np.eye(2)
    return (d * xi * (1 - xi) * (YX - YZ) / XZ + d * ((1 - xi) * T1 + xi * T2)
            + YX / YZ * B2 + YZ / YX * B1 - Q)
def Freg(x, y, z):
    s = lib.KVEC @ (x - z); YX = (y - x) @ (y - x); YZ = (y - z) @ (y - z)
    out = np.zeros((2, 2), complex)
    for i in range(2):
        for m in range(2):
            f = lambda xi, part: (np.cos if part == 0 else np.sin)(xi * s) * Nim(xi, x, y, z)[i, m] / (xi * YX + (1 - xi) * YZ)
            out[i, m] = quad(f, 0, 1, args=(0,), epsabs=1e-13, epsrel=1e-12)[0] - 1j * quad(f, 0, 1, args=(1,), epsabs=1e-13, epsrel=1e-12)[0]
    return out
def K2b(x, y, z, P):
    s = lib.KVEC @ (x - z)
    return P @ Freg(x, y, z) @ lib.Kv(z - x) + 11 / 6 * np.exp(-1j * s) * (P @ lib.Kv(y - x)) * (lib.Kv(y - z) @ (y - x)) / ((x - z) @ (x - z))
def K2a(x, y, z, P): return P @ Freg(x, y, z) @ lib.Kv(z - x)
x = np.array([0.2, -0.1]); y = np.array([-0.6, 0.7]); P = lib.Kv(np.array([0.5, 0.9]) - np.array([-0.8, -0.3]))
nang = 72; phis = np.linspace(0, 2 * np.pi, nang, endpoint=False)
for name, K in (('K2b (with 11/6)', K2b), ('K2a = P.K.F only', K2a)):
    print(name)
    for r in (1e-2, 1e-3, 1e-4):
        vals = np.array([r * r * K(x, y, x + r * np.array([np.cos(p), np.sin(p)]), P) for p in phis])
        c2 = np.mean(vals * np.exp(-2j * phis)); c1 = np.mean(vals * np.exp(-1j * phis))
        print('   r=%.0e  mean r^2K = %+.3e%+.3ei   |r^2K| max = %.4f   |2nd harmonic| = %.4f  |1st harmonic| = %.2e' % (
            r, vals.mean().real, vals.mean().imag, np.abs(vals).max(), abs(c2), abs(c1)))
