"""C rows of Group III after the w' integration (lines 72, 85, 105, 124) and with the plus prescription (76, 90, 109, 128).
Kernels only (the colour factors are those of the pre-integration lines, verified separately).
The common phase e^{-ik(y'-w)} is dropped; the p-dependent phase e^{-ik.(y'-z)p/k} is kept."""
import numpy as np
from fp import F, G, PI, I2, K, sq
import user23 as us

pi3 = us.pi3


def phase(c, k, p):
    return np.exp(-1j * (c.kv @ (c.pt['yp'] - c.pt['z'])) * p / k)


def tens3(c):
    """the three tensor structures K^m(y'-z) K^{k'}(x'-y') K^i(y-w) K^j(x-z) contracted with
    d_{k'm}d_{ij}, d_{jm}d_{ik'}, d_{im}d_{jk'}"""
    a, b, ci, dj = K(c.pt['yp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']), K(c.pt['y'] - c.pt['w']), K(c.pt['x'] - c.pt['z'])
    return (b @ a) * (ci @ dj), (dj @ a) * (ci @ b), (ci @ a) * (dj @ b)


def ker72(c, k, p, pref):
    """pref * (1/k) [t1/(p+k) - t2/p - t3/k]   (rows III(1): pref = -i g^4/(8 pi^4); V, VI: -i g^4/(4 pi^4))"""
    t1, t2, t3 = tens3(c)
    return pref / pi3 / k * (t1 / (p + k) - t2 / p - t3 / k)


def ker72_pieces(c, k, p, pref):
    """split: (regular part, coefficient g1 of 1/p) -- incl. the p-dependent phase"""
    t1, t2, t3 = tens3(c)
    ph = phase(c, k, p)
    reg = pref / pi3 / k * (t1 / (p + k) - t3 / k) * ph
    g1 = pref / pi3 / k * (-t2) * ph
    return reg, g1


def logterm72(c, k, Lam, vee, preflog):
    """the log term of lines 76 / 109 / 128: preflog/(2pi)^3 (1/k) log(...) (y'-z).(x-z) (x'-y').(y-w)/(...)"""
    t1, t2, t3 = tens3(c)
    return preflog / pi3 / k * (np.log(vee / Lam) + np.log(1 - k / vee)) * t2


def bracket85(c, k, p, d=2):
    x, z, w = c.pt['x'], c.pt['z'], c.pt['w']
    xz, xw, zw = x - z, x - w, z - w
    P = p + k
    S = sq(xz) - sq(xw)
    out = I2 * d * p * S / (2 * sq(zw) * P**2)
    out = out + p / (k * P) * I2 * (xz @ zw / sq(zw) - xz @ xw / (2 * sq(xw)) - S / (2 * sq(zw)))
    out = out + 1 / P * I2 * (xw @ zw / sq(zw) + xz @ xw / (2 * sq(xz)) - S / (2 * sq(zw)))
    sing = -(np.outer(xw, zw) / sq(zw) + np.outer(xw, xz) / (2 * sq(xz)))        # times 1/p   [k', m]
    out = out - p / k**2 * (np.outer(xz, zw) / sq(zw) - np.outer(xz, xw) / (2 * sq(xw)))
    out = out - 1 / k * (np.outer(zw, xz) / sq(zw) - np.outer(xw, xz) / (2 * sq(xw))
                         + np.outer(zw, xw) / sq(zw) + np.outer(xz, xw) / (2 * sq(xz)))
    return out, sing


def ker85_pieces(c, k, p):
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['w'])
    a, b = K(c.pt['yp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp'])
    reg, sing = bracket85(c, k, p)
    ph = phase(c, k, p)
    pre = -G**4 / (4 * PI**4) / pi3 / D
    return pre * np.einsum('m,K,Km->', a, b, reg) * ph, pre * np.einsum('m,K,Km->', a, b, sing) * ph


def logterm90(c, k, Lam, vee):
    x, z, w, yp, xp = (c.pt[n] for n in ('x', 'z', 'w', 'yp', 'xp'))
    xz, xw, zw = x - z, x - w, z - w
    val = ((xw @ (xp - yp)) * (zw @ (yp - z)) / (sq(zw) * sq(xw) * sq(yp - z) * sq(xp - yp))
           + (xz @ (yp - z)) * (xw @ (xp - yp)) / (2 * sq(xz) * sq(xw) * sq(yp - z) * sq(xp - yp)))
    return G**4 / (8 * PI**5) / pi3 / k * (np.log(vee / Lam) + np.log(1 - k / vee)) * val
