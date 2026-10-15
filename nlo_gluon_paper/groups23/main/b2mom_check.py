"""B_2 in momentum space (Calculation_main.tex, Omega section) vs 8 pi^3 x the vacuum two-gluon wave function
(part 1 + part 2 + contgg + contgg3, all verified from the matrix elements).  Gluon 1 = (a, i, k+, k), gluon 2 = (b, j, p+, p).
The wave functions are written in (P+, xi, P, ptilde); the measure dP+ dxi d2P d2pt = dk+ dp+ d2k d2p / P+."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'three_rho'))
import numerics as nm
_, F = nm.structure_constants(3)
PI = np.pi
g = 1.3
I2 = np.eye(2)
rng = np.random.default_rng(2)


def psi(kp, k, pp, p):
    """[c(rho), a, b, i, j] one-rho part, and [a, b, i, j] coefficient of {rho^a(-k), rho^b(-p)} (two-rho)"""
    Pp, P = kp + pp, k + p
    xi = kp / Pp                       # gluon 1 carries xi
    pt = k - xi * P
    den = P @ P * xi * (1 - xi) + pt @ pt
    s = np.sqrt(xi * (1 - xi))
    # contgg: gluon (c', l', xi) = gluon 1 (a, i); (b', j', 1-xi) = gluon 2 (b, j); rho^{a'} -> rho^c
    br = (np.einsum('m,ji->mij', pt, I2) - np.einsum('j,mi->mij', pt, I2) / (1 - xi) - np.einsum('i,mj->mij', pt, I2) / xi)
    t2 = -1j * g**2 * s / (8 * PI**3 * den * (P @ P)) * np.einsum('cba,m,mij->cabij', F, P, br)
    # contgg3: (c', j', xi) = gluon 1 (a, i); (b', i', 1-xi) = gluon 2 (b, j)
    t3 = 1j * g**2 * (1 - 2 * xi) * s / (16 * PI**3 * den) * np.einsum('cba,ij->cabij', F, I2)
    # part 1: gluon (a, i, xi) = gluon 1, (b, j, 1-xi, k-p) = gluon 2; notes' p -> k (gluon 1), notes' k-p -> p
    D = xi * (p @ p) + (1 - xi) * (k @ k)
    t1 = -1j * g**2 * s / (16 * PI**3 * D) * (1 / ((1 - xi) * (k @ k)) - 1 / (xi * (p @ p))) * np.einsum('abc,i,j->cabij', F, k, p)
    # part 2 (two-rho): {rho^a(-k), rho^b(-p)}
    two = g**2 / (16 * PI**3 * (k @ k) * (p @ p) * s) * np.einsum('i,j,ab->abij', k, p, np.ones((8, 8)))
    return (t1 + t2 + t3) / Pp, two / Pp, (t1 / Pp, (t2 + t3) / Pp)


def B2(kp, k, pp, p):
    """the notes' B_2^{ba}_{ji}(k+, k; p+, p): one-rho [c, a, b, i, j], two-rho [a, b, i, j]"""
    P = k + p
    s = np.sqrt(pp * kp)
    t1 = 1j * g**2 * s / (2 * (pp * (k @ k) + kp * (p @ p))) * (1 / (kp * (p @ p)) - 1 / ((k @ k) * pp)) * np.einsum('abc,j,i->cabij', F, p, k)
    two = g**2 / (2 * (p @ p) * (k @ k) * s) * np.einsum('j,i,ab->abij', p, k, np.ones((8, 8)))
    br = (np.einsum('lm,ji->lmij', I2, I2) / (pp + kp) - np.einsum('im,jl->lmij', I2, I2) / kp - np.einsum('jm,il->lmij', I2, I2) / pp)
    v = p * kp - pp * k
    t3 = -1j * g**2 / ((P @ P) * s * (pp + kp) * ((p @ p) / pp + (k @ k) / kp)) * np.einsum('abc,l,m,lmij->cabij', F, P, v, br)
    t4 = -1j * g**2 * (pp - kp) / (2 * s * (pp + kp)**2 * ((k @ k) / kp + (p @ p) / pp)) * np.einsum('abc,ij->cabij', F, I2)
    return t1 + t3 + t4, two, (t1, t3 + t4)


def ratio(a, b):
    m = np.abs(b) > 1e-9 * np.abs(b).max()
    r = a[m] / b[m]
    return np.round(r.mean(), 10), float(np.abs(r - r.mean()).max())


if __name__ == '__main__':
    for _ in range(3):
        kp, pp = rng.uniform(0.3, 2, size=2)
        k, p = rng.normal(size=2), rng.normal(size=2)
        one_w, two_w, (w1, w234) = psi(kp, k, pp, p)
        one_b, two_b, (b1, b34) = B2(kp, k, pp, p)
        print('B2/psi: total one-rho %s, two-rho %s, first term %s, terms 3+4 %s' % (
            ratio(one_b, one_w), ratio(two_b, two_w), ratio(b1, w1), ratio(b34, w234)))
