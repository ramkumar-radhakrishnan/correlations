"""Re-derivation of the non-instantaneous one-gluon -> one-gluon wave functions (Calculation_main.tex, Eqs. onegrho1/onegrho12
and onegrho2/onegrho22) from the matrix elements m2, m3 and hermiticity, compared with the momentum-space forms in the notes.
Colour: explicit SU(3) f^{abc}; the charge rho^c(+-r) is a random 8-vector v.  Momenta: random 2-vectors, plus momenta numbers.
Output: arrays over (a, b, i, j) of the coefficient of |g^b_j(p)> for incoming |g^a_i(k)>, delta functions stripped."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'three_rho'))
import numerics as nm
_, F = nm.structure_constants(3)
PI = np.pi
g = 1.0
I2 = np.eye(2)


def m3(bra1, bra2, ket):
    """<g^c_l(q) g^b_j(p)| H_ggg |g^a_i(k)> without delta: returns tensor T[c,b,a, l,j,i]; args (plus, vec)"""
    (qp, q), (pp, p), (kp, k) = bra1, bra2, ket
    pre = 1j * g / (8 * PI**1.5 * np.sqrt(kp * pp * qp))
    br = (np.einsum('i,jl->lji', p - q + (qp - pp) / kp * k, I2)
          + np.einsum('j,il->lji', k + q - (kp + qp) / pp * p, I2)
          + np.einsum('l,ij->lji', (kp + pp) / qp * q - p - k, I2))
    col = np.einsum('abc->cba', F)          # f^{abc} indexed [c, b, a]
    return pre * np.einsum('cba,lji->cbalji', col, br)


def Eg(pp, p):
    return p @ p / (2 * pp)


def psi1(kp, k, pp, p, v):
    """p+ > k+: valence adds gluon r = p - k (m2), then H_ggg merges (k, r) -> p.  Coefficient C[a,b,i,j]"""
    rp, r = pp - kp, p - k
    # <g^c_l(k) g^d_m(r)| H_g |g^a_i(k)> = delta^{ac} delta_{il} g r^m rho^d(-r) / (4 pi^1.5 r+^1.5)  (times 2 for the two terms, cancels the 1/2)
    em = g * r / (4 * PI**1.5 * rp**1.5)                  # [m], colour rho^d(-r) = v
    # <g^b_j(p)| H_ggg |g^c_l(k) g^d_m(r)> = conj(<g^c_l(k) g^d_m(r)| H_ggg |g^b_j(p)>)
    T = np.conj(m3((kp, k), (rp, r), (pp, p)))            # [c, d, b, l, m, j]
    den = (Eg(kp, k) - Eg(kp, k) - Eg(rp, r)) * (Eg(kp, k) - Eg(pp, p))
    # contract: c = a, l = i, d with v, m with em
    C = np.einsum('adbimj,d,m->abij', T, v, em) / den
    return C


def psi2(kp, k, pp, p, v):
    """p+ < k+: H_ggg splits k -> (p, r = k - p), then the valence absorbs r (m2 conjugated: rho^d(+r))"""
    rp, r = kp - pp, k - p
    T = m3((pp, p), (rp, r), (kp, k))                      # [c, d, a, l, m, i]: bra1 (c,l,p), bra2 (d,m,r)
    ab = g * r / (4 * PI**1.5 * rp**1.5)                   # conj of m2 entry; rho^d(+r) = v
    den = (Eg(kp, k) - Eg(pp, p) - Eg(rp, r)) * (Eg(kp, k) - Eg(pp, p))
    C = np.einsum('bdajmi,d,m->abij', T, v, ab) / den
    return C


def user_onegrho12(kp, k, pp, p, v):
    """Eq. onegrho12 as written (coefficient of |g^b_j(p)>, d^3p stripped)"""
    col = np.einsum('abc,c->ab', F, v)
    r = p - k
    den = (k @ k / kp - p @ p / pp) * ((k - p) @ (k - p)) / (kp - pp)
    br = (np.einsum('j,il->ijl', 2 * k - p + (pp - 2 * kp) / pp * p, I2) + np.einsum('i,jl->ijl', 2 * p - k - (2 * pp - kp) / kp * k, I2)
          + np.einsum('l,ij->ijl', (pp + kp) / (pp - kp) * (p - k) - k - p, I2))
    pre = -1j * g**2 / (8 * PI**3) / (np.sqrt(pp * (pp - kp) * kp) * (pp - kp)**1.5) / den
    return pre * np.einsum('ab,l,ijl->abij', col, r, br)


def user_onegrho1_xi(kp, k, pp, p, v):
    """the xi form after onegrho12 (coefficient per d xi d^2p; convert with dp+ = k+ dxi)"""
    xi = pp / kp
    col = np.einsum('abc,c->ab', F, v)
    vec = np.einsum('m,l->ml', p - k, p - xi * k) / (((p - k) @ (p - k)) * (p @ p - xi * k @ k))
    br = np.einsum('il,jm->ijlm', I2, I2) - 1 / xi * np.einsum('jl,im->ijlm', I2, I2) - 1 / (1 - xi) * np.einsum('ij,lm->ijlm', I2, I2)
    return -1j * g**2 / (4 * PI**3) * np.sqrt(xi) / (1 - xi) * np.einsum('ab,ml,ijlm->abij', col, vec, br) / kp


def user_onegrho22(kp, k, pp, p, v):
    """Eq. onegrho22 as written, with the obvious typo ((k++p+)/(k+-p+)) + (k-p)^l read as a product"""
    col = np.einsum('abc,c->ab', F, v)
    r = k - p
    d1 = k @ k / kp - r @ r / (kp - pp) - p @ p / pp
    d2 = k @ k / kp - p @ p / pp
    br = (np.einsum('i,jl->ijl', 2 * p - k + (kp - 2 * pp) / kp * k, I2) + np.einsum('j,il->ijl', 2 * k - p - (2 * kp - pp) / pp * p, I2)
          + np.einsum('l,ij->ijl', (kp + pp) / (kp - pp) * r - p - k, I2))
    pre = 1j * g**2 / (8 * PI**3) / (np.sqrt(kp * (kp - pp) * pp) * (kp - pp)**1.5) / (d1 * d2)
    return pre * np.einsum('ab,l,ijl->abij', col, r, br)


def user_onegrho2_xi(kp, k, pp, p, v):
    xi = pp / kp
    col = np.einsum('abc,c->ab', F, v)
    vec = np.einsum('m,l->ml', k - p, p - xi * k) / ((p @ p - xi * k @ k) * ((p - xi * k) @ (p - xi * k)))
    br = np.einsum('il,jm->ijlm', I2, I2) - 1 / xi * np.einsum('jl,im->ijlm', I2, I2) - 1 / (1 - xi) * np.einsum('lm,ij->ijlm', I2, I2)
    return 1j * g**2 / (4 * PI**3) * xi**1.5 / (1 - xi) * np.einsum('ab,ml,ijlm->abij', col, vec, br) / kp


if __name__ == '__main__':
    rng = np.random.default_rng(3)
    for t in range(3):
        k, p, v = rng.normal(size=2), rng.normal(size=2), rng.normal(size=8)
        kp = 1.3
        pp = kp * rng.uniform(1.2, 3.0)
        a = psi1(kp, k, pp, p, v)
        print('p+>k+: onegrho12 / derived = %s ;  xi-form / derived = %s' % (
            np.round((user_onegrho12(kp, k, pp, p, v) / a)[np.abs(a) > 1e-8 * np.abs(a).max()].mean(), 8),
            np.round((user_onegrho1_xi(kp, k, pp, p, v) / a)[np.abs(a) > 1e-8 * np.abs(a).max()].mean(), 8)))
        pp = kp * rng.uniform(0.1, 0.9)
        a = psi2(kp, k, pp, p, v)
        print('p+<k+: onegrho22 / derived = %s ;  xi-form / derived = %s' % (
            np.round((user_onegrho22(kp, k, pp, p, v) / a)[np.abs(a) > 1e-8 * np.abs(a).max()].mean(), 8),
            np.round((user_onegrho2_xi(kp, k, pp, p, v) / a)[np.abs(a) > 1e-8 * np.abs(a).max()].mean(), 8)))
