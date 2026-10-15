"""First-principles check of the vacuum one-gluon terms at order g^3 (Calculation_main.tex, Sec. 2.2):
Psi^{2a}, Psi^{2b} (instantaneous), Psi^3 (two gluons merge), Psi^4 (instantaneous pair, one absorbed),
Psi^5 (splitting, one absorbed), and the {rho,rho} / [rho,rho] split of Psi^4 + Psi^5.

Everything is pointwise in the momenta (before any p+ or transverse integration).
Matrix elements m1-m5 as in the appendix; Rayleigh-Schroedinger signs: order n carries (-1)^n / (E...E);
a two-gluon intermediate state carries 1/2.
A product of two charges rho^{c1}(u) rho^{c2}(v) is stored as a tensor T[c1, c2, (final colour), (final pol)]
together with its momentum arguments (u, v); colours are SU(3), f^{abc} explicit."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'three_rho'))
import numerics as nm
_, F = nm.structure_constants(3)
PI = np.pi
g = 1.3
I2 = np.eye(2)
NC = 3
D8 = np.eye(8)


# ------------------------------------------------------------------ matrix elements (appendix)
def em(P):
    """<g^a_i(P)|H_g|0> = g P^i rho^a(-P) / (4 pi^1.5 P+^1.5): vector [i] (colour = colour of rho)"""
    pp, p = P
    return g * p / (4 * PI**1.5 * pp**1.5)


def m3(bra1, bra2, ket):
    """<g^c_l(q) g^b_j(p)| H_ggg |g^a_i(k)>, delta stripped: T[c,b,a, l,j,i]"""
    (qp, q), (pp, p), (kp, k) = bra1, bra2, ket
    pre = 1j * g / (8 * PI**1.5 * np.sqrt(kp * pp * qp))
    br = (np.einsum('i,jl->lji', p - q + (qp - pp) / kp * k, I2)
          + np.einsum('j,il->lji', k + q - (kp + qp) / pp * p, I2)
          + np.einsum('l,ij->lji', (kp + pp) / qp * q - p - k, I2))
    return pre * np.einsum('abc,lji->cbalji', F, br)


def m4(Q, P):
    """<g^c_j(Q) g^b_i(P)|H_inst|0>: T[a(rho), c, b, j, i]; rho^a(-P-Q)"""
    (qp, q), (pp, p) = Q, P
    pre = 1j * g**2 * (qp - pp) / (2 * (2 * PI)**3 * np.sqrt(pp * qp) * (pp + qp)**2)
    return pre * np.einsum('abc,ji->acbji', F, I2)


def m5(Q, P):
    """<g^c_j(Q)|H_inst|g^b_i(P)>: T[a(rho), c, b, j, i]; rho^a(P-Q)"""
    (qp, q), (pp, p) = Q, P
    pre = 1j * g**2 * (pp + qp) / (2 * (2 * PI)**3 * np.sqrt(pp * qp) * (pp - qp)**2)
    return pre * np.einsum('abc,ji->acbji', F, I2)


def E(P):
    pp, p = P
    return p @ p / (2 * pp)


def sub(K, P):
    return (K[0] - P[0], K[1] - P[1])


# ------------------------------------------------------------------ first principles
def fp_psi2a(K, P):
    """|g^b_j(K)> <K|H_inst|P><P|H_g|0>/(E_K E_P);  rho^e(P-K) rho^a(-P): T[e, a, b, j]"""
    M = m5(K, P)                    # [e, b, a, j, i]
    return np.einsum('ebaji,i->eabj', M, em(P)) / (E(K) * E(P))


def fp_psi2b(K, P):
    """final gluon P (softer): |g^b_j(P)> <P|H_inst|K><K|H_g|0>/(E_P E_K);  rho^e(K-P) rho^a(-K): T[e,a,b,j]"""
    M = m5(P, K)
    return np.einsum('ebaji,i->eabj', M, em(K)) / (E(P) * E(K))


def fp_psi3(K, P):
    """-1/2 int dg1 dg2 over the pair state (g1, g2) = (P, K-P) (g1 is the integration variable, g2 fixed by the
    delta function), both m2 emission orders.
    Returns dict {('p','q'): T, ('q','p'): T}: rho(-P)rho(-(K-P)) and rho(-(K-P))rho(-P), T[c1, c2, d, l]."""
    Q = sub(K, P)
    out = {('p', 'q'): np.zeros((8, 8, 8, 2), complex), ('q', 'p'): np.zeros((8, 8, 8, 2), complex)}
    for g1, g2, n1, n2 in ((P, Q, 'p', 'q'),):
        # <K|H_ggg|g1 g2> = conj(<g1 g2|H_ggg|K>), T[c(g1), b(g2), d(K), l1, l2, l]
        V = np.conj(m3(g1, g2, K))
        E12 = E(g1) + E(g2)
        # m2: g1 emitted first (intermediate single gluon g1), then g2 emitted: rho(-g2) rho(-g1)
        out[(n2, n1)] += -0.5 * np.einsum('cbdxyl,y,x->bcdl', V, em(g2), em(g1)) / (E(K) * E12 * E(g1))
        # g2 emitted first, then g1: rho(-g1) rho(-g2)
        out[(n1, n2)] += -0.5 * np.einsum('cbdxyl,x,y->cbdl', V, em(g1), em(g2)) / (E(K) * E12 * E(g2))
    return out


def ab(Qg):
    """<0|H_g|g^b_j(Q)> = g Q^j rho^b(+Q)/(4 pi^1.5 Q+^1.5)"""
    return em(Qg)


def fp_psi4(K, P):
    """final gluon P; pair (P, Q=K-P) created instantaneously, Q absorbed: rho^x(Q) rho^a(-K): T[x, a, d, l]"""
    Q = sub(K, P)
    out = np.zeros((8, 8, 8, 2), complex)
    # (1/2) [ <P|H_g|g1 g2> ] with g1 = P, g2 = Q and g1 = Q, g2 = P; m4(Q', P') with bra <Q' P'|
    for g1, g2 in ((P, Q), (Q, P)):
        M = m4(g1, g2)                         # [a, c(g1), b(g2), j1, i2]
        if g1 is P:   # absorb g2 = Q: colour b, pol i2;  final = g1: colour c, pol j1
            out += 0.5 * np.einsum('acbji,i->bacj', M, ab(Q)) / (E(P) * (E(P) + E(Q)))
        else:         # absorb g1 = Q: colour c, pol j1; final = g2: colour b, pol i2
            out += 0.5 * np.einsum('acbji,j->cabi', M, ab(Q)) / (E(P) * (E(P) + E(Q)))
    return out


def fp_psi5(K, P):
    """final gluon P; K emitted, K -> (P, Q), Q absorbed: -1/2 * 2 terms; rho^x(Q) rho^a(-K): T[x, a, d, l]"""
    Q = sub(K, P)
    out = np.zeros((8, 8, 8, 2), complex)
    E12 = E(P) + E(Q)
    for g1, g2 in ((P, Q), (Q, P)):
        M = m3(g1, g2, K)                      # [c(g1), b(g2), a(K), l1, l2, i]
        if g1 is P:   # absorb g2
            out += -0.5 * np.einsum('cbaxyi,y,i->bacx', M, ab(Q), em(K)) / (E(P) * E12 * E(K))
        else:         # absorb g1 = Q, final g2 = P
            out += -0.5 * np.einsum('cbaxyi,x,i->caby', M, ab(Q), em(K)) / (E(P) * E12 * E(K))
    return out


# ------------------------------------------------------------------ the notes, as written
def w_instgg(K, P):
    """Eq. instgg: ig^3 f^{abc} rho^c(p-k) rho^a(-p) p^j (p+ + k+) sqrt(k+)/(16 pi^4.5 p+ (p+-k+)^2 k^2 p^2) |g^b_j(k)>"""
    (kp, k), (pp, p) = K, P
    s = 1j * g**3 * (pp + kp) * np.sqrt(kp) / (16 * PI**4.5 * pp * (pp - kp)**2 * (k @ k) * (p @ p))
    return s * np.einsum('abc,j->cabj', F, p)


def w_gguu(K, P):
    """Eq. eq:gg-uu: ig^3 f^{abc} rho^c(k-p) rho^a(-k) k^j (p+ + k+) sqrt(p+)/(16 pi^4.5 k+ (p+-k+)^2 k^2 p^2) |g^b_j(p)>"""
    (kp, k), (pp, p) = K, P
    s = 1j * g**3 * (pp + kp) * np.sqrt(pp) / (16 * PI**4.5 * kp * (pp - kp)**2 * (k @ k) * (p @ p))
    return s * np.einsum('abc,j->cabj', F, k)


def w_3grhorho(K, P):
    """Eq. 3grhorho (per dk+ dp+, i.e. divided by k+ for d xi): ig^3 f^{bad} rho^b(-p) rho^a(-k+p) [...]^l /(16 pi^4.5 xi (k-p)^2 D k^2 sqrt(k+))"""
    (kp, k), (pp, p) = K, P
    xi = pp / kp
    q = k - p
    D = (1 - xi) * (p @ p) + xi * (q @ q)
    B = ((2 * p - 2 * xi * k) * (p @ q) + (2 * (k @ p) - 2 / xi * (p @ p)) * q
         + ((1 + xi) / (1 - xi) * (q @ q) - k @ k + p @ p) * p)
    s = 1j * g**3 / (16 * PI**4.5 * xi * (q @ q) * D * (k @ k) * np.sqrt(kp)) / kp
    return np.einsum('bad,l->badl', F, B) * s           # T[b, a, d, l]: rho^b(-p) rho^a(-(k-p))


def w_psif(K, P):
    """Eq. psif (per dk+ dp+): -ig^3 f^{abc} rho^c(k-p) rho^a(-k) (k-p)^i (2xi-1) xi^1.5/(16 pi^4.5 D p^2 (1-xi)) /sqrt(k+) |g^b_i(p)>"""
    (kp, k), (pp, p) = K, P
    xi = pp / kp
    q = k - p
    D = xi * (q @ q) + (1 - xi) * (p @ p)
    s = -1j * g**3 * (2 * xi - 1) * xi**1.5 / (16 * PI**4.5 * D * (p @ p) * (1 - xi)) / np.sqrt(kp) / kp
    return s * np.einsum('abc,i->cabi', F, q)


def w_psiff(K, P):
    """Eq. psiff (per dk+ dp+)"""
    (kp, k), (pp, p) = K, P
    xi = pp / kp
    q = k - p
    D = xi * (q @ q) + (1 - xi) * (p @ p)
    vec = (((2 * xi - 1) * (k @ k) - p @ p - (1 + xi) / (1 - xi) * (q @ q)) * k
           + ((xi - 1) / xi * (k @ p) + (1 - xi**2) / xi * (k @ k)) * 2 * p)
    s = 1j * g**3 * xi**1.5 / (16 * PI**4.5 * (k @ k) * (p @ p) * D * (1 - xi)) / np.sqrt(kp) / kp
    return s * np.einsum('abc,i->cabi', F, vec)


def w_45tworho(K, P):
    """Eq. 4+5tworho, mapped back (p <-> k): coefficient of {rho^c(k-p), rho^a(-k)} |g^b_j(xi k+, p)>, per dk+ dp+"""
    (kp, k), (pp, p) = K, P
    xi = pp / kp
    # notes' variables: final gluon transverse momentum k_n = p, other p_n = k
    kn, pn = p, k
    qn = kn - pn
    D = xi * (qn @ qn) + (1 - xi) * (kn @ kn)
    vec = ((2 / xi * (pn @ pn - kn @ pn) - pn @ pn + 2 * (kn @ pn)) * kn
           - (kn @ kn + (1 + xi) / (1 - xi) * (qn @ qn)) * pn)
    s = 1j * g**3 * xi**1.5 / (32 * PI**4.5 * (kn @ kn) * (pn @ pn) * D * (1 - xi)) / np.sqrt(kp) / kp
    return s * np.einsum('abc,j->cabj', F, vec)        # S[c, a, b, j] for {rho^c(k-p), rho^a(-k)}


def w_45onerho(K, P):
    """Eq. 4+5onerho integrand before the transverse integral is not written; we compare the one-rho part of
    psif+psiff with -(N_c/2) ... directly below."""
    raise NotImplementedError


# ------------------------------------------------------------------ comparisons
def ratio(a, b):
    m = np.abs(b) > 1e-9 * np.abs(b).max()
    r = a[m] / b[m]
    return r.mean(), np.abs(r - r.mean()).max()


def split(T):
    """T[c1,c2,...] for rho^{c1}(u) rho^{c2}(v) -> (S, C): S[c1,c2,...] coefficient of {rho^{c1}(u), rho^{c2}(v)}
    and C[e,...] coefficient of rho^e(u+v)"""
    S = 0.5 * T
    C = 0.5j * np.einsum('xye,xy...->e...', F, T)
    return S, C


if __name__ == '__main__':
    rng = np.random.default_rng(11)
    res = {}
    for t in range(3):
        kp = 1.7
        k, p = rng.normal(size=2), rng.normal(size=2)
        # Psi^{2a}: final K, intermediate P with p+ < k+
        P = (kp * rng.uniform(0.1, 0.9), p)
        K = (kp, k)
        res.setdefault('instgg', []).append(ratio(w_instgg(K, P), fp_psi2a(K, P)))
        res.setdefault('gg-uu', []).append(ratio(w_gguu(K, P), fp_psi2b(K, P)))
        # Psi^3: symmetrised comparison
        Q = sub(K, P)
        # int dp FP(p) = int dp [A(p) + B(k-p)] rho(-p) rho(-(k-p)), A = fp[(p,q)], B = fp[(q,p)]
        tot = fp_psi3(K, P)[('p', 'q')] + fp_psi3(K, Q)[('q', 'p')]
        res.setdefault('3grhorho', []).append(ratio(w_3grhorho(K, P), tot))
        # Psi^4, Psi^5
        f4, f5 = fp_psi4(K, P), fp_psi5(K, P)
        res.setdefault('psif', []).append(ratio(w_psif(K, P), f4))
        res.setdefault('psiff', []).append(ratio(w_psiff(K, P), f5))
        S, C = split(f4 + f5)
        res.setdefault('4+5tworho', []).append(ratio(w_45tworho(K, P), S))
    for k_, v in res.items():
        print('%-28s' % k_, ['%.10f%+.10fi (spread %.1e)' % (r.real, r.imag, s) for r, s in v])
