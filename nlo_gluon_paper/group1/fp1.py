"""First-principles evaluation of the Group I integrands (one gluon in the final state).

Reference: Eq. (CS1) of Calculation_main.tex with the coefficients of the same file (as corrected in
../groups23/main): A^(1), A^(3) (coordinate form), B_3 (coordinate form, sign E3), B_2, C_1.
Everything is evaluated pointwise: random transverse points, random classical charges, random SU(3) adjoint
Wilson lines (general, or symmetric U = U^T when sym=True), fixed p+.

Cross section (CS1):  1/(2pi)^3 int_{w,w'} e^{-ik(w'-w)} sum_{a,i} Conj^a_i(w') V^a_i(w)  + c.c.
  Conj^a_i = -Abar^{dag a}_i(k,-w') + U^{ab'}(w') A^{dag b'}_i(k,-w')    (conjugate LO amplitude, charge at x')
           = i alpha_k [U(x') - U(w')]^{ab'} rho^{b'}(x') K^i(x'-w')
  Row II : V = U^{ab}(w) A^{(3)b}_i(k,-w) - Abar^{(3)a}_i(k,-w)
  Row III: V = int dp/2pi int_z Bbar3^{dag ba}_{ji}(k,w;p,-z) [-Abar^b_j(p,-z) + U^{bc}(z) A^c_j(p,-z)]
  Row IV : V = 2 int dp/2pi int_z [U^{ac}(w) U^{bd}(z) Abar^{dag b}_j(p,-z) B2^{dc}_{ji}(k,-w;p,-z)
                                    - Abar^{dag b}_j(p,-z) Bbar2^{ba}_{ji}(k,-w;p,-z)]
  Row V  : V = 2 int dp/2pi int_{x,z} [-Cbar1^{dag bca}_{jki}(k,w-z;p,z-x) Bbar2^{bc}_{jk}(k-p,-z;p,-x)
                                        + Cbar1^{dag dca}_{jki} U^{db}(x) U^{ce}(z) B2^{be}_{jk}(k-p,-z;p,-x)]
Positions (as in Group_I.tex): w measured gluon (amplitude), w' (conjugate), x' charge of the conjugate amplitude,
z the unobserved gluon p+, x and y the charges of the amplitude (x: charge of B_3 / B_2, y: charge of A(p))."""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'groups23'))
from fp import Cfg, F, NC, G, PI, I2, K, sq, alpha, A, Abar, B2one, B2two  # noqa: E402,F401

PI3 = (2 * PI)**3


def T4(kp, pp):
    """T_{ijkm} = d_ik d_jm - (k+/p+) d_jk d_im - (k+/(k+-p+)) d_ij d_km  (vertex bracket of B_3 and C_1)"""
    return (np.einsum('ik,jm->ijkm', I2, I2) - kp / pp * np.einsum('jk,im->ijkm', I2, I2)
            - kp / (kp - pp) * np.einsum('ij,km->ijkm', I2, I2))


def conj(c, kp):
    """Conj^a_i: conjugate of the LO amplitude [a, i]; general U placement U^{ab'}(w')"""
    r = (c.U['xp'] - c.U['wp']) @ c.rho['xp']
    return 1j * alpha(kp) * np.outer(r, K(c.pt['xp'] - c.pt['wp']))


# ------------------------------------------------------------------ B_3 in coordinate space
def B3one(c, kp, Z, pp, Uu, X, Lam, bar=False, sgt=-1, parts=('inst', 'ni', 't4')):
    """B^{ba}_{3ji}(k+, Z; p+, -Uu), one-rho part with charge at X: array [b, a, j, i].
    Z: position of the gluon k+ (annihilated by B_3), Uu: position of the gluon p+ (created).
    sgt: sign of the non-instantaneous Theta(p+ - k+ - Lambda) entry; -1 = corrected (E3), +1 = as written."""
    x, z, u = c.pt[X], c.pt[Z], c.pt[Uu]
    r = c.U[X] @ c.rho[X] if bar else c.rho[X]
    col = np.einsum('abc,c->ab', F, r)                      # f^{abc} rho^c   [a, b]
    Aq, Bq = sq(x - u), sq(x - z)                           # (x-u)^2, (x-z)^2
    gt, lt = pp > kp + Lam, pp < kp - Lam
    S = np.zeros((2, 2), complex)                           # [i, j]
    if 'inst' in parts and (gt or lt):
        S += 1j * G**2 / (2 * PI)**2 * np.sqrt(pp * kp) * (pp + kp) / (kp - pp)**2 / (pp * Aq - kp * Bq) * I2
    if 'ni' in parts and (gt or lt):
        V = kp * (x - z) - pp * (x - u)
        W = (x - z) - (x - u)
        pre = 1j * G**2 / (2 * PI**2 * kp) * np.sqrt(kp * pp) / (kp - pp) / (kp * Bq - pp * Aq)
        fac = sgt / sq(W) if gt else -pp * kp / sq(kp * (x - z) - pp * (x - u))
        S += pre * fac * np.einsum('ijkm,k,m->ij', T4(kp, pp), W, V)
    if 't4' in parts:
        S += (-1j * G**2 / (4 * PI**2 * np.sqrt(pp * kp)) * np.outer(K(x - z), K(x - u))
              * (kp * Bq + pp * Aq) / (pp * Aq - kp * Bq))
    return np.einsum('ab,ij->baji', col, S)


def B3two(c, kp, Z, pp, Uu, X, Y, bar=False):
    """two-rho part -g^2/(4 pi^2 sqrt(pk)) K^i(X-Z) K^j(Y-Uu) {rho^a(X), rho^b(Y)}: [b, a, j, i]"""
    rX = c.U[X] @ c.rho[X] if bar else c.rho[X]
    rY = c.U[Y] @ c.rho[Y] if bar else c.rho[Y]
    pre = -G**2 / (4 * PI**2 * np.sqrt(pp * kp)) * 2
    return pre * np.einsum('a,b,i,j->baji', rX, rY, K(c.pt[X] - c.pt[Z]), K(c.pt[Y] - c.pt[Uu]))


# ------------------------------------------------------------------ the rows (integrand at fixed p+)
def phase(c, w='w'):
    return np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt[w]))


def row3(c, kp, pp, Lam, sgt=-1, parts=('inst', 'ni', 't4')):
    """Group I Row III, three-rho part: integrand of int dp+/2pi int_{x,y,z,x',w,w'} (no c.c.)"""
    B = np.conj(B3one(c, kp, 'w', pp, 'z', 'x', Lam, bar=True, sgt=sgt, parts=parts))   # Bbar3^dag [b, a, j, i]
    amp = -Abar(c, pp, 'z', 'y') + np.einsum('bc,cj->bj', c.U['z'], A(c, pp, 'z', 'y'))  # [b, j]
    V = np.einsum('baji,bj->ai', B, amp)
    return phase(c) / PI3 * np.sum(conj(c, kp) * V)


def row3_4rho(c, kp, pp):
    """Group I Row III, four-rho part (the {rho,rho} entry of B_3): charges x, y2 in B_3, y in A(p)"""
    B = np.conj(B3two(c, kp, 'w', pp, 'z', 'x', 'yp', bar=True))
    amp = -Abar(c, pp, 'z', 'y') + np.einsum('bc,cj->bj', c.U['z'], A(c, pp, 'z', 'y'))
    return phase(c) / PI3 * np.sum(conj(c, kp) * np.einsum('baji,bj->ai', B, amp))


def row4(c, kp, pp):
    """Group I Row IV, three-rho part (one-rho B_2): integrand of int dp+/2pi int_{x,y,z,x',w,w'}"""
    Ad = np.conj(Abar(c, pp, 'z', 'y'))                     # Abar^dag(p,-z), charge y   [b, j]
    B = B2one(c, kp, 'w', pp, 'z', 'x')                     # B2^{bc}_{ji}(k,-w;p,-z): [b(p), c(k), j, i]
    Bb = B2one(c, kp, 'w', pp, 'z', 'x', bar=True)
    V = 2 * (np.einsum('ac,bd,bj,dcji->ai', c.U['w'], c.U['z'], Ad, B) - np.einsum('bj,baji->ai', Ad, Bb))
    return phase(c) / PI3 * np.sum(conj(c, kp) * V)


def row4_4rho(c, kp, pp):
    """Group I Row IV, four-rho part (two-rho B_2): charges y2 (gluon k), x (gluon p) in B_2, y in Abar"""
    Ad = np.conj(Abar(c, pp, 'z', 'y'))
    B = B2two(c, kp, 'w', pp, 'z', 'yp', 'x')
    Bb = B2two(c, kp, 'w', pp, 'z', 'yp', 'x', bar=True)
    V = 2 * (np.einsum('ac,bd,bj,dcji->ai', c.U['w'], c.U['z'], Ad, B) - np.einsum('bj,baji->ai', Ad, Bb))
    return phase(c) / PI3 * np.sum(conj(c, kp) * V)


def TC(kp, pp):
    """vertex bracket of C_1: d_im d_jk - (k+/p+) d_jm d_ik - (k+/(k+-p+)) d_km d_ij
    (i: gluon k+, j: gluon p+, k: gluon k+-p+)"""
    return (np.einsum('im,jk->ijkm', I2, I2) - kp / pp * np.einsum('jm,ik->ijkm', I2, I2)
            - kp / (kp - pp) * np.einsum('km,ij->ijkm', I2, I2))


def C1strip(c, kp, pp):
    """C^{bca}_{1jki}(k+, w-z; p+, z-x) without delta^2(w - z + (p/k)(z-x)): [b, c, a, j, k, i]
    b,j: gluon p+ (at x); c,k: gluon k+-p+ (at z); a,i: gluon k+ (at w)."""
    zx = c.pt['z'] - c.pt['x']
    pre = -G / (2 * PI * np.sqrt(2 * kp)) * np.sqrt(pp * (kp - pp)) / kp
    br = np.einsum('ijkm,m->ijk', TC(kp, pp), K(zx))          # [i, j, k]
    return pre * np.einsum('abc,ijk->bcajki', F, br)


def row5_part1(c, kp, pp):
    """Group I Row V, three-rho part (two-rho B_2): coefficient of delta^2(w - z + (p/k)(z-x)) in the integrand of
    int dp+/2pi int_{x,y,y',z,x',w,w'}.  w must be placed at z - (p/k)(z-x) by the caller."""
    C = C1strip(c, kp, pp)                                     # real: C^dag = C
    B = B2two(c, kp - pp, 'z', pp, 'x', 'y', 'yp')             # [b(p at x), c(k-p at z), j, k]
    Bb = B2two(c, kp - pp, 'z', pp, 'x', 'y', 'yp', bar=True)
    t1 = -np.einsum('bcajki,bcjk->ai', C, Bb)
    t2 = np.einsum('dcajki,db,ce,bejk->ai', C, c.U['x'], c.U['z'], B)
    V = 2 * (t1 + t2)
    return phase(c) / PI3 * np.sum(conj(c, kp) * V)


def row5_part2(c, kp, pp):
    """Group I Row V, two-rho part (one-rho B_2, charge y): coefficient of delta^2, as row5_part1"""
    C = C1strip(c, kp, pp)
    B = B2one(c, kp - pp, 'z', pp, 'x', 'y')
    Bb = B2one(c, kp - pp, 'z', pp, 'x', 'y', bar=True)
    t1 = -np.einsum('bcajki,bcjk->ai', C, Bb)
    t2 = np.einsum('dcajki,db,ce,bejk->ai', C, c.U['x'], c.U['z'], B)
    V = 2 * (t1 + t2)
    return phase(c) / PI3 * np.sum(conj(c, kp) * V)


# ------------------------------------------------------------------ A^(3), one-rho part (Row II, part 1)
def A3one(c, kp, W, X, bracket):
    """A^{(3)a}_{2i}(k+, -W), one-rho part with charge X; `bracket` is the number in the big square bracket"""
    pre = -1j * G**3 * NC / (16 * np.sqrt(2) * PI**3 * np.sqrt(kp))
    return pre * bracket * np.outer(c.rho[X], K(c.pt[X] - c.pt[W]))


def row2_part1(c, kp, bracket):
    A3 = A3one(c, kp, 'w', 'x', bracket)
    A3b = pre_bar = A3one(c, kp, 'w', 'x', bracket)
    A3b = np.einsum('ab,bi->ai', c.U['x'], pre_bar)
    V = np.einsum('ab,bi->ai', c.U['w'], A3) - A3b
    return phase(c) / PI3 * np.sum(conj(c, kp) * V)


# ------------------------------------------------------------------ A^(3), two-rho part (Row II, part 2), momentum form
def row2_part2(c, kp, kv, pv_, Bvec):
    """integrand of int d^2p d^2k int_{x,y,x',w,w'} for U(w) A^(3)(k,-w) - Abar^(3)(k,-w), two-rho part of the
    coordinate-space A^(3) of main.tex: e^{ik(w-y)} e^{ip(y-x)} i g^3 f^{abc}{rho^b(x),rho^c(y)}/(32 sqrt(2k) pi^5 p^2 k^2) B^i"""
    U, r, pt = c.U, c.rho, c.pt
    ph = np.exp(1j * kv @ (pt['w'] - pt['y'])) * np.exp(1j * pv_ @ (pt['y'] - pt['x']))
    pre = 1j * G**3 / (32 * np.sqrt(2 * kp) * PI**5 * (pv_ @ pv_) * (kv @ kv)) * ph
    A3 = pre * 2 * np.einsum('abc,b,c->a', F, r['x'], r['y'])
    A3b = pre * 2 * np.einsum('abc,b,c->a', F, U['x'] @ r['x'], U['y'] @ r['y'])
    V = np.outer(U['w'] @ A3 - A3b, Bvec)
    return phase(c) / PI3 * np.sum(conj(c, kp) * V)
