# Ref.3 / Li-Skokov comparison: the two cubic functionals evaluated for explicit classical configurations
# (random smooth SU(3) target Wilson lines, random smooth classical projectile charges) on a periodic lattice.
#   n4_LS,even(k)  : Li-Skokov three-Wilson-line term with the real part of E_2 (k-even part; E_1 not included)
#   n4_LS,odd(k)   : the same with i Im E_2
#   n4_DN(k)       : the paper's out-state term Delta N (eta=+1), in LS normalisation
# Canonical form (ls1_kernels.py): n4 = -(1/pi) int_p Gamma^j(p;u1,u2;k) f^{cde} M^{c*}_j(k) X^d(p) Y^e(k-p) + c.c.,
#   u -> Z^m (rotated WW-type block), 1 -> rho_bar (rotated charge).  Continuum FT f(p) = int d^2x e^{-ipx} f(x).
# usage: python3 -I ls3_lattice.py
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ls1_kernels import G2, GDN, comps          # noqa: E402

def su3():
    l = np.zeros((8, 3, 3), complex)
    l[0][0, 1] = l[0][1, 0] = 1; l[1][0, 1] = -1j; l[1][1, 0] = 1j; l[2][0, 0] = 1; l[2][1, 1] = -1
    l[3][0, 2] = l[3][2, 0] = 1; l[4][0, 2] = -1j; l[4][2, 0] = 1j; l[5][1, 2] = l[5][2, 1] = 1
    l[6][1, 2] = -1j; l[6][2, 1] = 1j; l[7] = np.diag([1, 1, -2]) / np.sqrt(3)
    t = l / 2
    f = np.real(-2j * (np.einsum('aij,bjk,cki->abc', t, t, t) - np.einsum('bij,ajk,cki->abc', t, t, t)))
    return t, f

def main(seed=31):
    print('=== configuration seed', seed)
    rng = np.random.default_rng(seed)
    N, Lbox = 48, 12.0
    a = Lbox / N
    m = np.fft.fftfreq(N, d=1.0 / N)
    P = 2 * np.pi * m / Lbox
    PX, PY = np.meshgrid(P, P, indexing='ij')
    P2 = PX ** 2 + PY ** 2; P2s = np.where(P2 == 0, 1.0, P2)
    t, f = su3()
    def smooth(field, s):
        return np.real(np.fft.ifft2(np.fft.fft2(field, axes=(-2, -1)) * np.exp(-0.5 * s * s * P2), axes=(-2, -1)))
    # target: V(x) = exp(i alpha^a t^a), smooth, dense
    alpha = smooth(rng.normal(size=(8, N, N)), 0.6); alpha *= 3.0 / alpha.std()
    H = np.einsum('axy,aij->xyij', alpha, t)
    w, vecs = np.linalg.eigh(H)
    V = np.einsum('xyij,xyj,xykj->xyik', vecs, np.exp(1j * w), vecs.conj())
    U = 2 * np.real(np.einsum('aij,xyjk,bkl,xyil->xyab', t, V, t, V.conj()))       # U^{ab} = 2 tr(t^a V t^b V^dag)
    print('unitarity of adjoint U: %.1e' % np.abs(np.einsum('xyab,xycb->xyac', U, U) - np.eye(8)).max())
    # projectile charges, smooth, neutral
    rho = smooth(rng.normal(size=(8, N, N)), 0.5); rho -= rho.mean(axis=(1, 2), keepdims=True); rho /= rho.std()
    FT = lambda fx: a * a * np.fft.fft2(fx, axes=(-2, -1)) * 1.0      # e^{-ipx} convention of numpy
    iFT = lambda fp: np.fft.ifft2(fp, axes=(-2, -1)) / (a * a)
    rp = FT(rho)
    # blocks: rho_bar^d(p) = FT[rho^e U^{ed}], Z^{m,d}(p) = FT[w^{m,e} U^{ed}], w^m(p1) = p1^m/p1^2 rho(p1)
    rot = lambda fx: np.einsum('exy,xyed->dxy', fx, U)
    rb = FT(rot(rho))
    wx = [iFT(np.where(P2 == 0, 0, comp / P2s) * rp) for comp in (PX, PY)]
    Z = [FT(rot(wm)) for wm in wx]                                            # w^m(x) = i x real (odd kernel)
    print('w(x) purely imaginary (odd kernel): max |Re w| / max |Im w| = %.1e' % max(np.abs(np.real(wm)).max() / np.abs(np.imag(wm)).max() for wm in wx))
    blocks = [rb, Z[0], Z[1]]
    idx = lambda v: (int(round(v[0] * Lbox / (2 * np.pi))) % N, int(round(v[1] * Lbox / (2 * np.pi))) % N)
    def M(k):
        i0 = idx(k)
        return np.array([Z[0][:, i0[0], i0[1]] - k[0] / (k @ k) * rb[:, i0[0], i0[1]],
                         Z[1][:, i0[0], i0[1]] - k[1] / (k @ k) * rb[:, i0[0], i0[1]]])          # [j, c]
    # kernel component index -> (A,B) block labels: '00'->(0,0), '1m'->(m,0), '2n'->(0,n), '1m2n'->(m,n)
    AB = [(0, 0), (1, 0), (2, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
    def n4(k, kern):
        Mk = M(k); T = 0.0 + 0j
        for ix in range(N):
            for iy in range(N):
                p = np.array([P[ix], P[iy]]); q = k - p
                if p @ p == 0 or q @ q == 0 or abs(p[0] * k[1] - p[1] * k[0]) < 1e-12:
                    continue                        # singular kernel points / collinear sign line (measure zero)
                C = comps(kern, p, k)               # (9, 2)
                iq = idx(q)
                for c_i, (A, B) in enumerate(AB):
                    X = blocks[A][:, ix, iy]; Y = blocks[B][:, iq[0], iq[1]]
                    col = np.einsum('cde,jc,d,e->j', f, Mk.conj(), X, Y)
                    T += C[c_i] @ col
        T *= -(1 / np.pi) / Lbox ** 2
        return 2 * T.real
    kre = lambda p, u1, u2, k: G2(p, u1, u2, k, part='re').real
    kim = lambda p, u1, u2, k: 1j * G2(p, u1, u2, k, part='im').imag
    out = []
    for km in [(2, 1), (1, 3), (3, -2)]:
        k = np.array(km, float) * 2 * np.pi / Lbox
        row = []
        for s in (1, -1):
            kk = s * k
            row.append((n4(kk, kre), n4(kk, kim), n4(kk, GDN)))
        (er, eo, dn), (er_, eo_, dn_) = row
        lo = (np.abs(M(k)) ** 2).sum() / np.pi
        print('k = %s (2pi/L):  LO n2 = %.4e' % (km, lo))
        print('   LS Re E_2 part : n4(k) = %+.4e, n4(-k) = %+.4e   (k-even check)' % (er, er_))
        print('   LS Im E_2 part : n4(k) = %+.4e, n4(-k) = %+.4e   (k-odd check)' % (eo, eo_))
        print('   Delta N (eta=1): n4(k) = %+.4e, n4(-k) = %+.4e' % (dn, dn_))
        print('   ratio Delta N / LS_even = %+.3f' % (dn / er))
        out.append(dn / er)
    print('ratios over the three k:', np.round(out, 3))

if __name__ == '__main__':
    for sd in (31, 97):
        main(sd)
