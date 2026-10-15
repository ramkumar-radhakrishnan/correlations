"""Transcriptions of the displayed equations of Group_II.tex and Group_III.tex (integrands).

Each function returns the integrand at fixed transverse points and fixed p+, WITHOUT int dp+/2pi and WITHOUT
the phase e^{-ik(w'-w)} (or its replacement after the w' integration), but WITH the prefactor 1/(2pi)^3 etc.
Charges are classical; {rho^a(x), rho^b(y)} = 2 rho^a(x) rho^b(y) at x != y.
Index placements are exactly as written in the notes."""
import numpy as np
from fp import F, NC, G, PI, I2, K, sq

pi3 = (2 * PI)**3


def P_(c):
    return {n: c.pt[n] for n in c.pt}


def Phit(c, k, p, W, Z, X):
    """the bracket of G2-I(1) / G3-I / G3-II:  delta/(2(p+k)(z-w)^2)[..] + 1/p[..] + 1/k[..]   [i, j]"""
    x, z, w = c.pt[X], c.pt[Z], c.pt[W]
    xz, xw, zw = x - z, x - w, z - w
    return (I2 * (sq(xz) - sq(xw)) / (2 * (p + k) * sq(zw))
            + 1 / p * (np.outer(xw, zw) / sq(zw) + np.outer(xw, xz) / (2 * sq(xz)))
            + 1 / k * (np.outer(zw, xz) / sq(zw) - np.outer(xw, xz) / (2 * sq(xw))))


def Rone(c, W, Z, X):
    """[f^{cbd} U^{ac}(w) rho^d(x) - f^{acd} U^{bc}(z) U^{de}(x) rho^e(x)]  ->  [a, b]"""
    U, r = c.U, c.rho
    return (np.einsum('cbd,ac,d->ab', F, U[W], r[X])
            - np.einsum('acd,bc,de,e->ab', F, U[Z], U[X], r[X]))


# ------------------------------------------------------------------ Group II, Row I, part 1 (line 46)
def g2_1_termA(c, k, p):
    U, r = c.U, c.rho
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['w'])
    kern = np.einsum('i,j,ij->', K(c.pt['xp'] - c.pt['wp']), K(c.pt['y'] - c.pt['z']), Phit(c, k, p, 'w', 'z', 'x')) / D
    L = (2 * np.einsum('ac,c,b->ab', U['wp'], r['xp'], r['y'])
         - 2 * np.einsum('bc,ad,ce,d,e->ab', U['z'], U['xp'], U['y'], r['xp'], r['y']))
    col = np.sum(L * Rone(c, 'w', 'z', 'x'))
    return -1 / pi3 * 1j * G**4 / (8 * PI**4) * kern * col


def g2_1_termB(c, k, p):
    U, r = c.U, c.rho
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['wp'])
    kern = np.einsum('i,j,ij->', K(c.pt['xp'] - c.pt['w']), K(c.pt['y'] - c.pt['z']), Phit(c, k, p, 'wp', 'z', 'x')) / D
    Lc = Rone(c, 'wp', 'z', 'x')          # [f^{c'bd'} U^{ac'}(w') rho^{d'}(x) - f^{ac'd'} U^{bc'}(z) U^{d'e'}(x) rho^{e'}(x)]
    R = (2 * np.einsum('ac,c,b->ab', U['w'], r['xp'], r['y'])
         - 2 * np.einsum('bc,ad,ce,d,e->ab', U['z'], U['xp'], U['y'], r['xp'], r['y']))
    return 1 / pi3 * 1j * G**4 / (8 * PI**4) * kern * np.sum(Lc * R)


# ------------------------------------------------------------------ Group II, Row I, part 2 (lines 54, 58, 62)
def PhiU(c, k, p, W, Z, X):
    from fp import Phi
    return Phi(c, k, W, p, Z, X)


def g2_2_colour54(c):
    U, r = c.U, c.rho
    L = (np.einsum('cbd,ac,d->ab', F, U['wp'], r['xp'])
         - np.einsum('acd,bc,ed,e->ab', F, U['z'], U['xp'], r['xp']))     # U^{e'd'}(x') rho^{e'}(x')
    return np.sum(L * Rone(c, 'w', 'z', 'x'))


def g2_2_colour58(c):
    U, r = c.U, c.rho
    t1 = np.einsum('Cbd,cbe,aC,ac,d,e->', F, F, U['wp'], U['w'], r['xp'], r['x'])
    t2 = np.einsum('Cbd,ace,aC,bc,eE,d,E->', F, F, U['wp'], U['z'], U['x'], r['xp'], r['x'])
    t3 = np.einsum('aCD,cbd,ac,bC,ED,E,d->', F, F, U['w'], U['z'], U['xp'], r['xp'], r['x'])
    t4 = NC * np.einsum('Ed,de,E,e->', U['xp'], U['x'], r['xp'], r['x'])
    return t1 - t2 - t3 + t4


def g2_2_pref(c, k, p):
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['w'])
    Dp = p * sq(c.pt['xp'] - c.pt['z']) + k * sq(c.pt['xp'] - c.pt['wp'])
    return 1 / pi3 * G**4 / (4 * PI**4) * p * k / (p + k)**2 / (D * Dp)


def g2_2_line54(c, k, p):
    kern = np.sum(PhiU(c, k, p, 'w', 'z', 'x') * PhiU(c, k, p, 'wp', 'z', 'xp'))
    return g2_2_pref(c, k, p) * kern * g2_2_colour54(c)


def g2_2_line58(c, k, p):
    kern = np.sum(PhiU(c, k, p, 'w', 'z', 'x') * PhiU(c, k, p, 'wp', 'z', 'xp'))
    return g2_2_pref(c, k, p) * kern * g2_2_colour58(c)


def g2_2_groups62(c, k, p, d=2):
    """the six groups of the bracket of line 62 (coefficients without the p+ factors): returns dict"""
    x, xp, z, w, wp = (c.pt[n] for n in ('x', 'xp', 'z', 'w', 'wp'))
    xz, xw, zw, xpz, xpwp, zwp = x - z, x - w, z - w, xp - z, xp - wp, z - wp
    S = sq(xz) - sq(xw)
    Sp = sq(xpz) - sq(xpwp)
    gr = {}
    gr['d'] = d * S * Sp / (4 * sq(zw) * sq(zwp))
    gr['k'] = (xz @ zw * Sp / (2 * sq(zw) * sq(zwp)) - xz @ xw * Sp / (4 * sq(xw) * sq(zwp))
               + xpz @ zwp * S / (2 * sq(zw) * sq(zwp)) - xpz @ xpwp * S / (4 * sq(zw) * sq(xpwp)))
    gr['p'] = (xw @ zw * Sp / (2 * sq(zw) * sq(zwp)) + xz @ xw * Sp / (4 * sq(xz) * sq(zwp))
               + xpwp @ zwp * S / (2 * sq(zw) * sq(zwp)) + xpz @ xpwp * S / (4 * sq(zw) * sq(xpz)))
    gr['kk'] = (xpz @ xz * (zw @ zwp) / (sq(zw) * sq(zwp)) - xpz @ xz * (xpwp @ zw) / (2 * sq(zw) * sq(xpwp))
                - xw @ zwp * (xz @ xpz) / (2 * sq(xw) * sq(zwp)) + xz @ xpz * (xpwp @ xw) / (4 * sq(xw) * sq(xpwp)))
    gr['pp'] = (xpwp @ xw * (zw @ zwp) / (sq(zw) * sq(zwp)) + xpwp @ xw * (xpz @ zw) / (2 * sq(zw) * sq(xpz))
                + xw @ xpwp * (xz @ zwp) / (2 * sq(xz) * sq(zwp)) + xz @ xpz * (xpwp @ xw) / (4 * sq(xz) * sq(xpz)))
    gr['pk'] = (xw @ zwp * (xpz @ zw) / (sq(zw) * sq(zwp)) - xw @ xpwp * (zw @ xpz) / (2 * sq(zw) * sq(xpwp))
                + xz @ xpz * (xw @ zwp) / (2 * sq(xz) * sq(zwp)) - xpz @ xz * (xpwp @ xw) / (4 * sq(xz) * sq(xpwp))
                + xpwp @ zw * (xz @ zwp) / (sq(zw) * sq(zwp)) + xpwp @ zw * (xz @ xpz) / (2 * sq(zw) * sq(xpz))
                - xz @ zwp * (xpwp @ xw) / (2 * sq(xw) * sq(zwp)) - xpwp @ xw * (xpz @ xz) / (4 * sq(xw) * sq(xpz)))
    return gr


def g2_2_line62(c, k, p):
    P = p + k
    gr = g2_2_groups62(c, k, p)
    kern = (gr['d'] + P / k * gr['k'] + P / p * gr['p'] + P**2 / k**2 * gr['kk'] + P**2 / p**2 * gr['pp']
            + P**2 / (p * k) * gr['pk'])
    return g2_2_pref(c, k, p) * kern * g2_2_colour58(c)


# ------------------------------------------------------------------ Group III, Rows I and II (lines 44, 55)
def g3_pref(c, k, p):
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['w'])
    kern = np.einsum('i,j,ij->', K(c.pt['xp'] - c.pt['wp']), K(c.pt['y'] - c.pt['z']), Phit(c, k, p, 'w', 'z', 'x')) / D
    return -1 / pi3 * 1j * G**4 / (4 * PI**4) * kern


def g3_1_line44(c, k, p):
    U, r = c.U, c.rho
    L = (np.einsum('aB,CD,Cb,B,D->ab', U['xp'], U['y'], U['z'], r['xp'], r['y'])
         - np.einsum('bD,DE,aC,C,E->ab', U['z'], U['y'], U['wp'], r['xp'], r['y']))
    return g3_pref(c, k, p) * np.sum(L * Rone(c, 'w', 'z', 'x'))


def g3_2_line55(c, k, p):
    U, r = c.U, c.rho
    L = (np.einsum('CD,aB,bC,D,B->ab', U['y'], U['xp'], U['z'], r['y'], r['xp'])
         - np.einsum('aB,b,B->ab', U['xp'], r['y'], r['xp']))
    return g3_pref(c, k, p) * np.sum(L * Rone(c, 'w', 'z', 'x'))


# ------------------------------------------------------------------ C rows: colour of the conj C' (as written)
def Lc_written(c):
    """[U^{e'c'}(y') U^{d'b}(z) rho^{e'}(x') - U^{bd'}(z) U^{c'e'}(x') rho^{e'}(x')]  -> [c', d', b]"""
    U, r = c.U, c.rho
    return (np.einsum('EC,Db,E->CDb', U['yp'], U['z'], r['xp'])
            - np.einsum('bD,CE,E->CDb', U['z'], U['xp'], r['xp']))


def Cbr(k, p):
    """[d_{k'm} d_{ij} - (P/p) d_{jm} d_{ik'} - (P/k) d_{im} d_{jk'}]  -> [k', m, i, j]"""
    P = p + k
    return (np.einsum('Km,ij->Kmij', I2, I2) - P / p * np.einsum('jm,iK->Kmij', I2, I2)
            - P / k * np.einsum('im,jK->Kmij', I2, I2))


def g3_3_part1_line68(c, k, p):
    """coefficient of delta^2[w' - ((p+k) y' - p z)/k]"""
    U, r = c.U, c.rho
    P = p + k
    kv = np.einsum('m,K,i,j,Kmij->', K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']),
                   K(c.pt['y'] - c.pt['w']), K(c.pt['x'] - c.pt['z']), Cbr(k, p))
    R2 = (2 * np.einsum('ac,c,b->ab', U['w'], r['y'], r['x'])
          - 2 * np.einsum('bc,ad,ce,d,e->ab', U['z'], U['y'], U['x'], r['y'], r['x']))
    col = np.einsum('CDa,CDb,ab->', F, Lc_written(c), R2)
    return -1 / pi3 * 1j * G**4 / (8 * PI**4) / k**2 * kv * col


def bigbracket80(c, k, p, d=2):
    """the bracket of line 80  ->  [k', m]"""
    x, z, w = c.pt['x'], c.pt['z'], c.pt['w']
    xz, xw, zw = x - z, x - w, z - w
    P = p + k
    S = sq(xz) - sq(xw)
    out = I2 * d * S / (2 * sq(zw))
    out = out + P / k * I2 * (xz @ zw / sq(zw) - xz @ xw / (2 * sq(xw)) - S / (2 * sq(zw)))
    out = out + P / p * I2 * (xw @ zw / sq(zw) + xz @ xw / (2 * sq(xz)) - S / (2 * sq(zw)))
    out = out - P**2 / p**2 * (np.outer(xw, zw) / sq(zw) + np.outer(xw, xz) / (2 * sq(xz)))      # [k', m]
    out = out - P**2 / k**2 * (np.outer(xz, zw) / sq(zw) - np.outer(xz, xw) / (2 * sq(xw)))
    out = out - P**2 / (p * k) * (np.outer(zw, xz) / sq(zw) - np.outer(xw, xz) / (2 * sq(xw))
                                  + np.outer(zw, xw) / sq(zw) + np.outer(xz, xw) / (2 * sq(xz)))
    return out


def g3_3_part2_line80(c, k, p):
    U, r = c.U, c.rho
    P = p + k
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['w'])
    kv = np.einsum('m,K,Km->', K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']), bigbracket80(c, k, p))
    col = np.einsum('CDa,CDb,ab->', F, Lc_written(c), Rone(c, 'w', 'z', 'x'))
    return -1 / pi3 * G**4 / (4 * PI**4) / k * p / (P * D) * kv * col


def g3_5_line101(c, k, p):
    U, r = c.U, c.rho
    kv = np.einsum('m,K,i,j,Kmij->', K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']),
                   K(c.pt['y'] - c.pt['w']), K(c.pt['x'] - c.pt['z']), Cbr(k, p))
    R = (np.einsum('bc,ce,ad,e,d->ab', U['z'], U['x'], U['y'], r['x'], r['y'])
         - np.einsum('ad,ce,cb,e,d->ab', U['w'], U['x'], U['z'], r['x'], r['y']))
    col = np.einsum('CDa,CDb,ab->', F, Lc_written(c), R)
    return -1 / pi3 * 1j * G**4 / (4 * PI**4) / k**2 * kv * col


def g3_6_line120(c, k, p):
    U, r = c.U, c.rho
    kv = np.einsum('m,K,i,j,Kmij->', K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']),
                   K(c.pt['y'] - c.pt['w']), K(c.pt['x'] - c.pt['z']), Cbr(k, p))
    R = (np.einsum('bc,ad,ce,d,e->ab', U['z'], U['y'], U['x'], r['y'], r['x'])
         - np.einsum('ad,d,b->ab', U['y'], r['y'], r['x']))
    col = np.einsum('CDa,CDb,ab->', F, Lc_written(c), R)
    return -1 / pi3 * 1j * G**4 / (4 * PI**4) / k**2 * kv * col


# ------------------------------------------------------------------ Group II, Row IV (line 92)
def g2_4_line92(c, k, p):
    """coefficient of delta^2[y-w+p/(p+k)(w-z)] delta^2[y'-w'+p/(p+k)(w'-z)]"""
    U, r = c.U, c.rho
    P = p + k
    kv = np.einsum('m,K,n,L,Kmij,Lnij->', K(c.pt['w'] - c.pt['z']), K(c.pt['x'] - c.pt['y']),
                   K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']), Cbr(k, p), Cbr(k, p))
    Lp = (np.einsum('CDa,EC,Db,E->ab', F, U['yp'], U['z'], r['xp'])
          - np.einsum('CDa,bD,CE,E->ab', F, U['z'], U['xp'], r['xp']))
    Ra = (np.einsum('cda,ec,db,e->ab', F, U['y'], U['z'], r['x'])
          - np.einsum('cda,bd,ce,e->ab', F, U['z'], U['x'], r['x']))
    return 1 / pi3 * G**4 / (4 * PI**4) * p * k / P**4 * kv * np.sum(Lp * Ra)
