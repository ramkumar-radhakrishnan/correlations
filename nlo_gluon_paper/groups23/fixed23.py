"""Corrected index placements for a general adjoint U (U^dag = U^T, U^{bc} != U^{cb}).
Only the colour factors change; kernels and prefactors are those of user23."""
import numpy as np
from fp import F, NC, G, PI, K, sq
import user23 as us

pi3 = us.pi3


def Rone_fix(c, W, Z, X):
    """[f^{cbd} U^{ac}(w) rho^d(x) - f^{acd} U^{cb}(z) U^{de}(x) rho^e(x)]"""
    U, r = c.U, c.rho
    return (np.einsum('cbd,ac,d->ab', F, U[W], r[X]) - np.einsum('acd,cb,de,e->ab', F, U[Z], U[X], r[X]))


def two_fix(c, W, Yk, Xp):
    """[U^{ac}(w){rho^c(Yk),rho^b(Xp)} - U^{cb}(z) U^{ad}(Yk) U^{ce}(Xp) {rho^d(Yk),rho^e(Xp)}]"""
    U, r = c.U, c.rho
    return (2 * np.einsum('ac,c,b->ab', U[W], r[Yk], r[Xp])
            - 2 * np.einsum('cb,ad,ce,d,e->ab', U['z'], U[Yk], U[Xp], r[Yk], r[Xp]))


def g2_1_termA(c, k, p):
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['w'])
    kern = np.einsum('i,j,ij->', K(c.pt['xp'] - c.pt['wp']), K(c.pt['y'] - c.pt['z']), us.Phit(c, k, p, 'w', 'z', 'x')) / D
    return -1 / pi3 * 1j * G**4 / (8 * PI**4) * kern * np.sum(two_fix(c, 'wp', 'xp', 'y') * Rone_fix(c, 'w', 'z', 'x'))


def g2_1_termB(c, k, p):
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['wp'])
    kern = np.einsum('i,j,ij->', K(c.pt['xp'] - c.pt['w']), K(c.pt['y'] - c.pt['z']), us.Phit(c, k, p, 'wp', 'z', 'x')) / D
    return 1 / pi3 * 1j * G**4 / (8 * PI**4) * kern * np.sum(Rone_fix(c, 'wp', 'z', 'x') * two_fix(c, 'w', 'xp', 'y'))


def g2_2_colour(c):
    """[f^{c'bd'} U^{ac'}(w') rho^{d'}(x') - f^{ac'd'} U^{c'b}(z) U^{d'e'}(x') rho^{e'}(x')] x Rone_fix"""
    return np.sum(Rone_fix(c, 'wp', 'z', 'xp') * Rone_fix(c, 'w', 'z', 'x'))


def g2_2_colour58(c):
    U, r = c.U, c.rho
    t1 = np.einsum('Cbd,cbe,aC,ac,d,e->', F, F, U['wp'], U['w'], r['xp'], r['x'])
    t2 = np.einsum('Cbd,ace,aC,cb,eE,d,E->', F, F, U['wp'], U['z'], U['x'], r['xp'], r['x'])
    t3 = np.einsum('aCD,cbd,ac,Cb,DE,E,d->', F, F, U['w'], U['z'], U['xp'], r['xp'], r['x'])
    t4 = NC * np.einsum('dE,de,E,e->', U['xp'], U['x'], r['xp'], r['x'])
    return t1 - t2 - t3 + t4


def g2_2_line54(c, k, p):
    kern = np.sum(us.PhiU(c, k, p, 'w', 'z', 'x') * us.PhiU(c, k, p, 'wp', 'z', 'xp'))
    return us.g2_2_pref(c, k, p) * kern * g2_2_colour(c)


def g3_1_line44(c, k, p):
    U, r = c.U, c.rho
    L = (np.einsum('aB,CD,Cb,B,D->ab', U['xp'], U['y'], U['z'], r['xp'], r['y'])
         - np.einsum('Db,DE,aC,C,E->ab', U['z'], U['y'], U['wp'], r['xp'], r['y']))
    return us.g3_pref(c, k, p) * np.sum(L * Rone_fix(c, 'w', 'z', 'x'))


def g3_2_line55(c, k, p):
    U, r = c.U, c.rho
    L = (np.einsum('CD,aB,Cb,D,B->ab', U['y'], U['xp'], U['z'], r['y'], r['xp'])
         - np.einsum('aB,b,B->ab', U['xp'], r['y'], r['xp']))
    return us.g3_pref(c, k, p) * np.sum(L * Rone_fix(c, 'w', 'z', 'x'))


def Lc_fix(c):
    """[U^{c'e'}(y') U^{d'b}(z) rho^{e'}(x') - U^{d'b}(z) U^{c'e'}(x') rho^{e'}(x')]"""
    U, r = c.U, c.rho
    return (np.einsum('CE,Db,E->CDb', U['yp'], U['z'], r['xp']) - np.einsum('Db,CE,E->CDb', U['z'], U['xp'], r['xp']))


def g3_3_part1_line68(c, k, p, sign=-1):
    kv = np.einsum('m,K,i,j,Kmij->', K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']),
                   K(c.pt['y'] - c.pt['w']), K(c.pt['x'] - c.pt['z']), us.Cbr(k, p))
    col = np.einsum('CDa,CDb,ab->', F, Lc_fix(c), two_fix(c, 'w', 'y', 'x'))
    return sign * (-1) / pi3 * 1j * G**4 / (8 * PI**4) / k**2 * kv * col


def g3_3_part2_line80(c, k, p, sign=-1):
    P = p + k
    D = p * sq(c.pt['x'] - c.pt['z']) + k * sq(c.pt['x'] - c.pt['w'])
    kv = np.einsum('m,K,Km->', K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']), us.bigbracket80(c, k, p))
    col = np.einsum('CDa,CDb,ab->', F, Lc_fix(c), Rone_fix(c, 'w', 'z', 'x'))
    return sign * (-1) / pi3 * G**4 / (4 * PI**4) / k * p / (P * D) * kv * col


def X1amp_fix(c):
    """[U^{cb}(z) U^{ce}(x) U^{ad}(y) rho^e(x) rho^d(y) - U^{ad}(w) U^{ce}(x) U^{cb}(z) rho^e(x) rho^d(y)]"""
    U, r = c.U, c.rho
    return (np.einsum('cb,ce,ad,e,d->ab', U['z'], U['x'], U['y'], r['x'], r['y'])
            - np.einsum('ad,ce,cb,e,d->ab', U['w'], U['x'], U['z'], r['x'], r['y']))


def X2amp_fix(c):
    U, r = c.U, c.rho
    return (np.einsum('cb,ad,ce,d,e->ab', U['z'], U['y'], U['x'], r['y'], r['x']) - np.einsum('ad,d,b->ab', U['y'], r['y'], r['x']))


def g3_5_line101(c, k, p, sign=-1):
    kv = np.einsum('m,K,i,j,Kmij->', K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']),
                   K(c.pt['y'] - c.pt['w']), K(c.pt['x'] - c.pt['z']), us.Cbr(k, p))
    col = np.einsum('CDa,CDb,ab->', F, Lc_fix(c), X1amp_fix(c))
    return sign * (-1) / pi3 * 1j * G**4 / (4 * PI**4) / k**2 * kv * col


def g3_6_line120(c, k, p, sign=-1):
    kv = np.einsum('m,K,i,j,Kmij->', K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']),
                   K(c.pt['y'] - c.pt['w']), K(c.pt['x'] - c.pt['z']), us.Cbr(k, p))
    col = np.einsum('CDa,CDb,ab->', F, Lc_fix(c), X2amp_fix(c))
    return sign * (-1) / pi3 * 1j * G**4 / (4 * PI**4) / k**2 * kv * col


def g2_4_line92(c, k, p):
    U, r = c.U, c.rho
    P = p + k
    kv = np.einsum('m,K,n,L,Kmij,Lnij->', K(c.pt['w'] - c.pt['z']), K(c.pt['x'] - c.pt['y']),
                   K(c.pt['wp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']), us.Cbr(k, p), us.Cbr(k, p))
    Lp = np.einsum('CDa,CDb->ab', F, Lc_fix(c))
    Ra = (np.einsum('cda,ce,db,e->ab', F, U['y'], U['z'], r['x']) - np.einsum('cda,db,ce,e->ab', F, U['z'], U['x'], r['x']))
    return 1 / pi3 * G**4 / (4 * PI**4) * p * k / P**4 * kv * np.sum(Lp * Ra)
