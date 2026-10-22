"""Transcription of Finite_terms.tex (as received), term by term.

Every function returns the integrand at fixed transverse points of a configuration c (fp.Cfg): random
points, classical charges (so {rho^a(u), rho^b(v)} = 2 rho^a(u) rho^b(v) at u != v), adjoint SU(3) Wilson
lines (symmetric U = U^T unless stated; Finite_terms.tex uses U^{bc} = U^{cb}).  Included: 1/(2pi)^3, the
prefactors, the phase.  Not included: the transverse integrals.  Where Finite_terms.tex keeps a p+ (or q+)
integral, the function `..._p(c, kp, pp)` returns the integrand of int dp+/2pi and `..._rest` the remaining
terms; the plus prescriptions are taken as written ([k/p]_+ f = (k/p)(f(p) - f(0)) on 0 < p < upper limit).
k+ = kp, vee = V, Lambda = L.  The logs are taken as log|.| (principal value, see the note)."""
import numpy as np
import env  # noqa: F401
from fp import F, NC, G, PI, I2, sq

PI3 = (2 * PI)**3
LOG = lambda t: np.log(np.abs(t))


def vec(c, a, b):
    return c.pt[a] - c.pt[b]


def Kv(c, a, b):
    d = vec(c, a, b)
    return d / sq(d)


def phase(c, a, b):
    """e^{-ik.(b-a)}"""
    return np.exp(-1j * c.kv @ (c.pt[b] - c.pt[a]))


# ------------------------------------------------------------------ colour factors (classical charges)
def CL(c):
    """[U^{ab'}(x') rho^{b'}(x') - U^{ab'}(w') rho^{b'}(x')]  -> vector a"""
    return (c.U['xp'] - c.U['wp']) @ c.rho['xp']


def col_G1_III(c):
    """f^{abd} CL^a [U^{de}(x)U^{bc}(z) rho^e(x) rho^c(y) - U^{de}(x) U^{bc}(y) rho^e(x) rho^c(y)]"""
    U, r = c.U, c.rho
    return np.einsum('abd,a,b,d->', F, CL(c), U['z'] @ r['y'] - U['y'] @ r['y'], U['x'] @ r['x'])


def col_G1_IV(c):
    """CL^a [f^{cdc'} U^{ac}(w) U^{bd}(z) U^{be}(y) rho^e(y) rho^{c'}(x) - f^{abc} U^{be}(y) U^{cd}(x) rho^e(y) rho^d(x)]"""
    U, r = c.U, c.rho
    t1 = np.einsum('cdC,c,d,C->', F, U['w'].T @ CL(c), U['z'].T @ (U['y'] @ r['y']), r['x'])
    t2 = np.einsum('abc,a,b,c->', F, CL(c), U['y'] @ r['y'], U['x'] @ r['x'])
    return t1 - t2


def col_G1_II_2a(c):
    """(CL^a U^{ah}(w) f^{hbc}{rho^b(x),rho^c(y)},  CL^a U^{bd}(x) U^{ce}(y) f^{abc}{rho^d(x),rho^e(y)})"""
    U, r = c.U, c.rho
    a = 2 * np.einsum('hbc,h,b,c->', F, U['w'].T @ CL(c), r['x'], r['y'])
    b = 2 * np.einsum('abc,a,b,c->', F, CL(c), U['x'] @ r['x'], U['y'] @ r['y'])
    return a, b


def col_G1_V1(c):
    """CL^a [f^{adc} U^{db}(x) U^{ce}(z){rho^e(y),rho^b(y')} - f^{abc} U^{cd}(y) U^{be}(y'){rho^d(y),rho^e(y')}]"""
    U, r = c.U, c.rho
    t1 = 2 * np.einsum('adc,a,d,c->', F, CL(c), U['x'] @ r['yp'], U['z'] @ r['y'])
    t2 = 2 * np.einsum('abc,a,b,c->', F, CL(c), U['yp'] @ r['yp'], U['y'] @ r['y'])
    return t1 - t2


def Rone(c, W, Z, X):
    """[f^{cbd} U^{ac}(W) rho^d(X) - f^{acd} U^{bc}(Z) U^{de}(X) rho^e(X)]  -> [a, b]"""
    U, r = c.U, c.rho
    return (np.einsum('cbd,ac,d->ab', F, U[W], r[X]) - np.einsum('acd,bc,de,e->ab', F, U[Z], U[X], r[X]))


def col_G2_I1(c, W, Wo):
    """[U^{ac'}(Wo){rho^{c'}(x'),rho^b(y)} - U^{bc'}(z) U^{ad'}(x') U^{c'e'}(y){rho^{d'}(x'),rho^{e'}(y)}] . Rone(W, z, x)
    (block A: W = w, Wo = w'; block B: W = w', Wo = w)"""
    U, r = c.U, c.rho
    L = (2 * np.einsum('ac,c,b->ab', U[Wo], r['xp'], r['y'])
         - 2 * np.einsum('bc,ad,ce,d,e->ab', U['z'], U['xp'], U['y'], r['xp'], r['y']))
    return np.sum(L * Rone(c, W, 'z', 'x'))


def col_G2_I2(c):
    """[f^{c'bd'} f^{cbd} U^{ac'}(w') U^{ac}(w) rho^{d'}(x') rho^d(x) - f^{c'bd'} f^{acd} U^{ac'}(w') U^{bc}(z) U^{de}(x) ..
        - f^{ac'd'} f^{cbd} U^{ac}(w) U^{bc'}(z) U^{e'd'}(x') .. + N_c U^{e'd}(x') U^{de}(x) rho^{e'}(x') rho^e(x)]"""
    U, r = c.U, c.rho
    t1 = np.einsum('Cbd,cbe,aC,ac,d,e->', F, F, U['wp'], U['w'], r['xp'], r['x'])
    t2 = np.einsum('Cbd,ace,aC,bc,eE,d,E->', F, F, U['wp'], U['z'], U['x'], r['xp'], r['x'])
    t3 = np.einsum('aCD,cbd,ac,bC,ED,E,d->', F, F, U['w'], U['z'], U['xp'], r['xp'], r['x'])
    t4 = NC * np.einsum('Ed,de,E,e->', U['xp'], U['x'], r['xp'], r['x'])
    return t1 - t2 - t3 + t4


def col_G3_I(c):
    U, r = c.U, c.rho
    L = (np.einsum('aB,CD,Cb,B,D->ab', U['xp'], U['y'], U['z'], r['xp'], r['y'])
         - np.einsum('bD,DE,aC,C,E->ab', U['z'], U['y'], U['wp'], r['xp'], r['y']))
    return np.sum(L * Rone(c, 'w', 'z', 'x'))


def col_G3_II(c):
    U, r = c.U, c.rho
    L = (np.einsum('CD,aB,bC,D,B->ab', U['y'], U['xp'], U['z'], r['y'], r['xp'])
         - np.einsum('aB,b,B->ab', U['xp'], r['y'], r['xp']))
    return np.sum(L * Rone(c, 'w', 'z', 'x'))


def Lc(c):
    """[U^{e'c'}(y') U^{d'b}(z) rho^{e'}(x') - U^{bd'}(z) U^{c'e'}(x') rho^{e'}(x')]  -> [c', d', b]"""
    U, r = c.U, c.rho
    return (np.einsum('EC,Db,E->CDb', U['yp'], U['z'], r['xp']) - np.einsum('bD,CE,E->CDb', U['z'], U['xp'], r['xp']))


def col_G3_III1(c):
    """f^{c'd'a} Lc x [U^{ac}(w){rho^c(y),rho^b(x)} - U^{bc}(z) U^{ad}(y) U^{ce}(x){rho^d(y),rho^e(x)}]"""
    U, r = c.U, c.rho
    R = (2 * np.einsum('ac,c,b->ab', U['w'], r['y'], r['x'])
         - 2 * np.einsum('bc,ad,ce,d,e->ab', U['z'], U['y'], U['x'], r['y'], r['x']))
    return np.einsum('CDa,CDb,ab->', F, Lc(c), R)


def col_G3_III2(c):
    return np.einsum('CDa,CDb,ab->', F, Lc(c), Rone(c, 'w', 'z', 'x'))


def col_G3_V(c):
    """f^{c'ca} [U^{e'c'}(y') rho^{e'}(x') - U^{c'e'}(x') rho^{e'}(x')] [U^{ce}(x) U^{ad}(y) rho^e(x) rho^d(y)
        - U^{ad}(w) U^{ce}(x) rho^e(x) rho^d(y)]"""
    U, r = c.U, c.rho
    L = U['yp'].T @ r['xp'] - U['xp'] @ r['xp']
    R = np.einsum('ce,ad,e,d->ca', U['x'], U['y'], r['x'], r['y']) - np.einsum('ad,ce,e,d->ca', U['w'], U['x'], r['x'], r['y'])
    return np.einsum('Cca,C,ca->', F, L, R)


def col_G3_VI(c):
    """f^{c'd'a} Lc x [U^{bc}(z) U^{ad}(y) U^{ce}(x) rho^d(y) rho^e(x) - U^{ad}(y) rho^d(y) rho^b(x)]"""
    U, r = c.U, c.rho
    R = (np.einsum('bc,ad,ce,d,e->ab', U['z'], U['y'], U['x'], r['y'], r['x']) - np.einsum('ad,d,b->ab', U['y'], r['y'], r['x']))
    return np.einsum('CDa,CDb,ab->', F, Lc(c), R)


# ------------------------------------------------------------------ transverse building blocks
def AB(c, X='x', W='w'):
    return sq(vec(c, X, 'z')), sq(vec(c, X, W))


def br4(c, W='w', X='x'):
    """[delta/(2(z-w)^2) + (x-z)^j(z-w)^i/((z-w)^2(x-z)^2) - (x-z)^j(x-w)^i/((x-w)^2(x-z)^2)
        - (x-w)^i(z-w)^j/((z-w)^2(x-w)^2)]  [i, j]"""
    xz, xw, zw = vec(c, X, 'z'), vec(c, X, W), vec(c, 'z', W)
    return (I2 / (2 * sq(zw)) + np.outer(zw, xz) / (sq(zw) * sq(xz)) - np.outer(xw, xz) / (sq(xw) * sq(xz))
            - np.outer(xw, zw) / (sq(zw) * sq(xw)))


def brlog(c, W='w', X='x'):
    """[(x-w)^i(z-w)^j/((z-w)^2(x-w)^2) + (x-z)^j(x-w)^i/(2(x-z)^2(x-w)^2)]  [i, j]"""
    xz, xw, zw = vec(c, X, 'z'), vec(c, X, W), vec(c, 'z', W)
    return np.outer(xw, zw) / (sq(zw) * sq(xw)) + np.outer(xw, xz) / (2 * sq(xz) * sq(xw))


def logIV(c, kp, V, L, X='x', W='w'):
    A, B = AB(c, X, W)
    return LOG(((V - kp) * A + kp * B) / (L * A + kp * B))


def ij(a, b, M):
    return np.einsum('i,j,ij->', a, b, M)


# ================================================================== GROUP I
def G1_II_2a(c, kp, V, L):
    """Row II, d^3N_2a (the spurious int_z dropped): returns (f^{hbc} term, f^{abc} term)"""
    x, y, w = c.pt['x'], c.pt['y'], c.pt['w']
    Lv = Kv(c, 'x', 'w') * np.log(sq(w - y) / sq(x - y)) - Kv(c, 'y', 'w') * np.log(sq(w - x) / sq(x - y))
    k = Kv(c, 'xp', 'wp') @ Lv
    ca, cb = col_G1_II_2a(c)
    pre = phase(c, 'w', 'wp') / PI3 * 1j * G**4 / (32 * PI**4) / kp * np.log(V / kp) * k
    return -pre * ca, pre * cb


def G1_III(c, kp, V, L):
    """Row III: N1, N2, N3a, N3b (finite) -> dict"""
    A, B = AB(c)
    xz, xw, zw, yz = vec(c, 'x', 'z'), vec(c, 'x', 'w'), vec(c, 'z', 'w'), vec(c, 'y', 'z')
    Xp, Yz, Zw, Xw, Xz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'z', 'w'), Kv(c, 'x', 'w'), Kv(c, 'x', 'z')
    p16, p8 = 1j * G**4 / (16 * PI**5) / kp, 1j * G**4 / (8 * PI**5) / kp
    XY = Xp @ Yz
    out = {}
    # N1
    n1 = (-p16 * np.log(1 - kp / V) * XY * (A + B) / (A - B)**2
          + p16 * XY * (B + A) / (A - B)**2 * LOG((V * A - kp * B) / ((L + kp) * A - kp * B))
          - p16 * np.log(V / (kp - L)) * XY * (A + B) / (A - B)**2
          + p16 * XY * (B + A) / (A - B)**2 * LOG(((kp - L) * A - kp * B) / (L * A - kp * B))
          - p8 * (Xp @ Xw) * (xz @ yz) / (A * sq(yz)) * LOG((kp * B - V * A) / (kp * B - L * A)))
    out['N1'] = n1
    # N2 (overall minus)
    L2 = LOG((B * kp - A * V) / (B * kp - A * (kp + L)))
    n2 = (-p8 * np.log(V / (L + kp)) * (Xp @ Xw) * (Yz @ Zw)
          + p8 * np.log(1 - kp / V) * XY * (zw @ xz) * B / (sq(zw) * (B - A)**2)
          - p8 * np.log(1 - kp / V) * XY * (zw @ xw) * A / (sq(zw) * (B - A)**2)
          + p8 * L2 * (Xp @ Zw) * (Yz @ xw) / (B - A)
          - p8 * L2 * (Xp @ Xw) * (Yz @ Zw) * A / (B - A)
          + p8 * L2 * XY * (zw @ xw) * A / (sq(zw) * (B - A)**2)
          - p8 * L2 * (Xp @ Zw) * (Yz @ Xz) * B / (B - A)
          + p8 * L2 * (Xp @ xz) * (Yz @ Zw) / (B - A)
          - p8 * L2 * XY * (zw @ xz) * B / (sq(zw) * (B - A)**2))
    out['N2'] = -n2
    # N3a
    L3 = LOG((B * kp - A * (kp - L)) / (B * kp - A * L))
    lkL = np.log((kp - L) / V)
    root = np.sqrt(A * B - (xw @ xz)**2)
    calA = (np.arctan(((kp - L) * A - kp * (xw @ xz)) / (kp * root)) - np.arctan((L * A - kp * (xw @ xz)) / (kp * root))) / root
    n3a = (-p8 * lkL * XY / sq(zw)
           - p8 * lkL * XY * (zw @ xw) * A / (sq(zw) * (B - A)**2)
           + p8 * lkL * XY * (zw @ xz) * B / (sq(zw) * (B - A)**2)
           - p8 * L3 * XY * (zw @ xw) * A / (sq(zw) * (B - A)**2)
           + p8 * L3 * XY * (zw @ xz) * B / (sq(zw) * (B - A)**2)
           - p8 * L3 * (Xp @ Zw) * (Yz @ xw) / (B - A)
           + p8 * L3 * (Xp @ Zw) * (Yz @ Xz) * B / (B - A)
           + p8 * L3 * (Yz @ Zw) * (Xp @ Xw) * A / (B - A)
           - p8 * L3 * (Yz @ Zw) * (Xp @ xz) / (B - A)
           + p8 * calA * (Xp @ Zw) * (Yz @ xw)
           + p8 * calA * (Yz @ Zw) * (Xp @ xz)
           - p8 * calA * (xw @ xz) * (Xp @ Zw) * (Yz @ xz) / A)
    out['N3a'] = n3a
    L4 = LOG(sq(kp * xw - (kp - L) * xz) / sq(kp * xw - L * xz))
    n3b = (-p8 * calA * (xw @ xz) * (Yz @ Zw) * (Xp @ xw) / B
           - p16 * L4 * (Xp @ Zw) * (Yz @ xz) / A
           + p16 * L4 * (Yz @ Zw) * (Xp @ xw) / B
           - p16 * L4 * XY / sq(zw))
    out['N3b'] = n3b
    pc = phase(c, 'w', 'wp') / PI3 * col_G1_III(c)
    return {k: v * pc for k, v in out.items()}


def G1_IV(c, kp, V, L):
    Xp, Yz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z')
    t = (1j * G**4 / (8 * PI**5) / kp * ij(Xp, Yz, br4(c)) * logIV(c, kp, V, L)
         - 1j * G**4 / (16 * PI**5) / kp * (Xp @ Yz) / sq(vec(c, 'z', 'w')) * np.log(V / (L + kp))
         + 1j * G**4 / (8 * PI**5) / kp * np.log(1 - kp / V) * ij(Xp, Yz, brlog(c)))
    return t * phase(c, 'w', 'wp') / PI3 * col_G1_IV(c)


def G1_V1_kernels(c):
    Xp, Yz, Zx, Ypx = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'z', 'x'), Kv(c, 'yp', 'x')
    return dict(a=(Xp @ Zx) * (Yz @ Ypx), b=(Ypx @ Zx) * (Xp @ Yz), q=(Xp @ Ypx) * (Zx @ Yz))


def G1_V1_p(c, kp, pp, which):
    """integrand of int_0^k dp/2pi (which='p') or int_0^k dq/2pi (which='q'), plus prescription included"""
    kk = G1_V1_kernels(c)
    col = col_G1_V1(c)
    if which == 'p':
        e = lambda p: np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z'])) * np.exp(-1j * (c.kv @ vec(c, 'z', 'x')) * p / kp)
        val = kk['a'] * e(pp) - kp / pp * kk['b'] * (e(pp) - e(0.0))
        return 1 / PI3 * 1j * G**4 / (16 * PI**4) / kp**2 * val * col
    e = lambda q: np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x'])) * np.exp(-1j * (c.kv @ vec(c, 'x', 'z')) * q / kp)
    return -1 / PI3 * 1j * G**4 / (16 * PI**4) / kp * kk['q'] * (e(pp) - e(0.0)) / pp * col


def G1_V1_rest(c, kp, V, L):
    kk = G1_V1_kernels(c)
    col = col_G1_V1(c)
    lg = np.log((kp - L) / V)
    t3 = -1 / PI3 * np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z'])) * 1j * G**4 / (32 * PI**5) / kp * lg * kk['b']
    t4 = -1 / PI3 * np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x'])) * 1j * G**4 / (32 * PI**5) / kp * lg * kk['q']
    return (t3 + t4) * col


# ================================================================== GROUP II
def G2_I1(c, kp, V, L):
    """Row I, first part: (block A, block B)"""
    Xp, Yz, Xpw = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'xp', 'w')
    p16, p32 = 1j * G**4 / (16 * PI**5) / kp, 1j * G**4 / (32 * PI**5) / kp
    lk, l1 = np.log(V / (L + kp)), np.log(1 - kp / V)
    kA = (-p16 * ij(Xp, Yz, br4(c, 'w')) * logIV(c, kp, V, L, 'x', 'w')
          + p32 * lk * (Xp @ Yz) / sq(vec(c, 'z', 'w'))
          - p16 * l1 * ij(Xp, Yz, brlog(c, 'w')))
    kB = (p16 * ij(Xpw, Yz, br4(c, 'wp')) * logIV(c, kp, V, L, 'x', 'wp')
          - p32 * lk * (Xpw @ Yz) / sq(vec(c, 'z', 'wp'))
          + p16 * l1 * ij(Xpw, Yz, brlog(c, 'wp')))
    ph = phase(c, 'w', 'wp') / PI3
    return ph * kA * col_G2_I1(c, 'w', 'wp'), ph * kB * col_G2_I1(c, 'wp', 'w')


def G2_I2(c, kp, V, L, d=2):
    """Row I, second part, 2a .. 2f -> dict"""
    x, xp, z, w, wp = (c.pt[n] for n in ('x', 'xp', 'z', 'w', 'wp'))
    xz, xw, zw, xpz, xpwp, zwp = x - z, x - w, z - w, xp - z, xp - wp, z - wp
    A, B, Ap, Bp, zw2, zwp2 = sq(xz), sq(xw), sq(xpz), sq(xpwp), sq(zw), sq(zwp)
    LX = LOG(((V - kp) * A + kp * B) / (L * A + kp * B))
    LXp = LOG(((V - kp) * Ap + kp * Bp) / (L * Ap + kp * Bp))
    l1, lk = np.log(1 - kp / V), np.log(V / (L + kp))
    D = Bp * A - B * Ap
    # zw' = z - w'
    br2a = ((xw @ zwp) * (xpz @ zw) / (zw2 * zwp2) - (xw @ xpwp) * (zw @ xpz) / (2 * zw2 * Bp)
            + (xz @ xpz) * (xw @ zwp) / (2 * A * zwp2) - (xpz @ xz) * (xpwp @ xw) / (4 * A * Bp)
            + (xpwp @ zw) * (xz @ zwp) / (zw2 * zwp2) + (xpwp @ zw) * (xz @ xpz) / (2 * zw2 * Ap)
            - (xz @ zwp) * (xpwp @ xw) / (2 * B * zwp2) - (xpwp @ xw) * (xpz @ xz) / (4 * B * Ap))
    br2b = ((xpwp @ xw) * (zw @ zwp) / (zw2 * zwp2) + (xpwp @ xw) * (xpz @ zw) / (2 * zw2 * Ap)
            + (xw @ xpwp) * (xz @ zwp) / (2 * A * zwp2) + (xz @ xpz) * (xpwp @ xw) / (4 * A * Ap))
    br2c = ((xpz @ xz) * (zw @ zwp) / (zw2 * zwp2) - (xpz @ xz) * (xpwp @ zw) / (2 * zw2 * Bp)
            - (xw @ zwp) * (xz @ xpz) / (2 * B * zwp2) + (xz @ xpz) * (xpwp @ xw) / (4 * B * Bp))
    br2d = ((xw @ zw) * (Ap - Bp) / (2 * zw2 * zwp2) + (xz @ xw) * (Ap - Bp) / (4 * A * zwp2)
            + (xpwp @ zwp) * (A - B) / (2 * zw2 * zwp2) + (xpz @ xpwp) * (A - B) / (4 * zw2 * Ap))
    br2e = ((xz @ zw) * (Ap - Bp) / (2 * zw2 * zwp2) - (xz @ xw) * (Ap - Bp) / (4 * B * zwp2)
            + (xpz @ zwp) * (A - B) / (2 * zw2 * zwp2) - (xpz @ xpwp) * (A - B) / (4 * zw2 * Bp))
    pre = G**4 / (8 * PI**5) / kp
    out = {}
    out['2a'] = pre / D * br2a * (LOG(((V - kp) * A + kp * B) / (kp * Bp + (V - kp) * Ap))
                                  - LOG((L * A + kp * B) / (kp * Bp + L * Ap)))
    out['2b'] = (pre * l1 / (Bp * A * B - B**2 * Ap) * br2b * (A - B * Ap / Bp)
                 - pre * A / (Bp * A * B - B**2 * Ap) * br2b * LX
                 + pre * Ap / (Bp**2 * A - B * Ap * Bp) * br2b * LXp)
    out['2c'] = pre * Bp / Ap / D * br2c * LXp - pre * B / A / D * br2c * LX
    out['2d'] = (-pre * lk / (A - B) / (Bp - Ap) * br2d + pre * A / (A - B) / D * br2d * LX
                 + pre * Ap / (Bp - Ap) / D * br2d * LXp)
    out['2e'] = (pre * lk / (A - B) / (Bp - Ap) * br2e - pre * B / (A - B) / D * br2e * LX
                 - pre * Bp / (Bp - Ap) / D * br2e * LXp)
    pf = G**4 * d / (32 * PI**5)
    out['2f'] = (pf * lk / kp * (Bp * B - A * Ap) / (zw2 * zwp2 * (A - B) * (Ap - Bp))
                 - pf / kp * A * B / (zw2 * zwp2 * (A - B)) * (Ap - Bp) / D * LX
                 + pf * (1 / V - 1 / (kp + L)) / (zw2 * zwp2)
                 + pf / kp * Bp * Ap / (zw2 * zwp2) * (A - B) / ((Ap - Bp) * D) * LXp)
    pc = phase(c, 'w', 'wp') / PI3 * col_G2_I2(c)
    return {k: v * pc for k, v in out.items()}


# ================================================================== GROUP III
def G3_I_II_kernel(c, kp, V, L):
    Xp, Yz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z')
    return (-1j * G**4 / (8 * PI**5) / kp * ij(Xp, Yz, br4(c)) * logIV(c, kp, V, L)
            + 1j * G**4 / (16 * PI**5) / kp * np.log(V / (kp + L)) * (Xp @ Yz) / sq(vec(c, 'z', 'w'))
            - 1j * G**4 / (8 * PI**5) / kp * np.log(1 - kp / V) * ij(Xp, Yz, brlog(c)))


def G3_I(c, kp, V, L):
    return phase(c, 'w', 'wp') / PI3 * G3_I_II_kernel(c, kp, V, L) * col_G3_I(c)


def G3_II(c, kp, V, L):
    return phase(c, 'w', 'wp') / PI3 * G3_I_II_kernel(c, kp, V, L) * col_G3_II(c)


def C_tens(c):
    """(y'-z)^m (x'-y')^{k'} (y-w)^i (x-z)^j / (...) contracted with d_{k'm}d_{ij}, d_{jm}d_{ik'}, d_{im}d_{jk'}"""
    a, b, ci, dj = Kv(c, 'yp', 'z'), Kv(c, 'xp', 'yp'), Kv(c, 'y', 'w'), Kv(c, 'x', 'z')
    return (b @ a) * (ci @ dj), (dj @ a) * (ci @ b), (ci @ a) * (dj @ b)


def C_phase(c, kp, pp):
    return np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * np.exp(-1j * (c.kv @ vec(c, 'yp', 'z')) * pp / kp)


def G3_C_p(c, kp, pp, pref, col):
    """integrand of int_0^{vee-k} dp/2pi of Rows III(1), V, VI as written (with [1/p]_+)"""
    t1, t2, t3 = C_tens(c)
    val = (t1 / (pp + kp) - t3 / kp) * C_phase(c, kp, pp) - t2 * (C_phase(c, kp, pp) - C_phase(c, kp, 0.0)) / pp
    return pref / PI3 / kp * val * col


def G3_C_rest(c, kp, V, L, preflog, col):
    t1, t2, t3 = C_tens(c)
    return -preflog / PI3 / kp * np.log(1 - kp / V) * t2 * C_phase(c, kp, 0.0) * col


def G3_III1_p(c, kp, pp):
    return G3_C_p(c, kp, pp, 1j * G**4 / (8 * PI**4), col_G3_III1(c))


def G3_III1_rest(c, kp, V, L):
    return G3_C_rest(c, kp, V, L, 1j * G**4 / (16 * PI**5), col_G3_III1(c))


def G3_V_p(c, kp, pp):
    return G3_C_p(c, kp, pp, 1j * G**4 / (4 * PI**4), col_G3_V(c))


def G3_V_rest(c, kp, V, L):
    return G3_C_rest(c, kp, V, L, 1j * G**4 / (8 * PI**5), col_G3_V(c))


def G3_VI_p(c, kp, pp):
    return G3_C_p(c, kp, pp, 1j * G**4 / (4 * PI**4), col_G3_VI(c))


def G3_VI_rest(c, kp, V, L):
    return G3_C_rest(c, kp, V, L, 1j * G**4 / (8 * PI**5), col_G3_VI(c))


def G3_III2_p(c, kp, pp, d=2):
    """Row III, second part: integrand of int_0^{vee-k} dp/2pi (as written: the coordinate x' does not appear in the
    measure; the bracket as written)"""
    x, z, w = c.pt['x'], c.pt['z'], c.pt['w']
    xz, xw, zw = x - z, x - w, z - w
    A, B = sq(xz), sq(xw)
    S = A - B
    P = pp + kp
    a, b = Kv(c, 'yp', 'z'), Kv(c, 'xp', 'yp')        # (y'-z)^m/(y'-z)^2, (x'-y')^{k'}/(x'-y')^2
    ab = a @ b
    reg = (ab * d * pp * S / (2 * sq(zw) * P**2)
           + pp / (kp * P) * ab * (xz @ zw / sq(zw) - xz @ xw / (2 * B) - S / (2 * sq(zw)))
           + 1 / P * ab * (xw @ zw / sq(zw) + xz @ xw / (2 * A) - S / (2 * sq(zw)))
           - pp / kp**2 * ((b @ xz) * (a @ zw) / sq(zw) - (b @ xz) * (a @ xw) / (2 * B))
           - 1 / kp * ((a @ xz) * (b @ zw) / sq(zw) - (a @ xz) * (b @ xw) / (2 * B) + (a @ xw) * (b @ zw) / sq(zw)
                       + (b @ xz) * (a @ xw) / (2 * A)))
    sing = -((b @ xw) * (a @ zw) / sq(zw) + (a @ xz) * (b @ xw) / (2 * A))      # times [1/p]_+
    Dn = lambda p: p * A + kp * B
    val = reg / Dn(pp) * C_phase(c, kp, pp) + (sing / Dn(pp) * C_phase(c, kp, pp) - sing / Dn(0.0) * C_phase(c, kp, 0.0)) / pp
    return G**4 / (4 * PI**4) / PI3 * val * col_G3_III2(c)


def G3_III2_rest(c, kp, V, L):
    x, z, w, yp, xp = (c.pt[n] for n in ('x', 'z', 'w', 'yp', 'xp'))
    xz, xw, zw = x - z, x - w, z - w
    val = ((xw @ (xp - yp)) * (zw @ (yp - z)) / (sq(zw) * sq(xw) * sq(yp - z) * sq(xp - yp))
           + (xz @ (yp - z)) * (xw @ (xp - yp)) / (2 * sq(xz) * sq(xw) * sq(yp - z) * sq(xp - yp)))
    return -G**4 / (8 * PI**5) / PI3 / kp * np.log(1 - kp / V) * val * C_phase(c, kp, 0.0) * col_G3_III2(c)
