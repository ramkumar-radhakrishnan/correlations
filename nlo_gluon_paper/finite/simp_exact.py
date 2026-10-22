"""Finite terms with the exact vee dependence (a = k+/vee kept, no expansion in k+/vee); only Lambda -> 0.

Every function returns the integrand at the transverse points of a configuration c (same normalization, colour
and phase as ft.py).  The p+ (or t = p+/k+) integrals that remain in Finite_terms.tex are done in closed form with
the exponential integral E1 (complex argument) and Ein(z) = E1(z) + log z + gamma_E.

check() compares with the expressions of the file (ft.py, plus the transcription of Group I Row V part 2 below):
  closed-form rows at Lambda = 1e-10 (difference O(Lambda/k+)), p+ integrals by direct quadrature,
  for vee/k+ = 1.3, 2.5, 10, 300.  Output: simp_exact.json."""
import json
import numpy as np
from scipy.special import exp1
import env  # noqa: F401
from fp import G, PI, sq, F, NC
import ft
from ft import vec, Kv, ij, br4, brlog, AB, PI3
from verify import cquad

EG = np.euler_gamma


def Ein(z):
    return exp1(z) + np.log(z) + EG


def av(kp, V):
    return kp / V, np.log(V / kp), np.log(1 - kp / V)


# ================================================================== GROUP I
def g1_iii(c, kp, V):
    """Row III, all four parts N1+N2+N3a+N3b, Lambda = 0, exact in a"""
    a, Lv, l1a = av(kp, V)
    A, B = AB(c)
    xz, xw, zw = vec(c, 'x', 'z'), vec(c, 'x', 'w'), vec(c, 'z', 'w')
    Xp, Yz, Zw, Xw, Xz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'z', 'w'), Kv(c, 'x', 'w'), Kv(c, 'x', 'z')
    p16, p8 = 1j * G**4 / (16 * PI**5) / kp, 1j * G**4 / (8 * PI**5) / kp
    XY, zw2 = Xp @ Yz, sq(zw)
    T1 = XY * (A + B) / (A - B)**2
    Q = XY * ((zw @ xz) * B - (zw @ xw) * A) / (zw2 * (B - A)**2)
    S4 = ((Xp @ Zw) * (Yz @ xw) - (Xp @ Xw) * (Yz @ Zw) * A - (Xp @ Zw) * (Yz @ Xz) * B + (Xp @ xz) * (Yz @ Zw)) / (B - A)
    lam2 = np.log(np.abs(A / B - a))
    r = np.sqrt(A * B - (xw @ xz)**2)
    A0 = (np.arctan((A - xw @ xz) / r) + np.arctan((xw @ xz) / r)) / r
    M0 = (Xp @ Zw) * (Yz @ xw) + (Yz @ Zw) * (Xp @ xz) - (xw @ xz) * (Xp @ Zw) * (Yz @ xz) / A
    L5 = np.log(zw2 / B)
    N3b = (-p8 * A0 * (xw @ xz) * (Yz @ Zw) * (Xp @ xw) / B
           - p16 * L5 * (Xp @ Zw) * (Yz @ xz) / A + p16 * L5 * (Yz @ Zw) * (Xp @ xw) / B - p16 * L5 * XY / zw2)
    tot = (p8 * Lv * ((Xp @ Xw) * (Yz @ Zw - Yz @ Xz) + XY / zw2 - S4)
           + lam2 * (p16 * T1 - p8 * (Xp @ Xw) * (Xz @ Yz) + p8 * (Q - S4))
           - l1a * (p16 * T1 + p8 * Q)
           + p8 * A0 * M0 + N3b)
    return tot * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G1_III(c)


def kernel_iv(c, kp, V, W='w', X='x', Xp_from='wp'):
    """[X'^i Y^j T^{ij} (Lv + lam) - X'.Y Lv / (2(z-w)^2) + log(1-a) X'^i Y^j Bl^{ij}],  lam = log((1-a)A/B + a)"""
    a, Lv, l1a = av(kp, V)
    A, B = AB(c, X, W)
    Xp, Yz = Kv(c, 'xp', Xp_from), Kv(c, 'y', 'z')
    lam = np.log((1 - a) * A / B + a)
    return (ij(Xp, Yz, br4(c, W, X)) * (Lv + lam) - 0.5 * (Xp @ Yz) / sq(vec(c, 'z', W)) * Lv
            + l1a * ij(Xp, Yz, brlog(c, W, X)))


def g1_iv(c, kp, V):
    return 1j * G**4 / (8 * PI**5) / kp * kernel_iv(c, kp, V) * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G1_IV(c)


def g1_v1(c, kp, V):
    """Row V part 1: the int_0^k dp, dq done.  E(b) = Ein(i b) + Lv = gamma + log(i b vee/k) + E1(i b)."""
    a, Lv, l1a = av(kp, V)
    kk = ft.G1_V1_kernels(c)
    beta = c.kv @ vec(c, 'z', 'x')
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    ex = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x']))
    val = (ez * (kk['a'] * (1 - np.exp(-1j * beta)) / (1j * beta) + kk['b'] * (Ein(1j * beta) + Lv))
           + ex * kk['q'] * (Ein(-1j * beta) + Lv))
    return 1j * G**4 / (32 * PI**5) / kp / PI3 * val * ft.col_G1_V1(c)


# ---- Group I Row V part 2 (2a, 2b finite): transcription of the file and master-integral form
def v2_scalars(c, d=2):
    x, y, z = c.pt['x'], c.pt['y'], c.pt['z']
    Xp, Zx = Kv(c, 'xp', 'wp'), Kv(c, 'z', 'x')
    A, Yx, Yz = sq(x - z), sq(y - x), sq(y - z)
    XZ = Xp @ Zx
    s = dict(
        s1=XZ * d / (2 * A) * (Yx - Yz),
        v2=XZ * ((y - z) @ (x - z) / A + (y - x) @ (y - z) / (2 * Yx) - (Yx - Yz) / (2 * A)),
        s3=(Xp @ (y - z)) * (Zx @ (x - z)) / A + (Zx @ (y - x)) * (Xp @ (y - z)) / (2 * Yx),
        v4=XZ * ((y - x) @ (x - z) / A - (y - x) @ (y - z) / (2 * Yz) - (Yx - Yz) / (2 * A)),
        s5=((Zx @ (y - x)) * (Xp @ (x - z)) / A - (Zx @ (y - x)) * (Xp @ (y - z)) / (2 * Yz)
            + (Zx @ (y - z)) * (Xp @ (x - z)) / A + (Xp @ (y - x)) * (Zx @ (y - z)) / (2 * Yx)),
        s6=(Xp @ (y - x)) * (Zx @ (x - z)) / A - (Xp @ (y - x)) * (Zx @ (y - z)) / (2 * Yz))
    return s, A, Yx, Yz


def col_V2a(c):
    """[U(x')-U(w')]^{ab'} [-f^{adc} f^{fbe} U^{db}(x) U^{ce}(z) + N_c U^{af}(x)] rho^{b'}(x') rho^f(y)"""
    U, r = c.U, c.rho
    L = (U['xp'] - U['wp']) @ r['xp']
    t1 = -np.einsum('adc,fbe,db,ce,a,f->', F, F, U['x'], U['z'], L, r['y'])
    return t1 + NC * L @ U['x'] @ r['y']


def col_V2b(c, Z='x'):
    """[U(x')-U(w')]^{ab'} N_c [U^{af}(y) - U^{af}(Z)] rho^{b'}(x') rho^f(y)"""
    U, r = c.U, c.rho
    L = (U['xp'] - U['wp']) @ r['xp']
    return NC * L @ (U['y'] - U[Z]) @ r['y']


def v2_direct(c, kp, V):
    """the braces of 2a / 2b finite as written (p+ = k+ t, q+ = k+ t), Lambda = 0, by quadrature"""
    a, Lv, l1a = av(kp, V)
    s, A, Yx, Yz = v2_scalars(c)
    beta = c.kv @ vec(c, 'z', 'x')
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    ex = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x']))
    Dh = lambda t: t * Yx + (1 - t) * Yz
    Dt = lambda t: (1 - t) * Yx + t * Yz
    f = lambda t: (np.exp(-1j * beta * t) / Dh(t) * (s['s1'] * t * (1 - t) + s['v2'] * (1 - t) + s['s3'] + s['v4'] * t - s['s5'])
                   - s['s3'] / t * (np.exp(-1j * beta * t) / Dh(t) - 1 / Yz))
    g = lambda t: -s['s6'] * ((np.exp(1j * beta * t) / Dt(t) - 1 / Yx) / t - np.exp(1j * beta * t) / Dt(t))
    return (ez * cquad(f, 0.0, 1.0) + ex * cquad(g, 0.0, 1.0)
            + Lv * (ez * s['s3'] / Yz + ex * s['s6'] / Yx))


def masters(beta, Y0, Y1):
    """M_n = int_0^1 t^n e^{-i beta t} / (Y0 + t (Y1 - Y0)) dt,  n = 0, 1, 2 (closed form)"""
    Dl = Y1 - Y0
    cc = Y0 / Dl
    M0 = np.exp(1j * beta * cc) * (exp1(1j * beta * cc) - exp1(1j * beta * (1 + cc))) / Dl
    m0 = (1 - np.exp(-1j * beta)) / (1j * beta)
    m1 = -np.exp(-1j * beta) / (1j * beta) + m0 / (1j * beta)
    M1 = m0 / Dl - cc * M0
    M2 = (m1 - cc * m0) / Dl + cc**2 * M0
    return M0, M1, M2


def v2_master(c, kp, V):
    """braces of 2a / 2b finite in closed form:
       e_z [s1 (M1-M2) + v2 (M0-M1) + v4 M1 + (s3 Yx/Yz - s5) M0 + s3 E(beta)/Yz] + e_x s6 [(Yz/Yx) N0 + E(-beta)/Yx]
       M_n: denominator t Yx + (1-t) Yz, phase e^{-i beta t};  N_n: (1-t) Yx + t Yz, phase e^{+i beta t}"""
    a, Lv, l1a = av(kp, V)
    s, A, Yx, Yz = v2_scalars(c)
    beta = c.kv @ vec(c, 'z', 'x')
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    ex = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x']))
    M0, M1, M2 = masters(beta, Yz, Yx)
    N0, _, _ = masters(-beta, Yx, Yz)
    E = lambda b: Ein(1j * b) + Lv
    return (ez * (s['s1'] * (M1 - M2) + s['v2'] * (M0 - M1) + s['v4'] * M1 + (s['s3'] * Yx / Yz - s['s5']) * M0
                  + s['s3'] * E(beta) / Yz)
            + ex * s['s6'] * (Yz / Yx * N0 + E(-beta) / Yx))


def g1_v2a(c, kp, V, braces=v2_master):
    pre = G**4 / (8 * PI**4) / kp / PI3 / (2 * PI)
    return pre * braces(c, kp, V) * col_V2a(c)


def g1_v2b_fin(c, kp, V, braces=v2_master):
    a, Lv, l1a = av(kp, V)
    pre = G**4 / (8 * PI**4) / kp / PI3 / (2 * PI)
    x, y, z = c.pt['x'], c.pt['y'], c.pt['z']
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    extra = ez * (11 / 6 + 2 * Lv) * (Kv(c, 'xp', 'wp') @ Kv(c, 'y', 'z')) * ((y - x) @ (y - z)) / (sq(y - x) * sq(x - z))
    return pre * (braces(c, kp, V) * col_V2b(c, 'x') + extra * col_V2b(c, 'z'))


# ================================================================== GROUP II
def g2_i1(c, kp, V):
    p16 = 1j * G**4 / (16 * PI**5) / kp
    kA = -p16 * kernel_iv(c, kp, V, 'w', 'x', 'wp')
    kB = p16 * kernel_iv(c, kp, V, 'wp', 'x', 'w')
    ph = ft.phase(c, 'w', 'wp') / PI3
    return ph * (kA * ft.col_G2_I1(c, 'w', 'wp') + kB * ft.col_G2_I1(c, 'wp', 'w'))


def g2_i2(c, kp, V, d=2):
    import simp
    a, Lv, l1a = av(kp, V)
    br, (A, B, Ap, Bp, zw2, zwp2) = simp.g2_i2_groups(c)
    D = Bp * A - B * Ap
    lam, lamp = np.log((1 - a) * A / B + a), np.log((1 - a) * Ap / Bp + a)
    pre = G**4 / (8 * PI**5) / kp
    pf = G**4 * d / (32 * PI**5) / kp
    R = {'2a': pre * br['a'] * (lam - lamp) / D,
         '2b': pre * br['b'] * ((l1a - Lv) / (B * Bp) - A * lam / (B * D) + Ap * lamp / (Bp * D)),
         '2c': pre * br['c'] * (Lv / (A * Ap) + (Bp / Ap * lamp - B / A * lam) / D),
         '2d': pre * br['d'] * (A * lam / ((A - B) * D) + Ap * lamp / ((Bp - Ap) * D)),
         '2e': -pre * br['e'] * (B * lam / ((A - B) * D) + Bp * lamp / ((Bp - Ap) * D)),
         '2f': pf / (zw2 * zwp2) * (-A * B * (Ap - Bp) / ((A - B) * D) * lam + Ap * Bp * (A - B) / ((Ap - Bp) * D) * lamp
                                    - (1 - a))}
    pc = ft.phase(c, 'w', 'wp') / PI3 * ft.col_G2_I2(c)
    return {k: v * pc for k, v in R.items()}


# ================================================================== GROUP III
def g3_i(c, kp, V):
    return -1j * G**4 / (8 * PI**5) / kp * kernel_iv(c, kp, V) * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G3_I(c)


def g3_ii(c, kp, V):
    return -1j * G**4 / (8 * PI**5) / kp * kernel_iv(c, kp, V) * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G3_II(c)


def crow(c, kp, V, pref, col):
    """Rows III(1), V, VI:  pref/(2 pi (2pi)^3 k+) e^{-ik.(y'-w)} col
         { t1 e^{i al k} [E1(i al k) - E1(i al vee)] - (t3/k) (1 - e^{-i al (vee-k)})/(i al)
           + t2 [gamma + log(i al vee) + E1(i al (vee-k))] },   al = k.(y'-z)/k+"""
    t1, t2, t3 = ft.C_tens(c)
    al = (c.kv @ vec(c, 'yp', 'z')) / kp
    Lr = V - kp
    J1 = np.exp(1j * al * kp) * (exp1(1j * al * kp) - exp1(1j * al * V))
    J3 = (1 - np.exp(-1j * al * Lr)) / (1j * al)
    K2 = EG + np.log(1j * al * V) + exp1(1j * al * Lr)
    val = t1 * J1 - t3 / kp * J3 + t2 * K2
    return pref / (2 * PI) / PI3 / kp * np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * val * col


def g3_iii1(c, kp, V):
    return crow(c, kp, V, 1j * G**4 / (8 * PI**4), ft.col_G3_III1(c))


def g3_v(c, kp, V):
    return crow(c, kp, V, 1j * G**4 / (4 * PI**4), ft.col_G3_V(c))


def g3_vi(c, kp, V):
    return crow(c, kp, V, 1j * G**4 / (4 * PI**4), ft.col_G3_VI(c))


def Ic(al, cc, Lr):
    """int_0^Lr e^{-i al p}/(p + cc) dp,  cc > 0"""
    return np.exp(1j * al * cc) * (exp1(1j * al * cc) - exp1(1j * al * (Lr + cc)))


def g3_iii2(c, kp, V, d=2):
    """Row III part 2: the p+ integral in closed form (partial fractions in p+, c_* = k+ B/A)"""
    x, z, w = c.pt['x'], c.pt['z'], c.pt['w']
    xz, xw, zw = x - z, x - w, z - w
    A, B = sq(xz), sq(xw)
    S = A - B
    av_ = Kv(c, 'yp', 'z')
    bv = Kv(c, 'xp', 'yp')
    ab = av_ @ bv
    c1 = ab * d * S / (2 * sq(zw))
    c2 = ab * (xz @ zw / sq(zw) - xz @ xw / (2 * B) - S / (2 * sq(zw)))
    c3 = ab * (xw @ zw / sq(zw) + xz @ xw / (2 * A) - S / (2 * sq(zw)))
    c4 = -((bv @ xz) * (av_ @ zw) / sq(zw) - (bv @ xz) * (av_ @ xw) / (2 * B))
    c5 = -((av_ @ xz) * (bv @ zw) / sq(zw) - (av_ @ xz) * (bv @ xw) / (2 * B) + (av_ @ xw) * (bv @ zw) / sq(zw)
           + (bv @ xz) * (av_ @ xw) / (2 * A))
    sing = -((bv @ xw) * (av_ @ zw) / sq(zw) + (av_ @ xz) * (bv @ xw) / (2 * A))
    al = (c.kv @ vec(c, 'yp', 'z')) / kp
    Lr = V - kp
    cs = kp * B / A
    Ik, Is = Ic(al, kp, Lr), Ic(al, cs, Lr)
    I0 = (1 - np.exp(-1j * al * Lr)) / (1j * al)
    I2k = 1 / kp - np.exp(-1j * al * Lr) / V - 1j * al * Ik
    # p/((p+k)^2 (p+c)) = a1/(p+c) + b1/(p+k) + g1/(p+k)^2
    a1, b1, g1 = -cs / (kp - cs)**2, cs / (kp - cs)**2, kp / (kp - cs)
    Ta = a1 * Is + b1 * Ik + g1 * I2k
    Tb = kp / (kp - cs) * Ik - cs / (kp - cs) * Is
    Tc = (Ik - Is) / (cs - kp)
    Td = I0 - cs * Is
    reg = (c1 * Ta + c2 / kp * Tb + c3 * Tc + c4 / kp**2 * Td + c5 / kp * Is) / A
    plus = sing / (kp * B) * (-(EG + np.log(1j * al * V) + exp1(1j * al * Lr)) - Is)
    val = reg + plus
    return G**4 / (4 * PI**4) / (2 * PI) / PI3 * np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * val * ft.col_G3_III2(c)


# ---- the same two closed forms written exactly as in the note (kappa = k.(y'-z), a = k+/vee, b = B/A)
def crow_note(c, kp, V, cpref, col):
    a = kp / V
    t1, t2, t3 = ft.C_tens(c)
    ka = c.kv @ vec(c, 'yp', 'z')
    br = (t1 * np.exp(1j * ka) * (exp1(1j * ka) - exp1(1j * ka / a)) - t3 * (1 - np.exp(-1j * ka * (1 - a) / a)) / (1j * ka)
          + t2 * (EG + np.log(1j * ka / a) + exp1(1j * ka * (1 - a) / a)))
    return cpref / kp / PI3 * np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * br * col


def g3_iii2_note(c, kp, V, d=2):
    a = kp / V
    x, z, w = c.pt['x'], c.pt['z'], c.pt['w']
    xz, xw, zw = x - z, x - w, z - w
    A, B = sq(xz), sq(xw)
    S = A - B
    av_, bv = Kv(c, 'yp', 'z'), Kv(c, 'xp', 'yp')
    ab = av_ @ bv
    c1 = ab * d * S / (2 * sq(zw))
    c2 = ab * (xz @ zw / sq(zw) - xz @ xw / (2 * B) - S / (2 * sq(zw)))
    c3 = ab * (xw @ zw / sq(zw) + xz @ xw / (2 * A) - S / (2 * sq(zw)))
    c4 = -((bv @ xz) * (av_ @ zw) / sq(zw) - (bv @ xz) * (av_ @ xw) / (2 * B))
    c5 = -((av_ @ xz) * (bv @ zw) / sq(zw) - (av_ @ xz) * (bv @ xw) / (2 * B) + (av_ @ xw) * (bv @ zw) / sq(zw)
           + (bv @ xz) * (av_ @ xw) / (2 * A))
    sing = -((bv @ xw) * (av_ @ zw) / sq(zw) + (av_ @ xz) * (bv @ xw) / (2 * A))
    ka, b, U = c.kv @ vec(c, 'yp', 'z'), B / A, (1 - a) / a
    J = lambda bb: np.exp(1j * ka * bb) * (exp1(1j * ka * bb) - exp1(1j * ka * (bb + U)))
    J1, Jb, J0 = J(1.0), J(b), (1 - np.exp(-1j * ka * U)) / (1j * ka)
    Gs = EG + np.log(1j * ka / a) + exp1(1j * ka * U)
    I1 = (b * (J1 - Jb) / (1 - b)**2 + (1 - a * np.exp(-1j * ka * U) - 1j * ka * J1) / (1 - b)) / (A * kp)
    I2 = (J1 - b * Jb) / (A * (1 - b))
    I3 = (J1 - Jb) / (kp * (b - 1) * A)
    I4 = kp * (J0 - b * Jb) / A
    I5 = Jb / A
    val = c1 * I1 + c2 / kp * I2 + c3 * I3 + c4 / kp**2 * I4 + c5 / kp * I5 + sing / (kp * B) * (-Gs - Jb)
    return G**4 / (8 * PI**5) / PI3 * np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * val * ft.col_G3_III2(c)


# ================================================================== checks
def fileC(pfun, rfun):
    def f(c, kp, V):
        return cquad(lambda p: pfun(c, kp, p), 0.0, V - kp) / (2 * PI) + rfun(c, kp, V, 0.0)
    return f


def file_V1(c, kp, V):
    a = cquad(lambda p: ft.G1_V1_p(c, kp, p, 'p'), 0.0, kp) / (2 * PI)
    b = cquad(lambda q: ft.G1_V1_p(c, kp, q, 'q'), 0.0, kp) / (2 * PI)
    return a + b + ft.G1_V1_rest(c, kp, V, 0.0)


def prep_C(c):
    """C rows: w' is fixed by the delta function; use the configuration's y' as is (as in ft.py)"""
    return c


LAM = 1e-10
ROWS = {
    'G1-III': (g1_iii, lambda c, k, V: sum(ft.G1_III(c, k, V, LAM).values())),
    'G1-IV': (g1_iv, lambda c, k, V: ft.G1_IV(c, k, V, LAM)),
    'G1-V(1)': (g1_v1, file_V1),
    'G1-V(2a)': (g1_v2a, lambda c, k, V: g1_v2a(c, k, V, v2_direct)),
    'G1-V(2b fin)': (g1_v2b_fin, lambda c, k, V: g1_v2b_fin(c, k, V, v2_direct)),
    'G2-I(1)': (g2_i1, lambda c, k, V: sum(ft.G2_I1(c, k, V, LAM))),
    'G2-I(2)': (lambda c, k, V: sum(g2_i2(c, k, V).values()), lambda c, k, V: sum(ft.G2_I2(c, k, V, LAM).values())),
    'G3-I': (g3_i, lambda c, k, V: ft.G3_I(c, k, V, LAM)),
    'G3-II': (g3_ii, lambda c, k, V: ft.G3_II(c, k, V, LAM)),
    'G3-III(1)': (g3_iii1, fileC(ft.G3_III1_p, ft.G3_III1_rest)),
    'G3-III(2)': (g3_iii2, fileC(ft.G3_III2_p, ft.G3_III2_rest)),
    'G3-V': (g3_v, fileC(ft.G3_V_p, ft.G3_V_rest)),
    'G3-VI': (g3_vi, fileC(ft.G3_VI_p, ft.G3_VI_rest)),
}
VEES = (1.3, 2.5, 10.0, 300.0)


def check(seeds=(1, 2, 3)):
    import fp
    out = {}
    for name, (s, f) in ROWS.items():
        worst = []
        for seed in seeds:
            c = fp.Cfg(seed, sym=True, names=('x', 'y', 'xp', 'yp', 'z', 'w', 'wp'))
            errs = []
            for V in VEES:
                a_, b_ = s(c, 1.0, V), f(c, 1.0, V)
                errs.append(float(abs(a_ - b_) / abs(b_)))
            worst.append(errs)
        out[name] = worst
        print('%-13s rel. diff. at vee/k+ = %s:' % (name, VEES), ['/'.join('%.0e' % e for e in w) for w in worst])
    # the 2b, 2c, 2d, 2e, 2f pieces separately
    import fp
    sep = []
    for seed in seeds:
        c = fp.Cfg(seed, sym=True, names=('x', 'y', 'xp', 'yp', 'z', 'w', 'wp'))
        for V in VEES:
            mine, fil = g2_i2(c, 1.0, V), ft.G2_I2(c, 1.0, V, LAM)
            sep.append(max(float(abs(mine[k] - fil[k]) / abs(fil[k])) for k in mine))
    out['G2-I(2) piecewise'] = sep
    print('G2-I(2) piece by piece, max rel. diff.: %.1e' % max(sep))
    # master integrals against quadrature
    rng = np.random.default_rng(7)
    mm = []
    for _ in range(20):
        beta, Y0, Y1 = rng.normal() * 3, rng.uniform(0.1, 3), rng.uniform(0.1, 3)
        M = masters(beta, Y0, Y1)
        for n in range(3):
            q = cquad(lambda t: t**n * np.exp(-1j * beta * t) / (Y0 + t * (Y1 - Y0)), 0.0, 1.0)
            mm.append(float(abs(M[n] - q) / abs(q)))
    out['masters'] = max(mm)
    print('master integrals M0, M1, M2 vs quadrature: %.1e' % max(mm))
    # the note's dimensionless forms against the code above
    nt = []
    for seed in seeds:
        c = fp.Cfg(seed, sym=True, names=('x', 'y', 'xp', 'yp', 'z', 'w', 'wp'))
        for V in VEES:
            for mine, note in ((g3_iii1(c, 1.0, V), crow_note(c, 1.0, V, 1j * G**4 / (16 * PI**5), ft.col_G3_III1(c))),
                               (g3_v(c, 1.0, V), crow_note(c, 1.0, V, 1j * G**4 / (8 * PI**5), ft.col_G3_V(c))),
                               (g3_iii2(c, 1.0, V), g3_iii2_note(c, 1.0, V))):
                nt.append(float(abs(mine - note) / abs(mine)))
    out['note_forms'] = max(nt)
    print('note forms (C rows, Row III part 2) vs code: %.1e' % max(nt))
    # the U(z) network of Group I Row III vs Group II Row I (A) + Group III Row II at finite vee (does not cancel)
    n2 = {}
    for V in (1.5, 3.0, 30.0):
        r = []
        for seed in (1, 2, 3, 4):
            c = fp.Cfg(seed, sym=True, names=('x', 'y', 'xp', 'yp', 'z', 'w', 'wp'))
            KG3 = -1j * G**4 / (8 * PI**5) * kernel_iv(c, 1.0, V)
            KIII = g1_iii(c, 1.0, V) / (ft.phase(c, 'w', 'wp') / PI3 * ft.col_G1_III(c))
            r.append(float(abs(KG3 + KIII) / abs(KG3)))
        n2[str(V)] = r
    out['N2_network_rel'] = n2
    print('U(z) network, |G1-III + (G2-I(1)A + G3-II)| / |G3 kernel| at vee/k+ = 1.5, 3, 30:', n2)
    return out


if __name__ == '__main__':
    json.dump(check(), open('simp_exact.json', 'w'), indent=1)
