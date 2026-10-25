"""Every row of Finite_terms_1.tex, written as   row = C_row * Lv + R_row(a),   Lv = log(vee/k+), a = k+/vee,
with C_row independent of a (the coefficient of the large log) and R_row(a) finite for a -> 0.

 1. C_row is read off from the closed forms (../finite/simp_exact.py, checked there against the file).
 2. Test: R_row(a) = row(a) - C_row Lv has a limit for a -> 0 (vee/k+ = 1e4, 1e6, 1e8).  A wrong C_row would make
    R_row grow like Lv.  Bounded oscillating terms (e^{-i kappa/a}, which vanish only after the y' integral) are
    removed explicitly in this test.
 3. c.c.: Group II rows are self-conjugate (row(c) = conj(row(c~)), c~ = the configuration with the amplitude and
    conjugate points exchanged), so no c.c. is added; Group I and Group III rows are not, so their c.c. is added.
Output: split.json"""
import json
import numpy as np
from scipy.special import exp1
import env  # noqa: F401
import fp
import ft
import simp_exact as SE
from ft import vec, Kv, ij, br4, AB, PI3
from fp import G, PI, sq

NAMES = ('x', 'y', 'xp', 'yp', 'z', 'w', 'wp')
EG = np.euler_gamma


def cfg(seed):
    return fp.Cfg(seed, sym=True, names=NAMES)


# ------------------------------------------------------------------ C_row (coefficient of Lv), same normalization as SE
def C_g1_iii(c, kp):
    A, B = AB(c)
    xz, xw, zw = vec(c, 'x', 'z'), vec(c, 'x', 'w'), vec(c, 'z', 'w')
    Xp, Yz, Zw, Xw, Xz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'z', 'w'), Kv(c, 'x', 'w'), Kv(c, 'x', 'z')
    p8 = 1j * G**4 / (8 * PI**5) / kp
    S4 = ((Xp @ Zw) * (Yz @ xw) - (Xp @ Xw) * (Yz @ Zw) * A - (Xp @ Zw) * (Yz @ Xz) * B + (Xp @ xz) * (Yz @ Zw)) / (B - A)
    k = p8 * ((Xp @ Xw) * (Yz @ Zw - Yz @ Xz) + (Xp @ Yz) / sq(zw) - S4)
    return k * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G1_III(c)


def Ck_iv(c, W='w', X='x', Xp_from='wp'):
    """Lv coefficient of the kernel K_IV: X'^i Y^j T^{ij} - X'.Y/(2(z-w)^2)"""
    Xp, Yz = Kv(c, 'xp', Xp_from), Kv(c, 'y', 'z')
    return ij(Xp, Yz, br4(c, W, X)) - 0.5 * (Xp @ Yz) / sq(vec(c, 'z', W))


def C_g1_iv(c, kp):
    return 1j * G**4 / (8 * PI**5) / kp * Ck_iv(c) * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G1_IV(c)


def C_g1_v1(c, kp):
    kk = ft.G1_V1_kernels(c)
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    ex = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x']))
    return 1j * G**4 / (32 * PI**5) / kp / PI3 * (ez * kk['b'] + ex * kk['q']) * ft.col_G1_V1(c)


def _v2_C(c):
    s, A, Yx, Yz = SE.v2_scalars(c)
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    ex = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x']))
    return ez * s['s3'] / Yz + ex * s['s6'] / Yx


def C_g1_v2a(c, kp):
    return G**4 / (8 * PI**4) / kp / PI3 / (2 * PI) * _v2_C(c) * SE.col_V2a(c)


def C_g1_v2b(c, kp):
    x, y, z = c.pt['x'], c.pt['y'], c.pt['z']
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    extra = ez * 2 * (Kv(c, 'xp', 'wp') @ Kv(c, 'y', 'z')) * ((y - x) @ (y - z)) / (sq(y - x) * sq(x - z))
    return G**4 / (8 * PI**4) / kp / PI3 / (2 * PI) * (_v2_C(c) * SE.col_V2b(c, 'x') + extra * SE.col_V2b(c, 'z'))


def C_g2_i1(c, kp):
    p16 = 1j * G**4 / (16 * PI**5) / kp
    ph = ft.phase(c, 'w', 'wp') / PI3
    return ph * (-p16 * Ck_iv(c, 'w', 'x', 'wp') * ft.col_G2_I1(c, 'w', 'wp')
                 + p16 * Ck_iv(c, 'wp', 'x', 'w') * ft.col_G2_I1(c, 'wp', 'w'))


def C_g2_i2(c, kp):
    import simp
    br, (A, B, Ap, Bp, zw2, zwp2) = simp.g2_i2_groups(c)
    pre = G**4 / (8 * PI**5) / kp
    return pre * (-br['b'] / (B * Bp) + br['c'] / (A * Ap)) * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G2_I2(c)


def C_g3_i(c, kp):
    return -1j * G**4 / (8 * PI**5) / kp * Ck_iv(c) * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G3_I(c)


def C_g3_ii(c, kp):
    return -1j * G**4 / (8 * PI**5) / kp * Ck_iv(c) * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G3_II(c)


def C_crow(pref, colf):
    def f(c, kp):
        t1, t2, t3 = ft.C_tens(c)
        return pref / (2 * PI) / PI3 / kp * np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * t2 * colf(c)
    return f


def C_g3_iii2(c, kp):
    x, z, w = c.pt['x'], c.pt['z'], c.pt['w']
    xz, xw, zw = x - z, x - w, z - w
    B = sq(xw)
    av_, bv = Kv(c, 'yp', 'z'), Kv(c, 'xp', 'yp')
    sing = -((bv @ xw) * (av_ @ zw) / sq(zw) + (av_ @ xz) * (bv @ xw) / (2 * sq(xz)))
    return G**4 / (4 * PI**4) / (2 * PI) / PI3 * np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * (-sing / (kp * B)) * ft.col_G3_III2(c)


# ------------------------------------------------------------------ bounded oscillating parts (vanish only after int_y')
def osc_crow(pref, colf):
    def f(c, kp, V):
        t1, t2, t3 = ft.C_tens(c)
        al = (c.kv @ vec(c, 'yp', 'z')) / kp
        J3osc = -np.exp(-1j * al * (V - kp)) / (1j * al)        # the e^{-i al (vee-k)} part of (1 - e^{...})/(i al)
        return pref / (2 * PI) / PI3 / kp * np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * (-t3 / kp * J3osc) * colf(c)
    return f


def osc_g3_iii2(c, kp, V, d=2):
    x, z, w = c.pt['x'], c.pt['z'], c.pt['w']
    xz, xw, zw = x - z, x - w, z - w
    A, B = sq(xz), sq(xw)
    av_, bv = Kv(c, 'yp', 'z'), Kv(c, 'xp', 'yp')
    c4 = -((bv @ xz) * (av_ @ zw) / sq(zw) - (bv @ xz) * (av_ @ xw) / (2 * B))
    al = (c.kv @ vec(c, 'yp', 'z')) / kp
    I0osc = -np.exp(-1j * al * (V - kp)) / (1j * al)
    return G**4 / (4 * PI**4) / (2 * PI) / PI3 * np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * (c4 / kp**2 * I0osc / A) * ft.col_G3_III2(c)


P_III1, P_VI = 1j * G**4 / (8 * PI**4), 1j * G**4 / (4 * PI**4)
ROWS = {
    'G1-III': (SE.g1_iii, C_g1_iii, None),
    'G1-IV': (SE.g1_iv, C_g1_iv, None),
    'G1-V(1)': (SE.g1_v1, C_g1_v1, None),
    'G1-V(2a)': (SE.g1_v2a, C_g1_v2a, None),
    'G1-V(2b fin)': (SE.g1_v2b_fin, C_g1_v2b, None),
    'G2-I(1)': (SE.g2_i1, C_g2_i1, None),
    'G2-I(2)': (lambda c, k, V: sum(SE.g2_i2(c, k, V).values()), C_g2_i2, None),
    'G3-I': (SE.g3_i, C_g3_i, None),
    'G3-II': (SE.g3_ii, C_g3_ii, None),
    'G3-III(1)': (SE.g3_iii1, C_crow(P_III1, ft.col_G3_III1), osc_crow(P_III1, ft.col_G3_III1)),
    'G3-III(2)': (SE.g3_iii2, C_g3_iii2, osc_g3_iii2),
    'G3-V': (SE.g3_v, C_crow(P_VI, ft.col_G3_V), osc_crow(P_VI, ft.col_G3_V)),
    'G3-VI': (SE.g3_vi, C_crow(P_VI, ft.col_G3_VI), osc_crow(P_VI, ft.col_G3_VI)),
}


def limit_test(seeds=(1, 2, 3), vees=(1e4, 1e6, 1e8)):
    out = {}
    for name, (row, C, osc) in ROWS.items():
        res = []
        for s in seeds:
            c = cfg(s)
            kp = 1.0
            Cv = C(c, kp)
            R = []
            for V in vees:
                r = row(c, kp, V) - Cv * np.log(V / kp)
                if osc is not None:
                    r -= osc(c, kp, V)
                R.append(r)
            spread = max(abs(R[i] - R[-1]) for i in range(len(R)))
            res.append(dict(spread=float(spread), R0=float(abs(R[-1])), C=float(abs(Cv)),
                            rel=float(spread / max(abs(R[-1]), abs(Cv)))))
        out[name] = res
        print('%-13s R_row(a) spread for a = 1e-4..1e-8, relative to |R|,|C|: %s' % (name, ' '.join('%.1e' % r['rel'] for r in res)))
    return out


# ------------------------------------------------------------------ c.c. test
def swapped(c, pairs):
    """the configuration with the points (and their U, rho) of each pair exchanged"""
    import copy
    d = copy.deepcopy(c)
    for p, q in pairs:
        for attr in ('pt', 'U', 'rho'):
            getattr(d, attr)[p], getattr(d, attr)[q] = getattr(c, attr)[q], getattr(c, attr)[p]
    return d


def cc_test(seeds=(1, 2, 3), V=7.0):
    """self-conjugacy: row(c) =? conj(row(c~)).  Group II: c~ exchanges the amplitude and conjugate points."""
    import simp_exact as S
    tests = {
        'G2-I(1)  [w<->w\']': (S.g2_i1, [('w', 'wp')]),
        'G2-I(2)  [w<->w\', x<->x\']': (lambda c, k, V: sum(S.g2_i2(c, k, V).values()), [('w', 'wp'), ('x', 'xp')]),
        'G1-IV    [w<->w\', x<->x\']': (S.g1_iv, [('w', 'wp'), ('x', 'xp')]),
        'G1-IV    [w<->w\']': (S.g1_iv, [('w', 'wp')]),
        'G3-I     [w<->w\', x<->x\']': (S.g3_i, [('w', 'wp'), ('x', 'xp')]),
        'G3-I     [w<->w\']': (S.g3_i, [('w', 'wp')]),
    }
    out = {}
    for name, (f, pairs) in tests.items():
        r = []
        for s in seeds:
            c = cfg(s)
            a_ = f(c, 1.0, V)
            b_ = np.conj(f(swapped(c, pairs), 1.0, V))
            r.append(float(abs(a_ - b_) / max(abs(a_), abs(b_))))
        out[name] = max(r)
        print('%-30s |row - conj(row~)| / |row| = %.1e' % (name, max(r)))
    return out


def g2_iv_cc(seeds=(1, 2, 3)):
    """Group II Row IV: the LO-type colour [U(y')-U(x')][U(y)-U(x)] rho(x') rho(x) and the kernel are symmetric under
    (x, y) <-> (x', y') with conj; check the colour and the kernel K(x-y).K(x'-y') e^{-ik(y'-y)} directly"""
    r = []
    for s in seeds:
        c = cfg(s)
        def val(c):
            U, rho = c.U, c.rho
            col = ((U['yp'] - U['xp']) @ rho['xp']) @ ((U['y'] - U['x']) @ rho['x'])
            return np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['y'])) * (Kv(c, 'x', 'y') @ Kv(c, 'xp', 'yp')) * col
        a_ = val(c)
        b_ = np.conj(val(swapped(c, [('y', 'yp'), ('x', 'xp')])))
        r.append(float(abs(a_ - b_) / abs(a_)))
    print('G2-IV  LO structure [y<->y\', x<->x\']: %.1e' % max(r))
    return max(r)


if __name__ == '__main__':
    out = dict(limit=limit_test(), cc=cc_test(), g2iv_cc=g2_iv_cc())
    json.dump(out, open('split.json', 'w'), indent=1)
