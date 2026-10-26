"""The terms of bfkl.tex (your BFKL list), transcribed term by term into the reordering engine of ../two_rho/core.py.

Units: N = (1/(2pi)^3) g^4/(16 pi^5) (1/k+) log(vee/k+).  So i g^4/(8 pi^5) -> 2i, i g^4/(16 pi^5) -> i,
i g^4/(32 pi^5) -> i/2, g^4/(8 pi^5) -> 2, i g^4/(32 pi^4) -> i pi/2, g^4 N_c/(8 pi^4) -> 2 pi N_c.
Kernels: ('KK', a, b, c, d) = K(a-b).K(c-d), K(X) = X/X^2;  ('ELL', a, b) = int_u K(a-u).K(b-u) in dim. reg.
(MS-bar, d = 2 - eps, Part2b_MSbar.pdf):  ELL(a,b) = -pi [2/eps + 2 gamma_E + log((a-b)^2 mu^2/4)].
Every row is mapped to the common external points: w (measured gluon, amplitude), w' (conjugate), z (soft gluon).

Two changes, flagged:
  G1 III: the user's N1 + N2 + N3a is written in the simplified form 2i[-K_b + K_h + K_f] O_III (identical, checked
          in three_rho_test.py to 1e-15); with fix F3 (as in Evolution_terms.tex) K_f O_III -> D_k O_III(z->w).
  G1 II 2a: X'.F(w;x,y) = (1/pi) int_z F_z (an identity, checked in checks.py)."""
import os
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
TWO = os.path.normpath(os.path.join(HERE, '..', 'two_rho'))
if TWO not in sys.path:
    sys.path.insert(0, TWO)
from core import term, rename, conj_swap  # noqa: E402

NC = 3
I = 1j
XP, WP, YP, B1, C1, D1, E1 = "x'", "w'", "y'", "b'", "c'", "d'", "e'"


def KK(a, b, c, d):
    return ('KK', a, b, c, d)


# ---------------------------------------------------------------- colour algebra: sums of (coef, Us, fs, rs)
def P(*brackets):
    out = [(1.0, [], [], [])]
    for br in brackets:
        out = [(c1 * c2, u1 + u2, f1 + f2, r1 + r2) for c1, u1, f1, r1 in out for c2, u2, f2, r2 in br]
    return out


def lin(*pairs):
    """lin((c1, bracket1), (c2, bracket2), ...)"""
    return [(c * cc, u, f, r) for c, br in pairs for cc, u, f, r in br]


def mono(c, Us, fs, rs):
    return [(c, list(Us), list(fs), list(rs))]


def anti(c, Us, fs, r1, r2):
    """c Us fs {r1, r2}"""
    return [(c, list(Us), list(fs), [r1, r2]), (c, list(Us), list(fs), [r2, r1])]


CL = lin((1, mono(1, [('a', B1, XP)], [], [(B1, XP)])), (-1, mono(1, [('a', B1, WP)], [], [(B1, XP)])))


def build(coef, kern, colour, tag, mp=None, cc=False):
    """coef x sum_k kc_k kern_k x colour; mp: renaming of points; cc: add the complex conjugate"""
    out = []
    for kc, kf in kern:
        for c, Us, fs, rs in colour:
            t = term(coef * kc * c, Us, fs, rs, kf, tag)
            if mp:
                t = rename(t, mp)
            out.append(t)
            if cc:
                out.append(conj_swap(t))
    return out


# ---------------------------------------------------------------- kernels
K_IV = [(1, [KK(XP, WP, 'z', 'w'), KK('y', 'z', 'x', 'z')]),
        (-1, [KK(XP, WP, 'x', 'w'), KK('y', 'z', 'x', 'z')]),
        (-1, [KK(XP, WP, 'x', 'w'), KK('y', 'z', 'z', 'w')])]
K_IVp = [(1, [KK(XP, 'w', 'z', WP), KK('y', 'z', 'x', 'z')]),       # block B of G2 I(1): w <-> w' in the gluon slots
         (-1, [KK(XP, 'w', 'x', WP), KK('y', 'z', 'x', 'z')]),
         (-1, [KK(XP, 'w', 'x', WP), KK('y', 'z', 'z', WP)])]
F_Z = [(1, [KK(XP, WP, 'x', 'w'), KK('y', 'z', 'x', 'z')]), (-1, [KK(XP, WP, 'x', 'w'), KK('y', 'z', 'w', 'z')]),
       (-1, [KK(XP, WP, 'y', 'w'), KK('x', 'z', 'y', 'z')]), (1, [KK(XP, WP, 'y', 'w'), KK('x', 'z', 'w', 'z')])]
K_b = [KK(XP, WP, 'x', 'w'), KK('x', 'z', 'y', 'z')]
K_h = [KK(XP, WP, 'z', 'w'), KK('y', 'z', 'x', 'z')]
K_f = [KK(XP, WP, 'y', 'z'), KK('z', 'w', 'z', 'w')]
D_k = [KK(XP, WP, 'y', 'w'), KK('x', 'z', 'w', 'z')]
T2 = [(1, [KK(YP, 'z', 'x', 'z'), KK(XP, YP, 'y', 'w')])]


# ---------------------------------------------------------------- Group I (with c.c.)
def g1_ii_1():
    col = P(CL, lin((1, mono(1, [('a', 'b', 'w')], [], [('b', 'x')])), (-1, mono(1, [('a', 'b', 'x')], [], [('b', 'x')]))))
    return build(NC, [(1, [KK(XP, WP, 'x', 'w'), ('ELL', 'x', 'w')])], col, 'G1 II(1)', cc=True)


def g1_ii_2a():
    col = lin((-1, P(CL, anti(1, [('a', 'h', 'w')], [('h', 'b', 'c')], ('b', 'x'), ('c', 'y')))),
              (1, P(CL, anti(1, [('b', 'd', 'x'), ('c', 'e', 'y')], [('a', 'b', 'c')], ('d', 'x'), ('e', 'y')))))
    return build(I / 2, F_Z, col, 'G1 II(2a)', cc=True)


def O_III(p='z'):
    return P(CL, lin((1, mono(1, [('d', 'e', 'x'), ('b', 'c', p)], [('a', 'b', 'd')], [('e', 'x'), ('c', 'y')])),
                     (-1, mono(1, [('d', 'e', 'x'), ('b', 'c', 'y')], [('a', 'b', 'd')], [('e', 'x'), ('c', 'y')]))))


def g1_iii(F3=True):
    out = build(2 * I, [(-1, K_b), (1, K_h)], O_III('z'), 'G1 III', cc=True)
    if F3:
        out += build(2 * I, [(1, D_k)], O_III('w'), 'G1 III (D)', cc=True)
    else:
        out += build(2 * I, [(1, K_f)], O_III('z'), 'G1 III (K_f)', cc=True)
    return out


def g1_iv():
    col = P(CL, lin((1, mono(1, [('a', 'c', 'w'), ('b', 'd', 'z'), ('b', 'e', 'y')], [('c', 'd', C1)], [('e', 'y'), (C1, 'x')])),
                    (-1, mono(1, [('b', 'e', 'y'), ('c', 'd', 'x')], [('a', 'b', 'c')], [('e', 'y'), ('d', 'x')]))))
    return build(2 * I, K_IV, col, 'G1 IV', cc=True)


def g1_v1():
    col = P(CL, lin((1, anti(1, [('d', 'b', 'x'), ('c', 'e', 'z')], [('a', 'd', 'c')], ('e', 'y'), ('b', YP))),
                    (-1, anti(1, [('c', 'd', 'y'), ('b', 'e', YP)], [('a', 'b', 'c')], ('d', 'y'), ('e', YP)))))
    tb = build(I / 2, [(1, [KK(XP, WP, 'y', 'z'), KK('z', 'x', YP, 'x')])], col, 'G1 V(1) b', mp={'z': 'w', 'x': 'z'}, cc=True)
    tq = build(I / 2, [(1, [KK(XP, WP, YP, 'x'), KK('z', 'x', 'y', 'z')])], col, 'G1 V(1) q', mp={'x': 'w'}, cc=True)
    return tb + tq


def _v_col(fpart=True):
    L = lin((1, mono(1, [('a', B1, XP)], [], [])), (-1, mono(1, [('a', B1, WP)], [], [])))
    R = lin((-1, mono(1, [('d', 'b', 'x'), ('c', 'e', 'z')], [('a', 'd', 'c'), ('f', 'b', 'e')], [])),
            (NC, mono(1, [('a', 'f', 'x')], [], [])))
    return P(L, R, mono(1, [], [], [(B1, XP), ('f', 'y')]))


def g1_v2a():
    col = _v_col()
    t1 = build(1, [(1, [KK(XP, WP, 'y', 'z'), KK('z', 'x', 'x', 'z')]), (0.5, [KK(XP, WP, 'y', 'z'), KK('z', 'x', 'y', 'x')])],
               col, 'G1 V(2a) z', mp={'z': 'w', 'x': 'z'}, cc=True)
    t2 = build(1, [(1, [KK(XP, WP, 'y', 'x'), KK('z', 'x', 'x', 'z')]), (-0.5, [KK(XP, WP, 'y', 'x'), KK('z', 'x', 'y', 'z')])],
               col, 'G1 V(2a) x', mp={'x': 'w'}, cc=True)
    return t1 + t2


def _lo_like(p, q):
    """[U^{ab'}(x') - U^{ab'}(w')][U^{af}(p) - U^{af}(q)] rho^{b'}(x') rho^f(y)"""
    L = lin((1, mono(1, [('a', B1, XP)], [], [])), (-1, mono(1, [('a', B1, WP)], [], [])))
    R = lin((1, mono(1, [('a', 'f', p)], [], [])), (-1, mono(1, [('a', 'f', q)], [], [])))
    return P(L, R, mono(1, [], [], [(B1, XP), ('f', 'y')]))


def g1_v2b_uv():
    return build(2 * NC, [(1, [KK(XP, WP, 'y', 'z'), ('ELL', 'y', 'z')])], _lo_like('y', 'z'), 'G1 V(2b UV)',
                 mp={'z': 'w'}, cc=True)


def g1_v2b_fin():
    tx = build(NC, [(1, [KK('y', 'x', XP, WP), KK('x', 'z', 'z', 'x')]), (-0.5, [KK('y', 'x', XP, WP), KK('y', 'z', 'z', 'x')])],
               _lo_like('y', 'x'), 'G1 V(2b fin) x', mp={'x': 'w'}, cc=True)
    tz = build(NC, [(1, [KK('y', 'z', XP, WP), KK('x', 'z', 'z', 'x')]), (0.5, [KK('y', 'x', 'z', 'x'), KK('y', 'z', XP, WP)])],
               _lo_like('y', 'x'), 'G1 V(2b fin) z', mp={'z': 'w', 'x': 'z'}, cc=True)
    # third term: (y-x).(y-z)/((y-x)^2 (x-z)^2) = K(z,z,x) - K(y,z,x)
    t3 = build(2 * NC, [(1, [KK(XP, WP, 'y', 'z'), KK('z', 'x', 'z', 'x')]), (-1, [KK(XP, WP, 'y', 'z'), KK('y', 'x', 'z', 'x')])],
               _lo_like('y', 'z'), 'G1 V(2b fin) 3', mp={'z': 'w', 'x': 'z'}, cc=True)
    return tx + tz + t3


# ---------------------------------------------------------------- Group II (no c.c.)
def _RA(W='w'):
    return lin((1, mono(1, [('a', 'c', W)], [('c', 'b', 'd')], [('d', 'x')])),
               (-1, mono(1, [('b', 'c', 'z'), ('d', 'e', 'x')], [('a', 'c', 'd')], [('e', 'x')])))


def g2_i1():
    LA = lin((1, anti(1, [('a', C1, WP)], [], (C1, XP), ('b', 'y'))),
             (-1, anti(1, [('b', C1, 'z'), ('a', D1, XP), (C1, E1, 'y')], [], (D1, XP), (E1, 'y'))))
    A = build(-I, K_IV, P(LA, _RA('w')), 'G2 I(1) A')
    RB = lin((1, mono(1, [('a', C1, WP)], [(C1, 'b', D1)], [(D1, 'x')])),
             (-1, mono(1, [('b', C1, 'z'), (D1, E1, 'x')], [('a', C1, D1)], [(E1, 'x')])))
    LB = lin((1, anti(1, [('a', 'c', 'w')], [], ('c', XP), ('b', 'y'))),
             (-1, anti(1, [('b', 'c', 'z'), ('a', 'd', XP), ('c', 'e', 'y')], [], ('d', XP), ('e', 'y'))))
    B = build(I, K_IVp, P(RB, LB), 'G2 I(1) B')
    return A + B


def O58():
    return lin((1, mono(1, [('a', C1, WP), ('a', 'c', 'w')], [(C1, 'b', D1), ('c', 'b', 'd')], [(D1, XP), ('d', 'x')])),
               (-1, mono(1, [('a', C1, WP), ('b', 'c', 'z'), ('d', 'e', 'x')], [(C1, 'b', D1), ('a', 'c', 'd')], [(D1, XP), ('e', 'x')])),
               (-1, mono(1, [('a', 'c', 'w'), ('b', C1, 'z'), (E1, D1, XP)], [('a', C1, D1), ('c', 'b', 'd')], [(E1, XP), ('d', 'x')])),
               (NC, mono(1, [(E1, 'd', XP), ('d', 'e', 'x')], [], [(E1, XP), ('e', 'x')])))


def g2_i2():
    k2b = [(1, [KK(XP, WP, 'x', 'w'), KK('z', 'w', 'z', WP)]), (0.5, [KK(XP, WP, 'x', 'w'), KK(XP, 'z', 'z', 'w')]),
           (0.5, [KK('x', 'w', XP, WP), KK('x', 'z', 'z', WP)]), (0.25, [KK('x', 'z', XP, 'z'), KK(XP, WP, 'x', 'w')])]
    k2c = [(1, [KK(XP, 'z', 'x', 'z'), KK('z', 'w', 'z', WP)]), (-0.5, [KK(XP, 'z', 'x', 'z'), KK(XP, WP, 'z', 'w')]),
           (-0.5, [KK('x', 'w', 'z', WP), KK('x', 'z', XP, 'z')]), (0.25, [KK('x', 'z', XP, 'z'), KK(XP, WP, 'x', 'w')])]
    return build(-2, k2b, O58(), 'G2 I(2b)') + build(2, k2c, O58(), 'G2 I(2c)')


def _lo4(p, pp):
    """[U^{e'c}(pp) - U^{e'c}(x')][U^{ec}(p) - U^{ce}(x)] rho^{e'}(x') rho^e(x)"""
    L = lin((1, mono(1, [(E1, 'c', pp)], [], [])), (-1, mono(1, [(E1, 'c', XP)], [], [])))
    R = lin((1, mono(1, [('e', 'c', p)], [], [])), (-1, mono(1, [('c', 'e', 'x')], [], [])))
    return P(L, R, mono(1, [], [], [(E1, XP), ('e', 'x')]))


def g2_iv():
    dr = build(-2 * NC, [(1, [KK('x', 'y', XP, YP), ('ELL', 'y', YP)])], _lo4('y', YP), 'G2 IV (eps)', mp={'y': 'w', YP: WP})
    h4 = build(2 * NC, [(1, [KK('w', 'z', WP, 'z'), KK('x', 'z', XP, 'z')])], _lo4('z', 'z'), 'G2 IV (H4)')
    return dr + h4


# ---------------------------------------------------------------- Group III (with c.c.)
def _L3(kind):
    if kind == 'I':
        return lin((1, mono(1, [('a', B1, XP), (C1, D1, 'y'), (C1, 'b', 'z')], [], [(B1, XP), (D1, 'y')])),
                   (-1, mono(1, [('b', D1, 'z'), (D1, E1, 'y'), ('a', C1, WP)], [], [(C1, XP), (E1, 'y')])))
    return lin((1, mono(1, [(C1, D1, 'y'), ('a', B1, XP), ('b', C1, 'z')], [], [(D1, 'y'), (B1, XP)])),
               (-1, mono(1, [('a', B1, XP)], [], [('b', 'y'), (B1, XP)])))


def g3_i():
    return build(-2 * I, K_IV, P(_L3('I'), _RA('w')), 'G3 I', cc=True)


def g3_ii():
    return build(-2 * I, K_IV, P(_L3('II'), _RA('w')), 'G3 II', cc=True)


def _Lc():
    return lin((1, mono(1, [(E1, C1, YP), (D1, 'b', 'z')], [], [(E1, XP)])),
               (-1, mono(1, [('b', D1, 'z'), (C1, E1, XP)], [], [(E1, XP)])))


def g3_iii1():
    R = lin((1, anti(1, [('a', 'c', 'w')], [], ('c', 'y'), ('b', 'x'))),
            (-1, anti(1, [('b', 'c', 'z'), ('a', 'd', 'y'), ('c', 'e', 'x')], [], ('d', 'y'), ('e', 'x'))))
    col = P(mono(1, [], [(C1, D1, 'a')], []), _Lc(), R)
    return build(I, T2, col, 'G3 III(1)', mp={YP: WP}, cc=True)


def g3_iii2():
    k = [(1, [KK('x', 'w', XP, YP), KK('z', 'w', YP, 'z')]), (0.5, [KK('x', 'z', YP, 'z'), KK('x', 'w', XP, YP)])]
    col = P(mono(1, [], [(C1, D1, 'a')], []), _Lc(), _RA('w'))
    return build(2, k, col, 'G3 III(2)', mp={YP: WP}, cc=True)


def g3_v():
    L = lin((1, mono(1, [(E1, C1, YP)], [], [(E1, XP)])), (-1, mono(1, [(C1, E1, XP)], [], [(E1, XP)])))
    R = lin((1, mono(1, [('c', 'e', 'x'), ('a', 'd', 'y')], [], [('e', 'x'), ('d', 'y')])),
            (-1, mono(1, [('a', 'd', 'w'), ('c', 'e', 'x')], [], [('e', 'x'), ('d', 'y')])))
    col = P(mono(1, [], [(C1, 'c', 'a')], []), L, R)
    return build(2 * I, T2, col, 'G3 V', mp={YP: WP}, cc=True)


def g3_vi():
    R = lin((1, mono(1, [('b', 'c', 'z'), ('a', 'd', 'y'), ('c', 'e', 'x')], [], [('d', 'y'), ('e', 'x')])),
            (-1, mono(1, [('a', 'd', 'y')], [], [('d', 'y'), ('b', 'x')])))
    col = P(mono(1, [], [(C1, D1, 'a')], []), _Lc(), R)
    return build(2 * I, T2, col, 'G3 VI', mp={YP: WP}, cc=True)


THREE = {'G1 II 2a': g1_ii_2a, 'G1 III': g1_iii, 'G1 IV': g1_iv, 'G1 V(1)': g1_v1, 'G2 I(1)': g2_i1,
         'G3 I': g3_i, 'G3 II': g3_ii, 'G3 III(1)': g3_iii1, 'G3 V': g3_v, 'G3 VI': g3_vi}
TWO_ROWS = {'G1 II(1)': g1_ii_1, 'G1 V(2a)': g1_v2a, 'G1 V(2b UV)': g1_v2b_uv, 'G1 V(2b fin)': g1_v2b_fin,
            'G2 I(2)': g2_i2, 'G2 IV': g2_iv, 'G3 III(2)': g3_iii2}


# ---------------------------------------------------------------- the y' -> z terms of Group III Rows III, V, VI
# At large p+ the phase e^{-ik.(y'-z) p+/k+} makes the y' integral of K(y'-z)^m G(y') equal to
#   e^{-ik.(z-w)} (-2 pi i k^m/k^2) G(z) k+/(p+ + k+) = int_{w'} e^{-ik.(w'-w)} K(w'-z)^m G(z) k+/(p+ + k+),
# so every p+ integrand that tends to a constant gives log(vee/k+) x [y' -> z in G, K(y'-z) -> K(w'-z)].
K_PC = [(1, [KK(WP, 'z', 'y', 'w'), KK(XP, 'z', 'x', 'z')])]


def _yz(colour):
    return [(c, [(i, j, 'z' if p == YP else p) for i, j, p in Us], fs, rs) for c, Us, fs, rs in colour]


def g3_iii1_P():
    R = lin((1, anti(1, [('a', 'c', 'w')], [], ('c', 'y'), ('b', 'x'))),
            (-1, anti(1, [('b', 'c', 'z'), ('a', 'd', 'y'), ('c', 'e', 'x')], [], ('d', 'y'), ('e', 'x'))))
    col = _yz(P(mono(1, [], [(C1, D1, 'a')], []), _Lc(), R))
    return build(-I, K_PC, col, 'G3 III(1) yz', cc=True)


def g3_v_P():
    L = lin((1, mono(1, [(E1, C1, YP)], [], [(E1, XP)])), (-1, mono(1, [(C1, E1, XP)], [], [(E1, XP)])))
    R = lin((1, mono(1, [('c', 'e', 'x'), ('a', 'd', 'y')], [], [('e', 'x'), ('d', 'y')])),
            (-1, mono(1, [('a', 'd', 'w'), ('c', 'e', 'x')], [], [('e', 'x'), ('d', 'y')])))
    col = _yz(P(mono(1, [], [(C1, 'c', 'a')], []), L, R))
    return build(-2 * I, K_PC, col, 'G3 V yz', cc=True)


def g3_vi_P():
    R = lin((1, mono(1, [('b', 'c', 'z'), ('a', 'd', 'y'), ('c', 'e', 'x')], [], [('d', 'y'), ('e', 'x')])),
            (-1, mono(1, [('a', 'd', 'y')], [], [('d', 'y'), ('b', 'x')])))
    col = _yz(P(mono(1, [], [(C1, D1, 'a')], []), _Lc(), R))
    return build(-2 * I, K_PC, col, 'G3 VI yz', cc=True)


def g3_iii2_P(sign=-1):
    k = [(1, [KK(XP, 'z', 'x', 'z'), KK(WP, 'z', 'z', 'w')]), (-0.5, [KK(XP, 'z', 'x', 'z'), KK(WP, 'z', 'x', 'w')])]
    col = _yz(P(mono(1, [], [(C1, D1, 'a')], []), _Lc(), _RA('w')))
    return build(2 * sign, k, col, 'G3 III(2) yz', cc=True)


THREE_YZ = {'G3 III(1) yz': g3_iii1_P, 'G3 V yz': g3_v_P, 'G3 VI yz': g3_vi_P}
TWO_YZ = {'G3 III(2) yz': g3_iii2_P}
