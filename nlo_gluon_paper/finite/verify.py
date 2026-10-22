"""Check Finite_terms.tex row by row against the exact p+ integrals of the first-principles (or previously verified)
integrands.  For every row R and several (vee, Lambda):

    D_R(vee, Lambda) = exact_R(vee, Lambda) - finite_R(vee, Lambda)

must be (coefficient) x log(vee/Lambda): the finite file must contain everything except the log(vee/Lambda) term.
So D_R / log(vee/Lambda) must not depend on vee or Lambda (up to O(Lambda/k+) where the file sets Lambda -> 0 in
a plus prescription).  Output: verify.json."""
import json
import numpy as np
from scipy.integrate import quad
import env  # noqa: F401
import fp
import fp1
import user23 as us
import fixed23 as fx
from integ import pv
import ft

PI = np.pi
NAMES = ('x', 'y', 'xp', 'yp', 'z', 'w', 'wp')
KP = 1.0


def cquad(f, a, b, pts=None):
    kw = dict(limit=800, epsabs=0, epsrel=1e-11)
    if pts:
        kw['points'] = pts
    return quad(lambda p: np.real(f(p)), a, b, **kw)[0] + 1j * quad(lambda p: np.imag(f(p)), a, b, **kw)[0]


def cfg(seed, sym=True):
    return fp.Cfg(seed, sym=sym, names=NAMES)


# ------------------------------------------------------------------ exact rows
def exact_G1_III(c, kp, V, L):
    A, B = ft.AB(c)
    ps = kp * B / A

    def g(parts):
        return lambda p: fp1.row3(c, kp, p, L, sgt=-1, parts=parts) * (p * A - kp * B)
    tot = pv(g(('inst', 'ni')), L, kp - L, ps, A) + pv(g(('inst', 'ni')), kp + L, V, ps, A) + pv(g(('t4',)), L, V, ps, A)
    return tot / (2 * PI)


def exact_G1_IV(c, kp, V, L):
    return cquad(lambda p: fp1.row4(c, kp, p), L, V - kp) / (2 * PI)


def exact_G1_V1(c, kp, V, L):
    z, x = c.pt['z'].copy(), c.pt['x'].copy()

    def f(p):
        c.pt['w'] = z - p / kp * (z - x)
        return fp1.row5_part1(c, kp, p)
    return cquad(f, L, kp - L) / (2 * PI)


def exact_noPhase(fun, c, kp, V, L):
    ph = ft.phase(c, 'w', 'wp')
    return ph * cquad(lambda p: fun(c, kp, p), L, V - kp) / (2 * PI)


def exact_C(fun, c, kp, V, L):
    """C rows: w' fixed by the delta function, w' = ((p+k) y' - p z)/k"""
    yp, z = c.pt['yp'].copy(), c.pt['z'].copy()

    def f(p):
        c.pt['wp'] = ((p + kp) * yp - p * z) / kp
        return fun(c, kp, p) * np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['w']))
    return cquad(f, L, V - kp) / (2 * PI)


# ------------------------------------------------------------------ finite file
def fin_G1_III(c, kp, V, L):
    return sum(ft.G1_III(c, kp, V, L).values())


def fin_G1_V1(c, kp, V, L):
    a = cquad(lambda p: ft.G1_V1_p(c, kp, p, 'p'), 0.0, kp) / (2 * PI)
    b = cquad(lambda q: ft.G1_V1_p(c, kp, q, 'q'), 0.0, kp) / (2 * PI)
    return a + b + ft.G1_V1_rest(c, kp, V, L)


def fin_C(pfun, rfun):
    def f(c, kp, V, L):
        return cquad(lambda p: pfun(c, kp, p), 0.0, V - kp) / (2 * PI) + rfun(c, kp, V, L)
    return f


ROWS = {
    'G1-III': (exact_G1_III, fin_G1_III),
    'G1-IV': (exact_G1_IV, ft.G1_IV),
    'G1-V(1)': (exact_G1_V1, fin_G1_V1),
    'G2-I(1)': (lambda c, k, V, L: exact_noPhase(lambda c, k, p: fx.g2_1_termA(c, k, p) + fx.g2_1_termB(c, k, p), c, k, V, L),
                lambda c, k, V, L: sum(ft.G2_I1(c, k, V, L))),
    'G2-I(2)': (lambda c, k, V, L: exact_noPhase(us.g2_2_line62, c, k, V, L),
                lambda c, k, V, L: sum(ft.G2_I2(c, k, V, L).values())),
    'G3-I': (lambda c, k, V, L: exact_noPhase(fx.g3_1_line44, c, k, V, L), ft.G3_I),
    'G3-II': (lambda c, k, V, L: exact_noPhase(fx.g3_2_line55, c, k, V, L), ft.G3_II),
    'G3-III(1)': (lambda c, k, V, L: exact_C(fx.g3_3_part1_line68, c, k, V, L), fin_C(ft.G3_III1_p, ft.G3_III1_rest)),
    'G3-III(2)': (lambda c, k, V, L: exact_C(fx.g3_3_part2_line80, c, k, V, L), fin_C(ft.G3_III2_p, ft.G3_III2_rest)),
    'G3-V': (lambda c, k, V, L: exact_C(fx.g3_5_line101, c, k, V, L), fin_C(ft.G3_V_p, ft.G3_V_rest)),
    'G3-VI': (lambda c, k, V, L: exact_C(fx.g3_6_line120, c, k, V, L), fin_C(ft.G3_VI_p, ft.G3_VI_rest)),
}
VL = [(6.0, 1e-3), (40.0, 1e-3), (6.0, 1e-5), (40.0, 1e-5), (300.0, 1e-4)]


def run(name, seeds=(1, 2, 3)):
    ex, fi = ROWS[name]
    res = []
    for seed in seeds:
        c = cfg(seed)
        rows = []
        for V, L in VL:
            e = ex(c, KP, V, L)
            f = fi(c, KP, V, L)
            rows.append(dict(V=V, L=L, exact=e, finite=f, ratio=(e - f) / np.log(V / L)))
        rats = np.array([r['ratio'] for r in rows])
        scale = max(abs(r['exact']) for r in rows)
        res.append(dict(seed=seed, spread=float(np.max(np.abs(rats - rats[-1])) / scale),
                        coef_over_scale=float(abs(rats[-1]) / scale), rows=rows))
    return res


def tojson(o):
    if isinstance(o, dict):
        return {k: tojson(v) for k, v in o.items()}
    if isinstance(o, list):
        return [tojson(v) for v in o]
    if isinstance(o, (complex, np.complexfloating)):
        return [float(o.real), float(o.imag)]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


if __name__ == '__main__':
    import sys
    names = sys.argv[1:] or list(ROWS)
    out = {}
    for n in names:
        out[n] = run(n)
        print('%-10s spread of D/log(vee/Lambda) (rel. to |exact|): %s' % (n, ['%.1e' % r['spread'] for r in out[n]]))
    try:
        old = json.load(open('verify.json'))
    except Exception:
        old = {}
    old.update(tojson(out))
    json.dump(old, open('verify.json', 'w'), indent=1)
