"""Compare the symmetrized three-rho content of the notes with the region-T
leading-log reference, kernel monomial by kernel monomial.

Reference: the region-T rows built from the factorized LL wave function
(rows/rowsT_sym.py of the earlier work, checked against an exact operator model).
Its rows are normalized as d3N = (1/(2pi)^3) g^4/(8 pi^5 k+) dY int [Row] = 2 N [Row];
the coefficient of log(vee/Lambda) is the region-T coefficient.

The three-rho (fully symmetrized) part of a row is
  three-rho words : the same word with commuting charges,
  four-rho words  : 1/2 sum_{i<j} S([X_i,X_j] * rest)   (check_pbw4.py).
"""
import itertools
import os
import sys
from collections import defaultdict
import numpy as np
import numerics as nm
from cancel3 import term_list, kernel_value, NC

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference'))
import rowsT_sym as RT  # noqa: E402

EXTMAP = {'w': 'w', 'wb': "w'", 'z': 'z'}
ORDER = ['w', "w'", 'z', 'p1', 'p2', 'p3']


# ---------------------------------------------------------------- generic classical 3-rho term
class C3:
    """coef * kernel(dot pairs) * colour(U, f) * rho rho rho (commuting)."""
    def __init__(self, coef, colour, rhos, kern, tag):
        self.coef, self.colour, self.rhos, self.kern, self.tag = coef, colour, rhos, kern, tag


def ref_terms():
    rows = {**RT.real_rows(), **RT.virt_rows()}
    out = []
    for name, L in rows.items():
        for t in L:
            t = RT.fix_kern(t)
            pos = lambda p: EXTMAP.get(p, p)
            Us = [('U', i, j, pos(p)) for i, j, p in t.Us]
            fs = [('f',) + tuple(f) for f in t.fs]
            rs = [(i, pos(p)) for i, p in t.rs]
            kern = [((pos(a), pos(b)), (pos(c), pos(d))) for (a, b), (c, d) in t.kern]
            c = 2 * t.c                     # 2 N per unit of [Row]
            if len(rs) == 3:
                out.append(C3(c, Us + fs, rs, kern, '3:' + name))
            elif len(rs) == 4:
                for i, j in itertools.combinations(range(4), 2):
                    (xi, pi), (yi, pj) = rs[i], rs[j]
                    g = 'g%d%d' % (i, j)
                    sub = lambda p: pi if p == pj else p
                    rest = [rs[k] for k in range(4) if k not in (i, j)]
                    newr = [(g, pi)] + [(a, sub(p)) for a, p in rest]
                    newc = [('U', u[1], u[2], sub(u[3])) for u in Us] + fs + [('f', xi, yi, g)]
                    newk = [((sub(a), sub(b)), (sub(cc), sub(d))) for (a, b), (cc, d) in kern]
                    out.append(C3(c * 0.5j, newc, newr, newk, '4:' + name))
    return out


def user_terms():
    out = []
    for key, cc, n, c, src, colour, word, cmap in term_list():
        out.append(('user', key, cc, n, c, src, colour, word, cmap))
    return out


def kk_key(a, b, c, d):
    def nv(p, q):
        return ((p, q), 1) if ORDER.index(p) < ORDER.index(q) else ((q, p), -1)
    v1, s1 = nv(a, b)
    v2, s2 = nv(c, d)
    if v1 == v2:
        return ('inv2', v1), s1 * s2
    return ('KK',) + tuple(sorted([v1, v2], key=lambda v: (ORDER.index(v[0]), ORDER.index(v[1])))), s1 * s2


GII_EXPANDED = [
    (1, [("x'", "w'", 'x', 'w'), ('w', 'z', 'y', 'z')]),
    (-1, [("x'", "w'", 'x', 'w'), ('x', 'z', 'y', 'z')]),
    (-1, [("x'", "w'", 'y', 'w'), ('w', 'z', 'x', 'z')]),
    (1, [("x'", "w'", 'y', 'w'), ('x', 'z', 'y', 'z')]),
]


def user_monomials(src):
    if src.key == 'GI.II':
        return [(c, [('KK',) + f for f in fs]) for c, fs in GII_EXPANDED]
    return [(float(c), fs) for c, fs in src.kexpl]


class Accum:
    def __init__(self, seed=17, symmetric=True):
        ts, self.f = nm.structure_constants(NC)
        rng = np.random.default_rng(seed)
        self.U = {l: nm.random_adjoint(NC, ts, rng, symmetric) for l in ORDER}
        self.r = {l: rng.normal(size=NC * NC - 1) for l in ORDER}
        self.data = defaultdict(lambda: defaultdict(complex))

    def col(self, facs):
        return nm.eval_colour(facs, self.f, self.U, self.r, NC * NC - 1)

    def add_ref(self, t, prefix='REF '):
        rl = []
        for _, p in t.rhos:
            if p not in rl:
                rl.append(p)
        assert len(rl) == 3, (t.tag, t.rhos)
        for perm in itertools.permutations(['p1', 'p2', 'p3']):
            m = dict(zip(rl, perm))
            g = lambda p: m.get(p, p)
            parts, sign = [], 1
            for (a, b), (c, d) in t.kern:
                k, s = kk_key(g(a), g(b), g(c), g(d))
                parts.append(k)
                sign *= s
            key = tuple(sorted(parts, key=str))
            facs = [('U', u[1], u[2], g(u[3])) if u[0] == 'U' else u for u in t.colour]
            facs += [('rho', i, g(p)) for i, p in t.rhos]
            self.data[key][prefix + t.tag] += t.coef * sign * self.col(facs) / 6

    def add_user(self, key_, cc, n, c, src, colour, word, cmap):
        rl = [p for _, p in word]
        name = 'USR ' + key_ + ('*' if cc else '')
        for perm in itertools.permutations(['p1', 'p2', 'p3']):
            m = dict(cmap)
            m.update(dict(zip(rl, perm)))
            facs = [('U', u[1], u[2], m[u[3]]) if u[0] == 'U' else u for u in colour]
            facs += [('rho', i, m[p]) for i, p in word]
            cv = self.col(facs)
            for kc, fs in user_monomials(src):
                parts, sign = [], 1
                for f in fs:
                    if f[0] == 'KK':
                        k, s = kk_key(*(m[l] for l in f[1:]))
                        parts.append(k)
                        sign *= s
                    else:
                        parts.append(('Phi', m[f[2][0]], m['w']))
                key = tuple(sorted(parts, key=str))
                self.data[key][name] += c * kc * sign * cv / 6


def run(seed=17):
    A = Accum(seed)
    for t in ref_terms():
        A.add_ref(t)
    for u in user_terms():
        A.add_user(*u[1:])
    return A


if __name__ == '__main__':
    A = run()
    tot_ref = {k: sum(v for n, v in d.items() if n.startswith('REF')) for k, d in A.data.items()}
    tot_usr = {k: sum(v for n, v in d.items() if n.startswith('USR')) for k, d in A.data.items()}
    print('reference: max |sum over rows| per kernel monomial =', max(abs(v) for v in tot_ref.values()),
          ' (scale', max(sum(abs(v) for n, v in d.items() if n.startswith('REF')) for d in A.data.values()), ')')
    print('notes    : max |sum over rows| per kernel monomial =', max(abs(v) for v in tot_usr.values()))
    print('number of kernel monomials:', len(A.data))
