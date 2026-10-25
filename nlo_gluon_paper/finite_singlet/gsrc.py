"""Colour operators of the three-rho finite rows with the index placements for a general adjoint U.

Group I rows: as written (they agree with first principles for a general U, ../group1).
Groups II and III: the general-U placements of ../groups23/fixed23.py (correction M4 of Groups_II_III_check.pdf).
In the file (convention U^{bc} = U^{cb}) these Wilson lines are written with the other index order:
  R^{ab}(w,z,x)         : U^{bc}(z) -> U^{cb}(z)
  two-charge bracket     : U^{bc'}(z) U^{c'e'}(y) -> U^{c'b}(z) U^{c'e'}(y)      (Group II Row I)
  Group III Row I        : U^{bd'}(z) U^{d'e'}(y) -> U^{d'b}(z) U^{d'e'}(y)
  Group III Row II       : U^{bc'}(z) -> U^{c'b}(z)
  Group III Rows III,V,VI: U^{e'c'}(y') -> U^{c'e'}(y'), U^{bd'}(z) -> U^{d'b}(z), U^{bc}(z) -> U^{cb}(z)
check(): the classical value of every operator equals the fixed23 colour function (general U)."""
import copy
import json
import numpy as np
import env  # noqa: F401
from sources import SOURCES
import fixed23 as FX
from fp import F
import numerics as nm

SRC = {s.key: s for s in SOURCES}
FLIPS = {
    'GII.I1': [('b', "c'", 'z'), ('b', 'c', 'z')],
    'GII.I2': [('b', "c'", 'z'), ('b', 'c', 'z')],
    'GIII.I': [('b', "d'", 'z'), ('b', 'c', 'z')],
    'GIII.II': [('b', "c'", 'z'), ('b', 'c', 'z')],
    'GIII.III': [("e'", "c'", "y'"), ('b', "d'", 'z'), ('b', 'c', 'z')],
    'GIII.V': [("e'", "c'", "y'")],
    'GIII.VI': [("e'", "c'", "y'"), ('b', "d'", 'z'), ('b', 'c', 'z')],
}


def general(key):
    src = SRC[key]
    if key not in FLIPS:
        return src
    g = copy.deepcopy(src)
    fl = set(FLIPS[key])
    sub = []
    for sign, colour, word in src.subterms:
        col = []
        for u in colour:
            if u[0] == 'U' and (u[1], u[2], u[3]) in fl:
                col.append(('U', u[2], u[1], u[3]))
            else:
                col.append(u)
        sub.append((sign, col, word))
    g.subterms = sub
    return g


# ---------------------------------------------------------------- numerical check against fixed23
NAME = {"x'": 'xp', "y'": 'yp', "w'": 'wp'}


class C:
    def __init__(self, seed, symmetric=False):
        rng = np.random.default_rng(seed)
        ts, _ = nm.structure_constants(3)
        self.U = {p: nm.random_adjoint(3, ts, rng, symmetric) for p in ('x', 'xp', 'y', 'yp', 'z', 'w', 'wp')}
        self.rho = {p: rng.normal(size=8) for p in self.U}


def classical(src, c):
    tot = 0.0
    for sign, colour, word in src.subterms:
        mult = 2 if src.form in ('A{BC}', '{BC}A') else 1
        fl = [('U', u[1], u[2], NAME.get(u[3], u[3])) if u[0] == 'U' else u for u in colour]
        fl += [('rho', i, NAME.get(p, p)) for i, p in word]
        tot += sign * mult * nm.eval_colour(fl, F, c.U, c.rho, 8)
    return tot


def reference(key, c):
    U, r = c.U, c.rho
    if key == 'GII.I1':
        return np.sum(FX.two_fix(c, 'wp', 'xp', 'y') * FX.Rone_fix(c, 'w', 'z', 'x'))
    if key == 'GII.I2':
        return np.sum(FX.Rone_fix(c, 'wp', 'z', 'x') * FX.two_fix(c, 'w', 'xp', 'y'))
    if key == 'GIII.I':
        L = (np.einsum('aB,CD,Cb,B,D->ab', U['xp'], U['y'], U['z'], r['xp'], r['y'])
             - np.einsum('Db,DE,aC,C,E->ab', U['z'], U['y'], U['wp'], r['xp'], r['y']))
        return np.sum(L * FX.Rone_fix(c, 'w', 'z', 'x'))
    if key == 'GIII.II':
        L = (np.einsum('CD,aB,Cb,D,B->ab', U['y'], U['xp'], U['z'], r['y'], r['xp'])
             - np.einsum('aB,b,B->ab', U['xp'], r['y'], r['xp']))
        return np.sum(L * FX.Rone_fix(c, 'w', 'z', 'x'))
    if key == 'GIII.III':
        return np.einsum('CDa,CDb,ab->', F, FX.Lc_fix(c), FX.two_fix(c, 'w', 'y', 'x'))
    if key == 'GIII.V':
        return np.einsum('CDa,CDb,ab->', F, FX.Lc_fix(c), FX.X1amp_fix(c))
    if key == 'GIII.VI':
        return np.einsum('CDa,CDb,ab->', F, FX.Lc_fix(c), FX.X2amp_fix(c))
    raise KeyError(key)


def check():
    out = {}
    for key in FLIPS:
        g = general(key)
        rel_gen, rel_sym = [], []
        for s in (1, 2, 3):
            c = C(s, symmetric=False)
            a, b = classical(g, c), reference(key, c)
            rel_gen.append(abs(a - b) / abs(b))
            cs = C(s, symmetric=True)
            a, b = classical(g, cs), classical(SRC[key], cs)
            rel_sym.append(abs(a - b) / abs(b))
        c = C(7, symmetric=False)
        as_written = abs(classical(SRC[key], c) - reference(key, c)) / abs(reference(key, c))
        out[key] = dict(general_vs_fixed23=float(max(rel_gen)), general_vs_file_for_symmetric_U=float(max(rel_sym)),
                        file_vs_fixed23_general_U=float(as_written))
    return out


if __name__ == '__main__':
    res = check()
    print(json.dumps(res, indent=1))
    json.dump(res, open('gsrc.json', 'w'), indent=1)
