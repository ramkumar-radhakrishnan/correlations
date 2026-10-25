"""Colour-singlet projection of every two-rho colour structure of the final finite result.

Projection: rho^a(u) rho^b(v) -> delta^{ab} <rho^f(u) rho^f(v)>/(N_c^2-1); the colour network is then written with
adjoint traces (../singlet/trace.py, exact for a general SU(N_c) adjoint U).  After the projection, monomials that
differ only in the point of one Wilson line are combined into differences D(p,q) = U(p) - U(q) (dcombine).

Pieces (general-U index placements, gsrc.py):
  LO-type                       : LO colour (all LO-like terms: RC, V_vee, L_* term, DGLAP, lin, Q1)
  two-rho rows                  : G1 V(2a), G1 V(2b) (both brackets), G2 I(2) (O_58), G3 III(2)
  two-rho parts of three-rho rows: G1 II, G1 III, G1 V(1), K_IV family, K_C family, merge by merge
Checks (random general SU(3) adjoint U):
  * trace form (before and after dcombine) = direct contraction with rho^a rho^b -> delta^{ab};
  * pair-rule two-rho parts = reordering engine (../two_rho/core.symmetrize) for the general-U operators;
  * two-rho row colours = first-principles colour functions (../group1/user1.py, ../groups23/fixed23.py).
Output: project.pkl, project.json"""
import itertools
import json
import os
import pickle
from collections import OrderedDict, defaultdict
import numpy as np
import sympy as sp
import env  # noqa: F401
import trace as TR
from trace import Nc
import numerics as nm
import gsrc
import compact as CP
from core import term as cterm, symmetrize

HERE = os.path.dirname(os.path.abspath(__file__))
TS, FF = nm.structure_constants(3)
PORDER = ['x', "x'", 'y', "y'", 'z', 'w', "w'"]


# ---------------------------------------------------------------- term lists (plain U) of the two-rho pieces
def T(c, Us, fs, rs, ncp=0):
    return dict(c=complex(c), ncp=ncp, Us=list(Us), fs=list(fs), rs=list(rs))


def CL(P):
    """conjugate LO amplitude [U^{ab'}(x') - U^{ab'}(w')] rho^{b'}(x'): list of (sign, U)"""
    return [(1, ('a', "b'", "x'")), (-1, ('a', "b'", "w'"))]


def lo_terms():
    out = []
    for s1, p1 in ((1, "y'"), (-1, "x'")):
        for s2, p2 in ((1, 'y'), (-1, 'x')):
            out.append(T(s1 * s2, [('c', "e'", p1), ('c', 'e', p2)], [], [("e'", "x'"), ('e', 'x')]))
    return out


def v2a_terms():
    out = []
    for s1, u in CL(None):
        out.append(T(-s1, [u, ('d', 'b', 'x'), ('c', 'e', 'z')], [('a', 'd', 'c'), ('f', 'b', 'e')], [("b'", "x'"), ('f', 'y')]))
        out.append(T(s1, [u, ('a', 'f', 'x')], [], [("b'", "x'"), ('f', 'y')], ncp=1))
    return out


def v2b_terms(q):
    """CL^a N_c [U^{af}(y) - U^{af}(q)] rho^f(y);  q = x (bracket {..}_V2) or z (the (11/6 + 2 L_*) bracket)"""
    out = []
    for s1, u in CL(None):
        for s2, p in ((1, 'y'), (-1, q)):
            out.append(T(s1 * s2, [u, ('a', 'f', p)], [], [("b'", "x'"), ('f', 'y')], ncp=1))
    return out


def R_terms(W, X, sfx=''):
    """R^{ab}(W, z, X) = f^{cbd} U^{ac}(W) rho^d(X) - f^{acd} U^{cb}(z) U^{de}(X) rho^e(X)  (general U)"""
    c, d, e = 'c' + sfx, 'd' + sfx, 'e' + sfx
    return [(1, [('a', c, W)], [(c, 'b', d)], [(d, X)]),
            (-1, [(c, 'b', 'z'), (d, e, X)], [('a', c, d)], [(e, X)])]


def o58_terms():
    """O_58 = R^{ab}(w', z, x') R^{ab}(w, z, x)"""
    out = []
    for s1, U1, f1, r1 in R_terms("w'", "x'", "'"):
        for s2, U2, f2, r2 in R_terms('w', 'x'):
            out.append(T(s1 * s2, U1 + U2, f1 + f2, r1 + r2))
    return out


def iii2_terms():
    """f^{c'd'a} [U^{c'e'}(y') U^{d'b}(z) - U^{d'b}(z) U^{c'e'}(x')] rho^{e'}(x') R^{ab}(w, z, x)"""
    out = []
    for s1, P in ((1, "y'"), (-1, "x'")):
        for s2, U2, f2, r2 in R_terms('w', 'x'):
            out.append(T(s1 * s2, [("c'", "e'", P), ("d'", 'b', 'z')] + U2, [("c'", "d'", 'a')] + f2, [("e'", "x'")] + r2))
    return out


TWO = OrderedDict([
    ('LO', lo_terms()),
    ('V2a', v2a_terms()),
    ('V2b1', v2b_terms('x')),
    ('V2b2', v2b_terms('z')),
    ('G2I2', o58_terms()),
    ('III2', iii2_terms()),
])

# three-rho rows: (source key, weight) as in Eqs. (famIV), (famC) of Finite_terms_simplified
THREE = OrderedDict([
    ('II', [('GI.II', 1)]),
    ('III', [('GI.III', 1)]),
    ('V1', [('GI.V1', 1)]),
    ('famIV', [('GI.IV', 1), ('GIII.I', -1), ('GIII.II', -1), ('GII.I1', -0.5)]),
    ('famC', [('GIII.III', 1), ('GIII.V', 2), ('GIII.VI', 2)]),
])


def three_groups(name):
    """merge (keep, gone) -> list of plain two-rho terms (pair rule), weights included"""
    acc = defaultdict(lambda: defaultdict(complex))
    for key, wgt in THREE[name]:
        for pair, d in CP.pair_terms(gsrc.general(key)).items():
            for k, c in d.items():
                acc[pair][k] += wgt * c
    out = OrderedDict()
    for pair in sorted(acc, key=lambda p: (PORDER.index(p[0]), PORDER.index(p[1]))):
        ts = []
        for (fs, Us, rs), c in acc[pair].items():
            if abs(c) < 1e-12:
                continue
            ts.append(T(c, [(u[1], u[2], u[3]) for u in Us], fs, rs))
        if ts:
            out[pair] = ts
    return out


# ---------------------------------------------------------------- projection and collection
def project(terms):
    acc = defaultdict(lambda: sp.Integer(0))
    for t in terms:
        c, ws = TR.singlet(t)
        if c == 0:
            continue
        acc[ws] += c
    return {ws: sp.expand(c) for ws, c in acc.items() if sp.expand(c) != 0}


# ---------------------------------------------------------------- differences after the projection
def _variants(w):
    """all rotations and reflections of a word: (sign, tuple)"""
    out = []
    for s, ww in ((1, list(w)), TR.transpose_word(w)):
        for k in range(len(ww)):
            out.append((s, tuple(ww[k:] + ww[:k])))
    return out


def _tr_letter(l):
    if l[0] == 'U':
        return 1, ('Ud', l[1])
    if l[0] == 'Ud':
        return 1, ('U', l[1])
    if l[0] == 'D':
        return 1, ('Dd', l[1], l[2])
    if l[0] == 'Dd':
        return 1, ('D', l[1], l[2])
    return -1, l


def transpose_word(w):
    s, out = 1, []
    for l in reversed(w):
        si, li = _tr_letter(l)
        s *= si
        out.append(li)
    return s, out


TR.transpose_word = transpose_word        # D letters: (U(p)-U(q))^T = U^dag(p) - U^dag(q)


def canon_mono(ws, c):
    ws2 = []
    for w in ws:
        s, k = TR.canon_word(w)
        c *= s
        ws2.append(k)
    return tuple(sorted(ws2)), c


def _orient(l):
    """D(p,q) with p before q in PORDER (sign)"""
    if l[0] in ('D', 'Dd') and PORDER.index(l[1]) > PORDER.index(l[2]):
        return -1, (l[0], l[2], l[1])
    return 1, l


def _try_merge(m1, m2):
    """m1, m2 tuples of words; returns (sign2, new key) if m1 and (sign2 * m2) differ in exactly one U/Ud letter"""
    if len(m1) != len(m2) or [len(w) for w in m1] != sorted([len(w) for w in m1]) and False:
        return None
    if sorted(len(w) for w in m1) != sorted(len(w) for w in m2):
        return None
    for perm in itertools.permutations(range(len(m2))):
        if any(len(m1[i]) != len(m2[perm[i]]) for i in range(len(m1))):
            continue
        options = [_variants(m2[perm[i]]) for i in range(len(m1))]
        for combo in itertools.product(*options):
            sign = 1
            diff = []
            for i, (s, w2) in enumerate(combo):
                sign *= s
                w1 = m1[i]
                for k in range(len(w1)):
                    if w1[k] != w2[k]:
                        diff.append((i, k, w1[k], w2[k]))
                        if len(diff) > 1:
                            break
                if len(diff) > 1:
                    break
            if len(diff) != 1:
                continue
            i, k, a, b = diff[0]
            if a[0] != b[0] or a[0] not in ('U', 'Ud'):
                continue
            lt = 'D' if a[0] == 'U' else 'Dd'
            so, newl = _orient((lt, a[1], b[1]))
            nw = list(m1)
            w = list(m1[i])
            w[k] = newl
            nw[i] = tuple(w)
            return sign, so, tuple(nw)
    return None


def dcombine(d):
    """combine c X[U(p)] - c X[U(q)] -> c X[U(p) - U(q)] repeatedly"""
    items = dict(d)
    changed = True
    while changed:
        changed = False
        keys = sorted(items, key=lambda k: (-sum(len(w) for w in k), str(k)))
        for n, k1 in enumerate(keys):
            for k2 in keys[n + 1:]:
                r = _try_merge(k1, k2)
                if r is None:
                    continue
                sign, so, nk = r
                c1, c2 = items[k1], items[k2]
                if sp.expand(c1 + sign * c2) != 0:
                    continue
                del items[k1], items[k2]
                nk, cn = canon_mono(nk, c1 * so)
                items[nk] = sp.expand(items.get(nk, 0) + cn)
                if items[nk] == 0:
                    del items[nk]
                changed = True
                break
            if changed:
                break
    return items


# ---------------------------------------------------------------- numerics
def word_matrix(w, U, Tg, e):
    M = np.eye(8, dtype=complex)
    for l in w:
        if l[0] == 'U':
            M = M @ U[l[1]]
        elif l[0] == 'Ud':
            M = M @ U[l[1]].T
        elif l[0] == 'D':
            M = M @ (U[l[1]] - U[l[2]])
        elif l[0] == 'Dd':
            M = M @ (U[l[1]] - U[l[2]]).T
        else:
            M = M @ Tg[e]
    return M


def eval_mono(ws, c, U):
    Tg = [-1j * FF[e] for e in range(8)]
    cc = complex(sp.sympify(c).subs(Nc, 3))
    hasT = any(l[0] == 'T' for w in ws for l in w)
    rng = range(8) if hasT else range(1)
    tot = 0
    for e in rng:
        v = cc
        for w in ws:
            v *= np.trace(word_matrix(w, U, Tg, e))
        tot += v
    return tot


def eval_dict(d, U):
    return sum(eval_mono(ws, c, U) for ws, c in d.items())


def direct(terms, U):
    """sum_a of the network with both rho's -> e_a (no 1/(N_c^2-1))"""
    tot = 0
    for t in terms:
        (i1, p1), (i2, p2) = t['rs']
        fl = [('U', i, j, p) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']]
        c = t['c'] * 3 ** t.get('ncp', 0)
        for a in range(8):
            e = np.zeros(8)
            e[a] = 1
            if p1 == p2:
                raise ValueError('equal points')
            tot += c * nm.eval_colour(fl + [('rho', i1, p1), ('rho', i2, p2)], FF, U, {p1: e, p2: e}, 8)
    return tot


def randU(seed):
    rng = np.random.default_rng(seed)
    return {p: nm.random_adjoint(3, TS, rng, False) for p in PORDER}


def check(terms, d, dc, seeds=(1, 2, 3)):
    r1, r2 = [], []
    for s in seeds:
        U = randU(s)
        a = direct(terms, U)
        scale = sum(abs(direct([t], U)) for t in terms) or 1
        r1.append(abs(a - eval_dict(d, U)) / scale)
        r2.append(abs(a - eval_dict(dc, U)) / scale)
    return float(max(r1)), float(max(r2))


# ---------------------------------------------------------------- engine check of the pair rule (general U)
def engine_check(name, seeds=(1, 2)):
    """pair rule vs ../two_rho/core.symmetrize, classical value of the whole two-rho part (general U, random rho)"""
    res = []
    for s in seeds:
        rng = np.random.default_rng(100 + s)
        U = {p: nm.random_adjoint(3, TS, rng, False) for p in PORDER}
        r = {p: rng.normal(size=8) for p in PORDER}
        mine = 0
        for pair, ts in three_groups(name).items():
            for t in ts:
                fl = [('U', i, j, p) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']] + [('rho', i, p) for i, p in t['rs']]
                mine += t['c'] * nm.eval_colour(fl, FF, U, r, 8)
        eng, scale = 0, 0
        for key, wgt in THREE[name]:
            src = gsrc.general(key)
            for sign, colour, word in src.subterms:
                Us = [(u[1], u[2], u[3]) for u in colour if u[0] == 'U']
                fs = [tuple(u[1:]) for u in colour if u[0] == 'f']
                labs = list(src.rho_labels)
                for w in CP.ordered(src.form, word):
                    for t in symmetrize(cterm(sign, Us, fs, w, [('TAG',) + tuple(labs)], key), (2,)):
                        sig = t['kern'][0][1:]
                        pair = [(labs[i], sig[i]) for i in range(3) if labs[i] != sig[i]][0]
                        keep = min(pair, key=labs.index)
                        mp = lambda p: keep if p in pair else p
                        fl = ([('U', i, j, mp(p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']]
                              + [('rho', i, mp(p)) for i, p in t['rs']])
                        v = wgt * t['c'] * 3 ** t.get('ncp', 0) * nm.eval_colour(fl, FF, U, r, 8)
                        eng += v
                        scale += abs(v)
        res.append(abs(mine - eng) / scale)
    return float(max(res))


# ---------------------------------------------------------------- first-principles colour functions for the two-rho rows
def two_row_reference(name, seed=5):
    import fixed23 as FX
    import user1 as U1
    c = gsrc.C(seed, symmetric=False)
    NM = {"x'": 'xp', "y'": 'yp', "w'": 'wp'}
    terms = TWO[name]
    mine = 0
    for t in terms:
        fl = [('U', i, j, NM.get(p, p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']]
        fl += [('rho', i, NM.get(p, p)) for i, p in t['rs']]
        mine += t['c'] * 3 ** t['ncp'] * nm.eval_colour(fl, FF, c.U, c.rho, 8)
    U, r = c.U, c.rho
    conj = (U['xp'] - U['wp']) @ r['xp']
    if name == 'LO':
        ref = ((U['yp'] - U['xp']) @ r['xp']) @ ((U['y'] - U['x']) @ r['x'])
    elif name == 'V2a':
        ref = conj @ (-np.einsum('adc,fbe,db,ce,f->a', FF, FF, U['x'], U['z'], r['y']) + 3 * U['x'] @ r['y'])
    elif name == 'V2b1':
        ref = conj @ (3 * (U['y'] - U['x']) @ r['y'])
    elif name == 'V2b2':
        ref = conj @ (3 * (U['y'] - U['z']) @ r['y'])
    elif name == 'G2I2':
        ref = FX.g2_2_colour58(c)
    elif name == 'III2':
        ref = np.einsum('CDa,CDb,ab->', FF, FX.Lc_fix(c), FX.Rone_fix(c, 'w', 'z', 'x'))
    return float(abs(mine - ref) / abs(ref))


def group1_reference(seed=6):
    """Group I three-rho operators (as written) = ../group1/user1.py colour functions, general U"""
    import user1 as U1
    c = gsrc.C(seed, symmetric=False)
    out = {}
    U, r = c.U, c.rho
    refs = {
        'GI.III': U1.col3(c),
        'GI.IV': U1.col4(c),
        'GI.V1': U1.col5a(c),
        'GI.II': ((U['xp'] - U['wp']) @ r['xp']) @ (2 * np.einsum('ah,hbc,b,c->a', U['w'], FF, r['x'], r['y'])
                                                     - 2 * np.einsum('abc,bd,ce,d,e->a', FF, U['x'], U['y'], r['x'], r['y'])),
    }
    for k, v in refs.items():
        out[k] = float(abs(gsrc.classical(gsrc.SRC[k], c) - v) / abs(v))
    return out


if __name__ == '__main__':
    store = OrderedDict(two=OrderedDict(), three=OrderedDict())
    rep = OrderedDict(two=OrderedDict(), three=OrderedDict())
    for name, terms in TWO.items():
        d = project(terms)
        dc = dcombine(d)
        chk = check(terms, d, dc)
        store['two'][name] = dict(raw=d, comb=dc)
        rep['two'][name] = dict(n_raw=len(d), n_comb=len(dc), traces_vs_direct=chk[0], combined_vs_direct=chk[1],
                                colour_vs_first_principles=two_row_reference(name))
        print(name, rep['two'][name])
    for name in THREE:
        groups = three_groups(name)
        g_store, g_rep = OrderedDict(), OrderedDict()
        for pair, terms in groups.items():
            d = project(terms)
            dc = dcombine(d)
            chk = check(terms, d, dc)
            g_store[pair] = dict(raw=d, comb=dc)
            g_rep['%s<-%s' % pair] = dict(n_terms=len(terms), n_raw=len(d), n_comb=len(dc), traces_vs_direct=chk[0],
                                          combined_vs_direct=chk[1])
        store['three'][name] = g_store
        rep['three'][name] = dict(groups=g_rep, pair_rule_vs_engine=engine_check(name))
        print(name, json.dumps(rep['three'][name]))
    rep['group1_as_written_vs_first_principles'] = group1_reference()
    rep['general_sources'] = json.load(open(os.path.join(HERE, 'gsrc.json')))
    pickle.dump(store, open(os.path.join(HERE, 'project.pkl'), 'wb'))
    json.dump(rep, open(os.path.join(HERE, 'project.json'), 'w'), indent=1)
