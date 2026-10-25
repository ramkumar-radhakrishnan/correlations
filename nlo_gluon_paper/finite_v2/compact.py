"""Two-rho remnants of the three-rho finite rows in compact form, by the pair rule.

For a colour operator O = sum T^{a1 a2 a3} rho^{a1}(u1) rho^{a2}(u2) rho^{a3}(u3) (charges in the written order;
{B,C} = BC + CB written out), the reordering into the fully symmetric part plus commutators gives the two-rho part
     O^(2) = 1/2 sum_{i<j}  O | rho^{ai}(ui) rho^{aj}(uj) -> i f^{ai aj g} rho^g(ui) delta(ui - uj).
For A{BC} this is [A,B]C + [A,C]B (the [B,C] terms cancel), for {BC}A it is [B,A]C + [C,A]B.
The merged point is always the one that comes first in the row's list of charge labels.

Terms are combined exactly; pairs that differ only in the point of one Wilson line are combined into
[U(p) - U(q)].  The result is checked numerically against the reordering engine (../two_rho/core.symmetrize,
which reproduced the Sec. 2 lists of Evolution_terms.tex exactly), merge pattern by merge pattern.
Output: compact.pkl, compact.json"""
import json
import os
import pickle
from collections import OrderedDict, defaultdict
import numpy as np
import env  # noqa: F401
import decomp as DC
import numerics as nm

TS, FF = nm.structure_constants(3)


# ---------------------------------------------------------------- generation
def ordered(form, word):
    if form == 'plain':
        return [list(word)]
    if form == 'A{BC}':
        A, B, C = word
        return [[A, B, C], [A, C, B]]
    if form == '{BC}A':
        B, C, A = word
        return [[B, C, A], [C, B, A]]
    raise ValueError(form)


def pair_terms(src):
    """-> dict merge(keep, gone) -> {key: coef};  key = (fs, Us, rs) with Us as ('U', i, j, p)"""
    order = list(src.rho_labels)
    out = defaultdict(lambda: defaultdict(complex))
    for sign, colour, word in src.subterms:
        Us = [('U', u[1], u[2], u[3]) for u in colour if u[0] == 'U']
        fs = [tuple(u[1:]) for u in colour if u[0] == 'f']
        for w in ordered(src.form, word):
            for i in range(3):
                for j in range(i + 1, 3):
                    (ai, pi), (aj, pj) = w[i], w[j]
                    # [rho^ai(pi), rho^aj(pj)] = i f^{ai aj g} rho^g(pi) delta(pi - pj)
                    keep, gone = (pi, pj) if order.index(pi) < order.index(pj) else (pj, pi)
                    mp = lambda p, gone=gone, keep=keep: keep if p == gone else p
                    rest = [(a, mp(p)) for k, (a, p) in enumerate(w) if k not in (i, j)]
                    rs = tuple(sorted(rest + [('g', keep)], key=lambda r: order.index(r[1]) if r[1] in order else 9))
                    fsg, sg = (ai, aj, 'g'), 1
                    if ai > aj:                                 # f antisymmetric: canonical order
                        fsg, sg = (aj, ai, 'g'), -1
                    key = (tuple(fs + [fsg]), tuple(('U', u[1], u[2], mp(u[3])) for u in Us), rs)
                    out[(keep, gone)][key] += 0.5j * sign * sg
    return out


def clean(d):
    return {k: v for k, v in d.items() if abs(v) > 1e-12}


# ---------------------------------------------------------------- combine into differences
def combine(d):
    """repeatedly combine c*X(..U^{ij}(p)..) - c*X(..U^{ij}(q)..) into c*X(..[U(p)-U(q)]^{ij}..)"""
    items = dict(d)
    changed = True
    while changed:
        changed = False
        keys = list(items)
        for n, k1 in enumerate(keys):
            if k1 not in items:
                continue
            for k2 in keys[n + 1:]:
                if k2 not in items or k1 == k2:
                    continue
                f1, U1, r1 = k1
                f2, U2, r2 = k2
                if f1 != f2 or r1 != r2 or len(U1) != len(U2):
                    continue
                diff = [m for m in range(len(U1)) if U1[m] != U2[m]]
                if len(diff) != 1:
                    continue
                m = diff[0]
                u1, u2 = U1[m], U2[m]
                if u1[0] != 'U' or u2[0] != 'U' or u1[1:3] != u2[1:3]:
                    continue
                c1, c2 = items[k1], items[k2]
                if abs(c1 + c2) > 1e-12:
                    continue
                newU = list(U1)
                newU[m] = ('D', u1[1], u1[2], u1[3:], u2[3:])
                nk = (f1, tuple(newU), r1)
                del items[k1], items[k2]
                items[nk] = items.get(nk, 0) + c1
                changed = True
                break
            if changed:
                break
    return items


# ---------------------------------------------------------------- numerics
def expand(key, c):
    """D factors back to U factors -> list of (coef, fs, Us, rs)"""
    fs, Us, rs = key
    terms = [(c, [])]
    for u in Us:
        if u[0] == 'U':
            terms = [(cc, L + [('U', u[1], u[2], u[3])]) for cc, L in terms]
        else:
            new = []
            for cc, L in terms:
                for s, pts in ((1, u[3]), (-1, u[4])):
                    new.append((s * cc, L + [('U', u[1], u[2], pts[0] if len(pts) == 1 else pts)]))
            terms = new
    return [(cc, fs, L, rs) for cc, L in terms]


def value(terms, U, r):
    tot = 0.0
    for c, fs, Us, rs in terms:
        fl = [u for u in Us] + [('f',) + tuple(f) for f in fs] + [('rho', i, p) for i, p in rs]
        tot += c * nm.eval_colour(fl, FF, U, r, 8)
    return tot


def engine_groups(row):
    """engine output (decomp.pkl) by unordered merged pair, mapped to the point kept by our rule"""
    store = pickle.load(open(os.path.join(DC.HERE, 'decomp.pkl'), 'rb'))[row]
    order = store['labels']
    out = defaultdict(list)
    for sig, name, ts in store['groups']:
        labs = order
        gone = [labs[i] for i in range(3) if labs[i] != sig[i]]
        kept = [sig[i] for i in range(3) if labs[i] != sig[i]]
        g, k = gone[0], kept[0]
        keep, gone_ = (g, k) if order.index(g) < order.index(k) else (k, g)
        for t in ts:
            mp = lambda p: keep if p in (g, k) else p
            out[(keep, gone_)].append((t['c'] * 3 ** t.get('ncp', 0), t['fs'],
                                       [('U', i, j, mp(p)) for i, j, p in t['Us']], [(i, mp(p)) for i, p in t['rs']]))
    return out


def check_row(row, mine, seeds=(1, 2, 3)):
    eng = engine_groups(row)
    res = {}
    for pair in set(eng) | set(mine):
        e = eng.get(pair, [])
        m = [x for k, c in mine.get(pair, {}).items() for x in expand(k, c)]
        rel = []
        for s in seeds:
            rng = np.random.default_rng(s)
            pts = {p for t in e + m for u in t[2] for p in [u[3]]} | {p for t in e + m for _, p in t[3]}
            U = {p: nm.random_adjoint(3, TS, rng, True) for p in sorted(pts)}
            r = {p: rng.normal(size=8) for p in sorted(pts)}
            a, b = value(e, U, r), value(m, U, r)
            scale = sum(abs(value([t], U, r)) for t in e + m) or 1.0
            rel.append(abs(a - b) / scale)
        res['%s<-%s' % pair] = float(max(rel))
    return res


if __name__ == '__main__':
    SRC = DC.SRC
    store, rep = OrderedDict(), OrderedDict()
    for row, key in DC.ROWS.items():
        pt = pair_terms(SRC[key])
        mine = OrderedDict()
        for pair, d in pt.items():
            d = clean(d)
            if not d:
                continue
            mine[pair] = combine(d)
        chk = check_row(row, mine)
        store[row] = dict(source=key, labels=SRC[key].rho_labels, cc=SRC[key].cc, form=SRC[key].form, groups=mine)
        rep[row] = dict(check=chk, n_terms={'%s<-%s' % p: len(v) for p, v in mine.items()})
        print('%-16s terms per merge: %-46s engine check: %s' % (row, rep[row]['n_terms'], {k: '%.0e' % v for k, v in chk.items()}))
    pickle.dump(store, open(os.path.join(DC.HERE, 'compact.pkl'), 'wb'))
    json.dump(rep, open(os.path.join(DC.HERE, 'compact.json'), 'w'), indent=1)


# ---------------------------------------------------------------- simplification: f f -> N_c delta, U U^T -> 1
def _perm_sign(seq, target):
    """sign of the permutation taking seq to target (both length 3, same elements)"""
    s = list(seq)
    sign = 1
    for i in range(3):
        j = s.index(target[i])
        if j != i:
            s[i], s[j] = s[j], s[i]
            sign = -sign
    return sign


def _relabel(key, old, new):
    fs, Us, rs = key
    m = lambda i: new if i == old else i
    fs = tuple(tuple(m(i) for i in f) for f in fs)
    Us2 = []
    for u in Us:
        if u[0] == 'U':
            Us2.append(('U', m(u[1]), m(u[2]), u[3]))
        else:
            Us2.append(('D', m(u[1]), m(u[2]), u[3], u[4]))
    rs = tuple((m(i), p) for i, p in rs)
    return (fs, tuple(Us2), rs)


def simplify_term(key, c, ncp=0):
    """apply f^{ab r}f^{ab s} = N_c delta^{rs} and U^{ij}(p)U^{ik}(p) = delta^{jk}, U^{ji}(p)U^{ki}(p) = delta^{jk}"""
    changed = True
    while changed:
        changed = False
        fs, Us, rs = key
        # f f
        for i in range(len(fs)):
            for j in range(i + 1, len(fs)):
                sh = [x for x in fs[i] if x in fs[j]]
                if len(sh) == 2:
                    r1 = [x for x in fs[i] if x not in sh][0]
                    r2 = [x for x in fs[j] if x not in sh][0]
                    e = _perm_sign(fs[i], (sh[0], sh[1], r1)) * _perm_sign(fs[j], (sh[0], sh[1], r2))
                    nfs = tuple(f for k, f in enumerate(fs) if k not in (i, j))
                    key = _relabel((nfs, Us, rs), r2, r1)
                    c, ncp = c * e, ncp + 1
                    changed = True
                    break
            if changed:
                break
        if changed:
            continue
        # U U^T
        for i in range(len(Us)):
            for j in range(i + 1, len(Us)):
                u, v = Us[i], Us[j]
                if u[0] != 'U' or v[0] != 'U' or u[3] != v[3]:
                    continue
                pair = None
                if u[1] == v[1] and u[1] not in (u[2], v[2]):
                    pair = (u[2], v[2])
                elif u[2] == v[2] and u[2] not in (u[1], v[1]):
                    pair = (u[1], v[1])
                if pair is None:
                    continue
                nUs = tuple(w for k, w in enumerate(Us) if k not in (i, j))
                key = _relabel((fs, nUs, rs), pair[1], pair[0])
                changed = True
                break
            if changed:
                break
    return key, c, ncp


def simplify_group(d):
    acc = defaultdict(complex)
    for k, c in d.items():
        k2, c2, n = simplify_term(k, c)
        acc[(k2, n)] += c2
    # recombine differences within each N_c power
    out = {}
    for n in sorted({n for _, n in acc}):
        sub = {k: v for (k, m), v in acc.items() if m == n and abs(v) > 1e-12}
        for k, v in combine(sub).items():
            out[(k, n)] = v
    return out


def expand_n(keyn, c):
    k, n = keyn
    return [(cc * 3 ** n, fs, Us, rs) for cc, fs, Us, rs in expand(k, c)]


def check_simplified(row, simp, seeds=(1, 2, 3)):
    eng = engine_groups(row)
    res = {}
    for pair in set(eng) | set(simp):
        e = eng.get(pair, [])
        m = [x for k, c in simp.get(pair, {}).items() for x in expand_n(k, c)]
        rel = []
        for s in seeds:
            rng = np.random.default_rng(s)
            pts = {u[3] for t in e + m for u in t[2]} | {p for t in e + m for _, p in t[3]}
            U = {p: nm.random_adjoint(3, TS, rng, True) for p in sorted(pts)}
            r = {p: rng.normal(size=8) for p in sorted(pts)}
            a, b = value(e, U, r), value(m, U, r)
            scale = sum(abs(value([t], U, r)) for t in e) or 1.0
            rel.append(abs(a - b) / scale)
        res['%s<-%s' % pair] = float(max(rel))
    return res


def run_simplify():
    store = pickle.load(open(os.path.join(DC.HERE, 'compact.pkl'), 'rb'))
    out, rep = OrderedDict(), OrderedDict()
    for row, d in store.items():
        g = OrderedDict((pair, simplify_group(terms)) for pair, terms in d['groups'].items())
        g = OrderedDict((p, t) for p, t in g.items() if t)
        chk = check_simplified(row, g)
        out[row] = dict(d, groups=g)
        rep[row] = dict(check=chk, n_terms={'%s<-%s' % p: len(v) for p, v in g.items()})
        print('%-16s simplified terms per merge: %-40s engine check max %.1e' % (row, rep[row]['n_terms'], max(chk.values())))
    pickle.dump(out, open(os.path.join(DC.HERE, 'compact_simpl.pkl'), 'wb'))
    json.dump(rep, open(os.path.join(DC.HERE, 'compact_simpl.json'), 'w'), indent=1)
