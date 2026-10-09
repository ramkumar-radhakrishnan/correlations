"""Exact class-by-class bookkeeping of the symmetrized three-rho terms.

Terms are collected from
  'notes'  : every three-rho term of the notes (incl. c.c. and the four-rho remainders),
  'fixA'   : remove the alpha part of the Row II kernel and add -[Nbar_2, Abar^(1)],
  'fixB'   : remove part 2 (the edge-kernel terms) of Row III,
  'Iprime' : the diagonal row I' of the region-T reference,
and sorted into kernel classes and colour networks with exact coefficients (units of N).
"""
import os
import sys
from collections import OrderedDict
import numpy as np
import sympy as sp
import numerics as nm
from engine import simplify, canonical_rename, tex_colour, Nc, I as SI
from cancel3 import term_list, NC
import classes3 as C3
import compare_ref as CR

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference'))
import rowsT_sym as RT  # noqa: E402

ALPHA_GII = [(c, fs) for c, fs in C3.GII_EXPANDED if ('w', 'z', 'y', 'z') not in fs and ('w', 'z', 'x', 'z') not in fs]


def to_sym(c):
    c = complex(c)
    return sp.nsimplify(round(c.real, 12)) + SI * sp.nsimplify(round(c.imag, 12))


def ct(name, cc, sub, coef, colour, word, monos, cmap):
    return dict(name=name, cc=cc, sub=sub, coef=coef, colour=colour, word=word, monos=monos, cmap=cmap)


def notes_terms():
    out = []
    for key, cc, n, c, src, colour, word, cmap in term_list():
        monos = [(sp.nsimplify(kc), fs) for kc, fs in C3.monomials(src)]
        out.append(ct(key, cc, n + 1, to_sym(c), colour, word, monos, cmap))
    return out


def fix_terms():
    out = []
    for key, cc, n, c, src, colour, word, cmap in term_list():
        if key == 'GI.II':
            monos = [(-sp.Integer(kc), [('KK',) + f for f in fs]) for kc, fs in ALPHA_GII]
            out.append(ct('fixA: -alpha(Row II)', cc, n + 1, to_sym(c), colour, word, monos, cmap))
        if key == 'GI.III':
            monos = [(-sp.nsimplify(kc), fs) for kc, fs in src.kexpl[1:]]
            out.append(ct('fixB: -part 2(Row III)', cc, n + 1, to_sym(c), colour, word, monos, cmap))
    # the commutator -[Nbar_2, Abar^(1)] of Row II (and c.c.)
    import comm as g1fix
    for t in g1fix.comm_terms():
        cc = t.tag.endswith('*')
        monos = [(sp.Integer(1), [('KK', a, b, c_, d) for (a, b), (c_, d) in t.kern])]
        ident = {l: l for l in ('w', "w'", 'z')}
        out.append(ct('fixA: commutator', cc, 0, to_sym(t.coef), t.colour, t.rhos, monos, ident))
    return out


def iprime_terms():
    out = []
    L = RT.virt_rows()['G1-Ia']
    for t in L:
        t = RT.fix_kern(t)
        if len(t.rs) != 3:
            continue
        pos = lambda p: {'wb': "w'"}.get(p, p)
        colour = [('U', i, j, pos(p)) for i, j, p in t.Us] + [('f',) + tuple(f) for f in t.fs]
        word = [(i, pos(p)) for i, p in t.rs]
        monos = [(sp.Integer(1), [('KK', pos(a), pos(b), pos(c_), pos(d)) for (a, b), (c_, d) in t.kern])]
        ident = {l: l for l in ('w', "w'", 'z')}
        out.append(ct("row I' (reference)", False, 0, to_sym(2 * t.c), colour, word, monos, ident))
    return out


def classify(terms):
    entries = []
    for t in terms:
        rl = []
        for _, p in t['word']:
            if p not in rl:
                rl.append(p)
        for kc, facs in t['monos']:
            r = C3.roles(facs, t['cmap'], rl)
            full = dict(t['cmap'])
            full.update(r)
            sign = 1
            parts = []
            for f in facs:
                if f[0] == 'KK':
                    k, s = C3.kkey(*(full.get(l, l) for l in f[1:]))
                    sign *= s
                    parts.append(k)
                else:
                    parts.append(('Phi', full[f[2][0]], full['w']))
            kk = tuple(sorted(parts, key=str))
            cls = C3.CLASS_OF.get(kk)
            assert cls is not None, (t['name'], kk)
            col = [('U', u[1], u[2], full.get(u[3], u[3])) if u[0] == 'U' else u for u in t['colour']]
            col += [('rho', i, full.get(p, p)) for i, p in t['word']]
            c1, col1, rules = simplify(t['coef'] * kc * sign, col)
            entries.append(dict(name=t['name'], cc=t['cc'], sub=t['sub'], cls=cls, coeff=sp.expand(c1),
                                col=canonical_rename(col1) if col1 else col1))
    return entries


def group(entries, seed=7):
    ts, f = nm.structure_constants(NC)
    rng = np.random.default_rng(seed)
    labs = ['w', "w'", 'z', 'x', "x'", 'y']
    Um = {l: nm.random_adjoint(NC, ts, rng, True) for l in labs}
    rv = {l: rng.normal(size=NC * NC - 1) for l in labs}
    res = OrderedDict((c, []) for c in ['alpha', "alpha'", 'beta', "beta'", 'gamma', "gamma'"])
    for e in entries:
        fp = nm.eval_colour(e['col'], f, Um, rv, NC * NC - 1) if e['col'] else 0.0
        if abs(fp) < 1e-12 or e['coeff'] == 0:
            continue
        nets = res[e['cls']]
        for net in nets:
            r = fp / net['fp']
            if abs(r - round(r.real)) < 1e-9 and abs(round(r.real)) == 1:
                net['members'].append((e, int(round(r.real))))
                break
        else:
            nets.append(dict(fp=fp, rep=e, members=[(e, 1)]))
    return res


def summarize(res):
    out = OrderedDict()
    for cls, nets in res.items():
        rows = []
        for net in nets:
            tot = sp.nsimplify(sp.expand(sum(e['coeff'] * s for e, s in net['members'])))
            rows.append((net, tot))
        out[cls] = rows
    return out


if __name__ == '__main__':
    for label, terms in [('notes as written', notes_terms()),
                         ('notes + fixA + fixB', notes_terms() + fix_terms()),
                         ("notes + fixA + fixB + row I'", notes_terms() + fix_terms() + iprime_terms())]:
        S = summarize(group(classify(terms)))
        print('=====', label)
        for cls, rows in S.items():
            bad = [(net, tot) for net, tot in rows if tot != 0]
            print(f'  {cls:7s}: {len(rows)} networks, {len(bad)} with nonzero sum')
            for net, tot in bad:
                print(f'       {tot}  x  {tex_colour(net["rep"]["col"])}')
