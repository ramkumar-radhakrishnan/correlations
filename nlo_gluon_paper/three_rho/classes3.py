"""Symbolic sorting of the symmetrized three-rho terms into kernel classes.

After the symmetrization the charges commute, and all positions except the
measured-gluon positions w, w' are dummy.  In every kernel monomial we rename
the three charge positions by their role:
    x  : the charge whose LO emission kernel ends at w   (KK(x-w)),
    x' : the charge whose LO emission kernel ends at w'  (KK(x'-w')),
    y  : the third charge,
and the soft-gluon position is z.  Every monomial is then one of
    alpha  : KK(x-w).KK(x'-w') KK(x-z).KK(y-z)      (soft gluon between x and y)
    alpha' : KK(x-w).KK(x'-w') KK(x'-z).KK(y-z)     (soft gluon between x' and y)
    beta   : KK(x-w).KK(x'-w') KK(w-z).KK(y-z)      (soft gluon between w and y)
    beta'  : KK(x-w).KK(x'-w') KK(w'-z).KK(y-z)     (soft gluon between w' and y)
    gamma  : KK(x'-w').KK(y-z) Phi(x)               (Row III of Group I)
    gamma' : KK(x-w).KK(y-z) Phi'(x')               (its c.c.)
Inside each class the colour networks are simplified and compared.
"""
from collections import OrderedDict
import numpy as np
import sympy as sp
import numerics as nm
from engine import simplify, canonical_rename, tex_colour, Nc, I as SI
from cancel3 import term_list, NC

GII_EXPANDED = [
    (1, [("x'", "w'", 'x', 'w'), ('w', 'z', 'y', 'z')]),
    (-1, [("x'", "w'", 'x', 'w'), ('x', 'z', 'y', 'z')]),
    (-1, [("x'", "w'", 'y', 'w'), ('w', 'z', 'x', 'z')]),
    (1, [("x'", "w'", 'y', 'w'), ('x', 'z', 'y', 'z')]),
]

ORD = ['x', "x'", 'y', 'w', "w'", 'z']
CLASS_OF = {
    (('KK', ('x', 'w'), ("x'", "w'")), ('KK', ('x', 'z'), ('y', 'z'))): 'alpha',
    (('KK', ('x', 'w'), ("x'", "w'")), ('KK', ("x'", 'z'), ('y', 'z'))): "alpha'",
    (('KK', ('x', 'w'), ("x'", "w'")), ('KK', ('y', 'z'), ('w', 'z'))): 'beta',
    (('KK', ('x', 'w'), ("x'", "w'")), ('KK', ('y', 'z'), ("w'", 'z'))): "beta'",
    (('KK', ("x'", "w'"), ('y', 'z')), ('Phi', 'x', 'w')): 'gamma',
    (('KK', ('x', 'w'), ('y', 'z')), ('Phi', "x'", "w'")): "gamma'",
}


CLASS_OF = {tuple(sorted(k, key=str)): v for k, v in CLASS_OF.items()}


def nvec(a, b):
    return ((a, b), 1) if ORD.index(a) < ORD.index(b) else ((b, a), -1)


def kkey(a, b, c, d):
    v1, s1 = nvec(a, b)
    v2, s2 = nvec(c, d)
    return ('KK',) + tuple(sorted([v1, v2], key=lambda v: (ORD.index(v[0]), ORD.index(v[1])))), s1 * s2


def monomials(src):
    if src.key == 'GI.II':
        return [(sp.Integer(c), [('KK',) + f for f in fs]) for c, fs in GII_EXPANDED]
    return [(sp.nsimplify(c), fs) for c, fs in src.kexpl]


def roles(facs, m, rl):
    """find the native labels playing x (at w), x' (at w') and y."""
    xw = xpw = None
    for f in facs:
        if f[0] == 'KK':
            p, q, r, s = (m.get(l, l) for l in f[1:])
            for (a, b) in ((p, q), (q, p), (r, s), (s, r)):
                if b == 'w' and a in rl:
                    xw = a
                if b == "w'" and a in rl:
                    xpw = a
        elif 'Phi' in f[1]:
            arg = f[2][0]
            wl = m['w']
            if wl == 'w':
                xw = arg
            else:
                xpw = arg
    third = [l for l in rl if l not in (xw, xpw)]
    assert xw and xpw and xw != xpw and len(third) == 1, (facs, m, xw, xpw)
    return {xw: 'x', xpw: "x'", third[0]: 'y'}


def classify_all():
    out = []
    for key, cc, n, c, src, colour, word, cmap in term_list():
        rl = [p for _, p in word]
        for kc, facs in monomials(src):
            r = roles(facs, cmap, rl)
            full = dict(cmap)
            full.update(r)
            sign = 1
            parts = []
            for f in facs:
                if f[0] == 'KK':
                    k, s = kkey(*(full[l] for l in f[1:]))
                    sign *= s
                    parts.append(k)
                else:
                    # Phi(x) of the native row depends on native z, w
                    parts.append(('Phi', full[f[2][0]], full['w']))
            kk = tuple(sorted(parts, key=str))
            cls = CLASS_OF.get(kk)
            assert cls is not None, (key, kk)
            col = [('U', u[1], u[2], full[u[3]]) if u[0] == 'U' else u for u in colour]
            col += [('rho', i, full[p]) for i, p in word]
            coef = sp.nsimplify(complex(c).imag) * SI * kc * sign
            c1, col1, rules = simplify(coef, col)
            out.append(dict(src=key, cc=cc, sub=n + 1, cls=cls, coeff=sp.expand(c1),
                            col=canonical_rename(col1) if col1 else col1, rules=rules,
                            raw=col, rawcoef=coef))
    return out


def fingerprints(entries, seed=7, symmetric=True):
    ts, f = nm.structure_constants(NC)
    rng = np.random.default_rng(seed)
    labs = ['w', "w'", 'z', 'x', "x'", 'y']
    Um = {l: nm.random_adjoint(NC, ts, rng, symmetric) for l in labs}
    rv = {l: rng.normal(size=NC * NC - 1) for l in labs}
    for e in entries:
        e['fp'] = nm.eval_colour(e['col'], f, Um, rv, NC * NC - 1) if e['col'] else 0.0
        e['val'] = complex(sp.N(e['coeff'].subs(Nc, NC))) * e['fp']
    return entries


def group(entries):
    """class -> list of networks; each network: representative + members with weights."""
    res = OrderedDict()
    for cls in ['alpha', "alpha'", 'beta', "beta'", 'gamma', "gamma'"]:
        nets = []
        for e in [e for e in entries if e['cls'] == cls]:
            if e['fp'] == 0 or abs(e['fp']) < 1e-12:
                continue
            for net in nets:
                r = e['fp'] / net['fp']
                if abs(r - round(r.real)) < 1e-9 and abs(round(r.real)) == 1:
                    net['members'].append((e, int(round(r.real))))
                    break
            else:
                nets.append(dict(fp=e['fp'], rep=e, members=[(e, 1)]))
        res[cls] = nets
    return res


if __name__ == '__main__':
    E = fingerprints(classify_all())
    G = group(E)
    for cls, nets in G.items():
        print(f'=== class {cls}: {sum(len(n["members"]) for n in nets)} contributions, {len(nets)} colour networks')
        for net in nets:
            tot = sum(e['coeff'] * s for e, s in net['members'])
            tot = sp.simplify(tot)
            srcs = ', '.join(f"{e['src']}{'*' if e['cc'] else ''}.{e['sub']}({sp.nsimplify(e['coeff'] * s)})" for e, s in net['members'])
            flag = 'CANCELS' if tot == 0 else f'RESIDUAL {tot}'
            print(f'  {tex_colour(net["rep"]["col"])}\n      {flag}   <- {srcs}')
