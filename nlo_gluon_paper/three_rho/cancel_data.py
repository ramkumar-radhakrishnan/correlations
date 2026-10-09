"""Data for the explicit cancellation check of the symmetrized three-rho terms.

Every symmetrized three-rho term is split into 'entries':
    one source  x  one sub-term  x  one kernel monomial  (x  amplitude or c.c.).
Each entry is brought to canonical variables
    x  : charge whose LO kernel ends at w,   x' : charge whose LO kernel ends at w',
    y  : third charge,                        z  : soft gluon,
and sorted into one of six kernel classes.  Inside a class the colour structures
are simplified and grouped into colour networks (equal up to a sign); the
coefficients of each network are summed exactly.
"""
from collections import OrderedDict
import numpy as np
import sympy as sp
import numerics as nm
from engine import simplify, canonical_rename, I as SI
from sources import SOURCES, Source, KK
from cancel3 import CANON, NC
import classes3 as C3

LETTERS = 'abcd'
CLASSES = ['alpha', "alpha'", 'beta', "beta'", 'gamma', "gamma'"]

# ------------------------------------------------------------------ the replacement D of Row III part 2
_III = [s for s in SOURCES if s.key == 'GI.III'][0]
D_SOURCE = Source(
    'D', r'Replacement $\mathcal D$ of Row III, part 2', 'correction', -2 * SI, ("w'", 'w'), r'K_{\mathcal D}',
    ["x'", 'x', 'y'], [(1, [KK("x'", "w'", 'y', 'w'), KK('x', 'z', 'w', 'z')])], 'plain', True,
    [(sg, [('U', u[1], u[2], 'w') if (u[0] == 'U' and u[3] == 'z') else u for u in col], word)
     for sg, col, word in _III.subterms],
    ["x'", 'x', 'y'], ["w'", 'w', 'z'])


def label(src, n, m, nmono, cc):
    return '%s.%d%s%s' % (src.key, n + 1, LETTERS[m] if nmono > 1 else '', '*' if cc else '')


def canon_map(src, cc):
    base = dict(CANON.get(src.key, {"w'": "w'", 'w': 'w', 'z': 'z'}))
    if cc:
        P1, P0 = src.phase
        base[P1], base[P0] = base[P0], base[P1]
    return base


def source_entries(src):
    """all entries of one source; also returns the renaming of every kernel monomial."""
    mult = 1 if src.form == 'plain' else 2
    monos = C3.monomials(src)
    entries, renamings = [], OrderedDict()
    for cc in ((False, True) if src.cc else (False,)):
        cmap = canon_map(src, cc)
        for n, (sign, colour, word) in enumerate(src.subterms):
            rl = [p for _, p in word]
            c0 = src.pref * mult * sign
            if cc:
                c0 = sp.conjugate(c0)
            for m, (kc, facs) in enumerate(monos):
                full = dict(cmap)
                full.update(C3.roles(facs, cmap, rl))
                ksign, parts, canon_facs = 1, [], []
                for f in facs:
                    if f[0] == 'KK':
                        k, s = C3.kkey(*(full[l] for l in f[1:]))
                        ksign *= s
                        parts.append(k)
                        canon_facs.append(('KK',) + tuple(full[l] for l in f[1:]))
                    else:
                        parts.append(('Phi', full[f[2][0]], full['w']))
                        canon_facs.append(('Phi', full[f[2][0]], full['w']))
                cls = C3.CLASS_OF[tuple(sorted(parts, key=str))]
                renamings[(cc, m)] = dict(full=full, cls=cls, ksign=ksign, kc=kc, native=facs, canon=canon_facs)
                col = [('U', u[1], u[2], full[u[3]]) if u[0] == 'U' else u for u in colour]
                col += [('rho', i, full[p]) for i, p in word]
                coeff, col1, rules = simplify(c0 * kc * ksign, col)
                entries.append(dict(label=label(src, n, m, len(monos), cc), src=src.key, sub=n + 1, mono=m, cc=cc,
                                    cls=cls, coeff=sp.nsimplify(sp.expand(coeff)),
                                    col=canonical_rename(col1) if col1 else col1, rules=rules))
    return entries, renamings


def all_entries(sources=SOURCES):
    out = []
    for src in sources:
        out += source_entries(src)[0]
    return out


# ------------------------------------------------------------------ the two corrections
def deleted_by_correction1(e):
    """Row II: alpha part of the kernel removed from the U(w) A^(3) sub-terms 1, 3 (monomials b, d)."""
    return e['src'] == 'GI.II' and e['sub'] in (1, 3) and e['mono'] in (1, 3)


def deleted_by_correction2(e):
    """Row III: part 2 (the Phi terms, monomial b) removed."""
    return e['src'] == 'GI.III' and e['mono'] == 1


def corrected_entries():
    keep = [e for e in all_entries() if not (deleted_by_correction1(e) or deleted_by_correction2(e))]
    return keep + source_entries(D_SOURCE)[0]


# ------------------------------------------------------------------ colour networks
def _setup(seed):
    ts, f = nm.structure_constants(NC)
    rng = np.random.default_rng(seed)
    labs = ['w', "w'", 'z', 'x', "x'", 'y']
    Um = {l: nm.random_adjoint(NC, ts, rng, True) for l in labs}
    rv = {l: rng.normal(size=NC * NC - 1) for l in labs}
    return f, Um, rv


def fingerprint(col, setup):
    f, Um, rv = setup
    return nm.eval_colour(col, f, Um, rv, NC * NC - 1)


def networks(entries, seeds=(7, 13)):
    """class -> list of networks {rep, members[(entry, sign)], total}.  The sign relating an entry
    to the representative is fixed with the first seed and confirmed with the second."""
    setups = [_setup(s) for s in seeds]
    res = OrderedDict((c, []) for c in CLASSES)
    for e in entries:
        fps = [fingerprint(e['col'], s) for s in setups]
        assert abs(fps[0]) > 1e-10, e['label']
        nets = res[e['cls']]
        for net in nets:
            r = fps[0] / net['fps'][0]
            if abs(abs(r) - 1) < 1e-9:
                s = int(round(r.real))
                assert all(abs(a - s * b) < 1e-9 * max(1, abs(a)) for a, b in zip(fps, net['fps'])), e['label']
                net['members'].append((e, s))
                break
        else:
            nets.append(dict(fps=fps, rep=e, members=[(e, 1)]))
    for nets in res.values():
        for k, net in enumerate(nets):
            net['total'] = sp.nsimplify(sp.expand(sum(e['coeff'] * s for e, s in net['members'])))
    return res


def numbered(res):
    """attach a global network number (class letter + index) to every network."""
    short = {'alpha': r'\alpha', "alpha'": r"\alpha'", 'beta': r'\beta', "beta'": r"\beta'",
             'gamma': r'\gamma', "gamma'": r"\gamma'"}
    for cls, nets in res.items():
        for k, net in enumerate(nets):
            net['name'] = '%s_{%d}' % (short[cls], k + 1)
    return res


if __name__ == '__main__':
    E = all_entries()
    print('entries as written:', len(E))
    for title, ents in (('as written', E), ('corrected', corrected_entries())):
        R = networks(ents)
        nz = sum(1 for nets in R.values() for n in nets if n['total'] != 0)
        print(title, {c: len(n) for c, n in R.items()}, 'networks with nonzero sum:', nz)
