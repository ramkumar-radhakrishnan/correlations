"""Three-rho finite rows: operator colour factor O = S(O) + O^(2).

For every three-rho row of Finite_terms_1.tex the colour operator O (the charges in the order written in the file;
{rho, rho} written out) is reordered with [rho^a(u), rho^b(v)] = i f^{abc} rho^c(u) delta(u - v):
  S(O)   fully symmetric (classical) three-rho part,
  O^(2)  two-rho part; each of its terms sets two charge points equal (merge pattern sigma).
The kernel of the row is a function of the charge points, so the two-rho part of the row is
   prefactor x sum_sigma K_row|_sigma x O^(2)_sigma     (+ c.c. for Groups I and III).
The colour words are those of ../three_rho/sources.py (identical to the finite file, see the note).
Output: decomp.pkl (per row: list of (sigma, terms)), decomp.json (sizes, checks)."""
import json
import os
import pickle
import sys
from collections import OrderedDict, defaultdict
import numpy as np
import env  # noqa: F401
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'two_rho')))
from core import term, symmetrize  # noqa: E402
from sources import SOURCES  # noqa: E402
from users import ordered_words  # noqa: E402

SRC = {s.key: s for s in SOURCES}
# finite row -> source key (same colour operator), and the charge labels of the row
ROWS = OrderedDict([
    ('G1-II (2a, 2b)', 'GI.II'),
    ('G1-III', 'GI.III'),
    ('G1-IV', 'GI.IV'),
    ('G1-V(1)', 'GI.V1'),
    ('G2-I(1) A', 'GII.I1'),
    ('G2-I(1) B', 'GII.I2'),
    ('G3-I', 'GIII.I'),
    ('G3-II', 'GIII.II'),
    ('G3-III(1)', 'GIII.III'),
    ('G3-V', 'GIII.V'),
    ('G3-VI', 'GIII.VI'),
])


def operator_terms(src):
    """O = sum_n sign_n colour_n word_n (form expanded), coefficient 1 per ordered word; kernel = TAG(charges)"""
    out = []
    tag = [('TAG',) + tuple(src.rho_labels)]
    for sign, colour, word in src.subterms:
        Us = [(u[1], u[2], u[3]) for u in colour if u[0] == 'U']
        fs = [tuple(u[1:]) for u in colour if u[0] == 'f']
        for w in ordered_words(src.form, word):
            out.append(term(sign, Us, fs, w, tag, src.key))
    return out


def decompose(src):
    O = operator_terms(src)
    two, three = [], []
    for t in O:
        two += symmetrize(t, (2,))
        three += symmetrize(t, (3,))
    groups = defaultdict(list)
    for t in two:
        groups[t['kern'][0][1:]].append(t)
    return O, three, OrderedDict(sorted(groups.items(), key=lambda kv: str(kv[0])))


def merge_name(src, sig):
    labs = src.rho_labels
    pairs = [(labs[i], sig[i]) for i in range(len(labs)) if labs[i] != sig[i]]
    return ', '.join('%s\\to %s' % (a, b) for a, b in pairs) if pairs else 'none'


# ---------------------------------------------------------------- numerical check of O = S(O) + O^(2)
def check_identity(src, seed=3):
    """evaluate the ordered operator O and S(O) + O^(2) with quantum charges: rho^a(u) -> T^a acting on a
    representation per point (here: tensor product of independent random 'charge' matrices is too large);
    instead use the engine's own consistency: the classical value of O equals that of S(O) (O^(2) has fewer rhos).
    """
    import numerics as nm
    ts, f = nm.structure_constants(3)
    rng = np.random.default_rng(seed)
    pts = sorted({p for t in operator_terms(src) for _, p in t['rs']} | {p for t in operator_terms(src) for *_, p in t['Us']})
    U = {p: nm.random_adjoint(3, ts, rng, True) for p in pts}
    r = {p: rng.normal(size=8) for p in pts}

    def val(ts_):
        tot = 0.0
        for t in ts_:
            fl = [('U', i, j, p) for i, j, p in t['Us']] + [('f',) + tuple(x) for x in t['fs']] + [('rho', i, p) for i, p in t['rs']]
            tot += complex(t['c']) * 3 ** t.get('ncp', 0) * nm.eval_colour(fl, f, U, r, 8)
        return tot
    O, three, groups = decompose(src)
    a, b = val(O), val(three)
    return float(abs(a - b) / max(abs(a), 1e-300))


if __name__ == '__main__':
    res = {}
    store = OrderedDict()
    for row, key in ROWS.items():
        src = SRC[key]
        O, three, groups = decompose(src)
        chk = check_identity(src)
        res[row] = dict(source=key, ordered_terms=len(O), symmetric_terms=len(three),
                        two_rho={merge_name(src, sig): len(ts) for sig, ts in groups.items()}, classical_O_vs_S=chk)
        store[row] = dict(source=key, labels=src.rho_labels, form=src.form, cc=src.cc,
                          groups=[(sig, merge_name(src, sig), ts) for sig, ts in groups.items()])
        print('%-16s %-9s O: %2d words, S(O) %3d, two-rho: %s   classical O vs S(O): %.1e'
              % (row, key, len(O), len(three), res[row]['two_rho'], chk))
    pickle.dump(store, open(os.path.join(HERE, 'decomp.pkl'), 'wb'))
    json.dump(res, open(os.path.join(HERE, 'decomp.json'), 'w'), indent=1)
