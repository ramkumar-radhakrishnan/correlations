"""Every term of the corrected list in operator order, with the user's own colour indices.
For each term the matching earlier entry gives its source, sub-term, kernel monomial and whether it is the c.c.;
the source fixes the operator form (ABC, A{B,C} or {B,C}A), and the c.c. reverses it."""
import numpy as np
import env
from core import term
import cancel_data as CD
from sources import SOURCES
import fixes as FX
import parse as P
from cmp_entries import fp_user, fp_entry

SRC = {s.key: s for s in SOURCES}
SRC['D'] = CD.D_SOURCE
KEYS = {'4rho.a': ['A1', 'A2'], '4rho.b': ['B1', 'B2'], 'I.II.1': ['GI.II'], 'I.II.2': ['GI.II'], 'I.III': ['GI.III'],
        'D': ['D'], 'I.IV': ['GI.IV'], 'I.V.1': ['GI.V1', 'GI.V2'], 'I.V.2': ['GI.V1', 'GI.V2'],
        'II.I': ['GII.I1', 'GII.I2'], 'III.I': ['GIII.I'], 'III.II': ['GIII.II'], 'III.III': ['GIII.III'],
        'III.V': ['GIII.V'], 'III.VI': ['GIII.VI']}
_CACHE = {}


def entries_with_maps(key):
    if key not in _CACHE:
        ents, ren = CD.source_entries(SRC[key])
        _CACHE[key] = [(en, ren[(en['cc'], en['mono'])]['full'], SRC[key]) for en in ents]
    return _CACHE[key]


def operator_form(e):
    """returns dict(form, order, coeff_factor, src, sub, cc): order = positions in operator order;
    form 'ABC' (one ordering, coefficient as written) or 'A{BC}' / '{BC}A' (coefficient halved)."""
    f = fp_user(e)
    for key in KEYS[e['block']]:
        for en, full, src in entries_with_maps(key):
            if en['cls'] != e['cls'] or not np.allclose(f, fp_entry(en), atol=1e-10):
                continue
            word = [full[p] for _, p in src.subterms[en['sub'] - 1][2]]
            if e['block'] in ('I.II.1', 'I.II.2'):
                # f^{abc} antisymmetry: the two written orderings of a pair are exactly rho {rho, rho};
                # keep them, and reverse each one for the c.c.
                wr = [p for _, p in e['rs']]
                order = wr[::-1] if en['cc'] else wr
                return dict(form='ABC', order=order, src=key, sub=en['sub'], mono=en['mono'], cc=en['cc'], label=en['label'])
            if src.form == 'plain':
                order, form = (word[::-1], 'ABC') if en['cc'] else (word, 'ABC')
            elif src.form == 'A{BC}':
                A, B, C = word
                order, form = ([B, C, A], '{BC}A') if en['cc'] else ([A, B, C], 'A{BC}')
            else:
                B, C, A = word
                order, form = ([A, B, C], 'A{BC}') if en['cc'] else ([B, C, A], '{BC}A')
            return dict(form=form, order=order, src=key, sub=en['sub'], mono=en['mono'], cc=en['cc'], label=en['label'])
    raise ValueError('no match %s %d' % (e['block'], e['n']))


def op_terms(e, of=None):
    """ordered terms (core format) of the operator form of entry e"""
    of = of or operator_form(e)
    byp = {p: (i, p) for i, p in e['rs']}
    c = e['pref'] * e['sign'] * e['mult']
    kern = P.CLASS_KERN[e['cls']]
    o = [byp[p] for p in of['order']]
    if of['form'] == 'ABC':
        words = [(c, o)]
    elif of['form'] == 'A{BC}':
        words = [(c / 2, o), (c / 2, [o[0], o[2], o[1]])]
    else:
        words = [(c / 2, o), (c / 2, [o[1], o[0], o[2]])]
    return [term(cc, e['Us'], e['fs'], w, kern, '%s.%d' % (e['block'], e['n'])) for cc, w in words]


if __name__ == '__main__':
    import sys
    from collections import OrderedDict, Counter
    from core import symmetrize, Config, pointwise_total
    sys.path.insert(0, env.ALL2)
    import sets as S2
    from two_cmp import ROW
    E = FX.corrected()
    stat = Counter()
    O = OrderedDict()
    for e in E:
        of = operator_form(e)
        same = [p for _, p in e['rs']] == of['order']
        stat[(e['block'], of['form'], of['cc'], 'order as written' if same else 'order changed')] += 1
        for t in op_terms(e, of):
            O.setdefault(ROW[e['block']], []).extend(symmetrize(t, (2,)))
    for k, v in sorted(stat.items()):
        print(k, v)
    A = S2.A_sets()
    for seed in (3, 7):
        cfg = Config(seed, True)
        for k in O:
            a = pointwise_total(O[k], cfg)[0].real
            b = pointwise_total(A[k], cfg)[0].real
            print('seed %d %-10s operator-ordered list %10.5f   earlier %10.5f   %s' % (seed, k, a, b, 'same' if abs(a - b) < 1e-9 * (1 + abs(b)) else 'DIFFERENT'))
