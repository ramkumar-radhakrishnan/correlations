"""Total of all two-rho terms proportional to log(vee/Lambda) vs JIMWLK x LO:
   Sec. 2 of Evolution_terms.tex (as written)  +  genuine terms (Sec. 1)  [+ the two-rho lists of 4rho.tex]."""
import json
from collections import defaultdict
import env  # noqa
from core import Config, pointwise_total, accumulate
import sets as S2
import target as TG
from parse_user import sec2_terms
from compare import report

GEN = ['g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III']


def nonzero_keys(sets, cfg):
    acc = defaultdict(lambda: defaultdict(complex))
    for lab, ts in sets.items():
        accumulate(ts, cfg, acc, False, label=lab)
    accumulate([dict(t, c=-t['c']) for t in TG.jimwlk()], cfg, acc, False, label='-JIMWLK')
    bad = []
    for k, d in acc.items():
        s = sum(d.values())
        sc = sum(abs(v) for v in d.values())
        if abs(s) > 1e-9 * (1 + sc):
            bad.append((k, s, {l: v for l, v in d.items() if abs(v) > 1e-12}))
    return acc, bad


if __name__ == '__main__':
    A = S2.A_sets()
    U2 = [t for _, t in sec2_terms('input/Evolution_terms.tex')]
    G = [t for k in GEN for t in A[k]]
    F = A['4rho']
    J = TG.jimwlk()
    out = {}
    for name, ts in (('sec2+genuine (as in Evolution_terms.tex)', U2 + G), ('sec2+genuine+4rho.tex lists', U2 + G + F)):
        r = []
        for seed in (3, 7, 11, 19):
            cfg = Config(seed, True)
            a, sa = pointwise_total(ts, cfg)
            j, sj = pointwise_total(J, cfg)
            r.append(abs(a - j) / sj)
        acc, bad = nonzero_keys({'all': ts}, Config(5, True))
        out[name] = dict(max_rel_pointwise=max(r), kernel_monomials=len(acc), differing=len(bad))
        print('%-45s  |sum - JIMWLK|/scale = %.2e   kernel monomials %d, differing %d' % (name, max(r), len(acc), len(bad)))
    cfg = Config(3, True)
    print('4rho.tex two-rho lists alone: value %.4f, scale %.4f, terms %d' % (pointwise_total(F, cfg)[0].real, pointwise_total(F, cfg)[1], len(F)))
    json.dump(out, open('total_check.json', 'w'), indent=1)
