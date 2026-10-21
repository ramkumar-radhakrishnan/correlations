"""Sec. 2 of Evolution_terms.tex (two-rho terms from reordering four and three rho's) against the earlier,
checked sets of ../two_rho_all (A_sets, after F3, F4) and the earlier note TwoRho_from_final.tex."""
import json
from collections import defaultdict
import env  # noqa
from core import Config, pointwise_total, accumulate
import sets as S2
from parse_user import sec2_terms

ROWMAP = [('From $4\\, \\rho\'s$|', 'A,B'), ('Group I (with c.c.)|Row II', '3:I.II'),
          ('Group I (with c.c.)|Row III|N1', '3:I.III'), ('Group I (with c.c.)|Row III|N2', 'D'),
          ('Group I (with c.c.)|Row IV', '3:I.IV'), ('Group I (with c.c.)|Row V', '3:I.V'),
          ('Group II|Row I', '3:II.I'), ('Group III\\, (with c.c.)|Row I', '3:III.I'),
          ('Group III\\, (with c.c.)|Row II', '3:III.II'), ('Group III\\, (with c.c.)|Row III', '3:III.III'),
          ('Group III\\, (with c.c.)|Row V', '3:III.V'), ('Group III\\, (with c.c.)|Row VI', '3:III.VI')]


def rows_user(T):
    out = defaultdict(list)
    for lab, t in T:
        for pre, key in ROWMAP:
            if lab == pre or lab.startswith(pre + '|'):
                out[key].append(t)
                break
        else:
            raise KeyError(lab)
    return out


def keydiff(a, b, cfg):
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate(a, cfg, acc, False, label='user')
    accumulate([dict(t, c=-t['c']) for t in b], cfg, acc, False, label='earlier')
    bad = []
    for k, d in acc.items():
        s = sum(d.values())
        sc = sum(abs(v) for v in d.values())
        if abs(s) > 1e-9 * (1 + sc):
            bad.append((str(k), abs(s)))
    return bad


if __name__ == '__main__':
    T = sec2_terms('input/Evolution_terms.tex')
    U = rows_user(T)
    A = S2.A_sets()
    res = {}
    for _, key in ROWMAP:
        if key in res:
            continue
        diffs = []
        for seed in (3, 7, 11):
            cfg = Config(seed, True)
            a = pointwise_total(U[key], cfg)[0].real
            b = pointwise_total(A[key], cfg)[0].real
            diffs.append(abs(a - b) / (1 + abs(b)))
        nk = len(keydiff(U[key], A[key], Config(5, True)))
        res[key] = dict(n_user=len(U[key]), max_rel_pointwise=max(diffs), differing_kernel_keys=nk)
        print('%-10s n=%3d  pointwise %.2e  differing kernel monomials %d' % (key, len(U[key]), max(diffs), nk))
    json.dump(res, open('compare_sec2.json', 'w'), indent=1)
