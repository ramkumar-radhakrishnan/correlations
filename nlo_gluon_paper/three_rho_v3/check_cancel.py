"""Do the symmetrised three-rho terms cancel?  Pointwise numerical test.
Every term is reduced to its fully symmetrised (classical) charge part; the three charges are put on three points
p1, p2, p3 in all 3! ways; terms are collected by kernel monomial (K.K products in p1, p2, p3, w, w', z) and the
colour factors are evaluated with random SU(3) adjoint Wilson lines and random charges.
symmetric=True: U = U^T (your convention U^{bc} = U^{cb}); symmetric=False: a general adjoint U.
Output: per class and in total, the number of kernel monomials and the largest residual."""
import json
import sys
from collections import defaultdict
import env  # noqa
from core import symmetrize
import acc3 as A3
import parse3 as P


def run(entries, seed, symmetric, by_label=True):
    cfg = A3.Config3(seed, symmetric)
    acc = defaultdict(lambda: defaultdict(complex))
    for e in entries:
        A3.accumulate3(symmetrize(P.as_term(e), (3,)), cfg, acc, label='%s.%d' % (e['block'], e['n']))
    out = []
    for k, d in acc.items():
        tot = sum(d.values())
        sc = max(abs(v) for v in d.values())
        out.append((k, abs(tot), sc, d))
    return out


def summary(entries, symmetric, seeds=(5, 17, 23)):
    res = []
    for s in seeds:
        r = run(entries, s, symmetric)
        bad = [x for x in r if x[1] > 1e-10 * max(1.0, x[2])]
        res.append(dict(seed=s, monomials=len(r), nonzero=len(bad), max_residual=max([x[1] for x in r] + [0.0]),
                        max_single=max(x[2] for x in r)))
    return res


if __name__ == '__main__':
    E = P.parse()
    R = {}
    for sym in (True, False):
        name = 'U symmetric (U^{bc}=U^{cb})' if sym else 'general U'
        R[name] = {'all terms': summary(E, sym)}
        for cls in ('alpha', "alpha'", 'beta', "beta'"):
            R[name][cls] = summary([e for e in E if e['cls'] == cls], sym)
        for k, v in R[name].items():
            print('%-28s %-9s monomials %3d  non-cancelling %s  largest residual %.2e  (largest single term %.2e)' % (
                name, k, v[0]['monomials'], [x['nonzero'] for x in v], max(x['max_residual'] for x in v), v[0]['max_single']))
    json.dump(R, open('cancel_check.json', 'w'), indent=1)
