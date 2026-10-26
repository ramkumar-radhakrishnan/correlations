"""Two-rho terms of bfkl.tex (genuine two-rho rows + the commutator terms of the three-rho rows, all c.c.'s) against
B - J, B = Sec. 4 (BFKL contribution) and J = Sec. 3 (JIMWLK contribution) of Evolution_terms.tex (parsed).

Comparison key by key (../two_rho/core.py): every kernel monomial at random points, colour numerically.  Terms whose
colour does not depend on U(z) are integrated over z in dim. reg. (int_z K(a,b,z) = pi[P - log(a-b)^2], int_z K(a,a,z) = 0,
the MS-bar convention of Part2b_MSbar.pdf); the rows with 1/eps (ELL factors) enter in the same form."""
import json
import os
import re
from collections import defaultdict
import numpy as np
import env  # noqa: F401
from core import symmetrize, accumulate_N, accumulate_dr, config_pair, pointwise_total
import rows as R

HERE = os.path.dirname(os.path.abspath(__file__))
EVOL = os.path.normpath(os.path.join(HERE, '..', 'bfkl_split', 'input', 'Evolution_terms.tex'))


def targets():
    """J: Sec. 3 of Evolution_terms.tex (collected form).  B: BFKL x LO of App. E (../bfkl/refP.py, units of N),
    and the (commented-out) Sec. 4 of Evolution_terms.tex for comparison."""
    import sys
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'check_v4')))
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'bfkl')))
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'bfkl', 'reference')))
    from parse_v4 import equations
    E = equations(EVOL)
    J = [e for e in E if e['sec'].startswith('JIMWLK') and e['kind'] == 'terms']
    B4 = [e for e in E if e['sec'].startswith('BFKL') and e['kind'] == 'terms']
    import refP
    return J, B4, refP.bfkl()


def file_terms(F3=True):
    out = {}
    for name, f in R.THREE.items():
        ts = f(F3=F3) if name == 'G1 III' else f()
        out[name] = [s for t in ts for s in symmetrize(t, keep=(2,))]
    for name, f in R.TWO_ROWS.items():
        out[name] = [s for t in f() for s in symmetrize(t, keep=(2,))]
    return out


def has_ell(t):
    return any(k[0] == 'ELL' for k in t['kern'])


def keyed(groups, seed):
    c1, c2 = config_pair(seed, True)
    acc = defaultdict(lambda: defaultdict(complex))
    for lab, ts in groups.items():
        pt = [t for t in ts if not has_ell(t)]
        dr = [t for t in ts if has_ell(t)]
        accumulate_N(pt, c1, c2, acc, label=lab)
        accumulate_dr(dr, c1, acc, label=lab)
    return acc


def residual(acc, tol=1e-9):
    bad, scale_all = [], 0.0
    for key, d in acc.items():
        tot = sum(d.values())
        sc = sum(abs(v) for v in d.values())
        scale_all += sc
        if key[0] == 'ell0':
            continue
        if abs(tot) > tol * max(1.0, sc):
            bad.append((key, abs(tot), sc, {k: complex(v) for k, v in d.items() if abs(v) > 1e-12}))
    return bad


def neg(ts):
    return [dict(t, c=-t['c']) for t in ts]


if __name__ == '__main__':
    J, B4, Bref = targets()
    print('Evolution file: JIMWLK eqs', [len(e['terms']) for e in J], ' BFKL eqs (commented out)', [len(e['terms']) for e in B4])
    Jc, Bc = J[-1]['terms'], Bref
    out = {}
    import target as TG
    for lab, a, b in (('Sec. 3 (collected) vs JIMWLKxLO', Jc, TG.jimwlk()),
                      ('Sec. 4 eq. 1 vs BFKLxLO (App. E)', B4[0]['terms'], Bref),
                      ('Sec. 4 eq. 2 vs BFKLxLO (App. E)', B4[1]['terms'], Bref)):
        r = []
        for seed in (1, 2, 3):
            acc = keyed({'a': a, 'b': neg(b)}, seed)
            bad = residual(acc)
            r.append(len(bad))
        out[lab] = r
        print('%-36s residual keys: %s' % (lab, r))
    for label, F3 in (('with F3', True), ('as written', False)):
        FT = file_terms(F3)
        groups = dict(FT)
        groups['-(B - J)'] = neg(Bc) + Jc
        res = []
        for seed in (1, 2, 3):
            acc = keyed(groups, seed)
            bad = residual(acc)
            res.append(dict(n_keys=len(acc), n_bad=len(bad),
                            worst=max([b[1] / max(1.0, b[2]) for b in bad], default=0.0)))
            if seed == 1 and bad:
                print(label, ': first residual keys')
                for k, a, s, d in sorted(bad, key=lambda b: -b[1])[:12]:
                    print('   %.3e  %s' % (a, k))
                    print('        ', {l: np.round(v, 4) for l, v in d.items()})
        out[label] = res
        print('%-10s' % label, res)
    json.dump(out, open('two_rho_engine.json', 'w'), indent=1)
