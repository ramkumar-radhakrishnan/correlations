"""Check of the new Evolution_terms.tex, equation by equation.

Sec. 1 (genuine) and Sec. 2 (from three and four rho's): against the earlier checked encodings
   (final_check: genuine terms, the old Sec. 2, and the two-rho lists of 4rho.tex).
Sum of Secs. 1, 2: against JIMWLK x LO (operator target).
Sec. 3: dN1, dN2, their sum, the collected form, against JIMWLK x LO.
Sec. 4: the two forms against each other and against BFKL x LO.
Pointwise comparisons at random points with random classical charges (no projection), for U = U^T (your convention)
and for a general adjoint U.  Output: check_all.json"""
import importlib.util
import json
import os
import sys
from collections import OrderedDict, defaultdict
import env  # noqa: F401
from core import Config, pointwise_total, conj_swap
import sets as S2
import target as TG
import refP
from parse_v4 import equations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'final_check'))
from parse_user import sec2_terms  # noqa: E402
import genuine_new as GN  # noqa: E402

OLD = os.path.join(HERE, '..', 'final_check', 'input', 'Evolution_terms.tex')
SEEDS = (3, 7, 11)


def rel(a, b, sym, seeds=SEEDS):
    """max over seeds of |a-b| / (scale of b + scale of a)"""
    r = []
    for s in seeds:
        cfg = Config(s, sym)
        va, sa = pointwise_total(a, cfg)
        vb, sb = pointwise_total(b, cfg)
        r.append(abs(va - vb) / max(sa + sb, 1e-300))
    return max(r)


def both(a, b):
    return dict(sym=rel(a, b, True), general=rel(a, b, False))


def neg(ts):
    return [dict(t, c=-t['c']) for t in ts]


def with_cc(ts):
    return ts + [conj_swap(t) for t in ts]


def by_tag(ts, pre):
    return [t for t in ts if t['tag'].startswith(pre)]


if __name__ == '__main__':
    E = equations()
    A = S2.A_sets()
    out = OrderedDict()

    def eq(sec, sub=None, subsub=None, k=0):
        L = [e for e in E if e['sec'].startswith(sec) and (sub is None or e['sub'].startswith(sub))
             and (subsub is None or e['subsub'].startswith(subsub))]
        return L[k]['terms']

    # ---------------------------------------------------------------- Sec. 1: genuine
    g = OrderedDict()
    g['Group I, Row II'] = both(eq('Genuine', 'Group I', 'Row II'), A['g:I.II'])
    v = GN.G1_V()
    v2a = by_tag(v, 'n:G1.V.2a')
    v2bf = by_tag(v, 'n:G1.V.2b.fin')
    vuv = by_tag(v, 'n:G1.V.2b.UV')
    vct = by_tag(v, 'n:G1.V.2b.CT')
    e2a, euv, efin, ect = (eq('Genuine', 'Group I', 'Row V', k) for k in range(4))
    g['Group I, Row V, 2a'] = both(e2a, v2a)
    g['Group I, Row V, 2b UV'] = both(euv, vuv)
    g['Group I, Row V, 2b finite (first)'] = both(efin, v2bf)
    g['Group I, Row V, 2b finite (second) = counterterm'] = both(ect, vct)
    g['Group I, Row V, 2b UV + second "2b finite" (should vanish)'] = dict(
        sym=max(abs(pointwise_total(euv + ect, Config(s, True))[0]) / pointwise_total(euv, Config(s, True))[1] for s in SEEDS),
        general=max(abs(pointwise_total(euv + ect, Config(s, False))[0]) / pointwise_total(euv, Config(s, False))[1] for s in SEEDS))
    g['Group I, Row V, all four'] = both(e2a + euv + efin + ect, A['g:I.V'])
    g['Group II, Row I'] = both(eq('Genuine', 'Group II', 'Row I'), A['g:II.I'])
    g['Group II, Row IV'] = both(eq('Genuine', 'Group II', 'Row IV'), A['g:II.IV'])
    g3 = eq('Genuine', 'Group III', 'Second')
    g['Group III, Row III (second part), as written'] = both(g3, A['g:III.III'])
    g['Group III, Row III (second part), written + c.c.'] = both(with_cc(g3), A['g:III.III'])
    out['genuine'] = g
    GEN_NEW = [t for e in E if e['sec'].startswith('Genuine') for t in e['terms']]
    GEN_OLD = [t for k in ('g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III') for t in A[k]]
    out['genuine, all'] = both(GEN_NEW, GEN_OLD)

    # ---------------------------------------------------------------- Sec. 2 vs the old Sec. 2
    old = defaultdict(list)
    for lab, t in sec2_terms(OLD):
        old[lab].append(t)
    s2 = OrderedDict()
    new2 = [e for e in E if e['sec'].startswith('Symmetrizing')]
    for e in new2:
        part = 'N1' if e['name'].startswith('1') else ('N2' if e['name'].startswith('2') else '')
        lab = '%s|%s|%s' % (e['sub'], e['subsub'], part)
        key = lab if lab in old else None
        if key is None:
            raise KeyError(lab)
        r = both(e['terms'], old[key])
        r.update(n_new=len(e['terms']), n_old=len(old[key]))
        s2['%s / %s / %s' % (e['sub'], e['subsub'], part)] = r
    out['sec2_vs_old'] = s2
    SEC2_NEW = [t for e in new2 for t in e['terms']]
    four_new = [t for e in new2 if e['sub'].startswith('From') for t in e['terms']]
    out['From 4 rho: new vs old Sec.2 + 4rho.tex lists (M1)'] = both(four_new, old["From $4\\, \\rho's$||"] + A['4rho'])

    # ---------------------------------------------------------------- totals vs JIMWLK
    J = TG.jimwlk()
    tot = GEN_NEW + SEC2_NEW
    out['Secs 1+2 vs JIMWLKxLO'] = both(tot, J)
    out['Secs 1+2 + 4rho.tex lists (M1) vs JIMWLKxLO'] = both(tot + A['4rho'], J)

    # ---------------------------------------------------------------- Sec. 3
    j1, j2, jc = (eq('JIMWLK', k=k) for k in range(3))
    s3 = OrderedDict()
    s3['dN1 + dN2 vs collected form'] = both(j1 + j2, jc)
    s3['collected form vs JIMWLKxLO'] = both(jc, J)
    s3['dN1 + dN2 vs JIMWLKxLO'] = both(j1 + j2, J)
    s3['collected form vs Secs 1+2'] = both(jc, tot)
    out['sec3'] = s3

    # ---------------------------------------------------------------- Sec. 4
    B = refP.bfkl()
    b1, b2 = (eq('BFKL', k=k) for k in range(2))
    s4 = OrderedDict()
    s4['first form vs collected form'] = both(b1, b2)
    s4['first form vs BFKLxLO'] = both(b1, B)
    s4['collected form vs BFKLxLO'] = both(b2, B)
    out['sec4'] = s4

    def show(d, ind=''):
        for k, v in d.items():
            if isinstance(v, dict) and 'sym' in v:
                extra = '  (n %d / %d)' % (v['n_new'], v['n_old']) if 'n_new' in v else ''
                print('%s%-62s U=U^T %.1e   general U %.1e%s' % (ind, k[:62], v['sym'], v['general'], extra))
            else:
                print(ind + k)
                show(v, ind + '   ')
    show(out)
    json.dump(out, open(os.path.join(HERE, 'check_all.json'), 'w'), indent=1)
