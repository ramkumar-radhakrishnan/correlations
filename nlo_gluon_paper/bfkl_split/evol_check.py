"""Part 1: does Sec. 1 + Sec. 2 of the new Evolution_terms.tex add up to Sec. 3, with no left-over?

Every two-rho term of the file is parsed (../check_v4/parse_v4.py; the \\begin{comment} block is removed first) and
evaluated pointwise at random transverse points with random classical charges, for U = U^T (your convention) and for a
general adjoint U.  Compared: Sec. 1 + Sec. 2, the two equations of Sec. 3 (dN_1 + dN_2) and the collected form of
Sec. 3, and JIMWLK x LO (../two_rho/target.py).  The new list of "From 4 rho's" (dN_2) is compared with the two-rho
lists of 4rho.tex (M1 of the v4 check).  Output: evol_check.json"""
import json
import os
import re
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'check_v4')))
import env  # noqa: F401,E402  (../check_v4/env.py)
from core import Config, pointwise_total  # noqa: E402
import sets as S2  # noqa: E402
import target as TG  # noqa: E402
from parse_v4 import equations  # noqa: E402

SRC = os.path.join(HERE, 'input', 'Evolution_terms.tex')
SEEDS = (3, 7, 11, 19, 23)


def rel(a, b, sym):
    r = []
    for s in SEEDS:
        cfg = Config(s, sym)
        va, sa = pointwise_total(a, cfg)
        vb, sb = pointwise_total(b, cfg)
        r.append(abs(va - vb) / max(sa + sb, 1e-300))
    return max(r)


def both(a, b):
    return dict(sym=rel(a, b, True), general=rel(a, b, False))


def neg(ts):
    return [dict(t, c=-t['c']) for t in ts]


if __name__ == '__main__':
    tex = open(SRC).read()
    tex = re.sub(r'\\begin\{comment\}.*?\\end\{comment\}', '', tex, flags=re.S)
    tmp = os.path.join(HERE, 'input', '_evol_nocomment.tex')
    open(tmp, 'w').write(tex)
    E = equations(tmp)
    os.remove(tmp)
    out = OrderedDict()
    sec = lambda p: [e for e in E if e['sec'].startswith(p) and e['kind'] == 'terms']
    G, S, J3 = sec('Genuine'), sec('Symmetrizing'), sec('JIMWLK')
    out['counts'] = dict(sec1_eqs=len(G), sec1_terms=sum(len(e['terms']) for e in G),
                         sec2_eqs=len(S), sec2_terms=sum(len(e['terms']) for e in S),
                         sec3_eqs=len(J3), sec3_terms=[len(e['terms']) for e in J3])
    sec1 = [t for e in G for t in e['terms']]
    sec2 = [t for e in S for t in e['terms']]
    four = [e for e in S if e['sub'].startswith('From')]
    out['four_rho_equations'] = [dict(name=e['name'], n=len(e['terms'])) for e in four]
    A = S2.A_sets()
    J = TG.jimwlk()
    jsum = [t for e in J3[:2] for t in e['terms']]
    jcol = J3[2]['terms']
    out['new 4rho list (dN_2) vs 4rho.tex two-rho lists (M1)'] = both(four[1]['terms'], A['4rho']) if len(four) > 1 else None
    out['Sec1 + Sec2 vs Sec3 (dN_1 + dN_2)'] = both(sec1 + sec2, jsum)
    out['Sec1 + Sec2 vs Sec3 (collected form)'] = both(sec1 + sec2, jcol)
    out['Sec3: dN_1 + dN_2 vs collected form'] = both(jsum, jcol)
    out['Sec3 collected form vs JIMWLKxLO'] = both(jcol, J)
    out['Sec1 + Sec2 vs JIMWLKxLO'] = both(sec1 + sec2, J)
    out['Sec1 + Sec2 without the new list vs JIMWLKxLO'] = both(sec1 + [t for e in S if e is not four[-1] for t in e['terms']], J)
    # size of the left-over relative to the terms themselves
    lo = []
    for s in SEEDS:
        cfg = Config(s, True)
        v, sc = pointwise_total(sec1 + sec2 + neg(jcol), cfg)
        lo.append(abs(v) / sc)
    out['left-over |Sec1+Sec2-Sec3| / sum|terms| (U=U^T)'] = max(lo)
    for k, v in out.items():
        print(k, v)
    json.dump(out, open(os.path.join(HERE, 'evol_check.json'), 'w'), indent=1)
