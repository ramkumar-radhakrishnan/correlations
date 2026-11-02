"""Sum of every two-rho list of BFKL.tex (written forms, 4-rho lists) against the 'BFKL contribution' (both forms),
in the texeval evaluator at the configurations of dpairs.py; with the per-row differences D this gives
sum of the commutator terms of the first forms - BFKL contribution = -sum D."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reorder as RO
import cmp_first as C
sys.path.insert(0, os.path.join(HERE, '..', 'evol_check'))
from cmp_three import kdot
from rows_user import rows
T3 = C.T3
R = rows('BFKL.tex')
two_lists = [kdot(C.prep(e[4])) for e in R if r'\mathcal{K}' in e[4] and r'\Bigg[' in e[4] and not e[0].startswith('BFKL')
             and r'\rho^{h}(y)' not in e[4] and 'rho^{a}(y)' not in e[4]]
final = [kdot(C.prep(e[4])) for e in R if e[0].startswith('BFKL')][:2]
D = json.load(open(os.path.join(HERE, 'dpairs.json')))
out = {}
for s in (5, 6):
    c = T3.configs(s)
    tw = sum(RO.two_value(t, c) for t in two_lists)
    fb = [RO.two_value(t, c) for t in final]
    sD = sum(complex(*v[str(s)]) for v in D.values())
    out[s] = dict(n_lists=len(two_lists), written_total=str(tw), bfkl=[str(x) for x in fb],
                  rel_written_vs_bfkl=abs(tw - fb[1]) / abs(fb[1]), sum_D=str(sD), rel_sumD=abs(sD) / abs(fb[1]))
    print(s, out[s], flush=True)
json.dump(out, open(os.path.join(HERE, 'total_tex.json'), 'w'), indent=1)
