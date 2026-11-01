"""reorder.py on the J (log k+/Lambda) three-rho terms of Group II Row I and Group III Row II against the checked
Sec. 2 lists of Evolution_terms.tex (validation of the tool on these rows)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reorder as RO
import cmp_first as C
sys.path.insert(0, os.path.join(HERE, '..', 'evol_check'))
from cmp_three import kdot
from rows_user import rows
T3 = C.T3
J = {lab: t for lab, t, cc in T3.row_terms('J', True)}
E = rows(os.path.join(HERE, '..', 'evol_check', 'input', 'Evolution_terms.tex'))
def sec2(group, row):
    return [kdot(C.prep(e[4])) for e in E if e[0].startswith('Symmetrizing') and e[1].split('\\')[0].strip() == group and (e[2] or '').strip() == row]
for name, labs, g, r, cc in (('G2 I', ['G2 I(1) A #2', 'G2 I(1) B #2'], 'Group II', 'Row I', False),
                             ('G3 II', ['G3 II #2'], 'Group III', 'Row II', True)):
    W = sec2(g, r)
    for s in (1, 2):
        c = T3.configs(s)
        a = sum(RO.remainder_cc(J[l], c, cc) for l in labs)
        b = sum(RO.two_value(t, c) for t in W)
        print(name, len(W), b / a, flush=True)
