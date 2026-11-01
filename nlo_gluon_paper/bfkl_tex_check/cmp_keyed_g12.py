r"""Group I (all rows) + Group II Row I (both parts): your 'written as' lists against the reference
B = B_old + J for the same rows, including the genuine two-rho rows (Group I Row II part 1, Row V 2a, 2b; Group II
Row I part 2): B_old from ../bfkl_final/rows.py (THREE and TWO_ROWS), J from Secs. 1 and 2 of Evolution_terms.tex.
Key by key (../two_rho/core.py, z integrated in dim. reg. where the colour does not depend on U(z)). -> cmp_keyed_g12.json"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cmp_written as W
TT = W.TT
import parse_v4 as P
from core import symmetrize
sys.path.insert(0, os.path.join(HERE, '..', 'bfkl_final'))
from two_rho_engine import keyed, residual, neg
import rows as R

E = TT.equations()
EVall = P.equations(W.EVOL)
EV = [e for e in EVall if e['sec'].startswith('Symmetrizing')]
GEN = [e for e in EVall if e['sec'].startswith('Genuine')]


def J2(sub, subsub):
    return [t for e in EV if e['kind'] == 'terms' and e['sub'].startswith(sub) and e['subsub'].strip() == subsub for t in e['terms']]


def J1(sub, subsub):
    return [t for e in GEN if e['kind'] == 'terms' and e['sub'].startswith(sub) and e['subsub'].strip() == subsub for t in e['terms']]


def two_rows(names):
    return [s for n in names for t in R.TWO_ROWS[n]() for s in symmetrize(t, keep=(2,))]


user = (W.pick(E, 'Group I', 'Row II', '') + W.pick(E, 'Group I', 'Row III', '') + W.pick(E, 'Group I', 'Row IV', '')
        + W.pick(E, 'Group II', 'Row I', 'First') + W.pick(E, 'Group II', 'Row I', 'Second'))
ref = (W.engine(['G1 II 2a', 'G1 III', 'G1 IV', 'G1 V(1)', 'G2 I(1)'])
       + two_rows(['G1 II(1)', 'G1 V(2a)', 'G1 V(2b UV)', 'G1 V(2b fin)', 'G2 I(2)'])
       + J2('Group I', 'Row II') + J2('Group I', 'Row III') + J2('Group I', 'Row IV') + J2('Group I', 'Row V')
       + J2('Group II', 'Row I') + J1('Group I', 'Row II') + J1('Group I', 'Row V') + J1('Group II', 'Row I'))
out = {}
for seed in (1, 2):
    bad = residual(keyed({'user': user, 'ref': neg(ref)}, seed))
    out[seed] = dict(n_bad=len(bad), worst=[str(b[0]) + ' %.3e / %.3e' % (b[1], b[2]) for b in sorted(bad, key=lambda b: -b[1])[:6]])
    print(seed, out[seed], flush=True)
json.dump(out, open(os.path.join(HERE, 'cmp_keyed_g12.json'), 'w'), indent=1)
