r"""The 'written as' two-rho lists of the three-rho rows of BFKL.tex, row by row.
Reference for row R:  B(R) = B_old(R) + J(R), because log(vee/Lambda) = log(k+/Lambda) + log(vee/k+):
  B_old(R): the two-rho (commutator) terms of the earlier, checked log(vee/k+) rows (../bfkl_final/rows.py, the
            reordering engine of ../two_rho/core.py; with F3 and the terms (B), y' -> z), units N log(vee/k+);
  J(R):     the two-rho list of the same row in Sec. 2 of Evolution_terms.tex (checked: Sec.1+Sec.2 = JIMWLK x LO).
Pointwise at random points with random classical charges, U = U^T and general U.  -> cmp_written.json"""
import json, os, sys
from collections import OrderedDict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import two_total as TT          # sets the paths (check_v4/env.py)
from core import symmetrize
sys.path.insert(0, os.path.join(HERE, '..', 'bfkl_final'))
import rows as R

EVOL = os.path.join(HERE, '..', 'evol_check', 'input', 'Evolution_terms.tex')

# user (section prefix, sub, subsub) -> (engine rows, Evolution_terms Sec. 2 (sub, subsub))
MAP = OrderedDict([
    ('G1 II', (('Group I', 'Row II', ''), ['G1 II 2a'], ('Group I', 'Row II'))),
    ('G1 III', (('Group I', 'Row III', ''), ['G1 III'], ('Group I', 'Row III'))),
    ('G1 IV', (('Group I', 'Row IV', ''), ['G1 IV'], ('Group I', 'Row IV'))),
    ('G1 V(1)', (None, ['G1 V(1)'], ('Group I', 'Row V'))),
    ('G2 I first part', (('Group II', 'Row I', 'First'), ['G2 I(1)'], ('Group II', 'Row I'))),
    ('G3 I', (('Group III', 'Row I', ''), ['G3 I'], ('Group III', 'Row I'))),
    ('G3 II', (('Group III', 'Row II', ''), ['G3 II'], ('Group III', 'Row II'))),
    ('G3 III first part', (('Group III', 'Row III', 'First'), ['G3 III(1)', 'G3 III(1) yz'], ('Group III', 'Row III'))),
    ('G3 V', (('Group III', 'Row V', ''), ['G3 V', 'G3 V yz'], ('Group III', 'Row V'))),
    ('G3 VI', (('Group III', 'Row VI', ''), ['G3 VI', 'G3 VI yz'], ('Group III', 'Row VI'))),
])


def engine(names):
    out = []
    for n in names:
        f = R.THREE.get(n) or R.THREE_YZ.get(n)
        ts = f(F3=True) if n == 'G1 III' else f()
        out += [s for t in ts for s in symmetrize(t, keep=(2,))]
    return out


def pick(E, sec, sub, subsub):
    return [t for e in E if e['kind'] == 'terms' and e['sec'].startswith(sec) and e['sub'].strip() == sub
            and e['subsub'].startswith(subsub) for t in e['terms']]


if __name__ == '__main__':
    E = TT.equations()
    import parse_v4 as P
    EV = [e for e in P.equations(EVOL) if e['sec'].startswith('Symmetrizing')]
    out = OrderedDict()
    for name, (user, eng, jrow) in MAP.items():
        J = [t for e in EV if e['kind'] == 'terms' and e['sub'].startswith(jrow[0]) and e['subsub'].strip() == jrow[1]
             for t in e['terms']]
        ref = engine(eng) + J
        if user is None:
            out[name] = dict(n_ref=len(ref), ref_is_zero=TT.both(ref, []) if ref else 'empty')
        else:
            u = pick(E, *user)
            out[name] = dict(n_user=len(u), n_J=len(J), vs_ref=TT.both(u, ref))
        print(name, out[name])
    json.dump(out, open(os.path.join(HERE, 'cmp_written.json'), 'w'), indent=1, default=float)
