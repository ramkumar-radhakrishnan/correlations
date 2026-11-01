"""D = written - remainder for Group II Row I (first part) and Group III Row II: dependence on U(z) and on z."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reorder as RO
import cmp_first as C
import cmp_reorder as CR
sys.path.insert(0, os.path.join(HERE, '..', 'evol_check'))
from cmp_three import kdot
from rows_user import rows
U = {r[3]: r[4] for r in rows('BFKL.tex')}
T3 = C.T3
def D(name, c):
    first, written, cc = CR.ROWS[name]
    F = [C.prep(C.plus(U[l])) for l in first]
    Wr = [kdot(C.prep(U[l])) for l in written]
    return sum(RO.two_value(t, c) for t in Wr) - sum(RO.remainder_cc(t, c, cc) for t in F)
out = {}
c = T3.configs(5)
c2 = T3.configs(5); c2.U['z'] = T3.configs(77).U['z']
c3 = T3.configs(5); c3.pt['z'] = T3.configs(78).pt['z']
for name in ('G2 I first part', 'G3 II'):
    out[name] = [complex(D(name, x)) for x in (c, c2, c3)]
    print(name, out[name], flush=True)
s = [out['G2 I first part'][i] + out['G3 II'][i] for i in range(3)]
print('sum', s)
json.dump({k: [str(v) for v in vs] for k, vs in out.items()}, open(os.path.join(HERE, 'zdep2.json'), 'w'), indent=1)
