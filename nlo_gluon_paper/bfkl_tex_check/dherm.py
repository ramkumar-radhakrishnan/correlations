"""D(row) = written - commutator terms of the first form for Group II Row I and Group III Row II:
D(c) against conj(D(c with w <-> w')) (hermitian or anti-hermitian?)"""
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
for name in ('G2 I first part', 'G3 II'):
    a = D(name, c)
    b = np.conj(D(name, T3.swapped_w(c)))
    out[name] = dict(D=str(a), Dcc=str(b), herm=abs(a - b) / abs(a), antiherm=abs(a + b) / abs(a))
    print(name, out[name], flush=True)
json.dump(out, open(os.path.join(HERE, 'dherm.json'), 'w'), indent=1)
