"""D = written - remainder for Group III Row I and Group I Row IV: does it depend on U(z)?"""
import os, sys
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
for name in ('G3 I', 'G1 IV'):
    first, written, cc = CR.ROWS[name]
    F = [C.prep(C.plus(U[l])) for l in first]
    Wr = [kdot(C.prep(U[l])) for l in written]
    c = T3.configs(5)
    def D(c):
        return sum(RO.two_value(t, c) for t in Wr) - sum(RO.remainder_cc(t, c, cc) for t in F)
    d1 = D(c)
    c2 = T3.configs(5)
    c2.U['z'] = T3.configs(77).U['z']
    d2 = D(c2)
    print(name, d1, d2, abs(d1 - d2) / abs(d1), flush=True)
