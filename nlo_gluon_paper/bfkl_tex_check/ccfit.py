"""Group II Row I (first part) and Group III Row II: written list W against r0 = remainder of the first form and
r1 = its c.c. (conj, w <-> w'); fit W = a r0 + b r1 on two seeds, test on a third."""
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
out = {}
for name in ('G2 I first part', 'G3 II'):
    first, written, cc = CR.ROWS[name]
    F = [C.prep(C.plus(U[l])) for l in first]
    Wr = [kdot(C.prep(U[l])) for l in written]
    rows_ = []
    for s in (5, 6, 7):
        c = T3.configs(s)
        r0 = sum(RO.remainder(t, c) for t in F)
        r1 = sum(np.conj(RO.remainder(t, T3.swapped_w(c))) for t in F)
        w = sum(RO.two_value(t, c) for t in Wr)
        rows_.append((r0, r1, w))
        print(name, s, r0, r1, w, flush=True)
    A = np.array([[r[0], r[1]] for r in rows_[:2]]); bvec = np.array([r[2] for r in rows_[:2]])
    a, b = np.linalg.solve(A, bvec)
    r0, r1, w = rows_[2]
    out[name] = dict(a=str(a), b=str(b), test=str(abs(a * r0 + b * r1 - w) / abs(w)))
    print(name, out[name], flush=True)
json.dump(out, open(os.path.join(HERE, 'ccfit.json'), 'w'), indent=1)
