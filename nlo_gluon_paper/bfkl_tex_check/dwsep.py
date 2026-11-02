"""Is D = written - commutator terms (Group II Row I first part, Group III Row II, and their sum) of the form
f(no w) + g(no w') (after the phase e^{-ik(w'-w)} is taken out)?  Then int_w e^{ik.w} (or int_w') kills it for
k != 0.  Test: mixed second difference E00 - E10 - E01 + E11, with w (point and U) and w' moved."""
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
def moved(c, which):
    new = T3.permuted(c, T3.CH)
    o = T3.configs(91)
    for n in which:
        new.pt[n], new.U[n], new.rho[n] = o.pt[n], o.U[n], o.rho[n]
    return new
def E(name, c):
    return D(name, c) / np.exp(-1j * (c.kv @ (c.pt['wp'] - c.pt['w'])))
out = {}
c = T3.configs(5)
cs = {'00': c, '10': moved(c, ['w']), '01': moved(c, ['wp']), '11': moved(c, ['w', 'wp'])}
tot = {k: 0 for k in cs}
for name in ('G2 I first part', 'G3 II'):
    e = {k: E(name, v) for k, v in cs.items()}
    for k in cs:
        tot[k] += e[k]
    mixed = e['00'] - e['10'] - e['01'] + e['11']
    out[name] = dict(mixed_rel=abs(mixed) / abs(e['00']), E={k: str(v) for k, v in e.items()})
    print(name, out[name], flush=True)
mixed = tot['00'] - tot['10'] - tot['01'] + tot['11']
out['sum'] = dict(mixed_rel=abs(mixed) / abs(tot['00']))
print('sum', out['sum'])
json.dump(out, open(os.path.join(HERE, 'dwsep.json'), 'w'), indent=1)
