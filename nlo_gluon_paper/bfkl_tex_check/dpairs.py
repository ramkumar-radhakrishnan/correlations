"""D(row) = written - remainder at one configuration, for every three-rho row; which ones cancel in pairs?"""
import os, sys, json
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
for seed in (5, 6):
    c = T3.configs(seed)
    for name, (first, written, cc) in CR.ROWS.items():
        F = [C.prep(C.plus(U[l])) for l in first] if name != 'G1 II' else [C.prep(C.ours('G1II2a'))]
        Wr = [kdot(C.prep(U[l])) for l in written]
        d = sum(RO.two_value(t, c) for t in Wr) - sum(RO.remainder_cc(t, c, cc) for t in F)
        out.setdefault(name, {})[seed] = [d.real, d.imag]
        print(seed, name, d, flush=True)
json.dump(out, open(os.path.join(HERE, 'dpairs.json'), 'w'), indent=1)
