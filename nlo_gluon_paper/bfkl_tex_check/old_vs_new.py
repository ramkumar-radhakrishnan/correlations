"""Group III Row II (and Group II Row I first part): our log(vee/k+) first form (BFKL.tex = Checks list) against the
earlier checked list B_old (../bfkl_final/input/bfkl.tex) + the J terms of the same row:
 (1) as classical three-rho integrands (symmetrized over the charges, c.c. as in the row);
 (2) their commutator (two-rho) terms against your 'written as' list."""
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
T3 = C.T3
U = {r[3]: r[4] for r in rows('BFKL.tex')}
OLD = {(r[0].split('\\')[0].strip(), r[1], r[2]): r[4] for r in rows(os.path.join(HERE, '..', 'bfkl_final', 'input', 'bfkl.tex'))}
J = {lab: t for lab, t, cc in T3.row_terms('J', True)}
SETS = {'G3 II': (('Group III', 'Row II', None), ['G3 II #2'], True),
        'G2 I first part': (('Group II', 'Row I', 'First part'), ['G2 I(1) A #2', 'G2 I(1) B #2'], False)}
out = {}
for name, (okey, jl, cc) in SETS.items():
    first, written, _ = CR.ROWS[name]
    F = [C.prep(C.plus(U[l])) for l in first]
    Bold = [C.prep(OLD[okey])] + [J[l] for l in jl]
    Wr = [kdot(C.prep(U[l])) for l in written]
    res = {}
    for s in (5, 6):
        c = T3.configs(s)
        a = sum(T3.sym_value(t, c, cc) for t in F)
        b = sum(T3.sym_value(t, c, cc) for t in Bold)
        ra = sum(RO.remainder_cc(t, c, cc) for t in F)
        rb = sum(RO.remainder_cc(t, c, cc) for t in Bold)
        w = sum(RO.two_value(t, c) for t in Wr)
        res[s] = dict(three_rho_rel=abs(a - b) / abs(b), rem_new_vs_written=abs(ra - w) / abs(w),
                      rem_old_vs_written=abs(rb - w) / abs(w))
        print(name, s, res[s], flush=True)
    out[name] = res
json.dump(out, open(os.path.join(HERE, 'old_vs_new.json'), 'w'), indent=1, default=float)
