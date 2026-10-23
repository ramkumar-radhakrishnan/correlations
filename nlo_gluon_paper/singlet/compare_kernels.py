"""Sum of all rows vs JIMWLK x LO, kernel structure by kernel structure (colour-singlet form, x <-> y symmetrized),
for U = U^T (your convention) and for a general U.  |rows - JIMWLK| relative to the sum of |single monomials| of the
rows in that structure.  Output: compare_kernels.json"""
import json
import os
import pickle
import numpy as np
import env  # noqa: F401
from core import Config
from render_singlet import _numeric

HERE = os.path.dirname(os.path.abspath(__file__))
D = pickle.load(open(os.path.join(HERE, 'singlet.pkl'), 'rb'))
R, J = D['total_written'], D['J_gen']


def val(coll, kk, cfg):
    sub = {kk: coll[kk]} if kk in coll else {}
    if not sub:
        return 0.0
    return 0.5 * (_numeric(sub, cfg, {'x': 'p1', 'y': 'p2'}) + _numeric(sub, cfg, {'x': 'p2', 'y': 'p1'}))


out = []
keys = list(R.keys()) + [k for k in J if k not in R]
for kk in keys:
    row = dict(kernel=str(kk), in_rows=kk in R, in_JIMWLK=kk in J)
    for sym in (True, False):
        rel = []
        for seed in (1, 2, 3):
            cfg = Config(seed, sym)
            a, b = val(R, kk, cfg), val(J, kk, cfg)
            sc = sum(abs(val({k: {ws: c}}, k, cfg)) for k in (kk,) for ws, c in R.get(kk, {}).items()) or 1.0
            rel.append(abs(a - b) / sc)
        row['sym' if sym else 'general'] = max(rel)
    out.append(row)
    print('%-95s rows:%d J:%d   U=U^T: %.1e   general U: %.1e' % (row['kernel'][:95], row['in_rows'], row['in_JIMWLK'],
                                                                   row['sym'], row['general']))
tot = {}
for sym in (True, False):
    rel = []
    for seed in (1, 2, 3):
        cfg = Config(seed, sym)
        a = sum(val(R, kk, cfg) for kk in R)
        b = sum(val(J, kk, cfg) for kk in J)
        rel.append(abs(a - b) / abs(b))
    tot['sym' if sym else 'general'] = max(rel)
print('total:', tot)
json.dump(dict(per_kernel=out, total=tot), open(os.path.join(HERE, 'compare_kernels.json'), 'w'), indent=1)
