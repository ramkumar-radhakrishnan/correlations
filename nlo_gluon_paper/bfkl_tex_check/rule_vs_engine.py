"""The ordering rule of reorder.py (every pair: (1/2) i f [earlier, later], the later point put on the earlier one)
applied to the engine's own three-rho terms (../bfkl_final/rows.py), against the engine's Weyl decomposition
(../two_rho/core.symmetrize, keep two-rho), pointwise in the engine's evaluator."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import two_total as TT
from core import symmetrize, term, rename, Config, pointwise_total
sys.path.insert(0, os.path.join(HERE, '..', 'bfkl_final'))
import rows as R


def rule(t):
    out = []
    rs = t['rs']
    assert len(rs) == 3
    for a in range(3):
        for b in range(a + 1, 3):
            (i, p), (j, q) = rs[a], rs[b]
            if p == q:
                continue
            d = '@%d%d' % (a, b)
            new_rs = [r for k, r in enumerate(rs) if k not in (a, b)] + [(d, p)]
            s = term(t['c'] * 0.5j, t['Us'], t['fs'] + [(i, j, d)], new_rs, t['kern'], t['tag'], t.get('ncp', 0))
            out.append(rename(s, {q: p}))
    return out


out = {}
for name in ('G3 II', 'G3 I', 'G3 V', 'G3 VI', 'G1 IV', 'G2 I(1)'):
    ts = R.THREE[name]()
    eng = [s for t in ts for s in symmetrize(t, keep=(2,))]
    mine = [s for t in ts for s in rule(t)]
    out[name] = TT.both(mine, eng)
    print(name, out[name], flush=True)
json.dump(out, open(os.path.join(HERE, 'rule_vs_engine.json'), 'w'), indent=1, default=float)
