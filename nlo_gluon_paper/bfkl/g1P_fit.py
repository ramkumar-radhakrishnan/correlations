import env  # noqa
import numpy as np
import sympy as sp
from collections import defaultdict
from core import Config, accumulate
import rowcmp
import userP as UP
import users as US
import refP
from sources import SOURCES
S = rowcmp.user_sets()
R = refP.rows()
ref = sum((R[k] for k in R if k.startswith('G1-')), [])
g1_4 = sum((S[k] for k in S if k.startswith('4:G1:')), [])
IV = UP.p_terms('GI.IV.P')
src = UP.PSRC['GI.III.P']
s = UP.mk('GI.III.P', 'GI.III', sp.I, 'K', [(-1, src.kexpl[0][1]), (2, src.kexpl[1][1])])
III = US.source_terms(s)
base = III + IV + g1_4
GII = [x for x in SOURCES if x.key == 'GI.II'][0]
cands = {
    'alpha_b_s13': US.source_terms(GII, monos=(1,), subs=(1, 3)),
    'alpha_d_s13': US.source_terms(GII, monos=(3,), subs=(1, 3)),
    'alpha_b_s24': US.source_terms(GII, monos=(1,), subs=(2, 4)),
    'alpha_d_s24': US.source_terms(GII, monos=(3,), subs=(2, 4)),
}
def vec(ts, seeds):
    out = {}
    for sd in seeds:
        cfg = Config(sd)
        acc = defaultdict(lambda: defaultdict(complex))
        accumulate(ts, cfg, acc, False, label='x')
        for k, d in acc.items():
            out[(sd, k)] = sum(d.values())
    return out
seeds = (5, 7)
r = vec(base + [dict(t, c=-t['c']) for t in ref], seeds)
V = {n: vec(ts, seeds) for n, ts in cands.items()}
keys = sorted(set(r) | set().union(*[set(v) for v in V.values()]), key=str)
A = np.array([[V[n].get(k, 0) for n in cands] for k in keys])
b = -np.array([r.get(k, 0) for k in keys])
sol, res, rk, sv = np.linalg.lstsq(A, b, rcond=None)
print('coefficients:', dict(zip(cands, np.round(sol, 6))))
print('residual after fit:', np.abs(A @ sol - b).max(), ' before:', np.abs(b).max(), 'rank', rk)
