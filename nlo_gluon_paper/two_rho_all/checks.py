"""Pointwise checks (no z integration): sum of the log(vee/Lambda) terms = JIMWLK x LO, sum of the log(vee) terms
= BFKL x LO, row by row values of A and P, and the comparison JIMWLK x LO vs BFKL x LO.  -> checks.json"""
import json
import env  # noqa
from core import Config, pointwise_total
import target as TG
import refP
import sets as S

A, P = S.A_sets(), S.P_sets()
J, B = TG.jimwlk(), refP.bfkl()
out = dict(points=[], rows={})
for seed in (3, 7, 11, 19):
    for sym in (True, False):
        cfg = Config(seed, sym)
        sa = sum(pointwise_total(ts, cfg)[0] for ts in A.values())
        sp_ = sum(pointwise_total(ts, cfg)[0] for ts in P.values())
        j, b = pointwise_total(J, cfg)[0], pointwise_total(B, cfg)[0]
        scale = max(abs(j), abs(b), 1e-300)
        out['points'].append(dict(seed=seed, sym=sym, sumA=sa.real, J=j.real, sumP=sp_.real, B=b.real,
                                  dA=abs(sa - j) / scale, dP=abs(sp_ - b) / scale))
        print(seed, sym, 'sumA %.6f J %.6f | sumP %.6f B %.6f' % (sa.real, j.real, sp_.real, b.real))
        if seed == 3 and sym:
            for k in S.ORDER:
                a = pointwise_total(A.get(k, []), cfg)[0].real if k in A else 0.0
                p = pointwise_total(P.get(k, []), cfg)[0].real if k in P else 0.0
                out['rows'][k] = dict(A=a, P=p)
json.dump(out, open('checks.json', 'w'), indent=1)
sym = [p for p in out["points"] if p["sym"]]
print("symmetric U: max rel. differences", max(p["dA"] for p in sym), max(p["dP"] for p in sym))
