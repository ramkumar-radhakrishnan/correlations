"""Split JIMWLK x LO and BFKL x LO by where the unobserved gluon z is attached:
   C : z attached only to charges (kernel K(a,b,z) with a,b charges)        -- present in both?
   J': z attached to the measured gluon (w or w') in JIMWLK x LO
   B': the measured gluon emitted by the gluon at z in BFKL x LO
Pointwise comparison of C_J and C_B. -> common.json"""
import json
import env  # noqa
from core import Config, pointwise_total
import target as TG
import refP

EXT = ('w', "w'")


def touches_w(t):
    for f in t['kern']:
        pts = f[1:]
        if 'z' in pts and any(p in EXT for p in pts):
            # a factor KK(a-b).KK(c-d) containing z and w/w' in the same difference or in the same factor
            a, b, c, d = pts
            if ('z' in (a, b) and any(p in EXT for p in (a, b))) or ('z' in (c, d) and any(p in EXT for p in (c, d))):
                return True
            if 'z' in (a, b) and 'z' in (c, d) and any(p in EXT for p in pts):
                return True
    return False


J, B = TG.jimwlk(), refP.bfkl()
CJ = [t for t in J if not touches_w(t)]
JW = [t for t in J if touches_w(t)]
CB = [t for t in B if not touches_w(t)]
BW = [t for t in B if touches_w(t)]
print('terms: J', len(J), 'C_J', len(CJ), "J'", len(JW), '| B', len(B), 'C_B', len(CB), "B'", len(BW))
out = []
for seed in (3, 7, 11, 19, 23):
    cfg = Config(seed, True)
    v = {k: pointwise_total(ts, cfg)[0].real for k, ts in (('C_J', CJ), ("J'", JW), ('C_B', CB), ("B'", BW))}
    out.append(dict(seed=seed, **v))
    print(seed, ' '.join('%s=%.6f' % kv for kv in v.items()))
json.dump(out, open('common.json', 'w'), indent=1)
