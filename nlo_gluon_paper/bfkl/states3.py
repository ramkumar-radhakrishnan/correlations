"""three-rho (symmetrized) log(vee/k+) terms: number of left-over structures, as written and after the fixes"""
import env  # noqa
import json
from collections import defaultdict
import users as US
import lsets
from sources import SOURCES
import cancel_data as CD
from acc3 import Config3, accumulate3, report3
SRC = {s.key: s for s in SOURCES}
L3w, L3c = lsets.sets((3,), False), lsets.sets((3,), True)


def total3(L, seed=7, extra=()):
    cfg = Config3(seed)
    acc = defaultdict(lambda: defaultdict(complex))
    for k, ts in L.items():
        accumulate3(ts, cfg, acc, label=k)
    for k, ts in extra:
        accumulate3(ts, cfg, acc, label=k)
    bad = report3(acc, show=False)
    return [len(acc), len(bad), max(abs(sum(d.values())) for d in acc.values())]


st3 = [('as written', total3(L3w)),
       ('+ P1', total3(L3w, extra=[('kh', lsets.fixKh((3,)))])),
       ('+ F3', total3(L3c)),
       ('+ F5', total3(L3c))]
print(st3)
json.dump(st3, open('states3.json', 'w'), default=float)
