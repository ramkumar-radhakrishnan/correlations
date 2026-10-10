import numpy as np
from collections import defaultdict
from core import Config, accumulate
import users as US
import user4
from compare import report
import rowcmp
S = rowcmp.user_sets()
W = []
for k in S:
    if k.startswith('4:'):
        W += S[k]
AB = S['A1.1'] + S['A1.2'] + S['A2.1'] + S['A2.2'] + S['B1.1'] + S['B1.2'] + S['B2.1'] + S['B2.2']
U4, bad = user4.terms()
for name, sets in [('4rho.tex lists  vs  exact 2rho of the 4rho terms', [(U4, 1), (W, -1)]),
                   ('4rho.tex lists  vs  exact 2rho of the 4rho terms - 2rho(A1..B2)', [(U4, 1), (W, -1), (AB, 1)]),
                   ('4rho.tex lists + 2rho(A1..B2)  vs  exact', [(U4, 1), (AB, 1), (W, -1)])]:
    for seed in (1, 2):
        cfg = Config(seed)
        acc = defaultdict(lambda: defaultdict(complex))
        for ts, s in sets:
            accumulate([dict(t, c=s * t['c']) for t in ts], cfg, acc, False, label='x')
        b = report(acc, show=False)
        tot = sum(abs(sum(d.values())) for d in acc.values())
        sc = sum(sum(abs(v) for v in d.values()) for d in acc.values())
        print('%-70s seed %d: keys %d nonzero %d  sum|res| %.3f (scale %.1f)' % (name, seed, len(acc), len(b), tot, sc))
