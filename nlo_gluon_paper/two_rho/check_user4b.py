import numpy as np
from collections import defaultdict
from core import Config, accumulate
import user4
import rowcmp
S = rowcmp.user_sets()
W = []
for k in S:
    if k.startswith('4:'):
        W += S[k]
AB = sum((S[k] for k in ('A1.1', 'A1.2', 'A2.1', 'A2.2', 'B1.1', 'B1.2', 'B2.1', 'B2.2')), [])
U4, bad = user4.terms()
cfg = Config(1)
acc = defaultdict(lambda: defaultdict(complex))
accumulate(U4, cfg, acc, False, label='4rho.tex')
accumulate(W, cfg, acc, False, label='exact W')
accumulate(AB, cfg, acc, False, label='2rho(A..B)')
for k, d in acc.items():
    print(k)
    print('      ', {l: round(v.real, 4) for l, v in d.items()})
