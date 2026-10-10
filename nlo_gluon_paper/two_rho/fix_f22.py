import itertools, numpy as np
from collections import defaultdict
from core import Config, accumulate, term
import user4
import rowcmp
S = rowcmp.user_sets()
W = sum((S[k] for k in S if k.startswith('4:')), [])
AB = sum((S[k] for k in ('A1.1', 'A1.2', 'A2.1', 'A2.2', 'B1.1', 'B1.2', 'B2.1', 'B2.2')), [])
U4, _ = user4.terms()
U4 = [t for t in U4 if t['tag'] != '4rho.tex F1.6']
LO, K = user4.LO, user4.K
U = user4.U
# candidates: the four occurrences of b' are in f^{deb'}, rho^{b'}(x'), U^{ab'}(y), U^{db'}(w'); rename two of them to 'g'
slots = ['f', 'rho', 'Uy', 'Uw']
cands = []
for pair in itertools.combinations(slots, 2):
    nm = {s: ('g' if s in pair else "b'") for s in slots}
    t = term(3, [U('a', nm['Uy'], 'y'), U('d', nm['Uw'], "w'"), U('e', 'b', 'z')], [('a', 'b', 'c'), ('d', 'e', nm['f'])],
             [(nm['rho'], "x'"), ('c', 'y')], [LO, K('y', 'y')], 'F2.2 cand %s' % (pair,))
    cands.append((pair, t))
for seed in (1, 2, 3):
    cfg = Config(seed)
    base = defaultdict(lambda: defaultdict(complex))
    accumulate(U4, cfg, base, False, label='u')
    accumulate([dict(t, c=-t['c']) for t in W], cfg, base, False, label='w')
    accumulate(AB, cfg, base, False, label='ab')
    res0 = {k: sum(d.values()) for k, d in base.items()}
    print('seed', seed, 'without F2.2: nonzero', sum(abs(v) > 1e-9 for v in res0.values()))
    for pair, t in cands:
        acc = defaultdict(lambda: defaultdict(complex))
        accumulate([t], cfg, acc, False, label='c')
        res = dict(res0)
        for k, d in acc.items():
            res[k] = res.get(k, 0) + sum(d.values())
        print('     ', pair, 'nonzero', sum(abs(v) > 1e-9 for v in res.values()))
