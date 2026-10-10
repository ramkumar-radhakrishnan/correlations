import itertools, numpy as np
from collections import defaultdict
from core import Config, accumulate
import users as US
import target as TG
from compare import report
from fullpt import g1_pointwise, g2iv_pointwise
import rowcmp
S = rowcmp.user_sets()
base = [k for k in S if not k.startswith(('A1', 'A2', 'B1', 'B2', 'D', '-', 'G:'))]


def run(extra, seed=11, d_terms=None, symmetric=True, g3sign=-1, g1fix=True):
    cfg = Config(seed, symmetric)
    acc = defaultdict(lambda: defaultdict(complex))
    for k in base:
        accumulate(S[k], cfg, acc, False, label=k)
    accumulate(S['G:G1_V'], cfg, acc, False, label='G1_V')
    accumulate(S['G:G2_I'], cfg, acc, False, label='G2_I')
    accumulate([dict(t, c=g3sign * t['c']) for t in S['G:G3_III']], cfg, acc, False, label='G3_III')
    accumulate(g1_pointwise(S['G:G1_II'], 1.0 if g1fix else 0.0), cfg, acc, False, label='G1_II')
    accumulate(g2iv_pointwise(S['G:G2_IV'], 0.5), cfg, acc, False, label='G2_IV')
    accumulate(S['-RowIII2'], cfg, acc, False, label='-RowIII2')
    accumulate(S['D'] if d_terms is None else d_terms, cfg, acc, False, label='D')
    for name, ts in extra:
        accumulate(ts, cfg, acc, False, label=name)
    accumulate([dict(t, c=-t['c']) for t in TG.jimwlk()], cfg, acc, False, label='-JIMWLK')
    bad = report(acc, show=False)
    mx = max([abs(sum(d.values())) for d in acc.values()])
    return len(acc), len(bad), mx


import json
res = {}
res['simple'] = run([('-RowIIaU', S['-RowIIaU'])]); print('correction 1, simple form      :', res['simple'])
res['commutator'] = run([('-RowIIa all', US.neg(US.RowII_alpha_all())), ('COMM', US.commutator_term())]); print('commutator', res['commutator'])
for order in itertools.permutations(range(3)):
    nm_ = ' '.join([["x'", 'x', 'y'][k] for k in order]); res['D ' + nm_] = run([('-RowIIaU', S['-RowIIaU'])], d_terms=US.D_terms_order(order)); print(nm_, res['D ' + nm_])
for seed in (21, 33):
    res['seed %d' % seed] = run([('-RowIIaU', S['-RowIIaU'])], seed=seed); print(seed, res['seed %d' % seed])
res['general U'] = run([('-RowIIaU', S['-RowIIaU'])], symmetric=False); print('general U', res['general U'])
json.dump({k: [int(v[0]), int(v[1]), float(v[2])] for k, v in res.items()}, open('variants.json', 'w'), indent=1)
