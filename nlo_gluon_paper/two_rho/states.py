"""Cumulative fixes: number of kernel structures (pointwise keys) where the two-rho total differs from JIMWLK x LO."""
import json
import numpy as np
from collections import defaultdict
from core import Config, accumulate
import users as US
import user4
import target as TG
from compare import report
from fullpt import g1_pointwise, g2iv_pointwise
import rowcmp

S = rowcmp.user_sets()
R3all = sum((S[k] for k in S if k.startswith(('GI.', 'GII.', 'GIII.', 'A1', 'A2', 'B1', 'B2'))), [])
U4w, _ = user4.terms()                      # as written (F2.2 cannot be evaluated: index b' four times)
U4c = user4.corrected_terms()
gen = S['G:G1_V'] + S['G:G2_I'] + g1_pointwise(S['G:G1_II']) + g2iv_pointwise(S['G:G2_IV'], 0.5)


def total(sets, seed=11):
    cfg = Config(seed)
    acc = defaultdict(lambda: defaultdict(complex))
    for name, ts in sets:
        accumulate(ts, cfg, acc, False, label=name)
    accumulate([dict(t, c=-t['c']) for t in TG.jimwlk()], cfg, acc, False, label='-JIMWLK')
    bad = report(acc, show=False)
    mx = max(abs(sum(d.values())) for d in acc.values())
    return len(acc), len(bad), mx


states = [
    ('as written', [('R3', R3all), ('4rho.tex', U4w), ('gen', gen), ('G3', S['G:G3_III'])]),
    ('+ 4rho.tex: drop the -8Nc term, fix the index of term 2', [('R3', R3all), ('4rho.tex', U4c), ('gen', gen), ('G3', S['G:G3_III'])]),
    ('+ Row III part 2 -> D', [('R3', R3all), ('4rho.tex', U4c), ('gen', gen), ('G3', S['G:G3_III']), ('-III2', S['-RowIII2']), ('D', S['D'])]),
    ('+ Row II: F_alpha removed from U(w)A3', [('R3', R3all), ('4rho.tex', U4c), ('gen', gen), ('G3', S['G:G3_III']), ('-III2', S['-RowIII2']), ('D', S['D']), ('-IIa', S['-RowIIaU'])]),
    ('+ Group III Row III part 2: sign', [('R3', R3all), ('4rho.tex', U4c), ('gen', gen), ('G3', US.neg(S['G:G3_III'])), ('-III2', S['-RowIII2']), ('D', S['D']), ('-IIa', S['-RowIIaU'])]),
]
out = []
for name, sets in states:
    n, nb, mx = total(sets)
    out.append(dict(state=name, keys=n, bad=nb, max=mx))
    print('%-60s structures %d  differing %d  max %.3e' % (name, n, nb, mx))
json.dump(out, open('states.json', 'w'), indent=1)
