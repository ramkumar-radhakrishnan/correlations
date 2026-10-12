"""Pointwise totals (no bookkeeping): coefficient of log(vee/Lambda) (A) vs JIMWLK x LO and coefficient of log(vee) (P) vs
BFKL x LO; two-rho terms and symmetrized three-rho terms, as written and corrected. -> allstates.json"""
import env  # noqa
import json
import numpy as np
from collections import defaultdict
from core import Config, accumulate
from compare import report
import target as TG
import refP
import allsets as AS
from acc3 import Config3, accumulate3, report3

J, B = TG.jimwlk(), refP.bfkl()
out = {}
for name, f, tgt in (('A', AS.A_sets, J), ('P', AS.P_sets, B)):
    for corr in (False, True):
        L = f((2,), corr)
        cfg = Config(11)
        acc = defaultdict(lambda: defaultdict(complex))
        for k, ts in L.items():
            accumulate(ts, cfg, acc, False, label=k)
        accumulate([dict(t, c=-t['c']) for t in tgt], cfg, acc, False, label='-T')
        bad = report(acc, show=False)
        mx = max(float(np.max(np.abs(sum(d.values())))) for d in acc.values())
        out['%s2_%s' % (name, 'corr' if corr else 'written')] = [len(bad), len(acc), float(mx)]
        L3 = f((3,), corr)
        cfg = Config3(7)
        acc = defaultdict(lambda: defaultdict(complex))
        for k, ts in L3.items():
            accumulate3(ts, cfg, acc, label=k)
        bad = report3(acc, show=False)
        mx = max(float(np.max(np.abs(sum(d.values())))) for d in acc.values())
        out['%s3_%s' % (name, 'corr' if corr else 'written')] = [len(bad), len(acc), float(mx)]
        print(name, corr, out['%s2_%s' % (name, 'corr' if corr else 'written')], out['%s3_%s' % (name, 'corr' if corr else 'written')])
json.dump(out, open('allstates.json', 'w'), indent=1)
