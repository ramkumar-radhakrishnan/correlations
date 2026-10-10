"""The two-rho terms of ThreeRho_to_TwoRho (results.pkl pieces) vs the exact Weyl two-rho part (users.R3)."""
import pickle, sys, os
import numpy as np
import sympy as sp
from collections import defaultdict
from core import Config, accumulate, term, rename, THREE
import users as US
from cancel3 import CANON
from engine import Nc
# results.pkl is written by ../three_rho/process.py
res = pickle.load(open(os.path.join(THREE, 'results.pkl'), 'rb'))


def prev_terms():
    out = []
    for sres in res:
        src = sres['src']
        kern_list = US.native_kernel(src)
        base = dict(CANON.get(src.key, {}))
        for n, sub in enumerate(sres['subs']):
            for pc in sub['pieces']:
                old, new = pc['subst']
                c0 = complex(sp.N((src.pref * pc['coeff_before']).subs(Nc, 3)))
                facs = pc['facs_before']
                Us = [(f[1], f[2], f[3]) for f in facs if f[0] == 'U']
                fs = [tuple(f[1:]) for f in facs if f[0] == 'f']
                rs = [(f[1], f[2]) for f in facs if f[0] == 'rho']
                for m, (kc, kern) in enumerate(kern_list):
                    kk = [(f[0],) + tuple(new if x == old else x for x in f[1:]) for f in kern]
                    t = term(2 * c0 * complex(sp.N(kc)), Us, fs, rs, kk, 'prev:%s.%d' % (src.key, n + 1))
                    t = rename(t, base)
                    out.append(t)
                    if src.cc:
                        out.append(rename(t, {"w'": 'w', 'w': "w'"}))
    return out


P = prev_terms()
R3 = US.R3()
for seed in (1, 2):
    cfg = Config(seed)
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate(P, cfg, acc, False, label='previous note')
    accumulate([dict(t, c=-t['c']) for t in R3], cfg, acc, False, label='-exact')
    nz = [(k, sum(d.values())) for k, d in acc.items() if abs(sum(d.values())) > 1e-9 * max(1, sum(abs(v) for v in d.values()))]
    print('seed', seed, 'keys', len(acc), 'nonzero', len(nz), 'max', max([abs(v) for _, v in nz] + [0]))
