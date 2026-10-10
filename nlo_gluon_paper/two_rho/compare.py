"""Compare the two-rho terms of the notes with JIMWLK x LO, key by key."""
import time
import pickle
import numpy as np
from collections import defaultdict
from core import Config, accumulate, accumulate_dr
import users as US
import target as TG


def group_of(tag):
    if tag.startswith('4rho'):
        return 'R4'
    if tag.startswith('G.'):
        return 'G'
    if tag == 'JIMWLK':
        return 'T'
    if tag.startswith('D.'):
        return 'D'
    return 'R3'


def build(seed=1, corrections=None, integrate=True, symmetric=True, extra=()):
    cfg = Config(seed, symmetric)
    acc = defaultdict(lambda: defaultdict(complex))
    sets = dict(R3=US.R3(), R4=US.R4())
    G = US.genuine()
    for k in ('G1_V', 'G2_I', 'G3_III'):
        accumulate(G[k], cfg, acc, integrate)
    accumulate_dr(G['G1_II'], cfg, acc)
    accumulate_dr(G['G2_IV'], cfg, acc)
    accumulate(sets['R3'], cfg, acc, integrate)
    accumulate(sets['R4'], cfg, acc, integrate)
    for name, ts in extra:
        accumulate(ts, cfg, acc, integrate, label=name)
    tg = TG.jimwlk()
    accumulate([dict(t, c=-t['c']) for t in tg], cfg, acc, integrate, label='-JIMWLK')
    return acc


def report(acc, tol=1e-9, top=60, show=True):
    bad = []
    for key, d in acc.items():
        tot = sum(d.values())
        scale = sum(abs(v) for v in d.values())
        if abs(tot) > tol * max(1.0, scale):
            bad.append((key, tot, scale, d))
    bad.sort(key=lambda x: -abs(x[1]))
    if show:
        print('keys:', len(acc), ' nonzero residual:', len(bad))
        for key, tot, scale, d in bad[:top]:
            print('%10.4f %10.4f  %s' % (abs(tot), scale, key))
            parts = defaultdict(complex)
            for lab, v in d.items():
                parts[lab] += v
            print('            ', {l: np.round(v, 4) for l, v in parts.items() if abs(v) > 1e-9})
    return bad


if __name__ == '__main__':
    t0 = time.time()
    acc = build()
    print('built in %.0fs' % (time.time() - t0))
    pickle.dump({k: dict(v) for k, v in acc.items()}, open('acc_written.pkl', 'wb'))
    report(acc)
