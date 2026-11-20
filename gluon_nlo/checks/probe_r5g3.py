import numpy as np, sys, time
from multiprocessing import Pool
import pnum, colour as C
EQ = {'6.9': ['C69'], '6.10': ['C610'], '6.11': ['C611'], '8.2': ['X1', 'X2'], '8.3': ['Y1', 'Y2'],
      '8.4': ['W1', 'W2'], '8.5': ['Z1', 'Z2']}
def eq_value(eq, P):
    return sum(pnum.term_value(lab, P) * C.colour(g, P) for g in EQ[eq] for lab in pnum.GROUPS[g])
LIM = {
 'rowV': [('z', 'x'), ('z', 'y'), ('y', 'x'), ('y', 'z'), ('xp', 'wp'), ('z', 'wp'), ('x', 'wp'), ('yp', 'x'), ('yp', 'y')],
 'g3':   [('z', 'yp'), ('z', 'w'), ('z', 'x'), ('z', 'xp'), ('z', 'y'), ('y', 'w'), ('x', 'z'), ('xp', 'yp'), ('yp', 'w'), ('x', 'y')],
}
NANG = 12
def job(args):
    eq, mover, target, seed, r = args
    rng = np.random.default_rng(seed)
    P0 = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    c = P0[target] if target != 'inf' else P0['w']
    v = 0
    for t in np.linspace(0, 2*np.pi, NANG, endpoint=False) + 0.1234:
        P = dict(P0); P[mover] = c + r*np.array([np.cos(t), np.sin(t)]); v += eq_value(eq, P)
    return (eq, mover, target, seed, r, r*r*v/NANG)
if __name__ == '__main__':
    t0 = time.time(); rng = np.random.default_rng(0)
    P = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    for eq in EQ: eq_value(eq, P)
    print('one evaluation of all 7 equations: %.2fs' % (time.time() - t0)); sys.stdout.flush()
    jobs = []
    for eq in EQ:
        lims = LIM['rowV'] if eq.startswith('6') else LIM['g3']
        if eq != '6.9': lims = [l for l in lims if 'yp' not in l] if eq.startswith('6') else lims
        for mover, target in lims:
            for seed in (1, 2):
                for r in (1e-1, 1e-2, 1e-3):
                    jobs.append((eq, mover, target, seed, r))
        for seed in (1, 2):
            for R in (1e1, 1e2, 1e3):
                jobs.append((eq, 'z', 'inf', seed, R))
    with Pool(4) as pool:
        res = pool.map(job, jobs, chunksize=4)
    out = {}
    for eq, mover, target, seed, r, v in res:
        out.setdefault((eq, mover, target, seed), []).append(v)
    for (eq, mover, target, seed), vs in out.items():
        print(f'({eq:4s}) {mover:2s}->{target:3s} seed{seed}: ' + ' '.join(f'{abs(v):.2e}' for v in vs))
