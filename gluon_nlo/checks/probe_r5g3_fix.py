"""(i) (6.11) coincidence limits after the index-sum fix; (ii) z -> infinity for every equation with a
physical far field: U(z) -> 1 (finite target) and U(z) -> fixed random matrix; colour included."""
import numpy as np
from multiprocessing import Pool
import pnum, colour as C
from probe_r5g3 import EQ
U_orig = C.U
def eq_value(eq, P):
    return sum(pnum.term_value(lab, P) * C.colour(g, P) for g in EQ[eq] for lab in pnum.GROUPS[g])
def job(a):
    eq, mover, target, seed, r, far = a
    rng = np.random.default_rng(seed)
    P0 = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    if far is not None:
        M = np.eye(8) if far == 'one' else U_orig(np.array([3.3, -1.7]))
        C.U = lambda x, M=M: M if np.linalg.norm(x) > 5 else U_orig(x)
    c = P0[target] if target != 'inf' else P0['w']
    v = 0
    for t in np.linspace(0, 2*np.pi, 12, endpoint=False) + 0.1234:
        P = dict(P0); P[mover] = c + r*np.array([np.cos(t), np.sin(t)]); v += eq_value(eq, P)
    C.U = U_orig
    return a, r*r*v/12
if __name__ == '__main__':
    jobs = [('6.11', m, t, s, r, None) for m, t in [('z','x'),('z','y'),('y','x'),('y','z'),('xp','wp'),('z','wp'),('x','wp')]
            for s in (1, 2) for r in (1e-1, 1e-2, 1e-3)]
    for eq in EQ:
        for far in ('one', 'const'):
            for s in (1, 2):
                for R in (1e1, 1e2, 1e3, 1e4):
                    jobs.append((eq, 'z', 'inf', s, R, far))
    with Pool(4) as pool:
        res = pool.map(job, jobs, chunksize=2)
    out = {}
    for (eq, m, t, s, r, far), v in res:
        out.setdefault((eq, m, t, s, far), []).append(v)
    for (eq, m, t, s, far), vs in out.items():
        print(f'({eq:4s}) {m:2s}->{t:3s} seed{s} {"" if far is None else "U(far)="+far:12s}: ' + ' '.join(f'{abs(v):.2e}' for v in vs))
