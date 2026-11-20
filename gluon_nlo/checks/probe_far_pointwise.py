"""z -> infinity, pointwise along fixed directions (no angular aliasing): R^2 |sum of terms x colour|,
U(z) -> 1 outside |z|>5 (finite target).  Integrable iff R^2|K| -> 0."""
import numpy as np
from multiprocessing import Pool
import pnum, colour as C
from probe_r5g3 import EQ
U_orig = C.U
def job(a):
    eq, seed, ang = a
    C.U = lambda x: np.eye(8) if np.linalg.norm(x) > 5 else U_orig(x)
    rng = np.random.default_rng(seed)
    P0 = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    n = np.array([np.cos(ang), np.sin(ang)]); out = []
    for R in (1e1, 1e2, 1e3, 1e4):
        P = dict(P0); P['z'] = P0['w'] + R*n
        v = sum(pnum.term_value(lab, P) * C.colour(g, P) for g in EQ[eq] for lab in pnum.GROUPS[g])
        out.append(R*R*abs(v))
    return a, out
if __name__ == '__main__':
    jobs = [(eq, s, ang) for eq in EQ for s in (1, 2) for ang in (0.3, 2.1, 4.0)]
    with Pool(4) as p:
        res = p.map(job, jobs)
    for (eq, s, ang), out in res:
        print(f'({eq:4s}) seed{s} dir {ang:.1f}: R^2|K| at R=1e1..1e4: ' + ' '.join(f'{v:.2e}' for v in out))
