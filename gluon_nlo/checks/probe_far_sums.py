"""Large-z tails: per-term origin, and the sums of equations that share kernels:
   S_III = (8.2)+(8.4)+(8.5)  (identical kernels a, b1, b2; different colour),  S_V = (6.10)+(6.11)."""
import numpy as np
from multiprocessing import Pool
import pnum, colour as C
U_orig = C.U
SUMS = {'(8.2)+(8.4)+(8.5)': ['X1', 'X2', 'W1', 'W2', 'Z1', 'Z2'], '(6.10)+(6.11)': ['C610', 'C611'],
        '(8.2) alone': ['X1', 'X2'], '(8.3) alone': ['Y1', 'Y2']}
def job(a):
    name, seed, ang = a
    C.U = lambda x: np.eye(8) if np.linalg.norm(x) > 5 else U_orig(x)
    rng = np.random.default_rng(seed); P0 = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    n = np.array([np.cos(ang), np.sin(ang)]); out = []
    for R in (1e1, 1e2, 1e3, 1e4):
        P = dict(P0); P['z'] = P0['w'] + R*n
        v = sum(pnum.term_value(l, P) * C.colour(g, P) for g in SUMS[name] for l in pnum.GROUPS[g])
        out.append(R*R*abs(v))
    return a, out
def per_term(seed=1, ang=0.3):
    C.U = lambda x: np.eye(8) if np.linalg.norm(x) > 5 else U_orig(x)
    rng = np.random.default_rng(seed); P0 = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    n = np.array([np.cos(ang), np.sin(ang)])
    for g in ['X1', 'X2', 'Y1', 'Y2']:
        for l in pnum.GROUPS[g]:
            # split the PLG-subtracted piece from the rest
            gg, c, meas, var, f = pnum.FUN[l]
            noPL = lambda **kw: f(**dict(kw, PLG=0.0))
            pnum.FUN['tmpA'] = (gg, c, meas, var, noPL)
            row = []
            for R in (1e2, 1e3, 1e4):
                P = dict(P0); P['z'] = P0['w'] + R*n; col = C.colour(g, P)
                full = pnum.term_value(l, P) * col; rest = pnum.term_value('tmpA', P) * col
                row.append((R*R*abs(full), R*R*abs(full - rest)))
            print(f'   {l:6s} R^2|term| (R=1e2,1e3,1e4): ' + ' '.join(f'{a:.2e}' for a, b in row) +
                  '   of which [ ]_+ piece: ' + ' '.join(f'{b:.2e}' for a, b in row))
if __name__ == '__main__':
    print('per-term large-z tails, seed 1, one direction:'); per_term()
    jobs = [(nm, s, ang) for nm in SUMS for s in (1, 2) for ang in (0.3, 2.1)]
    with Pool(4) as p: res = p.map(job, jobs)
    for (nm, s, ang), out in res:
        print(f'{nm:20s} seed{s} dir {ang}: R^2|K| R=1e1..1e4: ' + ' '.join(f'{v:.2e}' for v in out))
