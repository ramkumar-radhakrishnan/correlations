import numpy as np
from kernels_g1r3 import total
rng = np.random.default_rng(11)
order = ('x', 'y', 'z', 'xp', 'w', 'wp')
def pts(): return {k: rng.normal(size=2) for k in order}
T = lambda P: np.real(total(tuple(P[k] for k in order)))     # PV reading: log|.|

def ang_avg(P0, mover, target, r, n=64):
    vals = []
    for t in np.linspace(0, 2*np.pi, n, endpoint=False):
        P = {k: v.copy() for k, v in P0.items()}
        P[mover] = P0[target] + r * np.array([np.cos(t), np.sin(t)])
        vals.append(T(P))
    vals = np.array(vals)
    return vals.mean(), np.abs(vals).mean()

print('coincidence limits: r^2 * <K>_angle  and  r^2 * <|K|>_angle  (integrable in d^2r iff r^2<K> -> 0)')
for mover, target in [('z', 'w'), ('x', 'w'), ('x', 'z'), ('y', 'z'), ('y', 'x')]:
    print(f'  {mover} -> {target}:')
    for _ in range(3):
        P0 = pts()
        row = [ang_avg(P0, mover, target, r) for r in (1e-1, 1e-2, 1e-3, 1e-4)]
        rr = [1e-1, 1e-2, 1e-3, 1e-4]
        print('     r^2<K>:', ' '.join(f'{r*r*m:+.2e}' for r, (m, a) in zip(rr, row)),
              '   r^2<|K|>:', ' '.join(f'{r*r*a:.2e}' for r, (m, a) in zip(rr, row)))

print('large distances: R^2 * <K>_angle for one point sent to infinity (integrable iff -> 0)')
for mover in ('z', 'y', 'x'):
    print(f'  {mover} -> infinity:')
    for _ in range(2):
        P0 = pts(); c = P0['w']
        row = []
        for R in (1e1, 1e2, 1e3):
            vals = []
            for t in np.linspace(0, 2*np.pi, 64, endpoint=False):
                P = {k: v.copy() for k, v in P0.items()}; P[mover] = c + R*np.array([np.cos(t), np.sin(t)])
                vals.append(T(P))
            row.append(R*R*np.mean(vals))
        print('     ', ' '.join(f'{v:+.3e}' for v in row))
