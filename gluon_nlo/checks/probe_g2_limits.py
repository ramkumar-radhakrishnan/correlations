import numpy as np
from probe_g2 import total
rng = np.random.default_rng(5)
order = ('x', 'xp', 'z', 'w', 'wp')
def avg(P0, mover, c, r, n=96):
    v = 0
    for t in np.linspace(0, 2*np.pi, n, endpoint=False):
        P = dict(P0); P[mover] = c + r*np.array([np.cos(t), np.sin(t)]); v += total(tuple(P[k] for k in order))
    return v/n
print('Group II (7.3)-(7.8) summed kernel, colour factor stripped:')
for mover, target in [('z', 'w'), ('z', 'wp'), ('z', 'x'), ('z', 'xp'), ('x', 'w'), ('xp', 'wp'), ('w', 'wp')]:
    rows = []
    for _ in range(2):
        P0 = {k: rng.normal(size=2) for k in order}
        rows.append(' '.join(f'{r*r*avg(P0, mover, P0[target], r):+.1e}' for r in (1e-1, 1e-2, 1e-3, 1e-4)))
    print(f'  {mover:2s}->{target:2s} r^2<K>: ' + ' | '.join(rows))
for mover in ('z',):
    rows = []
    for _ in range(2):
        P0 = {k: rng.normal(size=2) for k in order}
        rows.append(' '.join(f'{R*R*avg(P0, mover, P0["w"], R):+.1e}' for R in (1e1, 1e2, 1e3)))
    print(f'  {mover}->inf  R^2<K>: ' + ' | '.join(rows))
