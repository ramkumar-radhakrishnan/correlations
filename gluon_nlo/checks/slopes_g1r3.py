"""Per-term coefficient of the non-integrable piece of Group I row III.
For each term: s_i = d[ r^2 <K_i>_angle ] / d ln(1/r) at z->w (r->0) and d[ R^2 <K_i> ]/d ln R at z->infinity."""
import numpy as np
from kernels_g1r3 import fns
rng = np.random.default_rng(11)
order = ('x', 'y', 'z', 'xp', 'w', 'wp')
labels = []
cnt = {}
for bi, s, k, f in fns:
    eq = {13: '6.4', 14: '6.5', 15: '6.6', 16: '6.7'}[bi]; cnt[eq] = cnt.get(eq, 0) + 1
    labels.append(f'({eq}) term {cnt[eq]}')
def avg(P0, r, center, n=128):
    out = np.zeros(len(fns))
    for t in np.linspace(0, 2*np.pi, n, endpoint=False):
        P = {k: v.copy() for k, v in P0.items()}; P['z'] = center + r*np.array([np.cos(t), np.sin(t)])
        out += np.array([np.real(s*f(*(P[k] for k in order))) for bi, s, k, f in fns])
    return out / n
for trial in range(2):
    P0 = {k: rng.normal(size=2) for k in order}
    a1, a2 = avg(P0, 1e-3, P0['w']), avg(P0, 1e-4, P0['w'])
    sw = (1e-8*a2 - 1e-6*a1) / np.log(10)
    b1, b2 = avg(P0, 1e2, P0['w']), avg(P0, 1e3, P0['w'])
    si = (1e6*b2 - 1e4*b1) / np.log(10)
    print(f'--- configuration {trial}:  per-term log-slopes  (z->w | z->inf)')
    for lab, a, b in zip(labels, sw, si):
        if abs(a) > 1e-6 or abs(b) > 1e-6:
            print(f'   {lab:16s}  {a:+.4f}   {b:+.4f}')
    print(f'   {"SUM":16s}  {sw.sum():+.4f}   {si.sum():+.4f}')
