"""Eq (7.1) kernel, transcribed from the typeset tensor form:
  K = (x'-w')^i (y-z)^j/((x'-w')^2 (y-z)^2) log[(x-z)^2/(x-w)^2] *
      [ d^{ij}/(2(z-w)^2) + (x-z)^j (z-w)^i/((z-w)^2 (x-z)^2) - (x-z)^j (x-w)^i/((x-w)^2 (x-z)^2) - (x-w)^i (z-w)^j/((z-w)^2 (x-w)^2) ]"""
import numpy as np
rng = np.random.default_rng(3)
def K(x, y, z, xp, w):
    A = (xp - w) / np.sum((xp - w)**2); Bv = (y - z) / np.sum((y - z)**2)
    s = lambda v: np.sum(v*v)
    L = np.log(s(x - z) / s(x - w))
    t = (A @ Bv) / (2*s(z - w)) + (A @ (z - w))*(Bv @ (x - z))/(s(z - w)*s(x - z)) \
        - (A @ (x - w))*(Bv @ (x - z))/(s(x - w)*s(x - z)) - (A @ (x - w))*(Bv @ (z - w))/(s(z - w)*s(x - w))
    return L * t
def avg(P, r, c, n=128):
    v = 0
    for t in np.linspace(0, 2*np.pi, n, endpoint=False):
        z = c + r*np.array([np.cos(t), np.sin(t)]); v += K(P['x'], P['y'], z, P['xp'], P['w'])
    return v/n
for _ in range(3):
    P = {k: rng.normal(size=2) for k in ('x', 'y', 'xp', 'w')}
    print('z->w   r^2<K>:', ' '.join(f'{r*r*avg(P, r, P["w"]):+.2e}' for r in (1e-1, 1e-2, 1e-3, 1e-4)),
          '   z->inf R^2<K>:', ' '.join(f'{R*R*avg(P, R, P["w"]):+.2e}' for R in (1e1, 1e2, 1e3)))
