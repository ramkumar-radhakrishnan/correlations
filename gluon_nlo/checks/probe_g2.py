import numpy as np
from scipy.optimize import brentq
from kernels_g2 import pieces, dot, pw, log
dperp = 2.0
fns = [(bi, s, eval('lambda x,xp,z,w,wp: ' + k)) for bi, s, k in pieces]
def total(P):
    return sum(s * f(*P) for bi, s, f in fns)
def parts(P):
    return [s * f(*P) for bi, s, f in fns]
def sq(a): return float(np.sum(a * a))
surf = {
 'D=(x\'-w\')^2(x-z)^2-(x-w)^2(x\'-z)^2': lambda x,xp,z,w,wp: sq(xp-wp)*sq(x-z)-sq(x-w)*sq(xp-z),
 '(x-z)^2-(x-w)^2':                       lambda x,xp,z,w,wp: sq(x-z)-sq(x-w),
 '(x\'-z)^2-(x\'-w\')^2':                 lambda x,xp,z,w,wp: sq(xp-z)-sq(xp-wp),
}
if __name__ == '__main__':
    rng = np.random.default_rng(2)
    for name, F in surf.items():
        print('\nsurface', name)
        found = 0
        while found < 4:
            x, xp, z, w, wp = rng.normal(size=(5, 2))
            # move w' (or w) along a line to hit the surface
            mover = 'w' if name.startswith('(x-z)') else 'wp'
            d = rng.normal(size=2)
            def g(t):
                P = dict(x=x, xp=xp, z=z, w=w, wp=wp); P[mover] = P[mover] + t * d
                return F(P['x'], P['xp'], P['z'], P['w'], P['wp'])
            ts = np.linspace(-3, 3, 601); vals = [g(t) for t in ts]
            idx = [i for i in range(600) if vals[i] * vals[i + 1] < 0]
            if not idx: continue
            t0 = brentq(g, ts[idx[0]], ts[idx[0] + 1], xtol=1e-15)
            found += 1
            row = []
            for eps in (1e-2, 1e-4, 1e-6):
                vs = []
                for sgn in (+1, -1):
                    P = dict(x=x, xp=xp, z=z, w=w, wp=wp); P[mover] = P[mover] + (t0 + sgn * eps) * d
                    vs.append(total(tuple(P[k] for k in ('x', 'xp', 'z', 'w', 'wp'))))
                row.append(vs)
            P = dict(x=x, xp=xp, z=z, w=w, wp=wp); P[mover] = P[mover] + (t0 + 1e-6) * d
            mx = max(abs(v) for v in parts(tuple(P[k] for k in ('x','xp','z','w','wp'))))
            print('  total at eps=1e-2,1e-4,1e-6 (+side,-side):',
                  ' | '.join(f'{a:+.3e},{b:+.3e}' for a, b in row), f'   largest single piece at 1e-6: {mx:.1e}')
