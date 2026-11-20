import numpy as np
from kernels_g1r3 import total, A0, A0_arctan, dot
rng = np.random.default_rng(4)
def pts():
    return dict(x=rng.normal(size=2), y=rng.normal(size=2), z=rng.normal(size=2),
                xp=rng.normal(size=2), w=rng.normal(size=2), wp=rng.normal(size=2))
order = ('x', 'y', 'z', 'xp', 'w', 'wp')
T = lambda P: total(tuple(P[k] for k in order))
parts = lambda P: total(tuple(P[k] for k in order), parts=True)

# 0. eq (6.8): arctan form == arccos form of A0
mx = 0
for _ in range(2000):
    P = pts(); a, b = A0(P['x'], P['z'], P['w']), A0_arctan(P['x'], P['z'], P['w'])
    mx = max(mx, abs(a - b) / abs(a))
print('(6.8) arctan form vs arccos form of A0: max rel diff', f'{mx:.1e}')

# 1. imaginary parts (logs of negative arguments)
im_pieces, im_tot = 0, 0
for _ in range(2000):
    P = pts(); v = parts(P)
    im_pieces = max(im_pieces, max(abs(np.imag(t)) for t in v)); im_tot = max(im_tot, abs(np.imag(sum(v))) / (abs(sum(v)) + 1e-300))
print(f'imaginary parts: largest single term |Im| = {im_pieces:.2e};  sum: max |Im|/|total| = {im_tot:.1e}')

def approach(name, setter, n=4):
    print(name)
    for _ in range(n):
        P0 = pts(); setter(P0)
        row = []
        for eps in (1e-2, 1e-4, 1e-6):
            vs = []
            for s in (1, -1):
                P = {k: v.copy() for k, v in P0.items()}; P['w'] = P0['w'] + s * eps * P0['_dir']
                vs.append(np.real(T(P)))
            row.append(vs)
        P = {k: v.copy() for k, v in P0.items()}; P['w'] = P0['w'] + 1e-6 * P0['_dir']
        big = max(abs(t) for t in parts(P))
        print('   sum at 1e-2,1e-4,1e-6 (+/-):', ' | '.join(f'{a:+.4e},{b:+.4e}' for a, b in row), f'  largest term {big:.1e}')

def on_bisector(P):        # put w on the surface |x-w| = |x-z| and move it radially
    r = P['z'] - P['x']; ang = rng.uniform(0, 2*np.pi)
    P['w'] = P['x'] + np.linalg.norm(r) * np.array([np.cos(ang), np.sin(ang)])
    P['_dir'] = (P['w'] - P['x']) / np.linalg.norm(P['w'] - P['x'])
def on_segment(P):         # put w on the segment x--z and move it perpendicular
    t = rng.uniform(0.2, 0.8); P['w'] = P['x'] + t * (P['z'] - P['x'])
    d = P['z'] - P['x']; P['_dir'] = np.array([-d[1], d[0]]) / np.linalg.norm(d)
def on_line_outside(P):    # w collinear with x,z but outside the segment
    t = rng.choice([-0.6, 1.7]); P['w'] = P['x'] + t * (P['z'] - P['x'])
    d = P['z'] - P['x']; P['_dir'] = np.array([-d[1], d[0]]) / np.linalg.norm(d)
approach('surface (x-z)^2 = (x-w)^2:', on_bisector)
approach('w on the segment x--z (A0 ~ pi/S):', on_segment)
approach('w collinear with x,z outside the segment:', on_line_outside)

# 2. z -> w : is the sum integrable in d^2(z-w)?  look at |sum| * |z-w|^2 (must ->0 for integrability)
print('z -> w: |z-w|^2 * |sum|  (-> 0 means integrable; const or growing means not):')
for _ in range(4):
    P0 = pts(); d = rng.normal(size=2); d /= np.linalg.norm(d)
    vals = []
    for r in (1e-1, 1e-2, 1e-3, 1e-4):
        P = {k: v.copy() for k, v in P0.items()}; P['z'] = P['w'] + r * d
        vals.append(r * r * abs(T(P)))
    print('   ', ' '.join(f'{v:.3e}' for v in vals))
print('z -> w: same for the single (6.7) term with log((z-w)^2)/(z-w)^2:')
from kernels_g1r3 import fns
f_last = fns[-1][3]
for _ in range(2):
    P0 = pts(); d = rng.normal(size=2); d /= np.linalg.norm(d)
    vals = []
    for r in (1e-1, 1e-2, 1e-3, 1e-4):
        P = {k: v.copy() for k, v in P0.items()}; P['z'] = P['w'] + r * d
        vals.append(r * r * abs(f_last(*(P[k] for k in order))))
    print('   ', ' '.join(f'{v:.3e}' for v in vals))
