"""(6.11) at z -> x: split into the main bracket (A) and the extra term (B) = -e^{-ik(x-w')} (1/k+)[p(k+-p)/k+^2 - 2](...)/(x-z)^2."""
import numpy as np, pnum
from pconv import indices
lab = '6.11'; g, c, meas, var, f = pnum.FUN[lab]
from terms_r5g3 import T
code = [t for t in T if t[0] == '6.11'][0][5]
i = code.index('-EXP(x,wp,1)')
codeA, codeB = code[1:i], code[i + 1:-1]          # strip outer parentheses
def mk(cd):
    idx = indices(cd); loops = ''.join(f' for {q} in (0,1)' for q in idx)
    return eval(f'lambda {pnum.ARGS}: sum(({cd}){loops})' if idx else f'lambda {pnum.ARGS}: ({cd})', vars(pnum))
fA, fB = mk(codeA), mk('-(' + codeB + ')')
def val(fn, P):
    pnum.FUN['6.tmp'] = (g, c, meas, var, fn); return pnum.term_value('6.tmp', P)
rng = np.random.default_rng(2)
for trial in range(3):
    P0 = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    rows = {'A': [], 'B': [], 'A+B': []}
    for r in (1e-1, 1e-2, 1e-3):
        a = b = 0
        for t in np.linspace(0, 2*np.pi, 16, endpoint=False) + 0.1:
            P = dict(P0); P['z'] = P0['x'] + r*np.array([np.cos(t), np.sin(t)])
            a += val(fA, P); b += val(fB, P)
        a, b = r*r*a/16, r*r*b/16
        rows['A'].append(a); rows['B'].append(b); rows['A+B'].append(a + b)
    print(f'config {trial}: ' + ' | '.join(f'{k}: ' + ' '.join(f'{v:+.3e}' for v in vs) for k, vs in rows.items()))

print('probe seed-1 configuration and a scan over random configurations:')
def resid(P0, r=1e-3, n=16):
    a = b = 0
    for t in np.linspace(0, 2*np.pi, n, endpoint=False) + 0.1:
        P = dict(P0); P['z'] = P0['x'] + r*np.array([np.cos(t), np.sin(t)])
        a += val(fA, P); b += val(fB, P)
    return r*r*a/n, r*r*b/n
for seed in range(1, 9):
    rng = np.random.default_rng(seed)
    P0 = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    a, b = resid(P0)
    print(f'  seed {seed}: A {a:+.3e}  B {b:+.3e}  A+B {a+b:+.3e}   |x-y|={np.linalg.norm(P0["x"]-P0["y"]):.2f}')
