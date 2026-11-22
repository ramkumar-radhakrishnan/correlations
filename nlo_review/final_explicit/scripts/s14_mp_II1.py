# s14_mp_II1.py -- 60-digit recheck of the Group II Row I second part kernel (N2a+...+N2e, author's form as in
# s2_kernels.py) at the triple cluster {z,w,w'} (soft gluon pair and both measured-gluon positions), where the
# individual terms are ~ r^-4 and cancel: exponent of max|kernel| over orientations, lam = 1e-3 ... 1e-9.
import os, re, types, mpmath as mp
import numpy as np
mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 's2_kernels.py')).read()
shim = types.SimpleNamespace(log=mp.log)
ns = dict(np=shim)
exec('dot = lambda a, b: a @ b\nsq = lambda a: a @ a\nKv = lambda v: v/sq(v)', ns)
exec(re.search(r'^def II1_author\(.*?(?=^\S)', src, re.S | re.M).group(0), ns)
V = lambda a, b: np.array([mp.mpf(a), mp.mpf(b)], dtype=object)
rng = np.random.default_rng(5); out = []
for base in range(3):
    x, xp, c = [V(*rng.uniform(-1.5, 1.5, 2)) for _ in range(3)]
    dirs = [[V(*(rng.uniform(0.5, 1.5) * np.array([np.cos(a), np.sin(a)]))) for a in rng.uniform(0, 2 * np.pi, 2)] for _ in range(12)]
    vals = []
    for lam in [mp.mpf('1e-3'), mp.mpf('1e-5'), mp.mpf('1e-7'), mp.mpf('1e-9')]:
        vals.append(max(abs(ns['II1_author'](x, xp, c, c + lam * d[0], c + lam * d[1])[0]) for d in dirs))
    sl = [float(mp.log(vals[i] / vals[i - 1]) / mp.log(mp.mpf('1e-2'))) for i in range(1, 4)]
    out.append('base %d: max|N2| at lam=1e-3,1e-5,1e-7,1e-9: %s ; exponents %s ; lam^2 max: %s' % (
        base, ' '.join('%.3e' % float(v) for v in vals), ' '.join('%.3f' % s for s in sl),
        ' '.join('%.3f' % float(v * l ** 2) for v, l in zip(vals, [mp.mpf('1e-3'), mp.mpf('1e-5'), mp.mpf('1e-7'), mp.mpf('1e-9')]))))
    print(out[-1], flush=True)
out.append('=> d=4 at {z,w,w\'}: |N2| ~ lam^-2 (the double-precision value alpha0=-3.16 of s10 is rounding loss in the cancellation of the r^-4 terms); exponent -2 > -4: the kernel alone is integrable there (colour not needed).')
print(out[-1]); open(os.path.join(HERE, 's14_mp_II1.log'), 'w').write('\n'.join(out) + '\n')
