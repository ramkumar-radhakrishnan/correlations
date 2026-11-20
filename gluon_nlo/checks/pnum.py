"""Numerical engine for the p+-integrated kernels of Group I row V and Group III."""
import numpy as np
from terms_r5g3 import T
from pconv import indices

KV = np.array([1.3, 0.4])      # transverse k (lattice-free, continuum)
Kp = 1.0                      # k+ (results are boost invariant)
ARGS = 'x,y,z,xp,yp,w,wp,p,K,PL0,PL1,PLQ,PLQK,PLG'

def V(a, b, i): return (a - b)[i]
def PW(a, b, n): return np.sum((a - b) ** 2) ** (n // 2)
def DOT(a, b, c, d): return np.sum((a - b) * (c - d))
def D(i, j): return 1.0 if i == j else 0.0
def EXP(a, b, f): return np.exp(-1j * (KV @ (a - b)) * f)
DPERP = 2.0

FUN = {}
for lab, g, c, meas, var, code in T:
    idx = indices(code)
    body = code.replace(' q/K', ' p/K').replace('q/K', 'p/K')
    body = body.replace('(K-q)', '(K-p)').replace('*q*', '*p*').replace('(q)', '(p)')
    import re
    body = re.sub(r'(?<![A-Za-z_])q(?![A-Za-z_])', 'p', body)          # integration variable -> p
    loops = ''.join(f' for {i} in (0,1)' for i in idx)
    src = f'lambda {ARGS}: sum(({body}){loops})' if idx else f'lambda {ARGS}: ({body})'
    FUN[lab] = (g, c, meas, var, eval(src))

# tanh-sinh nodes on (0,1)
_h = 0.04; _k = np.arange(-120, 121) * _h
_u = 0.5 * (1 + np.tanh(0.5 * np.pi * np.sinh(_k)))
_wu = 0.5 * _h * 0.5 * np.pi * np.cosh(_k) / np.cosh(0.5 * np.pi * np.sinh(_k)) ** 2
_m = (_u > 1e-14) & (_u < 1 - 1e-14); _u, _wu = _u[_m], _wu[_m]

PLW = {'PL0': (lambda p: Kp / p, 0.0), 'PL1': (lambda p: Kp / (Kp - p), Kp), 'PLQ': (lambda p: 1 / p, 0.0),
       'PLQK': (lambda p: Kp / p, 0.0), 'PLG': (lambda p: Kp / (p * (p + Kp)), 0.0)}
PLNAMES = ['PL0', 'PL1', 'PLQ', 'PLQK', 'PLG']

def integrand_parts(f, P, p):
    base = dict(P, p=p, K=Kp, **{n: 0.0 for n in PLNAMES})
    h = f(**base)
    parts = {}
    for n in PLNAMES:
        b2 = dict(base); b2[n] = 1.0
        g = f(**b2) - h
        if np.any(np.abs(g) > 0):
            parts[n] = g
    return h, parts

def term_value(lab, P):
    """value of one term at fixed transverse positions: coefficient * measure * int dp/(2pi) [...]"""
    g, c, meas, var, f = FUN[lab]
    mfac = {'1/K^2': 1 / Kp**2, '1/K': 1 / Kp, '1': 1.0}[meas]
    if lab.startswith('6.'):                                   # Row V: p in (0, k+)
        p = Kp * _u; w = Kp * _wu
        fv = np.vectorize(lambda pp: 0)                         # placeholder
        h = np.array([f(**dict(P, p=pp, K=Kp, **{n: 0.0 for n in PLNAMES})) for pp in p])
        tot = np.sum(w * h)
        for n in PLNAMES:
            gp = np.array([f(**dict(P, p=pp, K=Kp, **{m: (1.0 if m == n else 0.0) for m in PLNAMES})) for pp in p]) - h
            if np.max(np.abs(gp)) == 0: continue
            wf, p0 = PLW[n]
            g0 = f(**dict(P, p=p0, K=Kp, **{m: (1.0 if m == n else 0.0) for m in PLNAMES})) - \
                 f(**dict(P, p=p0, K=Kp, **{m: 0.0 for m in PLNAMES}))
            tot += np.sum(w * wf(p) * (gp - g0))
        return c * mfac * tot / (2 * np.pi)
    # Group III: p in (0, inf), contour rotated p = -i s t, s = sign(k.(z - y'))
    A = KV @ (P['z'] - P['yp']); s = 1.0 if A >= 0 else -1.0
    t = _u / (1 - _u); wt = _wu / (1 - _u) ** 2
    p = -1j * s * t * Kp; jac = -1j * s * Kp
    h = np.array([f(**dict(P, p=pp, K=Kp, **{n: 0.0 for n in PLNAMES})) for pp in p])
    tot = np.sum(wt * jac * h)
    if 'PLG' in FUN[lab][4].__code__.co_consts.__repr__() or True:
        gp = np.array([f(**dict(P, p=pp, K=Kp, **{m: (1.0 if m == 'PLG' else 0.0) for m in PLNAMES})) for pp in p]) - h
        if np.max(np.abs(gp)) > 0:
            g0 = f(**dict(P, p=0.0, K=Kp, **{m: (1.0 if m == 'PLG' else 0.0) for m in PLNAMES})) - \
                 f(**dict(P, p=0.0, K=Kp, **{m: 0.0 for m in PLNAMES}))
            tot += np.sum(wt * jac * PLW['PLG'][0](p) * (gp - g0))
    return c * mfac * tot / (2 * np.pi)

def group_value(group, P):
    return sum(term_value(lab, P) for lab, (g, *_r) in FUN.items() if g == group)

GROUPS = {}
for lab, (g, *_r) in FUN.items():
    GROUPS.setdefault(g, []).append(lab)
