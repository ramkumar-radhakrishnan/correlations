"""Configuration-by-configuration parity of the cubic (rho^3) terms: for term = i*K(k)*C + c.c. with real
colour factor C, oddness under k->-k <=> K(-k) = K(k)^*  (checked after the p+ integral)."""
import numpy as np, pnum
rng = np.random.default_rng(7)
mx = {}
for trial in range(6):
    P = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
    for lab in pnum.FUN:
        if lab.startswith('6.10') or lab.startswith('6.11') or lab.startswith('8.3'):
            continue                                   # rho^2 terms (no explicit i)
        pnum.KV = np.array([1.3, 0.4]);  a = pnum.term_value(lab, P) / pnum.FUN[lab][1]
        pnum.KV = -np.array([1.3, 0.4]); b = pnum.term_value(lab, P) / pnum.FUN[lab][1]
        mx[lab] = max(mx.get(lab, 0), abs(b - np.conj(a)) / (abs(a) + 1e-300))
pnum.KV = np.array([1.3, 0.4])
print('max |K(-k) - K(k)^*| / |K(k)| per cubic term:')
for k, v in mx.items(): print(f'   {k:6s} {v:.1e}   coefficient {pnum.FUN[k][1]}')
