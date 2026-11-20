# Remark (i): do the contact terms cancel between rows?  Lines with the same kernel (same kernel function after the
# contact substitution, same integration variables and phase) are grouped; the sum of their colour factors (with the
# printed numerical factors) is evaluated with random SU(3) Wilson lines.  A cancellation at the integrand level would
# require the sum to vanish identically.
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, lib, lines, sympy as sp
class RandField:
    def __init__(s, seed): s.rng = np.random.default_rng(seed); s.c = {}
    def U(s, x):
        k = tuple(np.round(x, 12))
        if k not in s.c:
            A = s.rng.normal(size=lib.NA) * 2.0
            M = np.einsum('a,abc->bc', A, np.real(1j * lib.TADJ)); lam, V = np.linalg.eigh(1j * M)
            s.c[k] = np.real((V * np.exp(-1j * lam)) @ V.conj().T)
        return s.c[k]
groups = {}
for u in lines.U:
    if u.name in ('deltaPV', 'delta(z) alone'): continue
    key = (u.name if u.row.startswith('III.') and u.row != 'III.II' else u.name, u.drop, u.keep)
    if u.row in ('II.I', 'III.II'): key = ('KT' + u.name[2:], u.drop, u.keep)
    groups.setdefault(key, []).append(u)
for key, us in groups.items():
    rows = [u.row for u in us]
    vals = []; parts = []
    for seed in range(3):
        col = lib.Colour(RandField(seed)); rng = np.random.default_rng(100 + seed)
        P = {p: rng.uniform(-1, 1, 2) for p in lines.G3 + ['w\'']}
        tot = 0; mx = 0
        for u in us:
            Q = u.full(P); c = u.coef * u.colour(col, Q); tot += c; mx = max(mx, abs(c))
            if u.colour2 is not None:
                c2 = u.coef * u.colour2(col, Q); mx = max(mx, abs(c2))
        vals.append(tot); parts.append(mx)
    cc = ' (+c.c. for III.II only)' if 'III.II' in rows and 'II.I' in rows else ''
    print('%-10s %-8s rows %-22s  |sum of colour factors| = %s   (largest single: %s)%s' % (
        key[0], '%s->%s' % (key[1], key[2]), '+'.join(rows), ' '.join('%.3f' % abs(v) for v in vals),
        ' '.join('%.2f' % p for p in parts), cc))
# symbolic sums of the printed colour factors for the classes with more than one row
print('\nsymbolic sums (printed colour factors, numerical factors included):')
for key, us in groups.items():
    if len(us) < 2: continue
    terms = {}
    for u in us:
        for c, kind, args in u.colour.terms:
            a = tuple(args)
            if kind == 'F': a = tuple(sorted(a))
            if kind == 'D': a = tuple(sorted(a))
            terms[(kind, a)] = terms.get((kind, a), 0) + u.coef * c
    s = ' '.join('%+g %s(%s)' % (v / (lib.NC if k[0] == 'D' else 1), ('N' if k[0] == 'D' else '') + k[0], ','.join(k[1])) for k, v in terms.items() if abs(v) > 1e-12)
    print('  %-8s %-8s %s' % (key[0], '%s->%s' % key[1:], s or '0'))
