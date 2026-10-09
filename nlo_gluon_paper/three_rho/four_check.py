"""(a) the four-rho terms of the notes cancel after symmetrization;
(b) their three-rho part (from reordering) equals the remainders A1, A2, B1, B2;
(c) row-by-row comparison of the full three-rho content with the reference."""
import itertools
import pickle
from collections import defaultdict
import numpy as np
import numerics as nm
import fourrho
import compare_ref as CR

NC = 3
ORD4 = ['w', "w'", 'z', 'p1', 'p2', 'p3', 'p4']


def kk4(a, b, c, d):
    def nv(p, q):
        return ((p, q), 1) if ORD4.index(p) < ORD4.index(q) else ((q, p), -1)
    v1, s1 = nv(a, b); v2, s2 = nv(c, d)
    if v1 == v2:
        return ('inv2', v1), s1 * s2
    return ('KK',) + tuple(sorted([v1, v2], key=lambda v: (ORD4.index(v[0]), ORD4.index(v[1])))), s1 * s2


def check_four(seed=3):
    ts, f = nm.structure_constants(NC)
    rng = np.random.default_rng(seed)
    Um = {l: nm.random_adjoint(NC, ts, rng, True) for l in ORD4}
    rv = {l: rng.normal(size=NC * NC - 1) for l in ORD4}
    acc = defaultdict(complex); scale = defaultdict(float)
    for c, row, lab, facs, word, kern, cc in fourrho.words():
        rl = [p for _, p in word]
        assert len(set(rl)) == 4
        for perm in itertools.permutations(['p1', 'p2', 'p3', 'p4']):
            m = dict(zip(rl, perm)); g = lambda p: m.get(p, p)
            parts, sgn = [], 1
            for (a, b), (cc_, d) in kern:
                k, s = kk4(g(a), g(b), g(cc_), g(d)); parts.append(k); sgn *= s
            key = tuple(sorted(parts, key=str))
            ff = [('U', u[1], u[2], g(u[3])) for u in facs] + [('rho', i, g(p)) for i, p in word]
            v = c * sgn * nm.eval_colour(ff, f, Um, rv, NC * NC - 1) / 24
            acc[key] += v; scale[key] += abs(v)
    return max(abs(v) for v in acc.values()), max(scale.values())


def pbw3_terms():
    out = []
    for c, row, lab, facs, word, kern, cc in fourrho.words():
        for i, j in itertools.combinations(range(4), 2):
            (xi, pi), (yi, pj) = word[i], word[j]
            g = 'g%d%d' % (i, j)
            sub = lambda p: pi if p == pj else p
            rest = [word[k] for k in range(4) if k not in (i, j)]
            newr = [(g, pi)] + [(a, sub(p)) for a, p in rest]
            newc = [('U', u[1], u[2], sub(u[3])) for u in facs] + [('f', xi, yi, g)]
            newk = [((sub(a), sub(b)), (sub(cc_), sub(d))) for (a, b), (cc_, d) in kern]
            out.append(CR.C3(c * 0.5j, newc, newr, newk, '4:' + row + ':' + lab))
    return out


if __name__ == '__main__':
    r, s = check_four()
    print(f'(a) four-rho terms of the notes, symmetrized: max |sum| per kernel monomial = {r:.2e} (scale {s:.2e})')
    A = CR.run()
    for t in pbw3_terms():
        A.add_ref(t, prefix='N4 ')
    keys = list(A.data.keys())
    def vec(pred):
        return np.array([sum(v for n, v in A.data[k].items() if pred(n)) for k in keys])
    ab = vec(lambda n: n.startswith('USR ') and n.split()[1][0] in 'AB')
    n4 = vec(lambda n: n.startswith('N4 '))
    print(f'(b) S3(A1+A2+B1+B2) vs three-rho part of the four-rho words: |A+B| = {np.abs(ab).max():.4f}, '
          f'|N4| = {np.abs(n4).max():.4f}, |difference| = {np.abs(ab-n4).max():.2e}')
    print('(c) full three-rho content per row: notes (genuine + four-rho part) vs reference (three-rho words + four-rho part)')
    rows = [
        ('Group I',      lambda n: n.startswith('USR GI.') or n.startswith('N4 4:G1:'), lambda n: n.startswith('REF 3:G1') or n.startswith('REF 4:G1')),
        ('G2-I',         lambda n: n.startswith('USR GII.') or n.startswith('N4 4:G2-I:'), lambda n: n in ('REF 3:G2-I', 'REF 4:G2-I')),
        ('G2-II',        lambda n: n.startswith('N4 4:G2-II:'), lambda n: n in ('REF 3:G2-II', 'REF 4:G2-II')),
        ('G2-III',       lambda n: n.startswith('N4 4:G2-III:'), lambda n: n in ('REF 3:G2-III', 'REF 4:G2-III')),
        ('G3-I',         lambda n: n.rstrip('*') == 'USR GIII.I' or n.startswith('N4 4:G3-I:'), lambda n: n in ('REF 3:G3-I', 'REF 4:G3-I')),
        ('G3-II',        lambda n: n.rstrip('*') == 'USR GIII.II' or n.startswith('N4 4:G3-II:'), lambda n: n in ('REF 3:G3-II', 'REF 4:G3-II')),
        ('G3-III',       lambda n: n.rstrip('*') == 'USR GIII.III', lambda n: n in ('REF 3:G3-III', 'REF 4:G3-III')),
        ('G3-IV',        lambda n: n.startswith('N4 4:G3-IV:'), lambda n: n in ('REF 3:G3-IV', 'REF 4:G3-IV')),
        ('G3-V',         lambda n: n.rstrip('*') == 'USR GIII.V', lambda n: n in ('REF 3:G3-V', 'REF 4:G3-V')),
        ('G3-VI',        lambda n: n.rstrip('*') == 'USR GIII.VI', lambda n: n in ('REF 3:G3-VI', 'REF 4:G3-VI')),
    ]
    for name, pu, pr in rows:
        u, r = vec(pu), vec(pr)
        print(f'   {name:10s} |notes| = {np.abs(u).max():8.4f}   |ref| = {np.abs(r).max():8.4f}   |notes - ref| = {np.abs(u - r).max():.3e}')
    tot_u = vec(lambda n: n.startswith('USR ') and n.split()[1][0] not in 'AB') + n4
    print(f'   total notes (with four-rho part from the words) = {np.abs(tot_u).max():.4f}')
    pickle.dump({'keys': keys, 'data': {k: dict(v) for k, v in A.data.items()}}, open('cmp_data.pkl', 'wb'))
