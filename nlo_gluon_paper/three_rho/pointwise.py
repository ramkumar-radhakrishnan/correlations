"""Pointwise numerical test (no symbolic simplification): sum of all symmetrized three-rho terms
at random points w, w', z, p1, p2, p3, random Wilson lines and random classical charges."""
import itertools
import numpy as np
import sympy as sp
import numerics as nm
import classes4 as C4
import classes6 as C6

NC = 3
def KKv(X, Y): return np.dot(X, Y) / (np.dot(X, X) * np.dot(Y, Y))

def mono_value(facs, P):
    v = 1.0
    for f in facs:
        if f[0] == 'KK':
            p, q, r, s = f[1:]
            v *= KKv(P[p] - P[q], P[r] - P[s])
        else:  # Phi(x) of Row III with native z, w
            x = P[f[2][0]]; z, w = P['z'], P['w']
            a, b = np.dot(x - w, x - w), np.dot(x - z, x - z)
            v *= (2 * (np.dot(z - w, x - z) * a - np.dot(z - w, x - w) * b) / (a - b) ** 2 - 1) / np.dot(z - w, z - w)
    return v

def total(terms, X, Um, rv, f):
    tot, scale = 0.0, 0.0
    for t in terms:
        rl = []
        for _, p in t['word']:
            if p not in rl: rl.append(p)
        c = complex(sp.N(t['coef']))
        for perm in itertools.permutations(['p1', 'p2', 'p3']):
            m = dict(t['cmap']); m.update(dict(zip(rl, perm)))
            P = {lab: X[m[lab]] for lab in m}
            k = sum(float(kc) * mono_value(fs, P) for kc, fs in t['monos'])
            facs = [('U', u[1], u[2], m.get(u[3], u[3])) if u[0] == 'U' else u for u in t['colour']]
            facs += [('rho', i, m[p]) for i, p in t['word']]
            v = c * k * nm.eval_colour(facs, f, Um, rv, NC * NC - 1) / 6
            tot += v; scale += abs(v)
    return tot, scale

if __name__ == '__main__':
    ts, f = nm.structure_constants(NC)
    sets = {'as written': C4.notes_terms(),
            'corrected': C4.notes_terms() + C4.fix_terms() + C6.row3_diag_terms()}
    for sym in (True, False):
        rng = np.random.default_rng(101)
        print('symmetric U (your convention)' if sym else 'general adjoint U')
        for trial in range(4):
            labs = ['w', "w'", 'z', 'p1', 'p2', 'p3']
            X = {l: rng.normal(size=2) for l in labs}
            Um = {l: nm.random_adjoint(NC, ts, rng, sym) for l in labs}
            rv = {l: rng.normal(size=NC * NC - 1) for l in labs}
            out = []
            for name, terms in sets.items():
                t, s = total(terms, X, Um, rv, f)
                out.append(f'{name}: |sum|/sum|terms| = {abs(t)/s:.1e}')
            print('   ', '   '.join(out))
