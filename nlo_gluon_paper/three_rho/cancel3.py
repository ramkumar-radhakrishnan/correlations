"""Do the symmetrized three-rho terms cancel?

Every three-rho term (four-rho remainders, Groups I-III, and the c.c. of
Groups I and III) is brought to the common form
    N int_{w,w'} e^{-ik(w'-w)} int_{z, p1, p2, p3} F(w, w', z; p1, p2, p3)
with classical (symmetrized) charges at p1, p2, p3 and the soft gluon at z.
The charges commute, so each term is summed over the 3! ways of placing
its three charges on p1, p2, p3 (divided by 6).  The terms cancel among each
other iff the sum vanishes pointwise (or after the z integration).
"""
import itertools
import numpy as np
import sympy as sp
from engine import I
from sources import SOURCES
import numerics as nm

NC = 3


def KKv(X, Y):
    return np.dot(X, Y) / (np.dot(X, X) * np.dot(Y, Y))


def Kv(X):
    return X / np.dot(X, X)


def kernel_value(src, P):
    """numerical kernel of a source; P maps native labels to 2-vectors."""
    tot = 0.0
    for c, facs in src.kexpl:
        val = float(c)
        for f in facs:
            if f[0] == 'KK':
                p, q, r, s = f[1:]
                val *= KKv(P[p] - P[q], P[r] - P[s])
            elif 'Phi' in f[1]:
                x = P[f[2][0]]
                z, w = P['z'], P['w']
                a, b = np.dot(x - w, x - w), np.dot(x - z, x - z)
                D = a - b
                val *= (2 * (np.dot(z - w, x - z) * a - np.dot(z - w, x - w) * b) / D**2 - 1) / np.dot(z - w, z - w)
            elif 'mathfrak F' in f[1]:
                xp, x, y = (P[l] for l in f[2])
                w, wp, z = P['w'], P["w'"], P['z']
                K = lambda a, b: KKv(a - z, b - z)
                Fi = Kv(x - w) * (K(w, y) - K(x, y)) - Kv(y - w) * (K(w, x) - K(x, y))
                val *= np.dot(Kv(xp - wp), Fi)
            else:
                raise ValueError(f)
        tot += val
    return tot


CANON = {
    'GI.V1': {"w'": "w'", 'z': 'w', 'x': 'z'},
    'GI.V2': {"w'": "w'", 'x': 'w', 'z': 'z'},
    'GIII.III': {"y'": "w'", 'w': 'w', 'z': 'z'},
    'GIII.V': {"y'": "w'", 'w': 'w', 'z': 'z'},
    'GIII.VI': {"y'": "w'", 'w': 'w', 'z': 'z'},
}


def rho_labels_of(src):
    return [p for _, p in src.subterms[0][2]]


def term_list():
    """(key, cc, coeff, src, sub, nativemap) for every symmetrized three-rho term."""
    out = []
    for src in SOURCES:
        mult = 1 if src.form == 'plain' else 2
        base = dict(CANON.get(src.key, {"w'": "w'", 'w': 'w', 'z': 'z'}))
        P1, P0 = src.phase
        for n, (sign, colour, word) in enumerate(src.subterms):
            c = complex(sp.N(src.pref)) * mult * sign
            out.append((src.key, False, n, c, src, colour, word, base))
            if src.cc:
                ccmap = dict(base)
                ccmap[P1], ccmap[P0] = base[P0], base[P1]
                out.append((src.key, True, n, np.conj(c), src, colour, word, ccmap))
    return out


def evaluate(terms, X, Um, rv, f, per_term=False):
    """X: canonical label -> vector; Um: canonical label -> U; rv: p1..p3 -> rho vector."""
    ncol = NC * NC - 1
    total = 0.0
    parts = []
    for key, cc, n, c, src, colour, word, cmap in terms:
        rl = [p for _, p in word]
        val = 0.0
        for perm in itertools.permutations(['p1', 'p2', 'p3']):
            m = dict(cmap)
            m.update(dict(zip(rl, perm)))
            P = {lab: X[m[lab]] for lab in m}
            k = kernel_value(src, P)
            facs = []
            for fc in colour:
                facs.append(('U', fc[1], fc[2], m[fc[3]]) if fc[0] == 'U' else fc)
            facs += [('rho', i, m[p]) for i, p in word]
            col = nm.eval_colour(facs, f, Um, rv, ncol)
            val += k * col
        val *= c / 6
        total += val
        parts.append(val)
    return (total, parts) if per_term else total


def random_config(rng, ts, symmetric=True):
    labs = ['w', "w'", 'z', 'p1', 'p2', 'p3']
    X = {l: rng.normal(size=2) for l in labs}
    Um = {l: nm.random_adjoint(NC, ts, rng, symmetric) for l in labs}
    rv = {l: rng.normal(size=NC * NC - 1) for l in labs}
    return X, Um, rv


if __name__ == '__main__':
    ts, f = nm.structure_constants(NC)
    rng = np.random.default_rng(3)
    terms = term_list()
    print('number of symmetrized three-rho terms (incl. c.c.):', len(terms))
    for trial in range(4):
        X, Um, rv = random_config(rng, ts)
        tot, parts = evaluate(terms, X, Um, rv, f, per_term=True)
        scale = sum(abs(p) for p in parts)
        print(f'trial {trial}: total = {tot:.6e}   sum |terms| = {scale:.6e}')
