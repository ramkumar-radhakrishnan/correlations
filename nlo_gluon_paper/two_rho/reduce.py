"""Reduce the colour networks of every kernel structure to an independent basis with exact
coefficients (polynomials in Nc, fixed from SU(3) and SU(4) fingerprints)."""
import pickle
from fractions import Fraction
import numpy as np
import sympy as sp
import core  # noqa: F401  (paths)
import numerics as nm
from engine import Nc
from tables import LABS


def fingerprints(cols, N, M=24, seed=0):
    ts, f = nm.structure_constants(N)
    rng = np.random.default_rng(seed + 100 * N)
    out = np.zeros((len(cols), M), complex)
    for m in range(M):
        U = {l: nm.random_adjoint(N, ts, rng, True) for l in LABS}
        r = {l: rng.normal(size=N * N - 1) for l in LABS}
        for k, col in enumerate(cols):
            out[k, m] = nm.eval_colour(col, f, U, r, N * N - 1)
    return out


def rat(x, tol=1e-7):
    fr = Fraction(x).limit_denominator(200)
    assert abs(float(fr) - x) < tol, x
    return sp.Rational(fr.numerator, fr.denominator)


def reduce_structure(nets):
    cols = [n['col'] for n in nets]
    F3, F4 = fingerprints(cols, 3), fingerprints(cols, 4)
    prio = sorted(range(len(nets)), key=lambda k: (0 if '-JIMWLK' in nets[k]['groups'] else 1, len(cols[k]), k))
    basis = []
    for k in prio:
        trial = basis + [k]
        for F in (F3, F4):
            s = np.linalg.svd(F[trial], compute_uv=False)
            if s[-1] < 1e-8 * s[0]:
                break
        else:
            basis = trial
    coeffs = []
    for k in range(len(nets)):
        a3 = np.linalg.lstsq(F3[basis].T, F3[k], rcond=None)[0]
        a4 = np.linalg.lstsq(F4[basis].T, F4[k], rcond=None)[0]
        assert np.abs(F3[basis].T @ a3 - F3[k]).max() < 1e-8 * (1 + np.abs(F3[k]).max())
        assert np.abs(F4[basis].T @ a4 - F4[k]).max() < 1e-8 * (1 + np.abs(F4[k]).max())
        row = []
        for x3, x4 in zip(a3, a4):
            assert abs(x3.imag) < 1e-8 and abs(x4.imag) < 1e-8
            r3, r4 = rat(x3.real), rat(x4.real)
            beta = r4 - r3
            alpha = r3 - 3 * beta
            row.append(sp.expand(alpha + beta * Nc))
        coeffs.append(row)
    groups = sorted({g for n in nets for g in n['groups']}, key=str)
    table = []
    for b, kb in enumerate(basis):
        contrib = {}
        for g in groups:
            v = sum(sp.sympify(n['groups'].get(g, 0)) * coeffs[k][b] for k, n in enumerate(nets))
            v = sp.factor(sp.expand(v))
            if v != 0:
                contrib[g] = v
        table.append(dict(col=cols[kb], contrib=contrib, total=sp.simplify(sum(contrib.values()))))
    relations = [(k, [(cols[kb], coeffs[k][b]) for b, kb in enumerate(basis) if coeffs[k][b] != 0])
                 for k in range(len(nets)) if k not in basis]
    return table, relations, [cols[k] for k in range(len(nets))]


if __name__ == '__main__':
    from tables import tex_kkey
    from engine import tex_colour
    Bk = pickle.load(open('book.pkl', 'rb'))
    out = {}
    nbad = 0
    for key, nets in Bk.items():
        table, rel, cols = reduce_structure(nets)
        out[key] = dict(table=table, relations=rel, cols=cols, nets=nets)
        print(tex_kkey(key), ': networks', len(nets), ' basis', len(table))
        for row in table:
            flag = 'ok' if row['total'] == 0 else 'NONZERO %s' % row['total']
            if row['total'] != 0:
                nbad += 1
            print('    ', flag, ' JIMWLK' if '-JIMWLK' in row['contrib'] else '', tex_colour(row['col'])[:90])
    print('basis elements with nonzero total:', nbad)
    pickle.dump(out, open('reduced.pkl', 'wb'))
