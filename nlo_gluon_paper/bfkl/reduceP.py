"""Basis reduction of the colour networks (exact coefficients, polynomials in Nc from SU(3), SU(4))."""
import env  # noqa
import pickle, sys
import sympy as sp
import reduce as RD          # two_rho/reduce.py
from book import tex_kkey
from engine import tex_colour


def reduce_structure(nets, prio_groups=('-BFKL', '+JIMWLK')):
    # same as two_rho/reduce.py but with a different priority for the basis
    import numpy as np
    cols = [n['col'] for n in nets]
    F3, F4 = RD.fingerprints(cols, 3), RD.fingerprints(cols, 4)
    prio = sorted(range(len(nets)), key=lambda k: (0 if any(g in nets[k]['groups'] for g in prio_groups) else 1, len(cols[k]), k))
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
            r3, r4 = RD.rat(x3.real), RD.rat(x4.real)
            beta = r4 - r3
            alpha = r3 - 3 * beta
            row.append(sp.expand(alpha + beta * RD.Nc))
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
    return table, relations, cols


if __name__ == '__main__':
    inp, outp = sys.argv[1], sys.argv[2]
    Bk = pickle.load(open(inp, 'rb'))
    out, nbad = {}, 0
    for key, nets in Bk.items():
        table, rel, cols = reduce_structure(nets)
        out[key] = dict(table=table, relations=rel, cols=cols, nets=nets)
        nb = sum(1 for r in table if r['total'] != 0)
        nbad += nb
        print('%-62s networks %3d basis %3d nonzero %d' % (tex_kkey(key)[:62], len(nets), len(table), nb))
    print('basis elements with nonzero total:', nbad)
    pickle.dump(out, open(outp, 'wb'))
