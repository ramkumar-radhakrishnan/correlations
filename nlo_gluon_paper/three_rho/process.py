"""Decompose every three-rho term into two-rho pieces and check the colour algebra."""
import pickle
import numpy as np
import sympy as sp
from engine import (Rho, pieces_for, make_piece, simplify, canonical_rename, check_contracted, Nc)
from sources import SOURCES
import numerics as nm

NC = 3


def process():
    results = []
    for src in SOURCES:
        sres = {'src': src, 'subs': []}
        for n, (sign, colour, word) in enumerate(src.subterms):
            facs_full = colour + [('rho', i, p) for i, p in word]
            check_contracted(facs_full, f'{src.key} sub {n}')
            rhos = [Rho(i, p) for i, p in word]
            pcs = []
            for coeff, X, Y, Z, tag in pieces_for(src.form, rhos):
                c0, facs0, (old, new) = make_piece(colour, sign * coeff, X, Y, Z)
                check_contracted(facs0, f'{src.key} sub {n} piece {tag}')
                c1, facs1, rules = simplify(c0, facs0)
                facs1c = canonical_rename(facs1) if facs1 else facs1
                pcs.append(dict(tag=tag, X=(X.idx, X.pos), Y=(Y.idx, Y.pos), Z=(Z.idx, Z.pos),
                                coeff_before=c0, facs_before=facs0, subst=(old, new),
                                coeff_after=c1, facs_after=facs1c, rules=rules))
            sres['subs'].append(dict(sign=sign, colour=colour, word=word, pieces=pcs))
        results.append(sres)
    return results


def colour_checks(results, nsamp=3, seed=1):
    ts, f = nm.structure_constants(NC)
    ncol = NC * NC - 1
    rng = np.random.default_rng(seed)
    labels = ["x'", 'x', 'y', 'z', 'w', "w'", "y'"]
    flags = []
    worst = 0.0
    for sym in (True, False):
        for _ in range(nsamp):
            Um = {p: nm.random_adjoint(NC, ts, rng, symmetric=sym) for p in labels}
            rv = {p: rng.normal(size=ncol) for p in labels}
            for sres in results:
                for n, sub in enumerate(sres['subs']):
                    for pc in sub['pieces']:
                        vb = complex(sp.N(pc['coeff_before'].subs(Nc, NC))) * nm.eval_colour(pc['facs_before'], f, Um, rv, ncol)
                        if pc['facs_after'] or pc['coeff_after'] != 0:
                            va = complex(sp.N(pc['coeff_after'].subs(Nc, NC))) * (nm.eval_colour(pc['facs_after'], f, Um, rv, ncol) if pc['facs_after'] else 1.0)
                        else:
                            va = 0.0
                        err = abs(vb - va) / (1 + abs(vb))
                        if sym:
                            worst = max(worst, err)
                        elif err > 1e-9:
                            flags.append((sres['src'].key, n, pc['tag'], pc['rules']))
    return worst, sorted(set((a, b, c, tuple(d)) for a, b, c, d in flags))


if __name__ == '__main__':
    res = process()
    worst, flags = colour_checks(res)
    print('symmetric-U worst relative error:', worst)
    print('pieces that need U^{bc}=U^{cb} (differ for general U):', len(flags))
    for fl in flags:
        print('   ', fl)
    npieces = sum(len(s['pieces']) for r in res for s in r['subs'])
    nzero = sum(1 for r in res for s in r['subs'] for p in s['pieces'] if p['coeff_after'] == 0)
    print('pieces:', npieces, 'vanishing:', nzero)
    from collections import Counter
    print(Counter(tuple(p['rules']) for r in res for s in r['subs'] for p in s['pieces']))
    pickle.dump(res, open('results.pkl', 'wb'))
