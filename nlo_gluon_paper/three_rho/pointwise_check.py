"""Pointwise numerical test, independent of the class/network bookkeeping: the sum of all
symmetrized three-rho terms at random points w, w', z, p1, p2, p3 with random Wilson lines and
random classical charges (summed over the 3! placements of the charges).  -> pointwise.json"""
import json
import numpy as np
import numerics as nm
import classes4 as C4
from pointwise import total, NC
from corrections import correction1_terms, correction2_terms, correction1_variant_terms


def main(trials=4):
    ts, f = nm.structure_constants(NC)
    sets = {'written': C4.notes_terms(),
            'corrected': C4.notes_terms() + correction1_terms() + correction2_terms(),
            'corrected_variant': C4.notes_terms() + correction1_variant_terms() + correction2_terms()}
    out = {}
    for sym in (True, False):
        rng = np.random.default_rng(101)
        rows = []
        for trial in range(trials):
            labs = ['w', "w'", 'z', 'p1', 'p2', 'p3']
            X = {l: rng.normal(size=2) for l in labs}
            Um = {l: nm.random_adjoint(NC, ts, rng, sym) for l in labs}
            rv = {l: rng.normal(size=NC * NC - 1) for l in labs}
            r = {}
            for name, terms in sets.items():
                t, s = total(terms, X, Um, rv, f)
                r[name] = (abs(t), s)
            rows.append(r)
        out['symmetric' if sym else 'general'] = rows
    json.dump(out, open('pointwise.json', 'w'), indent=1)
    return out


if __name__ == '__main__':
    o = main()
    for k, rows in o.items():
        print(k)
        for r in rows:
            print('   ', '   '.join(f'{n}: {a/s:.1e}' for n, (a, s) in r.items()))
