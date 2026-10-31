r"""Poles (1/eps) in the BFKL list, pointwise in the frame of ../llsplit/pole_check.py (charges symmetrized, c.c. of
Group I; Group II without c.c.).

The terms of the list with a pole are your Group II Row IV line 1 (-f(0) log(vee/k+), from the subtraction of
xi/(1-xi)_+) and the log(vee/k+) part of line 2 of Other_terms.tex: their poles (and finite parts) cancel.
The removed Group I Row V 2b UV line (E1) has a pole of its own, which nothing in the list cancels.
Output: pole_bfkl.json"""
import json
import os
import numpy as np
import rows as R
import check_rows as CR
import pole_check as PC
import texeval as TE

UB, UO, EX = CR.UB, CR.UO, CR.EX


def coeffs(tex, c, cc):
    es = (0.05, 0.1, 0.2)
    M = np.array([[1 / e, 1.0, e] for e in es])
    y = np.array([PC.value(tex, c, cc, e) for e in es])
    return np.linalg.solve(M, y)


def full(B, t):
    return B.prefix + ' (' + t + ') ' + B.suffix


def pole_items():
    """every term of the BFKL list that contains epsilon"""
    out = []
    for name, B in UB.items():
        for n, t in enumerate(B.terms):
            if (name, n) not in CR.REMOVED and r'\epsilon' in full(B, t):
                out.append(('BFKL.tex %s #%d' % (name, n + 1), full(B, t), not name.startswith('G2')))
    for name, L in EX.items():
        for n, s in enumerate(L):
            if s.get('B') is not None and r'\epsilon' in full(UO[name], s['B']):
                out.append(('extracted %s #%d' % (name, n + 1), full(UO[name], s['B']), not name.startswith('G2')))
    return out


if __name__ == '__main__':
    uv = full(UB['G1 V(2b UV)'], UB['G1 V(2b UV)'].terms[0])
    items = pole_items()
    out = dict(pole_items=[lab for lab, _, _ in items])
    for seed in (21, 22, 23):
        c = TE.cfg(seed)
        a = coeffs(uv, c, True)
        P = {lab: coeffs(t, c, cc) for lab, t, cc in items}
        tot = sum(P.values())
        big = max(abs(p[0]) for p in P.values())
        bigf = max(abs(p[1]) for p in P.values())
        out[str(seed)] = dict(pole_removed_G1V=complex(a[0]).real,
                              list_poles={lab: complex(p[0]).real for lab, p in P.items()},
                              list_pole_relative=float(abs(tot[0]) / big),
                              list_finite_relative=float(abs(tot[1]) / bigf))
        print(seed, out[str(seed)])
    json.dump(out, open(os.path.join(R.HERE, 'pole_bfkl.json'), 'w'), indent=1)
