r"""Poles (1/eps) in the BFKL list, pointwise in the frame of ../llsplit/pole_check.py (charges symmetrized, c.c. of
Group I).

 1. The two red lines of BFKL.tex with a pole (Group I Row V 2b UV, Group II Row IV line 1): their poles cancel each
    other, but their sum is not zero (log(vee/k+) times log[(w-w')^4/((x-w)^2 (x'-w')^2)] times the LO structure).
 2. The exact BFKL list (BFKL.tex without the red lines + the log(vee/k+) parts of Other_terms.tex): every term with a
    pole, and the pole of their sum.  The only one is the log(vee/k+) part of line 2 of Group II Row IV.
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
    """every term of the exact BFKL list that contains epsilon"""
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
    e4 = full(UB['G2 IV'], UB['G2 IV'].terms[0])
    items = pole_items()
    out = dict(pole_items=[lab for lab, _, _ in items])
    for seed in (21, 22, 23):
        c = TE.cfg(seed)
        a = coeffs(uv, c, True)
        b = coeffs(e4, c, False)
        P = {lab: coeffs(t, c, cc) for lab, t, cc in items}
        tot = sum(P.values())
        out[str(seed)] = dict(pole_G1V=complex(a[0]).real, pole_G2IV=complex(b[0]).real,
                              pole_sum_relative=float(abs(a[0] + b[0]) / max(abs(a[0]), abs(b[0]))),
                              finite_sum_relative=float(abs(a[1] + b[1]) / max(abs(a[1]), abs(b[1]))),
                              list_poles={lab: complex(p[0]).real for lab, p in P.items()},
                              list_pole_total=complex(tot[0]).real,
                              list_pole_over_red_G1V=float(np.real(tot[0] / a[0])))
        print(seed, out[str(seed)])
    json.dump(out, open(os.path.join(R.HERE, 'pole_bfkl.json'), 'w'), indent=1)
