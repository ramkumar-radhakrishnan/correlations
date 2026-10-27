r"""The two terms of BFKL.tex with a pole (Group I Row V 2b UV, Group II Row IV line 1): their poles cancel each other,
but their sum is not zero (it is log(vee/k+) times log[(w-w')^4/((x-w)^2 (x'-w')^2)] times the LO structure).
Pointwise in the frame of ../llsplit/pole_check.py (charges symmetrized, c.c. of Group I).  Output: pole_bfkl.json"""
import json
import os
import numpy as np
import rows as R
import check_rows as CR
import pole_check as PC
import texeval as TE

UB = CR.UB


def coeffs(tex, c, cc):
    es = (0.05, 0.1, 0.2)
    M = np.array([[1 / e, 1.0, e] for e in es])
    y = np.array([PC.value(tex, c, cc, e) for e in es])
    return np.linalg.solve(M, y)


if __name__ == '__main__':
    uv = UB['G1 V(2b UV)'].prefix + ' (' + UB['G1 V(2b UV)'].terms[0] + ') ' + UB['G1 V(2b UV)'].suffix
    e4 = UB['G2 IV'].prefix + ' (' + UB['G2 IV'].terms[0] + ') ' + UB['G2 IV'].suffix
    out = {}
    for seed in (21, 22, 23):
        c = TE.cfg(seed)
        a = coeffs(uv, c, True)
        b = coeffs(e4, c, False)
        out[str(seed)] = dict(pole_G1V=complex(a[0]).real, pole_G2IV=complex(b[0]).real,
                              pole_sum_relative=float(abs(a[0] + b[0]) / max(abs(a[0]), abs(b[0]))),
                              finite_sum_relative=float(abs(a[1] + b[1]) / max(abs(a[1]), abs(b[1]))))
        print(seed, out[str(seed)])
    json.dump(out, open(os.path.join(R.HERE, 'pole_bfkl.json'), 'w'), indent=1)
