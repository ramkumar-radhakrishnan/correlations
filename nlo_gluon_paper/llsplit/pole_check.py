r"""The 1/epsilon poles of the log(k+/Lambda) terms (Group I Row II part 1, Group I Row V 2b UV, Group II Row IV) cancel?
Pointwise, after renaming onto the frame w, w' (phase e^{-ik(w'-w)}), charges x, x' (symmetrized), with the c.c. of
Group I.  The pole coefficient a of a/eps + b + c eps is obtained from three values of eps.
Run with the Row V 2b UV coefficient as corrected (2 log(k+/Lambda), E1) and as written (1).  Output: pole_check.json"""
import itertools
import json
import os
import re
import numpy as np
import gen as GN
import texeval as TE
import three_check as T3

HERE = os.path.dirname(os.path.abspath(__file__))
E = np.e
CH2 = ('x', 'xp')
NAMES = ['G1 II(1)', 'G1 V(2b UV)', 'G2 IV']


def frame2(tex):
    code = TE.translate(tex)
    a, p = re.findall(r'KD\("([^"]*)","([^"]*)"\)', code)[0]
    charges = []
    for q in re.findall(r'CR\("[^"]*","([^"]*)"\)', code):
        if q not in charges:
            charges.append(q)
    assert len(charges) == 2, charges
    m = {a: 'wp', p: 'w'}
    m.update(dict(zip(charges, CH2)))
    return m


def perm2(c, perm):
    new = T3.permuted(c, T3.CH)
    for a, b in zip(CH2, perm):
        new.pt[a], new.U[a], new.rho[a] = c.pt[b], c.U[b], c.rho[b]
    return new


def value(tex, c, cc, eps):
    m = frame2(tex)
    v = 0.0
    for perm in itertools.permutations(CH2):
        cp = perm2(c, perm)
        v += TE.term_value(tex, TE.Env(T3.relabel(cp, m), KP=1.0, V=E, L=1.0 / E, EPS=eps))
        if cc:
            v += np.conj(TE.term_value(tex, TE.Env(T3.relabel(T3.swapped_w(cp), m), KP=1.0, V=E, L=1.0 / E, EPS=eps)))
    return v


def pole(tex, c, cc):
    es = (0.05, 0.1, 0.2)
    M = np.array([[1 / e, 1.0, e] for e in es])
    y = np.array([value(tex, c, cc, e) for e in es])
    return np.linalg.solve(M, y)[0]


if __name__ == '__main__':
    BLK, SPL = GN.all_splits()
    terms = {}
    for name in NAMES:
        B = BLK[name]
        n = 1 if name == 'G2 IV' else 0
        terms[name] = B.prefix + ' (' + SPL[name][n]['J'] + ') ' + B.suffix
    as_written = terms['G1 V(2b UV)'].replace(r'2\,' + GN.LKL, GN.LKL)
    out = {}
    for seed in (21, 22, 23):
        c = TE.cfg(seed)
        P = {name: pole(t, c, not name.startswith('G2')) for name, t in terms.items()}
        Pw = pole(as_written, c, True)
        scale = max(abs(v) for v in P.values())
        ratio = {k: float(np.real(v / P['G2 IV'])) for k, v in P.items()}
        out[str(seed)] = dict(poles_relative_to_G2_IV=ratio,
                              sum_corrected_relative=float(abs(sum(P.values())) / scale),
                              sum_as_written_relative=float(abs(P['G1 II(1)'] + Pw + P['G2 IV']) / scale))
        print(seed, out[str(seed)])
    json.dump(out, open(os.path.join(HERE, 'pole_check.json'), 'w'), indent=1)
