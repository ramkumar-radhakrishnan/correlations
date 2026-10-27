r"""Your two-rho four-rho lists of BFKL.tex (d^3N_1 'From 4 rho's' and d^3N_2) against the two four-rho lists that entered
the check of the BFKL sum (../final_check/Final_evolution_check.tex, Eqs. P12 and P13).

Pointwise, after summing over the two placements of the charges x, x' and adding the c.c. (conjugate, w <-> w').
Output: compare_4rho.json"""
import itertools
import json
import os
import re
import numpy as np
import rows as R
import texeval as TE
import three_check as T3

HERE = R.HERE
PRE = r'\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\log{\left(\frac{\vee}{k^{+}}\right)}'


def user(name):
    a, b = R.B_RANGES[name]
    t = R.align_text('BFKL', a, b)
    return t[re.search(r'&\s*=', t).end():]


def fec(label):
    S = open(os.path.join(HERE, '..', 'final_check', 'Final_evolution_check.tex')).read()
    j = S.index(r'\label{%s}' % label)
    i = S.rindex(r'\begin{align}', 0, j)
    t = S[i + len(r'\begin{align}'):j]
    t = t[t.index('={}&') + 4:]
    t = t.replace(r'\NNv', PRE).replace(r'\tfrac', r'\frac').replace('N_c', 'N_{c}')
    return TE.unwrap_colours(t, ('CC', 'JIM', 'BFK'))


def perm2(c, perm):
    new = T3.permuted(c, T3.CH)
    for a, b in zip(('x', 'xp'), perm):
        new.pt[a], new.U[a], new.rho[a] = c.pt[b], c.U[b], c.rho[b]
    return new


def sym(tex, c):
    v = 0.0
    for perm in itertools.permutations(('x', 'xp')):
        cp = perm2(c, perm)
        env = TE.Env(cp, KP=1.0, V=np.e, L=1e-9)
        v += TE.term_value(tex, env)
        v += np.conj(TE.term_value(tex, TE.Env(T3.swapped_w(cp), KP=1.0, V=np.e, L=1e-9)))
    return v


if __name__ == '__main__':
    mine = [fec('eq:P12'), fec('eq:P13')]
    yours = [user('4rho 1'), user('4rho 2')]
    out = {}
    for seed in (1, 2, 3):
        c = TE.cfg(seed)
        a = sum(sym(t, c) for t in yours)
        b = sum(sym(t, c) for t in mine)
        sep = [sym(t, c) for t in yours] + [sym(t, c) for t in mine]
        out[str(seed)] = dict(yours=[a.real, a.imag], checked=[b.real, b.imag], relative=float(abs(a - b) / max(abs(a), abs(b))),
                              separate=[[float(v.real), float(v.imag)] for v in sep])
        print(seed, out[str(seed)]['relative'], a, b)
    json.dump(out, open(os.path.join(HERE, 'compare_4rho.json'), 'w'), indent=1)
