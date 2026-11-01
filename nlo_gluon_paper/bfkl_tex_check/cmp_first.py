r"""The first forms (the log(vee/k+) terms of each row) of your BFKL.tex against the log(vee/k+) lists of
../full_check/Checks_JIMWLK_BFKL.tex (eq:B:...): text after normalization, and numerically (texeval, pointwise; three-rho
terms symmetrized over the charges, with c.c. for Groups I, III; log(vee/k+) = 1, epsilon = 0.37).  -> cmp_first.json"""
import json, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'full_check'))
import ll  # noqa
import three_check as T3
import texeval as TE
from rows_user import rows

CHK = open(os.path.join(HERE, '..', 'full_check', 'Checks_JIMWLK_BFKL.tex')).read()
MAP = {106: 'G1II2a', 126: 'G1IIIN1', 133: 'G1IIIN2', 188: 'G1IV', 247: 'G2I1A', 253: 'G2I1B', 295: 'G2I2b',
       300: 'G2I2c', 304: 'G2I2d', 309: 'G2I2e', 313: 'G2I2f', 345: 'G2IV', 359: 'G3I', 408: 'G3II', 462: 'G3III1',
       489: 'G3III2', 523: 'G3V', 549: 'G3VI'}
COLS = ('blue', 'red', 'green', 'orange', 'black', 'violet', 'gray', 'pink', 'purple', 'yellow')


def ours(lab):
    i = CHK.index(r'\label{eq:B:%s}' % lab)
    a = CHK.rindex(r'\begin{align}', 0, i)
    t = CHK[a + 13:CHK.index(r'\end{align}', i)]
    t = re.sub(r'\\label\{[^}]*\}', '', t).replace(r'\nc', '')
    return t[re.search(r'&\s*=', t).end():]


def norm(t):
    t = TE.unwrap_colours(t, COLS)
    t = re.sub(r'\\notag|\\nonumber|\\\\|&|\\,|\;|\\!|\s', '', t)
    return t.replace(r'\Bar', r'\bar').rstrip('.')


def prep(t):
    return T3.z_form(TE.unwrap_colours(t, COLS).replace(r'\Bar', r'\bar'))


def value(t, c):
    d = T3.degree(t, T3.configs(1))
    if d == 3:
        try:
            return T3.sym_value(t, c, False)
        except AssertionError:       # several phases (y' and w' forms): same labels on both sides, sum the placements
            import itertools
            return sum(TE.term_value(t, T3.env(T3.permuted(c, p))) for p in itertools.permutations(T3.CH))
    return TE.term_value(t, T3.env(c))


def plus(t):
    """a term that starts a line with no sign (G2 IV, third line) is read with +"""
    return re.sub(r'(\\\s*&\s*(?:\\nc)?\s*)(\\frac\{1\}\{\(2\\pi\)\^\{3\}\})', r'\1+\2', t)


if __name__ == '__main__':
    U = {r[3]: r[4] for r in rows('BFKL.tex')}
    out = {}
    for ln, lab in MAP.items():
        u, o = U[ln], ours(lab)
        same_text = norm(u) == norm(o)
        res = dict(ours='eq:B:' + lab, same_text=same_text)
        try:
            r = []
            for s in (1, 2):
                c = T3.configs(s)
                a, b = value(prep(plus(u)), c), value(prep(plus(o)), c)
                r.append(abs(a - b) / max(abs(b), 1e-300))
            res['rel'] = float(max(r))
        except Exception as e:  # noqa
            res['error'] = str(e)[:200]
        out[ln] = res
        print(ln, res)
    json.dump(out, open(os.path.join(HERE, 'cmp_first.json'), 'w'), indent=1)
