r"""Correct index sums for expressions in which some terms carry transverse indices (i, j, ...) and others do not.
texeval.term_value sums the WHOLE expression over every index assignment, so a term without a letter L is counted twice
for each such L.  Here every letter L gets a weight lam_L on each of its occurrences (COMP and DEL); a term containing L
has lam_L^2, a term without it lam_L^0.  The correct value is  prod_L (E_L[lam_L = 1] - E_L[lam_L = 0] / 2)  applied to
the (index-summed) expression, i.e. a sum over subsets S of letters with weight (-1/2)^{|S|} and lam = 0 on S."""
import itertools
import re
import numpy as np
import texeval as TE


def letters_of(code):
    a = set(re.findall(r'COMP\(D\("[^"]*","[^"]*"\),"([^"]*)"\)', code))
    b = set(sum((list(m) for m in re.findall(r'DEL\("([^"]*)","([^"]*)"\)', code)), []))
    return sorted(a | b)


def eval_code(code, names_of, post=None):
    """names_of(ix) -> dict of names for eval; post(v) -> number (e.g. colour contraction).  Returns the correctly
    index-summed value."""
    L = letters_of(code)
    tot = 0.0
    for S in itertools.product((0, 1), repeat=len(L)):          # S[k] = 1: letter k set to lam = 0
        lam = {l: (0.0 if s else 1.0) for l, s in zip(L, S)}
        w = (-0.5) ** sum(S)
        part = 0.0
        for vals in itertools.product((0, 1), repeat=len(L)):
            ix = dict(zip(L, vals))
            d = names_of(ix)
            d['COMP'] = lambda v, n, ix=ix: lam[n] * float(v[ix[n]])
            d['DEL'] = lambda a, b, ix=ix: lam[a] * lam[b] * (1.0 if ix[a] == ix[b] else 0.0)
            v = eval(code, {'np': np}, d)
            part += post(v) if post else v
        tot += w * part
    return tot


def term_value(tex, env):
    code = TE.translate(tex)
    return eval_code(code, env.names, post=lambda v: TE.scalar_or_col(v, env.cfg))
