r"""Two-rho remainder of a three-rho term, pointwise.
With the charges as operators, a product of three charges A B C is its symmetric (Weyl) part plus
   (1/2)([A,B] C + [B,C] A + [A,C] B)  + one-rho terms,
i.e. every pair contributes (1/2)[earlier, later];  [rho^i(p), rho^j(q)] = kappa f^{ijd} rho^d(p) delta(p-q).
The delta puts the two points of the pair on top of each other: the remainder is a two-rho function of the pair point and
the third point.  It is evaluated with the pair point and the third point at P1 = x and P2 = x' of the configuration and
summed over the two placements, i.e. exactly as a two-rho list with charges at x, x' (symmetrized).  {A,B} is kept as
AB + BA (ordered).  The c.c. of a term is conj(...) with w <-> w' (the order of the charges is reversed, kappa -> kappa*)."""
import itertools
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'full_check'))
import ll  # noqa
import three_check as T3
import texeval as TE
import texfix as TF
TE.term_value = TF.term_value          # correct index sums everywhere (see texfix.py)

CHG = ('x', 'y', 'xp')


class OEnv(TE.Env):
    def names(self, ix):
        d = super().names(ix)
        d['ACOM'] = lambda a, b: a * b + b * a
        return d


def colour_eval(code, env, ix):
    return eval(code, {'np': np}, env.names(ix))


def place(c, p, q, r, Pa, Pb):
    """pair (p, q) at the point Pa (with U, rho of Pa), third charge r at Pb"""
    new = T3.permuted(c, CHG)
    for n in (p, q):
        new.pt[n], new.U[n], new.rho[n] = c.pt[Pa], c.U[Pa], c.rho[Pa]
    new.pt[r], new.U[r], new.rho[r] = c.pt[Pb], c.U[Pb], c.rho[Pb]
    return new


def commute(col, p, q, kappa):
    out = []
    for coef, facs in col.m:
        idx = [k for k, f in enumerate(facs) if f[0] == 'r' and f[2] in (p, q)]
        assert len(idx) == 2, facs
        k1, k2 = idx
        i, j = facs[k1][1], facs[k2][1]
        d = '@%d' % len(out)
        new = [f for k, f in enumerate(facs) if k not in idx] + [('f', i, j, d), ('r', d, p)]
        out.append((coef * kappa / 2, new))
    return TE.Col(out)


def remainder(tex, c, kappa=1j):
    code = TE.translate(tex)
    letters = sorted(set(__import__('re').findall(r'COMP\(D\("[^"]*","[^"]*"\),"([^"]*)"\)', code)))
    charges = sorted(set(__import__('re').findall(r'CR\("[^"]*","([^"]*)"\)', code)))
    assert set(charges) == set(CHG), charges
    tot = 0.0
    for p, q in itertools.combinations(CHG, 2):
        r = [n for n in CHG if n not in (p, q)][0]
        for Pa, Pb in (('x', 'xp'), ('xp', 'x')):
            cfg = place(c, p, q, r, Pa, Pb)
            env = OEnv(cfg, KP=1.0, V=np.e, L=1.0 / np.e)
            # index sums done term by term (texfix: terms without a letter are not counted twice)
            tot += TF.eval_code(code, env.names, post=lambda v, cfg=cfg, p=p, q=q: commute(v, p, q, kappa).value(cfg))
    return tot


def remainder_cc(tex, c, cc, kappa=1j):
    v = remainder(tex, c, kappa)
    if cc:
        v += np.conj(remainder(tex, T3.swapped_w(c), kappa))
    return v


def two_value(tex, c, cc=False):
    """a two-rho list with charges at x, x': summed over the two placements"""
    v = 0
    for perm in (CHG, ('xp', 'y', 'x')):
        cp = T3.permuted(c, perm)
        v += TE.term_value(tex, T3.env(cp))
        if cc:
            v += np.conj(TE.term_value(tex, T3.env(T3.swapped_w(cp))))
    return v
