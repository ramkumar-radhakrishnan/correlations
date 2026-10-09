"""The two corrections of Group I as term lists (classes4 format), for the numerical checks.

correction 1 : remove the alpha part of the Row II kernel from the U(w) A^(3) sub-terms 1, 3 (and c.c.)
correction 2 : remove part 2 (the Phi terms) of Row III and add D (classes6.row3_diag_terms)
variant of 1 : remove the alpha part from all Row II sub-terms and add the commutator -[Nbar_2, Abar^(1)];
               it has the same symmetrized three-rho content as correction 1.
"""
import sympy as sp
import classes4 as C4
import classes6 as C6
from cancel3 import term_list


def correction1_terms():
    out = []
    for key, cc, n, c, src, colour, word, cmap in term_list():
        if key == 'GI.II' and n in (0, 2):
            monos = [(-sp.Integer(kc), [('KK',) + f for f in fs]) for kc, fs in C4.ALPHA_GII]
            out.append(C4.ct('correction 1', cc, n + 1, C4.to_sym(c), colour, word, monos, cmap))
    return out


def correction2_terms():
    return [t for t in C4.fix_terms() if t['name'].startswith('fixB')] + C6.row3_diag_terms()


def correction1_variant_terms():
    return [t for t in C4.fix_terms() if t['name'].startswith('fixA')]
