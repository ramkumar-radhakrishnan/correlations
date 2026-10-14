"""All two-rho terms, one label per source (after the fixes F1-F6, P1).

A_sets(): coefficient of log(vee/Lambda)            (units N  = N0 log(vee/Lambda))
P_sets(): coefficient of log(vee) = A + B           (units N0 log(vee));  B = P - A is the coefficient of log(vee/k+)
N0 = (1/(2pi)^3) g^4/(16 pi^5) (1/k+).

Labels:  3:<row>  two-rho terms from reordering the genuine three-rho terms of Evolution_three_rho.tex
         D        two-rho terms of the replacement D of Group I Row III part 2 (fix F3)
         g:<row>  genuine two-rho terms of Evolution_two_rho.tex
         4rho     two-rho lists of 4rho.tex (fixes F1, F2)
         A,B      two-rho terms of the three-rho remainders A1, A2, B1, B2 of the four-rho terms"""
from collections import OrderedDict
import env  # noqa
import users as US
import user4
import userP as UP
import lsets
import cancel_data as CD
from sources import SOURCES
from fullpt import g1_pointwise, g2iv_pointwise

SRC = {s.key: s for s in SOURCES}
THREE_LABELS = ['3:I.II', '3:I.III', 'D', '3:I.IV', '3:I.V', '3:II.I', '3:III.I', '3:III.II', '3:III.III', '3:III.V', '3:III.VI']
GEN_LABELS = ['g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III']
FOUR_LABELS = ['A,B', '4rho']
ORDER = THREE_LABELS + GEN_LABELS + FOUR_LABELS
SOURCES_OF = {'3:I.II': ['GI.II'], '3:I.III': ['GI.III'], 'D': [], '3:I.IV': ['GI.IV'], '3:I.V': ['GI.V1', 'GI.V2'],
              '3:II.I': ['GII.I1', 'GII.I2'], '3:III.I': ['GIII.I'], '3:III.II': ['GIII.II'], '3:III.III': ['GIII.III'],
              '3:III.V': ['GIII.V'], '3:III.VI': ['GIII.VI'], 'A,B': ['A1', 'A2', 'B1', 'B2']}


def src(*keys):
    out = []
    for k in keys:
        out += US.source_terms(SRC[k], (2,))
    return out


def alpha13():
    """the alpha part of Group I Row II in the U(w)A3 sub-terms 1, 3 (fix F4 moves it from log(vee/Lambda) to -log(vee/k+))"""
    return US.source_terms(SRC['GI.II'], (2,), monos=(1, 3), subs=(1, 3))


def four():
    return OrderedDict([('A,B', src('A1', 'A2', 'B1', 'B2')), ('4rho', user4.corrected_terms())])


def A_sets():
    out = OrderedDict()
    out['3:I.II'] = src('GI.II') + US.neg(alpha13())                                         # F4
    out['3:I.III'] = [t for t in src('GI.III') if not t['tag'].rstrip('*').endswith('b')]     # F3: part 2 removed
    out['D'] = US.source_terms(CD.D_SOURCE, (2,))                                             # F3: D added
    out['3:I.IV'] = src('GI.IV')
    out['3:I.V'] = src('GI.V1', 'GI.V2')
    out['3:II.I'] = src('GII.I1', 'GII.I2')
    for r in ('I', 'II', 'III', 'V', 'VI'):
        out['3:III.' + r] = src('GIII.' + r)
    G = US.genuine()
    out['g:I.II'] = g1_pointwise(G['G1_II'])                                                  # F6
    out['g:I.V'] = G['G1_V']
    out['g:II.I'] = G['G2_I']
    out['g:II.IV'] = g2iv_pointwise(G['G2_IV'], 0.5)
    out['g:III.III'] = US.neg(G['G3_III'])                                                    # F5
    out.update(four())
    return out


def P_sets():
    out = OrderedDict()
    out['3:I.II'] = US.neg(alpha13())                                                         # F4
    out['3:I.III'] = UP.p_terms('GI.III.P') + lsets.fixKh((2,))                               # P1
    out['3:I.IV'] = UP.p_terms('GI.IV.P')
    out['3:II.I'] = UP.p_terms('GII.I1.P') + UP.p_terms('GII.I2.P')
    for r in ('I', 'II', 'III', 'V', 'VI'):
        out['3:III.' + r] = UP.p_terms('GIII.%s.P' % r)
    out['g:II.I'] = UP.G2_I_P()
    out['g:II.IV'] = UP.G2_IV_P()
    out['g:III.III'] = UP.G3_III_P(-1.0)                                                      # F5
    out.update(four())
    return out
