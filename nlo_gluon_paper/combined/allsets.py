"""All terms by row, for the combined note.
A-side: coefficient of log(vee/Lambda) (target, JIMWLK);  P-side: coefficient of log(vee) (projectile, BFKL).
keep=(3,): symmetrized three-rho parts;  keep=(2,): two-rho parts."""
import env  # noqa
import sympy as sp
import users as US
import user4
import userP as UP
import lsets
import cancel_data as CD
from sources import SOURCES
from fullpt import g1_pointwise, g2iv_pointwise

SRC = {s.key: s for s in SOURCES}
ROW = {'A1': 'A,B', 'A2': 'A,B', 'B1': 'A,B', 'B2': 'A,B', 'GI.II': 'I.II', 'GI.III': 'I.III', 'GI.IV': 'I.IV',
       'GI.V1': 'I.V', 'GI.V2': 'I.V', 'GII.I1': 'II.I', 'GII.I2': 'II.I', 'GIII.I': 'III.I', 'GIII.II': 'III.II',
       'GIII.III': 'III.III', 'GIII.V': 'III.V', 'GIII.VI': 'III.VI'}
GII = SRC['GI.II']


def alpha13(keep):
    return US.source_terms(GII, keep, monos=(1, 3), subs=(1, 3))


def A_sets(keep, corrected=True):
    out = {}
    add = lambda k, ts: out.setdefault(k, []).extend(ts)
    for k, row in ROW.items():
        ts = US.source_terms(SRC[k], keep)
        if corrected and k == 'GI.III':
            ts = [t for t in ts if not t['tag'].rstrip('*').endswith('b')]
        if corrected and k == 'GI.II':
            bad = {id(t) for t in []}
            ts = ts + US.neg(alpha13(keep))
        add(row, ts)
    if corrected:
        add('D', US.source_terms(CD.D_SOURCE, keep))
    if keep == (2,):
        add('4rho', user4.corrected_terms() if corrected else user4.terms()[0])
        G = US.genuine()
        add('I.II', g1_pointwise(G['G1_II']))
        add('I.V', G['G1_V'])
        add('II.I', G['G2_I'])
        add('II.IV', g2iv_pointwise(G['G2_IV'], 0.5))
        add('III.III', US.neg(G['G3_III']) if corrected else G['G3_III'])
    return out


def P_sets(keep, corrected=True):
    """coefficient of log(vee): rows with P = A (four-rho rows and remainders) + region-P parts"""
    out = {}
    add = lambda k, ts: out.setdefault(k, []).extend(ts)
    for k in ('A1', 'A2', 'B1', 'B2'):
        add('A,B', US.source_terms(SRC[k], keep))
    for key in UP.PSRC:
        base = key[:-2]
        add(ROW[base], UP.p_terms(key, keep))
    if corrected:
        add('I.III', lsets.fixKh(keep))
        add('I.II', US.neg(alpha13(keep)))
    if keep == (2,):
        add('4rho', user4.corrected_terms() if corrected else user4.terms()[0])
        GP = UP.genuineP()
        add('II.I', GP['G2_I'])
        add('II.IV', GP['G2_IV'])
        add('III.III', UP.G3_III_P(-1.0 if corrected else 1.0))
    return out
