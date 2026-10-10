"""The log(vee/k+) terms, corrected, split by row and by origin:
   'P' = region-P limit of the row (from p+ >> k+),  'T' = minus the region-T limit (= minus the row's log(vee/Lambda) terms).
Returns {group_label: [terms]} for keep=(2,) (two-rho parts) or keep=(3,) (symmetrized three-rho parts)."""
import env  # noqa
import sympy as sp
import userP as UP
import users as US
from sources import SOURCES
import cancel_data as CD
from fullpt import g1_pointwise, g2iv_pointwise

SRC = {s.key: s for s in SOURCES}
ROW = {'GI.II': 'I.II', 'GI.III': 'I.III', 'GI.IV': 'I.IV', 'GI.V1': 'I.V', 'GI.V2': 'I.V', 'GII.I1': 'II.I',
       'GII.I2': 'II.I', 'GIII.I': 'III.I', 'GIII.II': 'III.II', 'GIII.III': 'III.III', 'GIII.V': 'III.V', 'GIII.VI': 'III.VI'}
PKEY = {'GI.III': 'GI.III.P', 'GI.IV': 'GI.IV.P', 'GII.I1': 'GII.I1.P', 'GII.I2': 'GII.I2.P', 'GIII.I': 'GIII.I.P',
        'GIII.II': 'GIII.II.P', 'GIII.III': 'GIII.III.P', 'GIII.V': 'GIII.V.P', 'GIII.VI': 'GIII.VI.P'}


def fixKh(keep):
    src = UP.PSRC['GI.III.P']
    return US.source_terms(UP.mk('GI.III.Kh', 'GI.III', sp.I, 'K', [(4, src.kexpl[1][1])]), keep)


def sets(keep=(2,), corrected=True):
    out = {}

    def add(lab, ts):
        out.setdefault(lab, []).extend(ts)
    for k, row in ROW.items():
        T = US.source_terms(SRC[k], keep)
        if corrected and k == 'GI.III':
            T = [t for t in T if not t['tag'].rstrip('*').endswith('b')] + US.source_terms(CD.D_SOURCE, keep)
        add(row + ' T', US.neg(T))
        if k in PKEY:
            add(row + ' P', UP.p_terms(PKEY[k], keep))
    if corrected:
        add('I.III P', fixKh(keep))
    if keep == (2,):
        G = US.genuine()
        add('I.II T', US.neg(g1_pointwise(G['G1_II'])))
        add('I.V T', US.neg(G['G1_V']))
        add('II.I P', UP.G2_I_P())
        add('II.I T', US.neg(G['G2_I']))
        add('II.IV P', UP.G2_IV_P())
        add('II.IV T', US.neg(g2iv_pointwise(G['G2_IV'], 0.5)))
        s = -1.0 if corrected else 1.0
        add('III.III P', UP.G3_III_P(s))
        add('III.III T', [dict(t, c=-s * t['c']) for t in G['G3_III']])
    return out
