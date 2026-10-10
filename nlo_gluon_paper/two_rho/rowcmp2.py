import re, numpy as np
from collections import defaultdict
from core import Config, accumulate
import users as US
import reference as RF
from compare import report
import rowcmp

S = rowcmp.user_sets()
R = RF.rows()


def sub(ts, ks):
    out = []
    for t in ts:
        m = re.search(r'\.(\d)[a-d]?\*?$', t['tag'])
        if m and int(m.group(1)) in ks:
            out.append(t)
    return out


def tagged(ts, pat):
    return [t for t in ts if re.search(pat, t['tag'])]


G = S
groups = [
    ('G2-I', S['GII.I1'] + S['GII.I2'] + S['4:G2-I:L5'] + S['G:G2_I'], R['G2-I']),
    ('G2-II', S['4:G2-II:L6'], R['G2-II']),
    ('G2-III', S['4:G2-III:L7'], R['G2-III']),
    ('G3-I', S['GIII.I'] + S['4:G3-I:L8'], R['G3-I']),
    ('G3-II', S['GIII.II'] + S['4:G3-II:L9'], R['G3-II']),
    ('G3-III', S['GIII.III'] + S['G:G3_III'], R['G3-III']),
    ('G3-IV', S['4:G3-IV:L10'], R['G3-IV']),
    ('G3-V', S['GIII.V'], R['G3-V']),
    ('G3-VI', S['GIII.VI'], R['G3-VI']),
    ('G1 Row I (normalization)', S['4:G1:L1.k2=1'] + S['4:G1:L1.k2=2'], R['G1-I_Nafter'] + R['G1-I_Nbarbefore']),
    ('G1 Row II direct (U A3), corrected', sub(S['GI.II'], (1, 3)) + S['-RowIIaU'] + S['4:G1:L2.k2=1'], R['G1-II_A3']),
    ('G1 Row III direct, corrected', sub(S['GI.III'], (1, 3)) + sub(S['-RowIII2'], (1, 3)) + S['4:G1:L4.k2=2'], R['G1-III']),
    ('G1 Row IV direct', sub(S['GI.IV'], (1, 3)) + S['4:G1:L3.k2=1'], R['G1-IV']),
    ('G1 Row V direct (3rho part)', sub(S['GI.V1'], (1, 3)) + sub(S['GI.V2'], (1, 3)), R['G1-V']),
    ('G1 Row V direct (3rho + genuine ff part)', sub(S['GI.V1'], (1, 3)) + sub(S['GI.V2'], (1, 3)) + tagged(S['G:G1_V'], r'2a\.dir'), R['G1-V']),
    ('D direct', sub(S['D'], (1, 3)), R['G1-Ia']),
    ('G1 unitarity combination, corrected', sub(S['GI.II'], (2, 4)) + sub(S['GI.III'], (2, 4)) + sub(S['-RowIII2'], (2, 4))
        + sub(S['GI.IV'], (2, 4)) + sub(S['GI.V1'], (2, 4)) + sub(S['GI.V2'], (2, 4)) + sub(S['D'], (2, 4))
        + S['4:G1:L2.k2=2'] + S['4:G1:L3.k2=2'] + S['4:G1:L4.k2=1'] + tagged(S['G:G1_V'], r'2b\.bar|\.x'), R['G1-II_A3bar']),
]
for name, us, rf in groups:
    cfg = Config(5)
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate(us, cfg, acc, False, label='notes')
    accumulate([dict(t, c=-t['c']) for t in rf], cfg, acc, False, label='-ref')
    bad = report(acc, show=False)
    mx = max([abs(sum(d.values())) for k, d in acc.items()] + [0])
    print('%-45s keys %3d  nonzero %3d   max|notes-ref| %.3e' % (name, len(acc), len(bad), mx))
    for key, tot, scale, d in bad[:8]:
        print('        %8.4f  %s   notes %8.4f ref %8.4f' % (abs(tot), key, d.get('notes', 0).real, -d.get('-ref', 0).real))
