"""Row-by-row two-rho comparison with the region-T reference, as written and corrected -> diag_rows.json"""
import json, re
import numpy as np
from collections import defaultdict
from core import Config, accumulate
import users as US
import reference as RF
from compare import report
from fullpt import g1_pointwise, g2iv_pointwise
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


def nophi(ts):
    return [t for t in ts if not any(f[0] == 'PHI' for f in t['kern'])]


neg = US.neg
G1II = g1_pointwise(S['G:G1_II'])
G3fix = neg(S['G:G3_III'])
rows = [
    ('Group II, Row I', S['GII.I1'] + S['GII.I2'] + S['4:G2-I:L5'] + S['G:G2_I'], None, R['G2-I']),
    ('Group II, Row II', S['4:G2-II:L6'], None, R['G2-II']),
    ('Group II, Row III', S['4:G2-III:L7'], None, R['G2-III']),
    ('Group II, Row IV', g2iv_pointwise(S['G:G2_IV'], 0.5), None, R['G2-IV']),
    ('Group III, Row I', S['GIII.I'] + S['4:G3-I:L8'], None, R['G3-I']),
    ('Group III, Row II', S['GIII.II'] + S['4:G3-II:L9'], None, R['G3-II']),
    ('Group III, Row III', S['GIII.III'] + S['G:G3_III'], S['GIII.III'] + G3fix, R['G3-III']),
    ('Group III, Row IV', S['4:G3-IV:L10'], None, R['G3-IV']),
    ('Group III, Row V', S['GIII.V'], None, R['G3-V']),
    ('Group III, Row VI', S['GIII.VI'], None, R['G3-VI']),
    ('Group I, all rows',
     nophi(S['GI.II'] + S['GI.III'] + S['GI.IV'] + S['GI.V1'] + S['GI.V2']) + sum((S[k] for k in S if k.startswith('4:G1:')), []) + G1II + S['G:G1_V'],
     S['GI.II'] + S['GI.III'] + S['GI.IV'] + S['GI.V1'] + S['GI.V2'] + S['-RowIII2'] + S['-RowIIaU'] + S['D']
     + sum((S[k] for k in S if k.startswith('4:G1:')), []) + G1II + S['G:G1_V'],
     sum((R[k] for k in R if k.startswith('G1-')), [])),
]


def dev(us, rf, seed=5):
    cfg = Config(seed)
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate(us, cfg, acc, False, label='notes')
    accumulate([dict(t, c=-t['c']) for t in rf], cfg, acc, False, label='-ref')
    bad = report(acc, show=False)
    mx = max([abs(sum(d.values())) for d in acc.values()] + [0])
    sc = max([abs(-d.get('-ref', 0)) for d in acc.values()] + [0])
    return mx, sc, len(bad), len(acc)


out = []
for name, us, usc, rf in rows:
    a = dev(us, rf)
    b = dev(usc, rf) if usc is not None else a
    out.append(dict(row=name, written=a, corrected=b))
    print('%-22s written %.2e (%d/%d keys bad)  corrected %.2e (%d bad)  size %.1f' % (name, a[0], a[2], a[3], b[0], b[2], a[1]))
json.dump(out, open('diag_rows.json', 'w'), indent=1)
