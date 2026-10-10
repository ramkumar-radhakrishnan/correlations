import numpy as np
from collections import defaultdict
from core import config_pair, accumulate_N, accumulate_dr
import reference as RF
from compare import report
import rowcmp
S = rowcmp.user_sets(); R = RF.rows()
g1_user = []
for k in S:
    if k.startswith(('GI.', '4:G1:')) or k in ('G:G1_V', 'D', '-RowIII2', '-RowIIaU'):
        g1_user += S[k]
g1_ref = []
for k in R:
    if k.startswith('G1-'):
        g1_ref += R[k]
cfg, cfg2 = config_pair(7)
acc = defaultdict(lambda: defaultdict(complex))
accumulate_N(g1_user, cfg, cfg2, acc, label='notes')
accumulate_dr(S['G:G1_II'], cfg, acc, label='notes G1-II dimreg')
accumulate_N([dict(t, c=-t['c']) for t in g1_ref], cfg, cfg2, acc, label='-ref')
bad = report(acc, show=False)
print('Group I corrected vs reference (z-free colour integrated, dim. reg.): keys', len(acc), 'nonzero', len(bad))
for key, tot, sc, d in bad:
    print('  %9.4f %s' % (abs(tot), key), {k: round(v.real, 4) for k, v in d.items()})
