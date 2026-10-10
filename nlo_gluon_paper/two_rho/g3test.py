import numpy as np
from collections import defaultdict
from core import Config, accumulate
import reference as RF
from compare import report
import rowcmp
S = rowcmp.user_sets(); R = RF.rows()
for name, sgn in (('as written', 1), ('genuine part 2 with sign flipped', -1)):
    cfg = Config(5)
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate(S['GIII.III'], cfg, acc, False, label='3rho part')
    accumulate([dict(t, c=sgn * t['c']) for t in S['G:G3_III']], cfg, acc, False, label='genuine part 2')
    accumulate([dict(t, c=-t['c']) for t in R['G3-III']], cfg, acc, False, label='-ref')
    bad = report(acc, show=False)
    print(name, 'nonzero keys', len(bad))
    for key, tot, sc, d in bad:
        print('   ', key, {k: round(v.real, 4) for k, v in d.items()})
