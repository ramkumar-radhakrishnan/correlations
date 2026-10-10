"""Reference check: the two-rho content of the region-T reference rows equals JIMWLK x LO pointwise in z."""
from collections import defaultdict
from core import Config, accumulate
import reference as RF
import target as TG
from compare import report
R = RF.rows()
print({k: len(v) for k, v in R.items()})
cfg = Config(3)
acc = defaultdict(lambda: defaultdict(complex))
for k, ts in R.items():
    accumulate(ts, cfg, acc, False, label=k)
accumulate([dict(t, c=-t['c']) for t in TG.jimwlk()], cfg, acc, False, label='-JIMWLK')
report(acc, top=5)
