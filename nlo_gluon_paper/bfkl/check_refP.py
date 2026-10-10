import env  # noqa
from collections import defaultdict
from core import Config, accumulate
from compare import report
import refP
R = refP.rows()
print({k: len(v) for k, v in R.items()})
B = refP.bfkl()
print('BFKL terms', len(B))
cfg = Config(3)
acc = defaultdict(lambda: defaultdict(complex))
for k, ts in R.items():
    accumulate(ts, cfg, acc, False, label=k)
accumulate([dict(t, c=-t['c']) for t in B], cfg, acc, False, label='-BFKL')
report(acc, top=5)
