import env  # noqa
from collections import defaultdict
from core import term, symmetrize
import fourrho
import refP
from acc4 import Config3, accumulate3, report3
K4 = (4,)
R = refP.rows(keep=K4)
W4 = []
for c, row, lab, Us, word, kern, cc in fourrho.words():
    kk = [('KK', a, b, c_, d) for (a, b), (c_, d) in kern]
    W4 += symmetrize(term(c, [(u[1], u[2], u[3]) for u in Us], [], word, kk, '4rho'), K4)
cfg = Config3(3)
acc = defaultdict(lambda: defaultdict(complex))
accumulate3(sum(R.values(), []), cfg, acc, label='refP')
print('region-P reference, four-rho total:', end=' '); report3(acc, top=0)
acc = defaultdict(lambda: defaultdict(complex))
accumulate3(W4, cfg, acc, label='notes')
accumulate3([dict(t, c=-t['c']) for t in sum(R.values(), [])], cfg, acc, label='-refP')
print('notes Eq.(4rho) minus region-P reference rows, four-rho:', end=' '); report3(acc, top=0)
