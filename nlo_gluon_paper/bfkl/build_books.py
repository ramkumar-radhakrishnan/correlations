import env  # noqa
import pickle, time, sys
import sympy as sp
import target as TG
import refP
import lsets
from book import Book, tex_kkey

t0 = time.time()
which = sys.argv[1] if len(sys.argv) > 1 else '2'
if which == '2':
    L = lsets.sets((2,), True)
    B = Book()
    for k, ts in L.items():
        for t in ts:
            B.add(t, k)
    print('rows added %.0fs' % (time.time() - t0)); sys.stdout.flush()
    for t in refP.bfkl():
        B.add(dict(t, c=-t['c']), '-BFKL')
    print('BFKL added %.0fs' % (time.time() - t0)); sys.stdout.flush()
    for t in TG.jimwlk():
        B.add(t, '+JIMWLK')
    B.finalize()
    out = 'book2.pkl'
else:
    L = lsets.sets((3,), True)
    B = Book()
    for k, ts in L.items():
        for t in ts:
            B.add(t, k)
    B.finalize()
    out = 'book3.pkl'
bad = 0
for key, nets in B.data.items():
    nb = sum(1 for n in nets if sp.nsimplify(sum(n['groups'].values())) != 0)
    bad += nb
    print('%-70s networks %3d  nonzero %d' % (tex_kkey(key)[:70], len(nets), nb))
print('structures', len(B.data), 'networks', sum(len(v) for v in B.data.values()), 'nonzero network totals', bad, ' %.0fs' % (time.time() - t0))
pickle.dump(B.dump(), open(out, 'wb'))
