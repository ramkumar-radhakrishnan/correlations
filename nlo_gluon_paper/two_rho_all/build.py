"""Exact tables by kernel structure and colour network, one group per source.
python3 build.py A   -> bookA.pkl, reducedA.pkl   (log(vee/Lambda) terms and -JIMWLK x LO)
python3 build.py P   -> bookP.pkl, reducedP.pkl   (log(vee) terms and -BFKL x LO)"""
import pickle
import sys
import time
import sympy as sp
import env  # noqa
import sets as S
import target as TG
import refP
import reduce as RD          # ../two_rho/reduce.py
import reduceP as RP         # ../bfkl/reduceP.py
from book import Book, tex_kkey

which = sys.argv[1] if len(sys.argv) > 1 else 'A'
t0 = time.time()
B = Book()
if which == 'A':
    tgt, L = '-JIMWLK', S.A_sets()
    for t in TG.jimwlk():
        B.add(dict(t, c=-t['c']), tgt)
else:
    tgt, L = '-BFKL', S.P_sets()
    for t in refP.bfkl():
        B.add(dict(t, c=-t['c']), tgt)
for k, ts in L.items():
    for t in ts:
        B.add(t, k)
B.finalize()
book = B.dump()
print('book built %.0fs: structures %d, networks %d' % (time.time() - t0, len(book), sum(len(v) for v in book.values())))
out, nbad = {}, 0
for key, nets in book.items():
    if which == 'A':
        table, rel, cols = RD.reduce_structure(nets)
    else:
        table, rel, cols = RP.reduce_structure(nets, (tgt,))
    out[key] = dict(table=table, relations=rel, cols=cols, nets=nets)
    nb = sum(1 for r in table if r['total'] != 0)
    nbad += nb
    print('%-66s networks %3d basis %3d nonzero %d' % (tex_kkey(key)[:66], len(nets), len(table), nb))
print('basis elements with nonzero total:', nbad, ' %.0fs' % (time.time() - t0))
pickle.dump(book, open('book%s.pkl' % which, 'wb'))
pickle.dump(out, open('reduced%s.pkl' % which, 'wb'))
