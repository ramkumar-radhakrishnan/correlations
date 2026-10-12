"""Exact tables for the coefficient of log(vee) (P, projectile): two-rho terms vs BFKL x LO, and symmetrized three-rho terms.
python3 build_P.py 2|3  -> bookP2.pkl / bookP3.pkl"""
import env  # noqa
import pickle, time, sys
import sympy as sp
import refP
import allsets as AS
from book import Book, tex_kkey

t0 = time.time()
which = sys.argv[1] if len(sys.argv) > 1 else '2'
keep = (2,) if which == '2' else (3,)
L = AS.P_sets(keep, True)
B = Book()
for k, ts in L.items():
    for t in ts:
        B.add(t, k)
print('rows added %.0fs' % (time.time() - t0)); sys.stdout.flush()
if which == '2':
    for t in refP.bfkl():
        B.add(dict(t, c=-t['c']), '-BFKL')
B.finalize()
bad = 0
for key, nets in B.data.items():
    nb = sum(1 for n in nets if sp.nsimplify(sum(n['groups'].values())) != 0)
    bad += nb
    print('%-70s networks %3d  nonzero %d' % (tex_kkey(key)[:70], len(nets), nb))
print('structures', len(B.data), 'networks', sum(len(v) for v in B.data.values()), 'nonzero network totals', bad, ' %.0fs' % (time.time() - t0))
pickle.dump(B.dump(), open('bookP%s.pkl' % which, 'wb'))
