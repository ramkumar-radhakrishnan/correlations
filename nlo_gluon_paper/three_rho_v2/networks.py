"""Exact colour-network tables of the symmetrized three-rho terms (list as sent with typos fixed, and corrected).
python3 networks.py -> networks.pkl"""
import pickle
import sympy as sp
import env
from core import symmetrize
from book import Book, tex_kkey
import parse as P
import fixes as FX

GROUP = {'4rho.a': 'A,B', '4rho.b': 'A,B', 'I.II.1': 'I.II', 'I.II.2': 'I.II', 'I.III': 'I.III', 'D': 'D', 'I.IV': 'I.IV',
         'I.V.1': 'I.V', 'I.V.2': 'I.V', 'II.I': 'II.I', 'III.I': 'III.I', 'III.II': 'III.II', 'III.III': 'III.III',
         'III.V': 'III.V', 'III.VI': 'III.VI'}


def build(entries, extra_group=None):
    B = Book()
    for e in entries:
        g = GROUP[e['block']]
        if extra_group and extra_group(e):
            g = extra_group(e)
        for t in symmetrize(P.as_term(e), (3,)):
            B.add(t, g)
    B.finalize()
    return B.dump()


def removed(g):
    return 'removed' in g


def sums(n):
    w = sp.nsimplify(sum(v for g, v in n['groups'].items() if g != 'D'))
    c = sp.nsimplify(sum(v for g, v in n['groups'].items() if not removed(g)))
    return w, c


if __name__ == '__main__':
    T = FX.typo_fixed()
    tag = lambda e: ('I.II (F4: removed)' if (e['block'], e['n']) in FX.F4_DROP else
                     ('I.III (F3: removed)' if e['block'] == 'I.III' and e['cls'] in ('gamma', "gamma'") else None))
    bk = build(T + FX.d_entries(), tag)
    nets = sum(len(v) for v in bk.values())
    nw = sum(1 for ns in bk.values() for n in ns if sums(n)[0] != 0)
    nc = sum(1 for ns in bk.values() for n in ns if sums(n)[1] != 0)
    ncor = sum(1 for ns in bk.values() for n in ns if any(not removed(g) for g in n['groups']))
    print('kernel structures %d, networks %d; as written (typos fixed) nonzero %d; after F3, F4: networks %d, nonzero %d' % (len(bk), nets, nw, ncor, nc))
    sent = build(P.parse())
    print('as sent: networks %d, nonzero %d' % (sum(len(v) for v in sent.values()), sum(1 for ns in sent.values() for n in ns if sp.nsimplify(sum(n['groups'].values())) != 0)))
    pickle.dump(dict(all=bk, sent=sent), open('networks.pkl', 'wb'))
