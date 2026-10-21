"""LaTeX fragments of the final note:
  - every two-rho term proportional to log(vee/Lambda), coloured by kernel class;
  - every two-rho term of the coefficient of log(vee) (BFKL), coloured by kernel class;
  - the sums by kernel structure and colour network against JIMWLK x LO and BFKL x LO.
Kernel classes (from the kernel of the term):
  J  the unobserved gluon z is attached to the measured gluon (w or w')        -> only in JIMWLK x LO
  C  the unobserved gluon z is exchanged between the two charges (x, x')        -> in both, with different coefficients
  B  the measured gluon is emitted by the gluon at z                            -> only in BFKL x LO"""
import json
import re
import sys
import env  # noqa
from collections import OrderedDict
import sympy as sp
from book import Book
from engine import tex_colour
import sets as S2
import target as TG
import refP
sys.path.insert(0, env.V2)
from two_from_final import coef_tex, kern_tex  # noqa: E402
import make_doc as MD  # noqa: E402
from genuine_new import genuine_new, G1_V  # noqa: E402
from parse_user import sec2_terms, parse_chunk  # noqa: E402

COL = {'J': 'JIM', 'C': 'CC', 'B': 'BFK'}


def klass_key(key):
    em, soft = MD.split_key(key)
    if any('z' in e for e in em):
        return 'B'
    if not soft:
        return 'C'
    m = re.match(r'K\(([^,]+),([^,]+),z\)', soft[0])
    a, b = m.group(1), m.group(2)
    return 'C' if a in ('x', "x'", 'y') and b in ('x', "x'", 'y') else 'J'


def klass_term(t):
    lo, soft = t['kern']
    if 'z' in lo[1:]:
        return 'B'
    a, b = soft[1], soft[3]
    return 'C' if a in ('x', "x'") and b in ('x', "x'") else 'J'


def book(groups):
    B = Book()
    for g, ts in groups.items():
        for t in ts:
            B.add(t, g)
    B.finalize()
    return B.dump()


def entries(bk, group):
    out = []
    for key in MD.ordered(bk):
        for net in bk[key]:
            c = sp.expand(net['groups'].get(group, 0))
            if c != 0:
                out.append((key, c, net['col']))
    return out


def lines(ents, count):
    L = []
    for key, c, col in ents:
        k = klass_key(key)
        count[k] += 1
        L.append(r'&\quad{\color{%s}%s\;%s\;%s}\notag\\' % (COL[k], coef_tex(c), kern_tex(key), tex_colour(col)))
    return L


def equation(title, label, ents, count):
    L = [r'\begin{align}', r'\frac{d^3N}{d^3k}\bigg|_{\rm %s}={}&\NN\int_{x,x\',z,w,w\'}e^{-ik\cdot(w\'-w)}\Bigg[\notag\\' % title]
    body = lines(ents, count)
    if not body:
        body = [r'&\quad 0\notag\\']
    body[-1] = body[-1][:-len(r'\notag\\')] + r'\Bigg]\label{%s}' % label
    return '\n'.join(L + body + [r'\end{align}'])


def check_sums(bk, groups, target, kind='A'):
    """per kernel structure: reduce the colour networks to an independent basis (../two_rho/reduce.py,
    ../bfkl/reduceP.py) and count the basis elements where sum over groups + target (target enters negated) != 0"""
    import reduce as RD
    import reduceP as RP
    bad, nets, keys, basis = 0, 0, 0, 0
    per = []
    for key in MD.ordered(bk):
        keys += 1
        sub = []
        for net in bk[key]:
            gr = {g: v for g, v in net['groups'].items() if g in groups or g == target}
            if gr:
                sub.append(dict(col=net['col'], groups=gr))
        nets += len(sub)
        if not sub:
            continue
        if kind == 'A':
            table, rel, cols = RD.reduce_structure(sub)
        else:
            table, rel, cols = RP.reduce_structure(sub, (target,))
        nb = sum(1 for r in table if r['total'] != 0)
        bad += nb
        basis += len(table)
        per.append((key, klass_key(key), len(table), nb))
    return dict(structures=keys, networks=nets, basis=basis, nonzero=bad), per


if __name__ == '__main__':
    A = S2.A_sets()
    N = genuine_new()
    gv = G1_V()
    gen = OrderedDict([
        ('G1-II', N['g:I.II']),
        ('G1-V 2a', [t for t in gv if t['tag'].startswith('n:G1.V.2a')]),
        ('G1-V 2b endpoints', [t for t in gv if t['tag'].startswith('n:G1.V.2b.fin')]),
        ('G1-V 2b UV', [t for t in gv if t['tag'].startswith('n:G1.V.2b.UV')]),
        ('G1-V 2b CT', [t for t in gv if t['tag'].startswith('n:G1.V.2b.CT')]),
        ('G2-I', N['g:II.I']), ('G2-IV', N['g:II.IV']), ('G3-III(2)', N['g:III.III'])])
    four = OrderedDict([('4rho', A['4rho'])])
    U2 = sec2_terms('input/Evolution_terms.tex')
    user = OrderedDict()
    for lab, t in U2:
        user.setdefault(lab, []).append(t)
    allA = OrderedDict(list(gen.items()) + list(four.items()) + list(user.items()))
    allA['-JIMWLK'] = [dict(x, c=-x['c']) for x in TG.jimwlk()]
    bkA = book(allA)
    resA, perA = check_sums(bkA, [g for g in allA if g != '-JIMWLK'], '-JIMWLK')
    print('log(vee/Lambda): sum vs JIMWLK x LO', resA)
    noFour, perNoFour = check_sums(bkA, [g for g in allA if g not in ('-JIMWLK', '4rho')], '-JIMWLK')
    print('without the 4rho.tex lists', noFour)

    P = S2.P_sets()
    PL = OrderedDict((k, v) for k, v in P.items())
    PL['-BFKL'] = [dict(x, c=-x['c']) for x in refP.bfkl()]
    bkP = book(PL)
    resP, perP = check_sums(bkP, [g for g in PL if g != '-BFKL'], '-BFKL', 'P')
    print('log(vee): sum vs BFKL x LO', resP)

    from collections import Counter
    frag = {}
    cnt = Counter()
    frag['gen'] = '\n\n'.join(equation(g, 'eq:g%d' % i, entries(bkA, g), cnt) for i, g in enumerate(gen))
    cntg = dict(cnt)
    cnt = Counter()
    frag['four'] = equation(r'4\rho\ {\rm lists}', 'eq:fourlist', entries(bkA, '4rho'), cnt)
    cnt4 = dict(cnt)
    cnt = Counter()
    frag['J'] = equation(r'JIMWLK\otimes LO', 'eq:Jsum', [(k, -c, col) for k, c, col in entries(bkA, '-JIMWLK')], cnt)
    cntJ = dict(cnt)
    cnt = Counter()
    PORDER = [k for k in S2.ORDER if k in P]
    frag['P'] = '\n\n'.join(equation(k.replace(':', '{:}'), 'eq:P%d' % i, entries(bkP, k), cnt) for i, k in enumerate(PORDER))
    cntP = dict(cnt)
    cnt = Counter()
    frag['B'] = equation(r'BFKL\otimes LO', 'eq:Bsum', [(k, -c, col) for k, c, col in entries(bkP, '-BFKL')], cnt)
    cntB = dict(cnt)
    # class counts of the user's Sec. 2 lines
    cu = Counter(klass_term(t) for _, t in U2)
    json.dump(dict(A=resA, A_without_4rho=noFour, P=resP, counts=dict(gen=cntg, four=cnt4, user_sec2=dict(cu), J=cntJ, P=cntP, B=cntB),
                   perA=[(str(k), c, n, r) for k, c, n, r in perA], perNoFour=[(str(k), c, n, r) for k, c, n, r in perNoFour], perP=[(str(k), c, n, r) for k, c, n, r in perP]),
              open('render.json', 'w'), indent=1)
    pass  # fragments are rendered by make_note.py
    import pickle
    pickle.dump(dict(bkA=bkA, bkP=bkP, gen=list(gen), user=list(user), P=PORDER,
                     perA=perA, perP=perP, perNoFour=perNoFour), open('books.pkl', 'wb'))
    print('counts', dict(gen=cntg, four=cnt4, user=dict(cu), J=cntJ, P=cntP, B=cntB))
