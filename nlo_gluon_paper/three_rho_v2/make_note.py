"""Write ThreeRho_check.tex from note_static.tex and the computed results.
Inputs: states3.json (check3b.py), networks.pkl (networks.py), two_cmp.json (two_cmp.py); builds the two-rho books itself."""
import json
import pickle
import re
import sys
from collections import OrderedDict
import sympy as sp
import env
from core import symmetrize
from book import Book
import target as TG
import reduce as RD
import parse as P
import fixes as FX
import oporder as OP
from engine import tex_colour
sys.path.insert(0, env.ALL2)
import sets as S2
import make_doc as MD   # ../two_rho_all/make_doc.py (kernel names, coefficient and colour formatting)

KTEX = {'alpha': r'K_{\alpha}', "alpha'": r"K_{\alpha'}", 'beta': r'K_{\beta}', "beta'": r"K_{\beta'}",
        'gamma': r'K_{\gamma}', "gamma'": r"K_{\gamma'}"}
BLOCK_PREF = {'4rho.a': r'-i\NN', '4rho.b': r'i\NN'}
ROWNAME = OrderedDict([('A,B', 'Four-$\\rho$ remainders'), ('I.II', 'Group I, Row II'), ('I.III', 'Group I, Row III'),
                       ('D', r'$\mathcal D$'), ('I.IV', 'Group I, Row IV'), ('I.V', 'Group I, Row V'), ('II.I', 'Group II, Row I'),
                       ('III.I', 'Group III, Row I'), ('III.II', 'Group III, Row II'), ('III.III', 'Group III, Row III'),
                       ('III.V', 'Group III, Row V'), ('III.VI', 'Group III, Row VI')])


def rho(i, p):
    return r'\rho^{%s}(%s)' % (i, p)


def U(i, j, p):
    return r'U^{%s%s}(%s)' % (i, j, p)


def coef(c):
    c = sp.nsimplify(c)
    if c == 1:
        return '+'
    if c == -1:
        return '-'
    s = sp.latex(c).replace(r'\frac{1}{2}', r'\tfrac12')
    return s if s.startswith('-') else '+' + s


def colour_tex(e, rho_part):
    f = r'f^{%s%s%s}' % e['fs'][0]
    return f + r'\,' + rho_part + r'\,' + r'\,'.join(U(*u) for u in e['Us'])


def line_sym(e, mark=None):
    """a term of the list in its written (symmetrized) order"""
    c = e['sign'] * e['mult']
    rp = r'\,'.join(rho(i, p) for i, p in e['rs'])
    body = r'%s\,%s\,%s' % (coef(c), KTEX[e['cls']], colour_tex(e, rp))
    if mark == 'red':
        return r'\textcolor{red}{%s}' % body
    if mark == 'blue':
        return r'\textcolor{blue}{%s}' % body
    return body


def line_op(e):
    """the same term in operator order"""
    of = OP.operator_form(e)
    byp = {p: (i, p) for i, p in e['rs']}
    o = [byp[p] for p in of['order']]
    c = e['sign'] * e['mult']
    if of['form'] == 'ABC':
        rp = r'\,'.join(rho(*r) for r in o)
    elif of['form'] == 'A{BC}':
        c = sp.Rational(c, 2)
        rp = r'%s\big\{%s,%s\big\}' % (rho(*o[0]), rho(*o[1]), rho(*o[2]))
    else:
        c = sp.Rational(c, 2)
        rp = r'\big\{%s,%s\big\}%s' % (rho(*o[0]), rho(*o[1]), rho(*o[2]))
    changed = (of['form'] != 'ABC') or ([p for _, p in e['rs']] != of['order'])
    body = r'%s\,%s\,%s' % (coef(c), KTEX[e['cls']], colour_tex(e, rp))
    return (r'\textcolor{blue}{%s}' % body) if changed else body, of


def block_eq(lines, pref, label):
    out = [r'\begin{align*}', r"&%s\int e^{-ik\cdot(w'-w)}\Big[\\" % pref]
    for n, (num, l) in enumerate(lines):
        end = r'\Big]' if n == len(lines) - 1 else ''
        out.append(r'&\quad\text{(%s)}\ \ %s%s\\' % (num, l, end))
    out[-1] = out[-1].rstrip('\\')
    out.append(r'\end{align*}')
    return '\n'.join(out)


BLOCK_TITLE = OrderedDict([('4rho.a', r'From four $\rho$\textquotesingle s, $d^3N_a$'), ('4rho.b', r'From four $\rho$\textquotesingle s, $d^3N_b$'),
                           ('I.II.1', 'Group I, Row II, part 1'), ('I.II.2', 'Group I, Row II, part 2 (c.c.)'),
                           ('I.III', 'Group I, Row III'), ('D', r"$\mathcal D$ (replaces the $\gamma,\gamma'$ terms of Group I, Row III)"),
                           ('I.IV', 'Group I, Row IV'), ('I.V.1', 'Group I, Row V, part 1'), ('I.V.2', 'Group I, Row V, part 2'),
                           ('II.I', 'Group II, Row I'), ('III.I', 'Group III, Row I'), ('III.II', 'Group III, Row II'),
                           ('III.III', 'Group III, Row III'), ('III.V', 'Group III, Row V'), ('III.VI', 'Group III, Row VI')])


def corrected_list():
    T = FX.typo_fixed() + FX.d_entries()
    out = []
    for blk, title in BLOCK_TITLE.items():
        es = [e for e in T if e['block'] == blk]
        lines = []
        for e in es:
            if (e['block'], e['n']) in FX.F4_DROP:
                mark, tag = 'red', r'\ \text{\textcolor{red}{[F4: removed]}}'
            elif e['block'] == 'I.III' and e['cls'] in ('gamma', "gamma'"):
                mark, tag = 'red', r'\ \text{\textcolor{red}{[F3: removed]}}'
            elif e.get('fix'):
                mark, tag = 'blue', r'\ \text{\textcolor{blue}{[%s]}}' % e['fix']
            elif e['block'] == 'D':
                mark, tag = 'blue', r'\ \text{\textcolor{blue}{[F3]}}'
            else:
                mark, tag = None, ''
            lines.append((str(e['n']), line_sym(e, mark) + tag))
        out.append(r'\paragraph{%s.}' % title)
        out.append(block_eq(lines, BLOCK_PREF.get(blk, r'i\NN'), blk))
    return '\n'.join(out)


def operator_list():
    C = FX.corrected()
    out = []
    for blk, title in BLOCK_TITLE.items():
        es = [e for e in C if e['block'] == blk]
        lines = []
        forms = set()
        for e in es:
            l, of = line_op(e)
            forms.add((of['form'], of['cc']))
            lines.append((str(e['n']), l))
        out.append(r'\paragraph{%s.}' % title)
        out.append(block_eq(lines, BLOCK_PREF.get(blk, r'i\NN'), blk))
    return '\n'.join(out)


def typo_lines():
    E0 = {(e['block'], e['n']): e for e in P.parse()}
    T = {(e['block'], e['n']): e for e in FX.typo_fixed()}
    out = []
    for key, lab, what in ((('III.I', 15), 'T1', r"Group III, Row I, term 15 ($\beta'$)"),
                           (('III.II', 13), 'T2', r"Group III, Row II, term 13 ($\beta'$)")):
        out.append(r'\paragraph{%s: %s.}' % (lab, what))
        out.append(r'\begin{align*}')
        out.append(r'\text{as sent:}\quad&%s\\' % line_sym(E0[key], 'red'))
        out.append(r'\text{corrected:}\quad&%s' % line_sym(T[key], 'blue'))
        out.append(r'\end{align*}')
    return '\n'.join(out)


def removed_lines():
    T = FX.typo_fixed()
    f4 = [(('%s.%d' % (e['block'], e['n'])), line_sym(e)) for e in T if (e['block'], e['n']) in FX.F4_DROP]
    f3 = [(('%s.%d' % (e['block'], e['n'])), line_sym(e)) for e in T if e['block'] == 'I.III' and e['cls'] in ('gamma', "gamma'")]
    D = [(str(e['n']), line_sym(e)) for e in FX.d_entries()]
    return block_eq(f4, r'i\NN', 'f4'), block_eq(f3, r'i\NN', 'f3'), block_eq(D, r'i\NN', 'D')


def state_table():
    st = json.load(open('states3.json'))
    return '\n'.join(r'%s & %d & %d & %s & $%s$\\' % (s['state'], s['terms'], s['keys'], ', '.join(map(str, s['bad'])), MD.sci(s['maxres']))
                     for s in st)


# ------------------------------------------------------------------ three-rho networks
def kname3(key):
    out = []
    for k in key:
        if k[0] == 'PHI':
            out.append(r'\phi_{%s}(%s)' % (k[2], k[1]))
        elif k[0] == 'inv2':
            out.append(r'\frac{1}{(%s-%s)^2}' % k[1])
        else:
            (a, b), (c, d) = k[1], k[2]
            out.append(r'\KK(%s-%s)\cdot\KK(%s-%s)' % (a, b, c, d))
    return r'\;'.join(out)


GL3 = {'A,B': 'A,B', 'I.II': 'I.II', 'I.III': 'I.III', 'D': r'\textcolor{blue}{$\mathcal D$}', 'I.IV': 'I.IV', 'I.V': 'I.V',
       'II.I': 'II.I', 'III.I': 'III.I', 'III.II': 'III.II', 'III.III': 'III.III', 'III.V': 'III.V', 'III.VI': 'III.VI',
       'I.II (F4: removed)': r'\textcolor{red}{I.II$_{\rm F4}$}', 'I.III (F3: removed)': r'\textcolor{red}{I.III$_{\rm F3}$}'}
ORD3 = list(GL3)


def network_tables():
    NW = pickle.load(open('networks.pkl', 'rb'))['all']
    out = []
    for n, (key, nets) in enumerate(NW.items()):
        out.append(r'\subsubsection*{Kernel structure %d: $%s$}' % (n + 1, kname3(key)))
        out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.37\textwidth}<{$}>{\raggedright\arraybackslash}p{0.39\textwidth}>{$}c<{$}>{$}c<{$}@{}}')
        out.append(r'\text{colour network} & contributions (row: coefficient) & \text{as written} & \text{fixed}\\ \hline\endhead')
        for net in nets:
            g = net['groups']
            items = r';\ '.join(r'%s:\,$%s$' % (GL3[k], MD.ctex(v, plus=False)) for k, v in sorted(g.items(), key=lambda kv: ORD3.index(kv[0])))
            w = sp.nsimplify(sum(v for k, v in g.items() if k != 'D'))
            c = sp.nsimplify(sum(v for k, v in g.items() if 'removed' not in k))
            fixed_only = all('removed' in k for k in g)
            ctx = (r'\text{--}' if fixed_only else MD.ctex(c, plus=False)) + ('' if fixed_only else (r'\ \text{\ok}' if c == 0 else r'\ \text{\bad}'))
            wtx = MD.ctex(w, plus=False) if 'D' not in g or len(g) > 1 else '0'
            if list(g) == ['D']:
                wtx = r'\text{--}'
            out.append(r'%s & %s & %s & %s\tabularnewline[2pt]' % (MD.colour(net['col']), items, wtx, ctx))
        out.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(out)


# ------------------------------------------------------------------ two-rho terms from the operator-ordered list
ROW2 = {'4rho.a': 'A,B', '4rho.b': 'A,B', 'I.II.1': '3:I.II', 'I.II.2': '3:I.II', 'I.III': '3:I.III', 'D': 'D',
        'I.IV': '3:I.IV', 'I.V.1': '3:I.V', 'I.V.2': '3:I.V', 'II.I': '3:II.I', 'III.I': '3:III.I', 'III.II': '3:III.II',
        'III.III': '3:III.III', 'III.V': '3:III.V', 'III.VI': '3:III.VI'}
TWO_ORDER = ['3:I.II', '3:I.III', 'D', '3:I.IV', '3:I.V', '3:II.I', '3:III.I', '3:III.II', '3:III.III', '3:III.V', '3:III.VI', 'A,B']
GL2 = dict(MD.GL)


def two_books():
    Bk = Book()
    for t in TG.jimwlk():
        Bk.add(dict(t, c=-t['c']), '-JIMWLK')
    for e in FX.corrected():
        for t in OP.op_terms(e):
            for s in symmetrize(t, (2,)):
                Bk.add(s, ROW2[e['block']])
    A = S2.A_sets()
    for k in ('g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III', '4rho'):
        for t in A[k]:
            Bk.add(t, k)
    Bk.finalize()
    return Bk.dump()


def merged(book):
    """groups merged into: three-rho list, genuine, 4rho lists, -JIMWLK"""
    def m(g):
        if g == '-JIMWLK':
            return g
        if g.startswith('g:'):
            return 'gen'
        if g == '4rho':
            return '4rho'
        return '3rho'
    out = {}
    for key, nets in book.items():
        nn = []
        for n in nets:
            gg = {}
            for g, v in n['groups'].items():
                gg[m(g)] = sp.nsimplify(gg.get(m(g), 0) + v)
            nn.append(dict(col=n['col'], groups={g: v for g, v in gg.items() if v != 0}))
        out[key] = nn
    return out


def sum_tables(book, keys):
    M = merged(book)
    out = []
    nb = 0
    for n, key in enumerate(keys):
        table, rel, cols = RD.reduce_structure(M[key])
        out.append(r'\subsubsection*{Structure $S_{%d}$: $%s$}' % (n + 1, MD.kname(key)))
        out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.40\textwidth}<{$}>{$}c<{$}>{$}c<{$}>{$}c<{$}>{$}c<{$}>{$}c<{$}@{}}')
        out.append(r'\text{colour structure} & \text{3$\rho$ list} & \text{genuine} & \text{\texttt{4rho.tex}} & \text{sum} & \text{JIMWLK}\\ \hline\endhead')
        for row in table:
            c = row['contrib']
            if not c:
                continue
            nb += 1
            a, g, f = (c.get(k, 0) for k in ('3rho', 'gen', '4rho'))
            j = -c.get('-JIMWLK', 0)
            s = sp.nsimplify(sp.expand(a + g + f))
            ok = r'\ \text{\ok}' if sp.simplify(s - j) == 0 else r'\ \text{\bad}'
            out.append(r'%s & %s & %s & %s & %s & %s%s\tabularnewline[2pt]' % (
                MD.colour(row['col']), MD.ctex(a, plus=False), MD.ctex(g, plus=False), MD.ctex(f, plus=False),
                MD.ctex(s, plus=False), MD.ctex(j, plus=False), ok))
        out.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(out), nb


def literal_table():
    d = json.load(open('two_cmp.json'))
    rows = []
    for k, v in d['rows'].items():
        same = abs(v['literal'] - v['operator']) < 1e-9 * (1 + abs(v['operator']))
        rows.append(r'%s & $%.5f$ & $%.5f$ & %s\\' % (GL2[k], v['literal'], v['operator'], 'same' if same else r'\textcolor{red}{different}'))
    tot = '\n'.join(r'%d & $%.6f$ & $%.6f$ & $%.6f$\\' % (t['seed'], t['literal'], t['operator'], t['J']) for t in d['totals'])
    return '\n'.join(rows), tot


def main():
    T = open('note_static.tex').read()
    f4, f3, D = removed_lines()
    book = two_books()
    keys = MD.ordered(book)
    sums, nb = sum_tables(book, keys)
    ent = {}
    for g in TWO_ORDER:
        tab, n = MD.entries_table(book, g, keys, 'S')
        ent[g] = tab
    lit, tot = literal_table()
    rep = {
        '@@STATES@@': state_table(), '@@TYPOS@@': typo_lines(), '@@F4LINES@@': f4, '@@F3LINES@@': f3, '@@DLINES@@': D,
        '@@CORRECTED@@': corrected_list(), '@@NETWORKS@@': network_tables(), '@@OPLIST@@': operator_list(),
        '@@LITERAL@@': lit, '@@LITTOT@@': tot, '@@SUMS@@': sums, '@@NB@@': str(nb), '@@NS@@': str(len(keys)),
        '@@STRUCT@@': '\n'.join(r'$S_{%d}$ & $%s$\\' % (n + 1, MD.kname(k)) for n, k in enumerate(keys)),
    }
    for g in TWO_ORDER:
        rep['@@TWO:%s@@' % g] = ent[g]
    for k, v in rep.items():
        assert k in T, k
        T = T.replace(k, v)
    left = re.findall(r'@@[^@]*@@', T)
    assert not left, left
    open('ThreeRho_check.tex', 'w').write(T)
    print('written ThreeRho_check.tex; two-rho basis structures', nb)


if __name__ == '__main__':
    main()
