"""Write TwoRho_all.tex from static.tex and the computed tables.
Inputs: bookA.pkl, reducedA.pkl, bookP.pkl, reducedP.pkl (build.py), checks.json (checks.py)."""
import json
import pickle
import re
import sympy as sp
import env  # noqa
from engine import tex_colour
from kernels import tex_kernel
import make_tex as MT          # ../three_rho/make_tex.py (input rendering of the three-rho sources)
import cancel_data as CD
import make_two_tex as M2      # ../two_rho/make_two_tex.py (lists of 4rho.tex)
import sets as S

GL = {'3:I.II': r'I.II$^{3\rho}$', '3:I.III': r'I.III$^{3\rho}$', 'D': r'$\mathcal D$', '3:I.IV': r'I.IV$^{3\rho}$',
      '3:I.V': r'I.V$^{3\rho}$', '3:II.I': r'II.I$^{3\rho}$', '3:III.I': r'III.I$^{3\rho}$', '3:III.II': r'III.II$^{3\rho}$',
      '3:III.III': r'III.III$^{3\rho}$', '3:III.V': r'III.V$^{3\rho}$', '3:III.VI': r'III.VI$^{3\rho}$',
      'g:I.II': r'I.II$^{\rm gen}$', 'g:I.V': r'I.V$^{\rm gen}$', 'g:II.I': r'II.I$^{\rm gen}$', 'g:II.IV': r'II.IV$^{\rm gen}$',
      'g:III.III': r'III.III$^{\rm gen}$', 'A,B': r'A,B', '4rho': r'4$\rho$'}
TITLE = {'3:I.II': 'Group I, Row II', '3:I.III': 'Group I, Row III (part 1)', 'D': r'$\mathcal D$, the replacement of Group I, Row III, part 2',
         '3:I.IV': 'Group I, Row IV', '3:I.V': 'Group I, Row V', '3:II.I': 'Group II, Row I', '3:III.I': 'Group III, Row I',
         '3:III.II': 'Group III, Row II', '3:III.III': 'Group III, Row III (part 1)', '3:III.V': 'Group III, Row V',
         '3:III.VI': 'Group III, Row VI', 'g:I.II': 'Group I, Row II', 'g:I.V': 'Group I, Row V', 'g:II.I': 'Group II, Row I',
         'g:II.IV': 'Group II, Row IV', 'g:III.III': 'Group III, Row III (part 2)', 'A,B': 'Remainders A1, A2, B1, B2',
         '4rho': r'The lists of \texttt{4rho.tex}'}
FIX = {
    '3:I.II': r"""\emph{Fix F4.} In sub-terms 1 and 3 (the $U(w)\mathbb A^{(3)}$ term) the kernel $\mathfrak F^i$ is replaced by $\mathfrak F^i_\beta$: the monomials b and d (the $\mathfrak F^i_\alpha$ part) are removed. $\mathfrak F^i_\alpha$ is the commutator $[\mathbb N^{(2)},\mathbb A^{(1)}]$; with the $p^+$ ordering of the virtual factor it multiplies $-\log(\vee/k^+)$ instead of $\log(k^+/\Lambda)$ (combined note, Part I). It reappears in the coefficient of $\log\vee$, Sec.~\ref{sec:P}.""",
    '3:I.III': r"""\emph{Fix F3.} Only part 1 (monomial a, the kernel $K^{(1)}$) is kept here. Part 2 (monomial b, the edge kernel $\Phi$) is replaced by $\mathcal D$, Sec.~\ref{src:D}.""",
    'D': r"""This is the term that replaces part 2 of Group I, Row III (fix F3). It was found in the three-$\rho$ cancellation: with it all symmetrized three-$\rho$ terms cancel. Its two-$\rho$ terms follow from the charge order $\rho(x')\rho(x)\rho(y)$ of Row III.""",
}
SECTION_REF = {'3:I.II': ['GI.II'], '3:I.III': ['GI.III'], '3:I.IV': ['GI.IV'], '3:I.V': ['GI.V1', 'GI.V2'],
               '3:II.I': ['GII.I1', 'GII.I2'], '3:III.I': ['GIII.I'], '3:III.II': ['GIII.II'], '3:III.III': ['GIII.III'],
               '3:III.V': ['GIII.V'], '3:III.VI': ['GIII.VI'], 'A,B': ['A1', 'A2', 'B1', 'B2']}


# ------------------------------------------------------------------ formatting
def ctex(c, plus=True):
    c = sp.nsimplify(sp.expand(c))
    s = sp.latex(c).replace('N_{c}', 'N_c')
    s = (s.replace(r'\frac{N_c}{2}', r'\tfrac{N_c}{2}').replace(r'\frac{1}{2}', r'\tfrac12').replace(r'\frac{1}{4}', r'\tfrac14')
         .replace(r'\frac{3}{2}', r'\tfrac32').replace(r'\frac{3 N_c}{2}', r'\tfrac{3N_c}{2}').replace(r'\frac{N_c}{4}', r'\tfrac{N_c}{4}'))
    if plus and not s.startswith('-'):
        s = '+' + s
    return s


def cpre(c, first=False):
    c = sp.nsimplify(sp.expand(c))
    if c == 1:
        return '' if first else '+'
    if c == -1:
        return '-'
    return ctex(c, plus=not first) + r'\,'


def colour(col):
    return tex_colour(col).replace(r'\,', r'\,\allowbreak ')


def split_key(key):
    """kernel key -> (emission part, soft part).  The soft part is the factor K(a,b,z) of the unobserved gluon:
    the z-factor between two charges if there is one, otherwise the only factor containing z twice."""
    facs = []
    for k in key:
        if k[0] == 'inv2':
            a, b = k[1]
            facs.append(('soft2' if b == 'z' else 'other', r'K(%s,%s,z)' % (a, a) if b == 'z' else r'\frac{1}{(%s-%s)^2}' % (a, b), a, a))
        else:
            (a, b), (c, d) = k[1], k[2]
            if b == 'z' and d == 'z':
                facs.append(('zz', r'K(%s,%s,z)' % (a, c), a, c))
            else:
                facs.append(('em', r'\KK(%s-%s)\cdot\KK(%s-%s)' % (a, b, c, d), a, c))
    ext = ('w', "w'")
    cand = [f for f in facs if f[0] in ('zz', 'soft2') and f[2] not in ext and f[3] not in ext]
    if len(cand) != 1:
        cand = [f for f in facs if f[0] in ('zz', 'soft2')][:1]
    soft = [cand[0][1]] if cand else []
    em = []
    for f in facs:
        if cand and f is cand[0]:
            continue
        if f[0] == 'zz':
            em.append(r'\KK(%s-z)\cdot\KK(%s-z)' % (f[2], f[3]))
        else:
            em.append(f[1])
    return em, soft


def kname(key):
    em, soft = split_key(key)
    return r'\;'.join(em + soft)


SOFT_ORDER = ['K(x,x,z)', "K(x',x',z)", "K(x,x',z)", 'K(w,w,z)', "K(w',w',z)", "K(w,w',z)", 'K(x,w,z)', "K(x,w',z)",
              "K(x',w,z)", "K(x',w',z)"]
EM_ORDER = [r"\KK(x-w)\cdot\KK(x'-w')", r"\KK(x'-w)\cdot\KK(x'-w')", r"\KK(w-z)\cdot\KK(w'-z)", r"\KK(x'-w')\cdot\KK(w-z)",
            r"\KK(x-w)\cdot\KK(w'-z)"]


def key_rank(key):
    em, soft = split_key(key)
    e = EM_ORDER.index(em[0]) if len(em) == 1 and em[0] in EM_ORDER else 99
    s = SOFT_ORDER.index(soft[0]) if len(soft) == 1 and soft[0] in SOFT_ORDER else 99
    return (e, s, str(key))


def ordered(R):
    keys = sorted(R, key=key_rank)
    return keys


def deref(t):
    """remove references to labels of other notes"""
    t = t.replace(r'App.~\ref{app:F}', r'Part I of the combined note')
    return t


# ------------------------------------------------------------------ inputs of the three-rho sources
def src_header(src):
    L = [r'\begin{align}']
    head = r'\frac{d^3N}{d^3k}\bigg|_{\rm %s}={}&%s\,\NN\int %s\,%s(%s)\,\Big\{' % (
        src.key, MT.tex_pref(src.pref), MT.phase_tex(src), src.kname, ','.join(src.kargs))
    L.append(head + r'\notag\\')
    for n, (sign, col, word) in enumerate(src.subterms):
        sg = '+' if sign > 0 else '-'
        line = r'&\quad (%d)\ %s\,%s\,%s' % (n + 1, sg, tex_colour(col, keep_order=True), MT.tex_word(word, src.form))
        if n == len(src.subterms) - 1:
            L.append(line + r'\Big\}' + (r'+{\rm c.c.}' if src.cc else '') + r'\label{eq:in-%s}' % src.key)
        else:
            L.append(line + r'\notag\\')
    L.append(r'\end{align}')
    L.append(r'with $%s(%s)=%s$.' % (src.kname, ','.join(src.kargs), tex_kernel(src.kexpl)))
    form = {'plain': r'Eq.~\eqref{eq:plain}', 'A{BC}': r'Eq.~\eqref{eq:Aanti}', '{BC}A': r'Eq.~\eqref{eq:antiA}'}[src.form]
    w0 = src.subterms[0][2]
    abc = (w0[2][1], w0[0][1], w0[1][1]) if src.form == '{BC}A' else (w0[0][1], w0[1][1], w0[2][1])
    L.append(r'Reordering: %s with $A$, $B$, $C$ the charges at $%s$, $%s$, $%s$.' % (form, *abc))
    if src.cc:
        L.append(r"The c.c.\ is included in the entries below (it is the same term with the phase reversed and the charges in reverse order; in common variables, $w\leftrightarrow w'$).")
    return '\n'.join(L)


def three_input(label):
    out = []
    if label == 'D':
        src = CD.D_SOURCE
        out.append(FIX['D'])
        out.append(src_header(src))
        return '\n'.join(out)
    for key in SECTION_REF[label]:
        src = S.SRC[key]
        if len(SECTION_REF[label]) > 1:
            out.append(r'\paragraph{\lab{%s}: %s.}' % (key, src.title))
        rem = MT.REMARKS.get(key, '')
        if rem:
            out.append(deref(rem))
        out.append(src_header(src))
    if label in FIX:
        out.append(FIX[label])
    return '\n'.join(out)


# ------------------------------------------------------------------ entries
def entries(book, group, keys):
    rows = []
    for n, key in enumerate(keys):
        for net in book[key]:
            c = net['groups'].get(group, 0)
            if sp.nsimplify(c) != 0:
                rows.append((n + 1, key, c, net['col']))
    return rows


def entries_table(book, group, keys, tag):
    rows = entries(book, group, keys)
    if not rows:
        return '', 0
    L = [r'\begin{longtable}{@{}>{$}r<{$} >{\raggedright\arraybackslash$}p{0.30\textwidth}<{$} >{$}r<{$} >{\raggedright\arraybackslash$}p{0.50\textwidth}<{$}@{}}',
         r'\text{S} & \text{kernel} & \text{coeff.} & \text{colour structure}\\ \hline\endhead']
    for n, key, c, col in rows:
        L.append(r'%s_{%d} & %s & %s & %s\tabularnewline' % (tag, n, kname(key), ctex(c), colour(col)))
    L.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(L), len(rows)


# ------------------------------------------------------------------ comparison tables
def match_tables(R, keys, target, tname, tag):
    out = []
    for n, key in enumerate(keys):
        d = R[key]
        has = any(target in r['contrib'] for r in d['table'])
        out.append(r'\subsubsection*{Structure $%s_{%d}$: $%s$}' % (tag, n + 1, kname(key)))
        if not has:
            out.append(r'%s has no term with this kernel: the terms cancel among themselves.' % tname)
        out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.34\textwidth}<{$}>{\raggedright\arraybackslash}p{0.42\textwidth}>{$}c<{$}>{$}c<{$}@{}}')
        out.append(r'\text{colour structure} & contributions (source: coefficient) & \text{sum} & \text{%s}\\ \hline\endhead' % tname.split('$')[0])
        for row in d['table']:
            if not row['contrib']:
                continue
            notes = [(g, v) for g, v in row['contrib'].items() if g != target]
            notes.sort(key=lambda gv: S.ORDER.index(gv[0]))
            tv = -row['contrib'].get(target, 0)
            s = sp.nsimplify(sp.expand(sum(v for _, v in notes)))
            items = r';\ '.join(r'%s:\,$%s$' % (GL[g], ctex(v, plus=False)) for g, v in notes) or r'--'
            ok = r'\ \text{\ok}' if sp.simplify(s - tv) == 0 else r'\ \text{\bad}'
            out.append(r'%s & %s & %s & %s%s\tabularnewline[2pt]' % (colour(row['col']), items, ctex(s, plus=False) if notes else '0',
                                                                     ctex(tv, plus=False), ok))
        out.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(out)


def decomp(R, keys, target, lhs, lab, tag):
    out = [r'\begin{align}', r"%s={}&\int e^{-ik\cdot(w'-w)}\int_{x,x',z}\Big\{\notag\\" % lhs]
    for n, key in enumerate(keys):
        d = R[key]
        terms = [(row['col'], -row['contrib'][target]) for row in d['table'] if target in row['contrib']]
        if not terms:
            continue
        out.append(r'&\quad +\underbrace{%s}_{%s_{%d}}\,\Big[\notag\\' % (kname(key), tag, n + 1))
        for j, (col, c) in enumerate(terms):
            end = r'\Big]' if j == len(terms) - 1 else ''
            out.append(r'&\qquad\qquad %s%s%s\notag\\' % (cpre(c), tex_colour(col), end))
    out[-1] = out[-1].rstrip('\\').rstrip().replace(r'\notag', r'\Big\}\label{%s}' % lab)
    out.append(r'\end{align}')
    return '\n'.join(out)


def relations(R, keys, tag):
    out = []
    for n, key in enumerate(keys):
        d = R[key]
        if not d['relations']:
            continue
        out.append(r'\paragraph{Structure $%s_{%d}$, $%s$.}' % (tag, n + 1, kname(key)))
        out.append(r'\begin{align*}')
        lines = []
        for k, combo in d['relations']:
            rhs = ''.join(r'%s%s' % (cpre(c, first=(j == 0)), tex_colour(col)) for j, (col, c) in enumerate(combo))
            lines.append(r'&%s\\&\qquad=%s' % (tex_colour(d['cols'][k]), rhs))
        out.append(r'\\'.join(lines))
        out.append(r'\end{align*}')
    return '\n'.join(out)


def summary_table(RA, keysA, RP, keysP):
    """which structure appears where"""
    out = []
    allk = list(keysA) + [k for k in keysP if k not in keysA]
    for key in sorted(allk, key=key_rank):
        a = 'S_{%d}' % (keysA.index(key) + 1) if key in keysA else '--'
        p = "S'_{%d}" % (keysP.index(key) + 1) if key in keysP else '--'
        j = r'\ok' if key in keysA and any('-JIMWLK' in r['contrib'] for r in RA[key]['table']) else '--'
        b = r'\ok' if key in keysP and any('-BFKL' in r['contrib'] for r in RP[key]['table']) else '--'
        nA = sum(len(r['contrib']) - ('-JIMWLK' in r['contrib']) for r in RA[key]['table']) if key in keysA else 0
        out.append(r'$%s$ & $%s$ & $%s$ & %s & %s\\' % (kname(key), a, p, j, b))
    return '\n'.join(out)


def sci(x):
    if x == 0:
        return '0'
    m, e = ('%.1e' % x).split('e')
    e = int(e)
    return m if e == 0 else r'%s\times10^{%d}' % (m, e)


def checks_tables(ch):
    pts = '\n'.join(r'%d & %s & $%.6f$ & $%.6f$ & $%s$ & $%.6f$ & $%.6f$ & $%s$\\' % (
        p['seed'], 'symmetric' if p['sym'] else 'general', p['sumA'], p['J'], sci(p['dA']), p['sumP'], p['B'], sci(p['dP']))
        for p in ch['points'])
    rows = []
    for k in S.ORDER:
        r = ch['rows'][k]
        a, p = r['A'], r['P']
        rel = 'P=A' if abs(p - a) < 1e-12 * (1 + abs(a)) else ('P=0' if abs(p) < 1e-14 else '')
        rows.append(r'%s & $%.6f$ & $%.6f$ & $%.6f$ & %s\\' % (GL[k], a, p, p - a, rel))
    sa = sum(ch['rows'][k]['A'] for k in S.ORDER)
    sp_ = sum(ch['rows'][k]['P'] for k in S.ORDER)
    p0 = [p for p in ch['points'] if p['seed'] == 3 and p['sym']][0]
    rows.append(r'\hline sum & $%.6f$ & $%.6f$ & $%.6f$ & \\' % (sa, sp_, sp_ - sa))
    rows.append(r'target & $\mathcal J=%.6f$ & $\mathcal B=%.6f$ & $\mathcal B-\mathcal J=%.6f$ & \\' % (p0['J'], p0['B'], p0['B'] - p0['J']))
    return pts, '\n'.join(rows)


# ------------------------------------------------------------------ main
def main():
    BA, RA = pickle.load(open('bookA.pkl', 'rb')), pickle.load(open('reducedA.pkl', 'rb'))
    BP, RP = pickle.load(open('bookP.pkl', 'rb')), pickle.load(open('reducedP.pkl', 'rb'))
    ch = json.load(open('checks.json'))
    keysA, keysP = ordered(RA), ordered(RP)
    T = open('static.tex').read()
    rep = {}
    nA = 0
    counts = {}
    for g in S.ORDER:
        tab, n = entries_table(BA, g, keysA, 'S')
        rep['@@TAB:%s@@' % g] = tab
        counts[g] = n
        nA += n
    nP = 0
    for g in S.ORDER:
        tab, n = entries_table(BP, g, keysP, "S'")
        rep['@@PTAB:%s@@' % g] = tab
        nP += n
    for g in S.THREE_LABELS + ['A,B']:
        rep['@@IN:%s@@' % g] = three_input(g)
    rep['@@USER4@@'] = M2.user4_list()
    rep['@@MATCHA@@'] = match_tables(RA, keysA, '-JIMWLK', r'JIMWLK$\otimes$LO', 'S')
    rep['@@MATCHP@@'] = match_tables(RP, keysP, '-BFKL', r'BFKL$\otimes$LO', "S'")
    rep['@@JDECOMP@@'] = decomp(RA, keysA, '-JIMWLK', r'\mathcal J', 'eq:Jdecomp', 'S')
    rep['@@BDECOMP@@'] = decomp(RP, keysP, '-BFKL', r'\mathcal B', 'eq:Bdecomp', "S'")
    rep['@@RELA@@'] = relations(RA, keysA, 'S')
    rep['@@RELP@@'] = relations(RP, keysP, "S'")
    rep['@@STRUCTTABLE@@'] = summary_table(RA, keysA, RP, keysP)
    pts, rows = checks_tables(ch)
    rep['@@POINTS@@'] = pts
    rep['@@ROWVALUES@@'] = rows
    rep['@@NA@@'] = str(nA)
    rep['@@NSA@@'] = str(len(keysA))
    rep['@@NSP@@'] = str(len(keysP))
    rep['@@NNETA@@'] = str(sum(len(BA[k]) for k in keysA))
    rep['@@NNETP@@'] = str(sum(len(BP[k]) for k in keysP))
    rep['@@NBASA@@'] = str(sum(len(RA[k]['table']) for k in keysA))
    rep['@@NBASP@@'] = str(sum(len(RP[k]['table']) for k in keysP))
    raw = {g: len(ts) for g, ts in S.A_sets().items()}
    rep['@@NRAW@@'] = str(sum(raw.values()))
    rep['@@COUNTTABLE@@'] = '\n'.join(r'%s & %s & %d & %d\\' % (GL[g], TITLE[g], raw[g], counts[g]) for g in S.ORDER)
    for k, v in rep.items():
        assert k in T, k
        T = T.replace(k, v)
    left = re.findall(r'@@[^@]*@@', T)
    assert not left, left
    open('TwoRho_all.tex', 'w').write(T)
    print('written TwoRho_all.tex; entries A', nA, 'P', nP)


if __name__ == '__main__':
    main()
