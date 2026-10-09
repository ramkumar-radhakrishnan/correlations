"""Generate ThreeRho_cancellation.tex: the explicit check that the symmetrized three-rho terms cancel.

Input:  cancel_data.py (entries, classes, colour networks), crosscheck.json (crosscheck.py),
        pointwise.json (pointwise_check.py).
"""
import json
import sympy as sp
import cancel_data as CD
from sources import SOURCES
from engine import tex_colour
from kernels import tex_factor

CLS_TEX = {'alpha': r'\alpha', "alpha'": r"\alpha'", 'beta': r'\beta', "beta'": r"\beta'",
           'gamma': r'\gamma', "gamma'": r"\gamma'"}
CLS_KERNEL = {
    'alpha': r"\KK(x-w)\cdot\KK(x'-w')\;\KK(x-z)\cdot\KK(y-z)",
    "alpha'": r"\KK(x-w)\cdot\KK(x'-w')\;\KK(x'-z)\cdot\KK(y-z)",
    'beta': r"\KK(x-w)\cdot\KK(x'-w')\;\KK(w-z)\cdot\KK(y-z)",
    "beta'": r"\KK(x-w)\cdot\KK(x'-w')\;\KK(w'-z)\cdot\KK(y-z)",
    'gamma': r"\KK(x'-w')\cdot\KK(y-z)\;\Phi_{w}(x)",
    "gamma'": r"\KK(x-w)\cdot\KK(y-z)\;\Phi_{w'}(x')",
}
CLS_WORDS = {
    'alpha': r'soft gluon between the charges $x$ and $y$',
    "alpha'": r"soft gluon between the charges $x'$ and $y$",
    'beta': r'soft gluon between the measured gluon at $w$ and the charge $y$',
    "beta'": r"soft gluon between the measured gluon at $w'$ and the charge $y$",
    'gamma': r'the edge kernel of Row III, part 2',
    "gamma'": r'the edge kernel of Row III, part 2 (c.c.)',
}
SEC_OF = [('Three-$\\rho$ remainders of the four-$\\rho$ terms (no c.c.)', ['A1', 'A2', 'B1', 'B2']),
          ('Group I (with c.c.)', ['GI.II', 'GI.III', 'GI.IV', 'GI.V1', 'GI.V2']),
          ('Group II (no c.c.)', ['GII.I1', 'GII.I2']),
          ('Group III (with c.c.)', ['GIII.I', 'GIII.II', 'GIII.III', 'GIII.V', 'GIII.VI'])]
FORM_TEX = {'plain': r'ordered product $\rho\rho\rho\to S(\rho\rho\rho)$, factor $1$',
            'A{BC}': r'$\rho\{\rho,\rho\}\to 2S(\rho\rho\rho)$, factor $2$',
            '{BC}A': r'$\{\rho,\rho\}\rho\to 2S(\rho\rho\rho)$, factor $2$'}
NATIVE_NOTE = {
    'GI.II': r"The kernel $\KK^i(x'-w')\,\mathfrak F^i(w;x,y)$ is written with the $z$-integral representation of $\mathfrak F^i$ (previous note, App.~A): four monomials, a--d.",
    'GI.III': r"$\Phi(x)$ is the edge kernel of the notes (it depends on $x$, $z$ and $w$); monomial b is part 2 of the row.",
    'GI.V1': r"The measured gluon of the amplitude sits at $z$ and the soft gluon at $x$; we first rename $z\to w$, $x\to z$.",
    'GI.V2': r"The measured gluon of the amplitude sits at $x$; we first rename $x\to w$.",
    'GIII.III': r"The measured gluon of the conjugate amplitude sits at $y'$; we first rename $y'\to w'$.",
    'GIII.V': r"Phase $e^{-ik\cdot(y'-w)}$; we first rename $y'\to w'$.",
    'GIII.VI': r"Phase $e^{-ik\cdot(y'-w)}$; we first rename $y'\to w'$.",
}


def ctex(c, plus=True):
    s = sp.latex(sp.nsimplify(c)).replace(' ', '')
    s = s.replace(r'\frac{i}{2}', r'\tfrac{i}{2}')
    if plus and not s.startswith('-'):
        s = '+' + s
    return s


def lab(e, colour=None):
    t = r'\lab{%s}' % e['label']
    return r'\textcolor{%s}{%s}' % (colour, t) if colour else t


def mono_native_tex(kc, facs):
    body = r'\,'.join(tex_factor(f) for f in facs)
    kc = sp.nsimplify(kc)
    if kc == 1:
        return body
    if kc == -1:
        return '-' + body
    return sp.latex(kc).replace(r'\frac{1}{2}', r'\tfrac12') + r'\,' + body


def mono_canon_tex(facs):
    out = []
    for f in facs:
        if f[0] == 'KK':
            out.append(r'\KK(%s-%s)\cdot\KK(%s-%s)' % f[1:])
        else:
            out.append(r'\Phi_{%s}(%s)' % (f[2], f[1]))
    return r'\,'.join(out)


def rename_tex(full, src):
    labels = list(dict.fromkeys([p for _, p in src.subterms[0][2]] + list(src.fixed_labels)))
    moves = [(a, full[a]) for a in labels if a in full and full[a] != a]
    if not moves:
        return r'\text{none}'
    return r',\ '.join('%s\\to %s' % mv for mv in moves)


def source_block(src):
    """LaTeX for one source: input, renaming of every kernel monomial, and its entries."""
    entries, ren = CD.source_entries(src)
    monos = __import__('classes3').monomials(src)
    L = []
    L.append(r'\subsection{\lab{%s}: %s}\label{src:%s}' % (src.key, src.title, src.key))
    P1, P0 = src.phase
    L.append(r'Prefactor $%s\NN$; %s; %s; phase $e^{-ik\cdot(%s-%s)}$.' % (
        ctex(src.pref, plus=False), FORM_TEX[src.form], 'with c.c.' if src.cc else 'no c.c.', P1, P0))
    if src.key in NATIVE_NOTE:
        L.append(NATIVE_NOTE[src.key])
    if len(monos) == 1:
        L.append(r'Kernel: $K=%s$.' % mono_native_tex(*monos[0]))
    else:
        L.append(r'Kernel $K=\sum$ of the monomials')
        L.append(r'\begin{center}\small\begin{tabular}{@{}l >{$}l<{$}@{}}')
        for m, (kc, f) in enumerate(monos):
            L.append(r'%s: & %s\\' % (CD.LETTERS[m], mono_native_tex(kc, f)))
        L.append(r'\end{tabular}\end{center}')
    L.append(r'The symmetrized term (commuting charges):')
    L.append(r'\begin{align*}')
    L.append(r'&%s\NN\int e^{-ik\cdot(%s-%s)}\,K\,\Big\{\\' % (
        ctex(src.pref * (1 if src.form == 'plain' else 2), plus=False), P1, P0))
    for n, (sg, col, word) in enumerate(src.subterms):
        last = n == len(src.subterms) - 1
        L.append(r'&\qquad (%d)\ \ %s\,%s%s' % (n + 1, '+' if sg > 0 else '-',
                                                 tex_colour(col + [('rho', i, p) for i, p in word]),
                                                 (r'\Big\}' + (r'+{\rm c.c.}' if src.cc else '')) if last else r'\\'))
    L.append(r'\end{align*}')
    L.append(r'\emph{Renaming to common variables} (all replacements simultaneous; $\varepsilon$ is the sign from $\KK(-X)=-\KK(X)$):')
    L.append(r'\begin{center}\small\begin{tabular}{@{}l l >{$}l<{$} >{$}l<{$} >{$}c<{$} >{$}c<{$}@{}}')
    L.append(r'mon. & term & \text{renaming} & \text{kernel in common variables} & \varepsilon & \text{class}\\ \hline')
    for (cc, m), r in ren.items():
        L.append(r'%s & %s & %s & %s & %s & %s\\' % (
            CD.LETTERS[m] if len(monos) > 1 else '', 'c.c.' if cc else 'amp.', rename_tex(r['full'], src),
            mono_canon_tex(r['canon']), '+' if r['ksign'] > 0 else '-', CLS_TEX[r['cls']]))
    L.append(r'\end{tabular}\end{center}')
    L.append(r'\emph{Entries} (coefficient in units of $\NN$; colour structure in common variables):')
    L.append(r'\begin{center}\small\begin{longtable}{@{}l >{$}c<{$} >{$}r<{$} >{$}l<{$} l@{}}')
    L.append(r'entry & \text{class} & \text{coeff.} & \text{colour structure} & \\ \hline\endhead')
    for e in entries:
        rule = {('UU', 'delta'): '(U1)', ('UU(sym)', 'delta'): r'(U1$^\star$)'}.get(tuple(e['rules']), '')
        L.append(r'%s & %s & %s & %s & %s\\' % (lab(e), CLS_TEX[e['cls']], ctex(e['coeff']), tex_colour(e['col']), rule))
    L.append(r'\end{longtable}\end{center}\addtocounter{table}{-1}')
    return '\n'.join(L)


def contrib(e, s, colour=None):
    c = ctex(e['coeff'])
    t = r'%s\,$(%s)%s$' % (lab(e), c, '' if s > 0 else '(-1)')
    return r'\textcolor{%s}{%s}' % (colour, t) if colour else t


def is_deleted(e):
    return CD.deleted_by_correction1(e) or CD.deleted_by_correction2(e)


def network_table(cls, nets, mode):
    """mode 'written': all networks, members as written (deleted ones in red);
       mode 'corrected': only networks changed by the corrections."""
    L = [r'\begin{longtable}{@{}>{\raggedright}p{0.05\textwidth}>{\raggedright$}p{0.36\textwidth}<{$}>{\raggedright}p{0.42\textwidth}>{$}r<{$}@{}}',
         r'& \text{colour structure } C & entries: coefficient $\times$ sign & \text{sum}\\ \hline\endhead']
    for net in nets:
        written = [(e, s) for e, s in net['members'] if e['src'] != 'D']
        corrected = [(e, s) for e, s in net['members'] if not is_deleted(e)]
        if mode == 'written':
            mem = written
            items = [contrib(e, s, 'red' if is_deleted(e) else None) for e, s in mem]
        else:
            if set(id(e) for e, _ in written) == set(id(e) for e, _ in corrected):
                continue
            mem = corrected
            items = [contrib(e, s, 'blue' if e['src'] == 'D' else None) for e, s in mem]
        tot = sp.nsimplify(sum(e['coeff'] * s for e, s in mem))
        if not mem:
            items = [r'\emph{no entries left}']
        tot_tex = r'0\ \text{\ok}' if tot == 0 else r'\textcolor{red}{%s\ \text{\bad}}' % ctex(tot)
        L.append(r'$%s$ & %s & %s & %s\tabularnewline[2pt]' % (net['name'], tex_colour(net['rep']['col']).replace(r'\,', r'\,\allowbreak '),
                                                              r';\ '.join(items), tot_tex))
    L.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(L)


def build_networks():
    written = CD.all_entries()
    dents = CD.source_entries(CD.D_SOURCE)[0]
    res = CD.numbered(CD.networks(written + dents))
    for nets in res.values():
        for net in nets:
            w = [(e, s) for e, s in net['members'] if e['src'] != 'D']
            c = [(e, s) for e, s in net['members'] if not is_deleted(e)]
            net['tot_written'] = sp.nsimplify(sum(e['coeff'] * s for e, s in w)) if w else None
            net['tot_corrected'] = sp.nsimplify(sum(e['coeff'] * s for e, s in c)) if c else None
            net['n_written'], net['n_corrected'] = len(w), len(c)
    return written, dents, res


def fmt(x):
    m, e = ('%.1e' % x).split('e')
    return r'%s\times10^{%d}' % (m, int(e))


def main():
    written, dents, res = build_networks()
    cc = json.load(open('crosscheck.json'))
    pw = json.load(open('pointwise.json'))
    n_entries = len(written)
    nets_w = [n for nets in res.values() for n in nets if n['n_written']]
    n_net = len(nets_w)
    n_bad = sum(1 for n in nets_w if n['tot_written'] != 0)
    nets_c = [n for nets in res.values() for n in nets if n['n_corrected']]
    n_net_c = len(nets_c)
    n_bad_c = sum(1 for n in nets_c if n['tot_corrected'] != 0)
    assert n_bad_c == 0
    T = open('cancel_static.tex').read()
    rep = {}
    rep['@@NENTRIES@@'] = str(n_entries)
    rep['@@NNET@@'] = str(n_net)
    rep['@@NOK@@'] = str(n_net - n_bad)
    rep['@@NBAD@@'] = str(n_bad)
    rep['@@NNETC@@'] = str(n_net_c)
    rep['@@NENTRIESC@@'] = str(len(written) - sum(1 for e in written if is_deleted(e)) + len(dents))
    rep['@@NDEL1@@'] = str(sum(1 for e in written if CD.deleted_by_correction1(e)))
    rep['@@NDEL2@@'] = str(sum(1 for e in written if CD.deleted_by_correction2(e)))
    rep['@@ND@@'] = str(len(dents))
    # summary table per class
    rows = []
    for cls, nets in res.items():
        w = [n for n in nets if n['n_written']]
        c = [n for n in nets if n['n_corrected']]
        rows.append(r'$%s$ & %s & %d & %d & %d & %d & %d\\' % (
            CLS_TEX[cls], CLS_KERNEL[cls], sum(n['n_written'] for n in w), len(w),
            sum(1 for n in w if n['tot_written'] == 0), sum(1 for n in w if n['tot_written'] != 0),
            len(c)))
    rep['@@SUMMARYROWS@@'] = '\n'.join(rows)
    # per-source blocks
    blocks = []
    srcs = {s.key: s for s in SOURCES}
    for title, keys in SEC_OF:
        blocks.append(r'\subsection*{%s}' % title)
        blocks.append(r'\addcontentsline{toc}{subsection}{%s}' % title)
        for k in keys:
            blocks.append(source_block(srcs[k]).replace(r'\subsection{', r'\subsubsection{'))
    rep['@@SOURCEBLOCKS@@'] = '\n\n'.join(blocks)
    # network tables as written
    tabs = []
    for cls, nets in res.items():
        nets = [n for n in nets if n['n_written']]
        nbad = sum(1 for n in nets if n['tot_written'] != 0)
        tabs.append(r'\subsection{Class $%s$: %s}' % (CLS_TEX[cls], CLS_WORDS[cls]))
        tabs.append(r'Kernel $K_{%s}=%s$. %d entries, %d colour networks; %d cancel, %d do not.' % (
            CLS_TEX[cls], CLS_KERNEL[cls], sum(n['n_written'] for n in nets), len(nets), len(nets) - nbad, nbad))
        tabs.append(network_table(cls, nets, 'written'))
    rep['@@NETWORKTABLES@@'] = '\n'.join(tabs)
    # diagnosis table: the networks that do not cancel
    drows = []
    for cls, nets in res.items():
        for n in nets:
            if not n['n_written'] or n['tot_written'] == 0:
                continue
            keep = [(e, s) for e, s in n['members'] if e['src'] != 'D' and not is_deleted(e)]
            d1 = [(e, s) for e, s in n['members'] if CD.deleted_by_correction1(e)]
            d2 = [(e, s) for e, s in n['members'] if CD.deleted_by_correction2(e)]
            dd = [(e, s) for e, s in n['members'] if e['src'] == 'D']
            sm = lambda L: ctex(sum(e['coeff'] * s for e, s in L), plus=False) if L else r'\text{--}'
            drows.append(r'$%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$\\' % (
                n['name'], ctex(n['tot_written'], plus=False), sm(keep), sm(d1), sm(d2), sm(dd),
                ctex(n['tot_corrected'], plus=False) if n['n_corrected'] else r'\text{empty}'))
    rep['@@DIAGROWS@@'] = '\n'.join(drows)
    # D entries
    rep['@@DBLOCK@@'] = source_block(CD.D_SOURCE).replace(r'\subsection{', r'\subsection*{')
    # deleted entries of correction 1
    del1 = [e for e in written if CD.deleted_by_correction1(e)]
    rep['@@DEL1LIST@@'] = r',\ '.join(r'%s\,$(%s)$' % (lab(e), ctex(e['coeff'])) for e in del1)
    del2 = [e for e in written if CD.deleted_by_correction2(e)]
    rep['@@DEL2LIST@@'] = r',\ '.join(r'%s\,$(%s)$' % (lab(e), ctex(e['coeff'])) for e in del2)
    # corrected tables
    tabs = []
    for cls, nets in res.items():
        changed = [n for n in nets if n['n_written'] != n['n_corrected'] or any(e['src'] == 'D' for e, _ in n['members'])]
        if not changed:
            tabs.append(r'\paragraph{Class $%s$.} No entry deleted or added; all %d networks cancel as in Sec.~\ref{sec:networks}.' % (
                CLS_TEX[cls], len(nets)))
            continue
        if all(n['n_corrected'] == 0 for n in changed):
            tabs.append(r'\paragraph{Class $%s$.} All %d entries are deleted; the class is empty.' % (
                CLS_TEX[cls], sum(n['n_written'] for n in changed)))
            continue
        tabs.append(r'\paragraph{Class $%s$.} %d of the %d networks change:' % (CLS_TEX[cls], len(changed), len([n for n in nets if n['n_written']])))
        tabs.append(network_table(cls, changed, 'corrected'))
    rep['@@CORRTABLES@@'] = '\n'.join(tabs)
    # numbers from the cross-checks
    mx = lambda L: max(a for a, b in L)
    sc = lambda L: max(b for a, b in L)
    rep['@@FOURRHO@@'] = fmt(max(a for a, b in cc['four_rho']))
    rep['@@FOURRHOSCALE@@'] = '%.0f' % max(b for a, b in cc['four_rho'])
    rep['@@ABN4@@'] = fmt(mx(cc['AB_vs_N4']))
    rep['@@ABN4SCALE@@'] = '%.1f' % sc(cc['AB_vs_N4'])
    rowtex = []
    names = {'G2-I': 'Group II, Row I (both terms + its four-$\\rho$ part)', 'G2-II': 'Group II, Row II (four-$\\rho$ part only)',
             'G2-III': 'Group II, Row III (four-$\\rho$ part only)', 'G3-I': 'Group III, Row I', 'G3-II': 'Group III, Row II',
             'G3-III': 'Group III, Row III', 'G3-IV': 'Group III, Row IV (four-$\\rho$ part only)', 'G3-V': 'Group III, Row V',
             'G3-VI': 'Group III, Row VI'}
    for k, v in cc['rows'].items():
        rowtex.append(r'%s & $%s$ & %.1f & %s\\' % (names[k], fmt(mx(v)), sc(v), 'agree' if sc(v) > 1e-10 else 'both vanish'))
    p = cc['group1_pieces']
    g1n = {'norm': 'Group I, Row I (normalization; four-$\\rho$ part)',
           'Abar^dag UU B2 (Row IV direct)': 'Group I, Row IV, sub-terms 1, 3',
           'Cbar^dag UU B2 (Row V direct)': 'Group I, Row V, sub-terms 1, 3',
           'U A3          (Row II direct)': 'Group I, Row II, sub-terms 1, 3 ($U(w)\\mathbb A^{(3)}$)',
           'B3bar^dag U A (Row III direct)': 'Group I, Row III, sub-terms 1, 3',
           'unitarity combination for Abar(3)dag': 'Group I, sub-terms 2, 4 of Rows II--V (together)'}
    for k in ['norm', 'Abar^dag UU B2 (Row IV direct)', 'Cbar^dag UU B2 (Row V direct)',
              'U A3          (Row II direct)', 'B3bar^dag U A (Row III direct)', 'unitarity combination for Abar(3)dag']:
        v = p[k]
        ok = mx(v) < 1e-10
        rowtex.append(r'%s & $%s$ & %.1f & %s\\' % (g1n[k], fmt(mx(v)), sc(v), 'agree' if ok else r'\textcolor{red}{differ}'))
    g = cc['group1']
    rowtex.append(r'Group I, all rows, as written & $%s$ & %.1f & \textcolor{red}{differ}\\' % (
        fmt(max(x['written'] for x in g)), max(x['scale'] for x in g)))
    rep['@@XROWS@@'] = '\n'.join(rowtex)
    a = cc['group1_after']
    after = []
    an = {'II_direct': r'Row II, sub-terms 1, 3, with correction 1 & $U(w)\mathbb A^{(3)}$',
          'III_direct': r'Row III, sub-terms 1, 3, part 1 only & $\bar{\mathbb B}_3^\dagger U\mathbb A$',
          'D_direct': r"$\mathcal D$, sub-terms 1, 3 & row I$'$ (soft gluon between $w$ and $U\rho$)",
          'unitarity': r'sub-terms 2, 4 of Rows II--V and of $\mathcal D$ & $\bar{\mathbb A}^{(3)\dagger}$ (unitarity)'}
    for k, v in a.items():
        after.append(r'%s & $%s$ & %.1f\\' % (an[k], fmt(mx(v)), sc(v)))
    after.append(r'Group I, all rows, corrected & all & $%s$ & %.1f\\' % (fmt(max(x['corrected'] for x in g)), max(x['scale'] for x in g)))
    rep['@@AFTERROWS@@'] = '\n'.join(after)
    rep['@@G1VARIANT@@'] = fmt(max(x['corrected_variant'] for x in g))
    # pointwise
    prow = []
    for k in range(len(pw['symmetric'])):
        s, gen = pw['symmetric'][k], pw['general'][k]
        r = lambda d, n: d[n][0] / d[n][1]
        prow.append(r'%d & $%s$ & $%s$ & $%s$ & $%s$ & $%s$\\' % (
            k + 1, fmt(r(s, 'written')), fmt(r(s, 'corrected')), fmt(r(s, 'corrected_variant')),
            fmt(r(gen, 'written')), fmt(r(gen, 'corrected'))))
    rep['@@POINTROWS@@'] = '\n'.join(prow)
    for k, v in rep.items():
        assert k in T, k
        T = T.replace(k, v)
    assert '@@' not in T
    open('ThreeRho_cancellation.tex', 'w').write(T)
    print('written ThreeRho_cancellation.tex:', n_entries, 'entries,', n_net, 'networks,', n_bad, 'not cancelling;',
          'corrected:', n_net_c, 'networks,', n_bad_c, 'not cancelling')


if __name__ == '__main__':
    main()
