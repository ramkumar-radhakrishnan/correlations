r"""Writes Rest_terms.tex: every term of the cross section that is neither JIMWLK (log(k+/Lambda)) nor BFKL (log(vee/k+)),
row by row in your notation, all black.

It is Remaining_rows.tex without the colour marks: removed lines (red) are left out, changed parts (blue) are kept as they
are, and the two subtractions of a BFKL term hidden in line 1 are written as integrals,
    -(B)                  = + ... \int_0^{vee-k+} dp+/(2 pi) 1/(p+ + k+) ...      (Group III Rows III, V, VI)
    -log(vee/k+) H4       = - ... \int_{k+/vee}^1 dxi/xi ...                       (Group II Row IV)
using \int_0^{vee-k+} dp+/(p+ + k+) = \int_{k+/vee}^1 dxi/xi = log(vee/k+) exactly (see Why_log_vee.tex)."""
import os
import re
import make_bfkl as MB
import rows as R
import assemble as AS

HERE = MB.HERE
UO, EX, HB = MB.UO, MB.EX, MB.HB
LOGV = r'\log{\left(\frac{\vee}{k^{+}}\right)}'
LOGV_H4 = r'\log{\left[\frac{\vee}{k^{+}}\right]}'


def _end(t, i):
    """index after the brace group that opens at t[i] == '{' (escaped braces \{ \} do not count)"""
    depth, j = 0, i
    while True:
        c = t[j]
        if c == '\\':
            j += 2
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return j + 1
        j += 1


def strip_marks(t):
    while r'\old{' in t:
        i = t.index(r'\old{')
        t = t[:i] + t[_end(t, i + 4):]
    while r'\new{' in t:
        i = t.index(r'\new{')
        j = _end(t, i + 4)
        t = t[:i] + t[i + 5:j - 1] + t[j:]
    return t


def b_integral(t):
    """-(B) with log(vee/k+) written as (1/2pi) \int_0^{vee-k+} dp+ / (p+ + k+)"""
    t = MB.neg(t)
    for a, b in ((r'{16\pi^{5}}', r'{8\pi^{4}}'), (r'{8\pi^{5}}', r'{4\pi^{4}}')):
        if a in t:
            assert t.count(a) == 1
            t = t.replace(a, b)
            break
    else:
        raise ValueError(t[:80])
    assert t.count(LOGV) == 1
    return t.replace(LOGV, r'\int_{0}^{\vee-k^{+}}\frac{dp^{+}}{2\pi}\,\frac{1}{(p^{+}+k^{+})}\,')


def h4_integral():
    t = re.sub(r'^\s*&?\s*\+?\s*', '', MB.UB['G2 IV'].terms[1])
    assert t.count(LOGV_H4) == 1
    return '-' + t.replace(LOGV_H4, r'\int_{\frac{k^{+}}{\vee}}^{1}\frac{d\xi}{\xi}')


def rest_items(name):
    items = []
    for (t, style) in MB.remaining_items(name):
        if style == 'red':
            continue
        if style == 'blue':          # the two subtractions (the remainders of the logarithms carry style None)
            continue
        items.append((strip_marks(t), None))
    if name in HB:
        items.append((b_integral(HB[name]), None))
    if name == 'G2 IV':
        items.append((h4_integral(), None))
    return items


NOTE = {
    'G1 II(1)': r'''As in your file. $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}\to0$ for $\vee\to\infty$. The pole bracket times $\frac{11}{3}$ is the $\beta_0$ (running-coupling) term.''',
    'G1 II(2b)': r'''As in your file, with $R^i$ below.''',
    'G1 III N1': r'''The logarithms are what is left after their $\log\frac{\vee}{k^+}$ parts are taken out (those are in \texttt{BFKL\_rows.pdf}), with $\Lambda\to0$.''',
    'G1 III N2': r'''As $d^3N_1$.''',
    'G1 III N3a': r'''$\Lambda\to0$; $\mathcal A_0$ is $\mathcal A$ at $\Lambda=0$, Eq.~\eqref{eq:A0}.''',
    'G1 III N3b': r'''$\Lambda\to0$. Without the last line of your $d^3N_{3b}$: for $\Lambda\to0$ it is $\log\frac{(z-w)^2}{(x-w)^2}$ times the edge kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$, which is not integrable at $z\to w$. Its partner $\propto\log\frac{k^+}{\Lambda}$ was replaced by $\mathcal D$; its finite part belongs to the $\mathcal D$ region, which is still open (C4).''',
    'G1 IV': r'''The two logarithms are the remainders.''',
    'G1 V(1)': 'As in your file.', 'G1 V(2a)': 'As in your file.',
    'G1 V(2b UV)': r'''$2\log(1-\frac{\Lambda}{k^+})\to0$. The $-\frac{11}{6}$ is the $\beta_0$ term that cancels the $\frac{11}{3}$ of Row II part 1 (running coupling).''',
    'G1 V(2b fin)': r'''$\Lambda\to0$ in the limits of the plus-prescribed integrals.''',
    'G2 I(1) A': r'''Remainders. Your line 2 is not here: $\log\frac{\vee}{\Lambda+k^+}\to\log\frac{\vee}{k^+}$, so it is entirely BFKL.''',
    'G2 I(1) B': r'''As the first line.''',
    'G2 I(2a)': r'''$\Lambda\to0$ in the second logarithm.''',
    'G2 I(2b)': 'Remainders.', 'G2 I(2c)': 'Remainders.', 'G2 I(2d)': 'Remainders.', 'G2 I(2e)': 'Remainders.',
    'G2 I(2f)': r'''Remainders; $\frac{1}{k^++\Lambda}\to\frac{1}{k^+}$.''',
    'G2 IV': r'''The DGLAP term. Line 2 keeps $\log(\frac{\vee}{k^+}-1)$: it cancels the $\log\frac{\vee}{k^+}$ of the subtraction in $\frac{\xi}{(1-\xi)_+}$ of line 1. The last line subtracts the BFKL term $H_4$ that sits in the $\frac{1-\xi}{\xi}$ term of line 1 ($\int_{k^+/\vee}^1\frac{d\xi}{\xi}=\log\frac{\vee}{k^+}$; \texttt{Why\_log\_vee.pdf}, Sec.~10). In line 3, $\frac{k^+}{k^++\Lambda}\to1$.''',
    'G3 I': 'Remainders.', 'G3 II': 'Remainders.',
    'G3 III(1)': r'''Line 2 keeps $\log(\frac{\vee}{k^+}-1)$: it cancels the $\log\frac{\vee}{k^+}$ of the subtraction in $\left[\frac{1}{p^+}\right]_+$ of line 1. The last line subtracts the BFKL term (B) that sits in the $-\frac{1}{k^+}$ term of line 1 ($\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$; \texttt{Why\_log\_vee.pdf}). Together the lines are finite for $\vee\to\infty$.''',
    'G3 III(2)': r'''As the first part; here (B) sits in the $-\frac{p^+}{k^{+2}}$ term of line 1. The measure is $\int_{w,x',y',z,x}$ (C5).''',
    'G3 V': r'''As Row III.''', 'G3 VI': r'''As Row III.''',
}


def rows_tex():
    out = []
    for gname, gl, rows in MB.GROUPS:
        out.append(r'\section{%s}\label{sec:R%s}' % (gname, gl))
        for rname, names, special in rows:
            out.append(r'\subsection{%s}' % rname)
            if special == 'four':
                out.append(r'Only four $\rho$\textquotesingle s: nothing is left.')
                continue
            for name in names:
                sub = MB.SUB[name]
                out.append(r'\paragraph{%s%s.}' % (rname, (', ' + sub) if sub else ''))
                if name not in UO:
                    out.append(r'Nothing is left: this part has only $\log\frac{k^+}{\Lambda}$ and $\log\frac{\vee}{k^+}$ terms.')
                    continue
                out.append(NOTE[name])
                out.append(MB.align(UO[name], rest_items(name), MB.lab('R', name)))
                if name == 'G1 II(2b)':
                    out.append('where\n\\begin{align}' + R.align_text('Other_terms', 50, 55) + '\\end{align}')
    return '\n'.join(out)


DOC = MB.PREAMBLE + r'''
\title{The rest of the terms: everything apart from JIMWLK and BFKL}
\date{}
\begin{document}
\maketitle
\noindent Every term of the cross section that is neither JIMWLK ($\propto\log\frac{k^+}{\Lambda}$, your JIMWLK file) nor BFKL ($\propto\log\frac{\vee}{k^+}$, \texttt{BFKL\_rows.pdf}), row by row in your notation, with $\Lambda\to0$ and $\vee$ exact. It is \texttt{Remaining\_rows.pdf} without the colour marks; row by row it equals \texttt{Finite\_terms\_rows.pdf} ($@@COK@@$). So
\begin{equation*}
\frac{d^3N_{\rm NLO}}{d^3k}=\log\frac{k^+}{\Lambda}\,\big[\text{JIMWLK}\otimes\text{LO}\big]+\log\frac{\vee}{k^+}\,\big[\text{BFKL}\otimes\text{LO}\big]+\big[\text{this file}\big].
\end{equation*}
\begin{itemize}
\item No term here grows with $\vee$. Two logarithms of $\vee$ are still written out, and both are cancelled inside their own row:
\begin{itemize}
\item $\log(\frac{\vee}{k^+}-1)$ in line 2 of Group III Rows III, V, VI and Group II Row IV. It cancels the $\log\frac{\vee}{k^+}$ hidden in the subtraction of $\left[\frac{1}{p^+}\right]_+$ (resp.\ $\frac{\xi}{(1-\xi)_+}$) in line 1 of the same row.
\item $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}$ in the last line of Group III Rows III, V, VI, and $\int_{k^+/\vee}^1\frac{d\xi}{\xi}$ in Group II Row IV. Each removes the BFKL term hidden in line 1 of its row (it is in the BFKL list).
\end{itemize}
Both are explained step by step in \texttt{Why\_log\_vee.pdf}.
\item $\log(1-\frac{k^+}{\vee})\to0$ for $\vee\to\infty$.
\item Running coupling: Group I Row II part 1 ($\frac{11}{3}$) and Row V 2b UV ($-\frac{11}{6}$). DGLAP: Group II Row IV. See \texttt{Finite\_terms\_final.pdf}, Secs.~4 and 5.
\item Open: the finite part of the $\mathcal D$ region (C4, Group I Row III $d^3N_{3b}$).
\end{itemize}
\begin{align}
@@A0@@
\label{eq:A0}
\end{align}

@@ROWS@@
\end{document}
'''

if __name__ == '__main__':
    doc = DOC.replace('@@COK@@', MB.tex_num(MB.C_ok)).replace('@@A0@@', MB.A0).replace('@@ROWS@@', rows_tex())
    assert r'\new{' not in doc.split(r'\begin{document}')[1] and r'\old{' not in doc.split(r'\begin{document}')[1]
    open(os.path.join(HERE, 'Rest_terms.tex'), 'w').write(doc)
    print('written Rest_terms.tex')
