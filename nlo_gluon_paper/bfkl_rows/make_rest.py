r"""Writes Rest_terms.tex: every term of the cross section that is neither JIMWLK (log(k+/Lambda)) nor BFKL (an explicit
log(vee/k+)), row by row in your notation, all black.

It is Remaining_rows.tex without the colour marks: removed lines (red) are left out, changed parts (blue) are kept as they
are.  Only vee > k+, p+ > Lambda are used; Lambda -> 0 where possible, vee exact; the p+ and xi integrals are kept as
written (no large-p+ step, nothing subtracted from them)."""
import os
import re
import make_bfkl as MB
import rows as R
import assemble as AS

HERE = MB.HERE
UO, EX = MB.UO, MB.EX


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


def rest_items(name):
    items = []
    for (t, style) in MB.remaining_items(name):
        if style == 'red':
            continue
        assert style is None, (name, style)
        items.append((strip_marks(t), None))
    return items


NOTE = {
    'G1 II(1)': r'''As in your file. $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ stays here: its coefficient depends on $\vee$, and $0\le-\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}\le\log^22$ for all $\vee>k^+$ (maximum at $\vee=2k^+$). The pole bracket times $\frac{11}{3}$ is the $\beta_0$ (running-coupling) term.''',
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
    'G2 IV': r'''Line 2: $\log(1-\frac{k^+}{\vee})$, the remainder of $\log(\frac{\vee}{k^+}-1)$ (its $\log\frac{\vee}{k^+}$ part is in the BFKL list). Line 1, the $\xi$ integral from $\frac{k^+}{\vee}$ with $\frac{\xi}{(1-\xi)_+}$, is as you wrote it. In line 3, $\frac{k^+}{k^++\Lambda}\to1$.''',
    'G3 I': 'Remainders.', 'G3 II': 'Remainders.',
    'G3 III(1)': r'''Line 2: $\log(1-\frac{k^+}{\vee})$, the remainder of $\log(\frac{\vee}{k^+}-1)$ (its $\log\frac{\vee}{k^+}$ part is in the BFKL list). Line 1, the $p^+$ integral up to $\vee-k^+$ with $\left[\frac{1}{p^+}\right]_+$, is as you wrote it.''',
    'G3 III(2)': r'''As the first part. The measure is $\int_{w,x',y',z,x}$ (C5).''',
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
\noindent Every term of the cross section that is neither JIMWLK ($\propto\log\frac{k^+}{\Lambda}$, your JIMWLK file) nor BFKL (an explicit $\log\frac{\vee}{k^+}$, \texttt{BFKL\_rows.pdf}), row by row in your notation. It is \texttt{Remaining\_rows.pdf} without the colour marks. Row by row it is the remainder of the exact split of your rows (check C, $@@COK@@$). So
\begin{equation*}
\frac{d^3N_{\rm NLO}}{d^3k}=\log\frac{k^+}{\Lambda}\,\big[\text{JIMWLK}\otimes\text{LO}\big]+\log\frac{\vee}{k^+}\,\big[\text{the list of \texttt{BFKL\_rows.pdf}}\big]+\big[\text{this file}\big].
\end{equation*}
\begin{itemize}
\item \textbf{Rule.} The only constraints are $\vee>k^+$ and $p^+>\Lambda$; $\Lambda\to0$ wherever possible, $\vee$ exact. Only explicit logarithms are split. Every $p^+$ and $\xi$ integral is kept as you wrote it, with no large-$p^+$ step and nothing subtracted from it.
\item \textbf{Withdrawn from the previous version:}
\begin{itemize}
\item the subtractions $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}\times$(B) (Group III Rows III, V, VI) and $\int_{k^+/\vee}^1\frac{d\xi}{\xi}\,H_4$ (Group II Row IV);
\item keeping $\log(\frac{\vee}{k^+}-1)$ whole in line 2 of these rows. It is now $\log(1-\frac{k^+}{\vee})$, as in every other row.
\end{itemize}
\item \textbf{$\vee$ dependence here:}
\begin{itemize}
\item the $p^+$ integrals up to $\vee-k^+$ (line 1 of Group III Rows III, V, VI);
\item the $\xi$ integrals from $\frac{k^+}{\vee}$ (lines 1 and 3 of Group II Row IV);
\item $\log(1-\frac{k^+}{\vee})$;
\item $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ (Group I Row II part 1; $\le\log^22$);
\item $R^i$ (Group I Row II part 2).
\end{itemize}
\item \textbf{Running coupling:} Group I Row II part 1 ($\frac{11}{3}$) and Row V 2b UV ($-\frac{11}{6}$). \textbf{Poles:} Group I Row II part 1 and Row V 2b UV, Group II Row IV lines 1 and 2. The pole $\times\log\frac{\vee}{k^+}$ of Group II Row IV line 2 is in the BFKL list (\texttt{BFKL\_rows.pdf}, section ``Poles'').
\item \textbf{Open:} the finite part of the $\mathcal D$ region (C4, Group I Row III $d^3N_{3b}$).
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
