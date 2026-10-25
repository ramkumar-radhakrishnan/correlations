"""Corrected copy of the new Evolution_terms.tex: removed text in red, new text in blue.
  E1  Group I Row V: the fourth equation is the counterterm (label).
  E2  'From 4 rho's': the two-rho lists of 4rho.tex (M1) added.
  E3  singlet forms: <rho^f(x) rho^f(y)>/(N_c^2-1), integral over y, K(a,b,z) defined.
  E4  Sec. 3 = sum of Secs. 1, 2; BFKL with (1/2) log(vee/Lambda); the total with the finite remainder added.
The result is parsed back: Secs. 1+2 must equal JIMWLK x LO (U = U^T) and Sec. 3.
Output: Evolution_terms_corrected.tex"""
import os
import re
import env  # noqa: F401
from core import Config, pointwise_total
import target as TG
from parse_v4 import equations
from check_all import both

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'input', 'Evolution_terms.tex')
DST = os.path.join(HERE, 'Evolution_terms_corrected.tex')


def once(s, old, new, count=1):
    assert s.count(old) == count, (old[:80], s.count(old))
    return s.replace(old, new)


if __name__ == '__main__':
    s = open(SRC).read()
    # macros
    s = once(s, r'\numberwithin{equation}{section}',
             '\\numberwithin{equation}{section}\n\\newcommand{\\new}[1]{{\\color{blue}#1}}\n\\newcommand{\\old}[1]{{\\color{red}#1}}')

    # E1: second 'd^3N_{2b finite}' is the counterterm
    lab = r'\frac{d^{3}N_{2b\,\rm finite\, NLO}^{\rm}}{d^{3}k}'
    i = s.index(lab)
    j = s.index(lab, i + 1)
    s = s[:j] + r'\frac{d^{3}N_{2b\,\old{\rm finite}\new{\rm CT}\, NLO}^{\rm}}{d^{3}k}' + s[j + len(lab):]
    k = s.index(r'\end{align}', j) + len(r'\end{align}')
    s = s[:k] + ('\n\\new{This is the counterterm of Part 2b: it is $-d^{3}N_{2b\\,\\rm UV}$ term by term, so the two cancel. '
                 'The label ``finite\'\' of the first $d^3N_{2b}$ list is the right one; this one is not a finite part.}') + s[k:]

    # E2: the two-rho lists of 4rho.tex in 'From 4 rho's'
    body = open(os.path.join(HERE, 'fourlists_body.tex')).read()
    body = re.sub(r'\\color\{\w+\}', r'\\color{blue}', body).replace(r'\tfrac', r'\frac').replace(r'&\quad', '&')
    four = ('\n\\new{Missing (M1): the two-$\\rho$ terms of \\texttt{4rho.tex} (``Two $\\rho$ contribution\'\', $dN^{(1)}$, $dN^{(2)}$ '
            'with c.c., with the fixes F1, F2 of \\texttt{Final\\_evolution\\_check.pdf}). They belong to ``From $4\\,\\rho$\'s\'\'. '
            'Without them Secs.~1 and 2 do not add up to the JIMWLK contribution of Sec.~3.}\n'
            '\\begin{align}\n{\\color{blue}\\frac{d^{3}N_{\\mathrm{NLO}}}{d^{3}k}} & {\\color{blue}= \\frac{1}{(2\\pi)^{3}}\\,\\frac{g^{4}}{16\\pi^{5}}'
            '\\,\\frac{1}{k^{+}}\\,\\log\\Big[\\frac{\\vee}{\\Lambda}\\Big]\\int_{x,x\',z,w,w\'}e^{-ik(w\'-w)}}\\,{\\color{blue}\\Bigg[}\\notag\\\\\n'
            + body.strip() + '{\\color{blue}\\Bigg]}\n\\end{align}\n')
    m = s.index(r'\textbf{From three $\rho$')
    s = s[:m] + four + s[m:]

    # E3: singlet forms
    for head in ('Color singlet form', 'In color singlet form, we get'):
        a = s.index(head)
        b = s.index(r'\end{align}', a)
        seg = s[a:b]
        seg = seg.replace(head, head + ' \\new{(projectile in a colour singlet: $\\rho^b(x\')\\,\\rho^c(x)\\to\\delta^{bc}\\,'
                          '\\langle\\rho^f(x)\\rho^f(y)\\rangle/(N_c^2-1)$ with $x\'\\to y$; $K(a,b,z)\\equiv\\mathcal K(a-z)\\cdot\\mathcal K(b-z)$, '
                          '$K(a,a,z)=\\frac{1}{(a-z)^2}$; $U^\\dagger=U^T$, any $U$)}', 1)
        seg = once(seg, r"\int_{x,x',z,w,w'}", r"\new{\frac{1}{N_c^2-1}}\int_{x,\old{x'}\new{y},z,w,w'}")
        seg = re.sub(r"(e\^\{-ik\\?c?d?o?t?\(?w'-w\)?\}\\?,?\s*)\\Bigg\[", lambda mm: mm.group(1) + r"\new{\langle\rho^f(x)\rho^f(y)\rangle}\Bigg[", seg, count=1)
        assert r'\langle\rho^f(x)\rho^f(y)\rangle}\Bigg[' in seg, head
        s = s[:a] + seg + s[b:]

    # E4: Sec. 3 is the sum of Secs. 1, 2; BFKL with (1/2) log(vee/Lambda); the total with the finite remainder
    a3 = s.index(r'\section{JIMWLK contribution}')
    s = s[:a3] + s[a3:].replace(r'\section{JIMWLK contribution}', r'\section{JIMWLK contribution}' + '\n' +
        '\\new{Eqs.~(3.1)--(3.3) are the sum of all $\\log\\frac{\\vee}{\\Lambda}$ terms of Secs.~1 and 2 (with the $4\\rho$ lists added): '
        'their bracket, $\\mathcal J$, is JIMWLK$\\otimes$LO, for any $U$. In the final result (Sec.~5) half of it, '
        '$\\frac12\\log\\frac{\\vee}{\\Lambda}\\,\\mathcal J$, is the JIMWLK evolution of the target; the other half combines with the '
        '$\\log\\frac{\\vee}{k^+}$ terms of the rows into the BFKL evolution of the projectile, Sec.~4.}\n', 1)
    a = s.index(r'\section{BFKL contribution}')
    b = s.index(r'\end{document}')
    seg = s[a:b]
    n = seg.count(r'\log\Big[\frac{\vee}{\Lambda}\Big]')
    assert n == 3, n
    seg = seg.replace(r'\log\Big[\frac{\vee}{\Lambda}\Big]', r'\new{\frac{1}{2}}\log\Big[\frac{\vee}{\Lambda}\Big]')
    assert seg.count(r'\frac{d^{3}N_{\mathrm{NLO}}}{d^{3}k}') == 3
    seg = seg.replace(r'\frac{d^{3}N_{\mathrm{NLO}}}{d^{3}k}', r'\frac{d^{3}N_{\new{\mathrm{BFKL}}}}{d^{3}k}')
    seg = seg.replace(r'\section{BFKL contribution}', r'\section{BFKL contribution}' + '\n' +
        '\\new{The bracket below, $\\mathcal B$, is BFKL$\\otimes$LO (Eq.~(4.2) for any $U$). It is not among the $\\log\\frac{\\vee}{\\Lambda}$ '
        'terms of Secs.~1, 2: in all of them the measured gluon is emitted by a valence charge, $\\mathcal K(x-w)\\cdot\\mathcal K(x\'-w\')$, '
        'while in 23 of the 35 terms of Eq.~(4.2) it is emitted by the gluon at $z$. These terms come from the upper limits of the $p^+$ '
        'integrals, $p^+\\to\\vee$, $p^++k^+\\to\\vee$, with $\\log\\frac{\\vee}{k^+}$: the $\\log\\frac{\\vee}{k^+}$ terms of all rows sum to '
        '$\\log\\frac{\\vee}{k^+}(\\mathcal B-\\mathcal J)$. With the other half of Sec.~3 they give the BFKL contribution '
        '$\\frac12\\log\\frac{\\vee}{\\Lambda}\\,\\mathcal B$ below, plus a finite remainder (Sec.~5). '
        '(Left-hand sides renamed from \\old{$d^3N_{\\rm NLO}$} to \\new{$d^3N_{\\rm BFKL}$}.)}\n', 1)
    total = r'''
\section{\new{Total: evolution and finite remainder}}
\new{Your rows, with rule A, $\log\frac{k^+}{\Lambda}=\log\frac{\vee}{\Lambda}-\log\frac{\vee}{k^+}$, and with $\mathcal J$, $\mathcal B$ the brackets $\Big[\dots\Big]$ of Secs.~3 and 4 (without the logarithm), give}
\begin{align}
{\color{blue}\frac{d^{3}N_{\rm LL}}{d^{3}k}}&{\color{blue}=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{x,x',z,w,w'}e^{-ik(w'-w)}\Big[\log\Big[\frac{\vee}{\Lambda}\Big]\,\mathcal J+\log\Big[\frac{\vee}{k^{+}}\Big]\big(\mathcal B-\mathcal J\big)\Big]}\notag\\
&{\color{blue}=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{x,x',z,w,w'}e^{-ik(w'-w)}\Big[\frac{1}{2}\log\Big[\frac{\vee}{\Lambda}\Big]\,\mathcal J+\frac{1}{2}\log\Big[\frac{\vee}{\Lambda}\Big]\,\mathcal B+\log\Big[\frac{k^{+}}{\sqrt{\vee\Lambda}}\Big]\big(\mathcal J-\mathcal B\big)\Big],}
\end{align}
\new{using $\log\frac{\vee}{k^+}=\frac12\log\frac{\vee}{\Lambda}-\log\frac{k^+}{\sqrt{\vee\Lambda}}$:}
\begin{itemize}
\item \new{$\frac12\log\frac{\vee}{\Lambda}\,\mathcal J$: JIMWLK evolution of the target (Sec.~3 with $\frac12$);}
\item \new{$\frac12\log\frac{\vee}{\Lambda}\,\mathcal B$: BFKL evolution of the projectile (Sec.~4);}
\item \new{$\log\frac{k^+}{\sqrt{\vee\Lambda}}\,(\mathcal J-\mathcal B)$: finite. It is the distance in rapidity of the measured gluon from the centre of the window and vanishes at mid rapidity, $k^+=\sqrt{\vee\Lambda}$.}
\item \new{In colour-singlet form, $\mathcal J$ and $\mathcal B$ are the brackets of Eqs.~(3.4) and (4.3).}
\item \new{Counting $\log\frac{\vee}{k^+}$ itself as finite would leave BFKL inside the finite part, $\log\frac{\vee}{k^+}(\mathcal B-\mathcal J)$, which at mid rapidity is $\frac12\log\frac{\vee}{\Lambda}(\mathcal B-\mathcal J)$, not finite.}
\end{itemize}
'''
    s = s[:a] + seg + total + s[b:]
    open(DST, 'w').write(s)

    # parse back
    E = equations(DST)
    S12 = [t for e in E if e['sec'].startswith(('Genuine', 'Symmetrizing')) for t in e['terms']]
    Jc = [e for e in E if e['sec'].startswith('JIMWLK') and e['kind'] == 'terms'][2]['terms']
    res = dict(S12_vs_J=both(S12, TG.jimwlk()), S12_vs_sec3=both(S12, Jc))
    print('corrected:', res)
    import json
    json.dump(res, open(os.path.join(HERE, 'corrected.json'), 'w'), indent=1, default=float)
    print('equations:', len(E), 'terms:', sum(len(e['terms']) for e in E))
