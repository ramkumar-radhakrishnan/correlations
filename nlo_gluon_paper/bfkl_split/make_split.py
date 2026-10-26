r"""Finite_terms_split.tex:
  Sec. 1  Evolution_terms.tex: Sec. 1 + Sec. 2 = Sec. 3? (numbers from evol_check.json)
  Secs. 2-4  the BFKL (log(vee/k+)) terms of every row of Finite_terms_1.tex, row by row
             (the corrected list of ../bfkl_check, checked there against d/dlog(vee) of every row)
  Sec. 5  their sum
  Secs. 6-8  the remaining terms of every row, in the format of Finite_terms_1.tex:
             every log containing vee replaced by its remainder (split_logs.py, exact identities),
             Group III Rows III, V, VI and Group II Row IV edited by hand (p+ and xi integrals, checked in split_logs.py).
Row = BFKL terms + remaining terms, exactly, for any vee and Lambda."""
import json
import numpy as np
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'bfkl_check')))
import make_corrected as MC  # noqa: E402
import split_logs as SL  # noqa: E402

SRC = os.path.join(HERE, 'input', 'Finite_terms_1.tex')
OUT = os.path.join(HERE, 'Finite_terms_split.tex')
LV = SL.LVEE
T2 = MC.T2


# ---------------------------------------------------------------- helpers
def match_brace(t, i):
    """t[i] == '{' -> index after the matching '}'"""
    d = 0
    for k in range(i, len(t)):
        d += {'{': 1, '}': -1}.get(t[k], 0)
        if d == 0:
            return k + 1
    raise ValueError


def accept(t):
    """accept the changes of the red/blue copy: drop red, keep blue, remove the colour commands"""
    out, i = [], 0
    while i < len(t):
        if t.startswith(r'{\color{red}', i):
            i = match_brace(t, i)
            continue
        if t.startswith(r'{\color{blue}', i):
            j = match_brace(t, i)
            out.append(accept(t[i + len(r'{\color{blue}'):j - 1]))
            i = j
            continue
        if t.startswith(r'\old{', i):
            i = match_brace(t, i + 4)
            continue
        if t.startswith(r'\new{', i):
            j = match_brace(t, i + 4)
            out.append(accept(t[i + 5:j - 1]))
            i = j
            continue
        if t.startswith(r'\nc ', i):
            i += 4
            continue
        out.append(t[i])
        i += 1
    return ''.join(out)


# ---------------------------------------------------------------- Sec. 1
def evol_section():
    r = json.load(open(os.path.join(HERE, 'evol_check.json')))
    e = lambda x: r'%.0e' % x
    def tex(x):
        m, p = ('%.0e' % x).split('e')
        return r'$%s\times10^{%d}$' % (m, int(p))
    c = r['counts']
    return r"""
\section{\texttt{Evolution\_terms.tex}: does Sec.~1 + Sec.~2 give Sec.~3?}
\textbf{Yes, exactly, with no left-over.} Every two-$\rho$ term of the file was parsed and evaluated at random transverse points with random charges, with your convention $U^{bc}=U^{cb}$:
\begin{center}\small
\begin{tabular}{@{}lc@{}}\hline
comparison & $|{\rm difference}|/|{\rm terms}|$ (5 configurations)\\\hline
Sec.~1 + Sec.~2 vs.\ Sec.~3, $d^3N_1+d^3N_2$ & @@A@@\\
Sec.~1 + Sec.~2 vs.\ Sec.~3, collected form & @@B@@\\
Sec.~3 vs.\ JIMWLK$\otimes$LO & @@C@@\\
new list of ``From $4\rho$'s'' ($d^3N_2$) vs.\ the two-$\rho$ lists of \texttt{4rho.tex} & @@D@@\\
Sec.~1 + Sec.~2 \emph{without} the new list vs.\ Sec.~3 & @@E@@\\\hline
\end{tabular}
\end{center}
\begin{itemize}
\item Sec.~1 has @@N1@@ terms (@@E1@@ equations), Sec.~2 has @@N2@@ terms (@@E2@@ equations). Their sum is Sec.~3 term by term, and Sec.~3 is JIMWLK$\otimes$LO.
\item The new list $d^3N_2$ of ``From $4\rho$'s'' was what was missing before: without it the difference is @@E@@.
\item Two small points: in the colour-singlet form of Sec.~3, $\int_{x,,z,w,w'}$ should read $\int_{x,y,z,w,w'}$. For a general (non-symmetric) $U$ only the collected form of Sec.~3 is exact; with $U^{bc}=U^{cb}$ all forms agree.
\item Nothing in Secs.~1, 2 gives the BFKL bracket of Sec.~4: all $\log\frac{\vee}{\Lambda}$ terms add up to JIMWLK$\otimes$LO alone. The BFKL bracket comes from \texttt{Finite\_terms\_1.tex}, with $\log\frac{\vee}{k^+}$ (Sec.~\ref{sec:sum}).
\end{itemize}
""".replace('@@A@@', tex(r['Sec1 + Sec2 vs Sec3 (dN_1 + dN_2)']['sym'])) \
        .replace('@@B@@', tex(r['Sec1 + Sec2 vs Sec3 (collected form)']['sym'])) \
        .replace('@@C@@', tex(r['Sec3 collected form vs JIMWLKxLO']['sym'])) \
        .replace('@@D@@', tex(r['new 4rho list (dN_2) vs 4rho.tex two-rho lists (M1)']['sym'])) \
        .replace('@@E@@', r'%.1f\%%' % (100 * r['Sec1 + Sec2 without the new list vs JIMWLKxLO']['sym'])) \
        .replace('@@N1@@', str(c['sec1_terms'])).replace('@@E1@@', str(c['sec1_eqs'])) \
        .replace('@@N2@@', str(c['sec2_terms'])).replace('@@E2@@', str(c['sec2_eqs']))


# ---------------------------------------------------------------- Secs. 2-4: BFKL terms, row by row
B_NOTES = dict(
    NOTE_II=r"""
\textbf{Row II.} Part 1: from $\log\big[\frac{\vee}{k^+}\big]+\log\big[\frac{\vee}{k^+}-1\big]=2\log\big[\frac{\vee}{k^+}\big]+\log\big[1-\frac{k^+}{\vee}\big]$. Part 2a: all of it (nothing in it depends on $z$, so $\int_z$ is dropped). Part 2b: none.
""",
    NOTE_III=r"""
\textbf{Row III.} $d^3N_1$: from $\log\Big|\frac{\vee(x-z)^{2}-k^{+}(x-w)^{2}}{(\Lambda+k^{+})(x-z)^{2}-k^{+}(x-w)^{2}}\Big|$, $\log\big[\frac{\vee}{k^+-\Lambda}\big]$ and $\log\Big|\frac{k^{+}(x-w)^{2}-\vee(x-z)^{2}}{k^{+}(x-w)^{2}-\Lambda(x-z)^{2}}\Big|$. The first two give $+\log\frac{\vee}{k^+}$ and $-\log\frac{\vee}{k^+}$ with the same kernel and cancel; the third gives the term above. $d^3N_2$: the first term from $\log\big[\frac{\vee}{\Lambda+k^+}\big]$, the six terms in the bracket from $\log\Big|\frac{(x-w)^{2}k^{+}-(x-z)^{2}\vee}{(x-w)^{2}k^{+}-(x-z)^{2}(k^{+}+\Lambda)}\Big|$. $d^3N_{3a}$: from $\log\big[\frac{k^+-\Lambda}{\vee}\big]=-\log\frac{\vee}{k^+}+\log(1-\frac{\Lambda}{k^+})$. $d^3N_{3b}$: none.
""",
    NOTE_IV=r"""
\textbf{Row IV.} The first term is from $\log\big[\frac{\vee}{\Lambda+k^+}\big]$, the second from
\begin{equation*}
\log\Big[\frac{(\vee-k^{+})(x-z)^{2}+k^{+}(x-w)^{2}}{\Lambda(x-z)^{2}+k^{+}(x-w)^{2}}\Big]=\log\frac{\vee}{k^{+}}+\log\Big[\frac{\big(1-\frac{k^+}{\vee}\big)(x-z)^2+\frac{k^+}{\vee}(x-w)^2}{\frac{\Lambda}{k^+}(x-z)^2+(x-w)^2}\Big].
\end{equation*}
Together (the $\delta^{ij}$ part cancels the $\frac{1}{(z-w)^2}$ of the first term) the kernel is
\begin{equation*}
\frac{ig^4}{8\pi^5}\frac{1}{k^+}\Bigg[\frac{(x'-w')\cdot(z-w)\,(y-z)\cdot(x-z)}{(x'-w')^2(z-w)^2(y-z)^2(x-z)^2}-\frac{(x'-w')\cdot(x-w)\,(y-z)\cdot(x-z)}{(x'-w')^2(x-w)^2(y-z)^2(x-z)^2}-\frac{(x'-w')\cdot(x-w)\,(y-z)\cdot(z-w)}{(x'-w')^2(x-w)^2(y-z)^2(z-w)^2}\Bigg],
\end{equation*}
which is integrable at $z\to w$; the first term alone is not, and with it alone the three-$\rho$ terms could not cancel. The same holds for Group II Row I part 1 and Group III Rows I, II.
""",
    NOTE_V=r"""
\textbf{Row V.} From $\log\big[\frac{k^+-\Lambda}{\vee}\big]=-\log\frac{\vee}{k^+}+\log(1-\frac{\Lambda}{k^+})$ and from the $-2\log\frac{\vee}{k^+}$ of 2b UV. The $p^+$, $q^+$ integrals (from $0$ to $k^+$) have none.
""",
    NOTE_G2I1=r"""
\textbf{Row I, part 1.} As Group I Row IV, in both blocks.
""",
    G2I2_INTRO=r"""
\textbf{Row I, part 2.} The logarithms $\log\Big[\frac{(\vee-k^{+})(x-z)^{2}+k^{+}(x-w)^{2}}{\Lambda(x-z)^{2}+k^{+}(x-w)^{2}}\Big]$ and $\log\Big[\frac{(\vee-k^{+})(x'-z)^{2}+k^{+}(x'-w')^{2}}{\Lambda(x'-z)^{2}+k^{+}(x'-w')^{2}}\Big]$ each contain $\log\frac{\vee}{k^+}$. In 2b and in 2c the two terms of the row combine into the factors $\frac{1}{(x-w)^2(x'-w')^2}$ and $\frac{1}{(x-z)^2(x'-z)^2}$. 2a: none ($\vee$ cancels between numerator and denominator). 2d, 2e, 2f: their $\log\big[\frac{\vee}{k^++\Lambda}\big]$ cancels exactly against the $\log\frac{\vee}{k^+}$ of the two logarithms in the same rows: none.
""",
    G2I2_NOTE=r"""
The first term of $d^3N_{2c}$ has the kernel $\mathcal K(z-w)\cdot\mathcal K(z-w')\,\frac{(x-z)\cdot(x'-z)}{(x-z)^2(x'-z)^2}$: the measured gluon is emitted by the soft gluon $z$ in the amplitude and in the conjugate amplitude (second line of your Eq.~(LOBFKL)).
""",
    NOTE_G3C=r"""
\textbf{Group III Rows III, V, VI.} From the $\left[\frac{1}{p^+}\right]_+\delta_{jm}\delta_{ik'}$ term:
\begin{equation*}
\int_0^{\vee-k^+}\!dp^+\left[\frac{1}{p^+}\right]_+e^{-ik\cdot(y'-z)\frac{p^+}{k^+}}=-\Big[\log\frac{\vee}{k^+}+\log\Big(1-\frac{k^+}{\vee}\Big)+\gamma_E+\log\big(ik\cdot(y'-z)\big)+E_1\Big(ik\cdot(y'-z)\frac{\vee-k^+}{k^+}\Big)\Big]
\end{equation*}
(in part 2 in the same way, with $\frac{1}{p^+(x-z)^2+k^+(x-w)^2}\to\frac{1}{k^+(x-w)^2}$ at $p^+=0$). The other terms of these rows have no $\log\frac{\vee}{k^+}$.

\paragraph{Why the $\left[\frac{1}{p^+}\right]_+$ term contains $\log\frac{\vee}{k^+}$.} $\left[\frac{1}{p^+}\right]_+f=\frac{f(p^+)-f(0)}{p^+}$ is finite at $p^+\to0$: this term has no $\log\Lambda$. But the $p^+$ integral runs up to $\vee-k^+$, and the subtraction $-\frac{f(0)}{p^+}$ falls only like $\frac{1}{p^+}$. Splitting at $p^+=k^+$, with $f(p^+)=e^{-ik\cdot(y'-z)p^+/k^+}$:
\begin{equation*}
\int_0^{\vee-k^+}\!dp^+\,\frac{f(p^+)-f(0)}{p^+}=\underbrace{\int_0^{k^+}\!dp^+\,\frac{f(p^+)-f(0)}{p^+}+\int_{k^+}^{\vee-k^+}\!dp^+\,\frac{f(p^+)}{p^+}}_{\text{finite for }\vee\to\infty\ (f\text{ oscillates})}\;-\;f(0)\log\Big(\frac{\vee}{k^+}-1\Big).
\end{equation*}
The $p^+$ integral done numerically (no closed form used; \texttt{pplus\_check.py}), $k\cdot(y'-z)=0.7$:
\begin{center}\small
\begin{tabular}{@{}ccc@{}}\hline
$\frac{\vee}{k^+}$ & $I=\int_0^{\vee-k^+}dp^+\left[\frac{1}{p^+}\right]_+f$ & $I+\log\frac{\vee}{k^+}$\\\hline
@@PPTAB@@
$\infty$ & & $-\gamma_E-\log\big(ik\cdot(y'-z)\big)=@@PPLIM@@$\\\hline
\end{tabular}
\end{center}
$I$ grows like $-\log\frac{\vee}{k^+}$; $I+\log\frac{\vee}{k^+}$ is finite. The same holds in part 2 (\texttt{pplus\_check.json}). This $\log\frac{\vee}{k^+}$ belongs where it is put: before the subtraction the same term was $f(0)\int_\Lambda^{\vee-k^+}\frac{dp^+}{p^+}+\int_0^{\vee-k^+}dp^+\left[\frac{1}{p^+}\right]_+f$. Its $f(0)\log\frac{\vee}{\Lambda}$ is in \texttt{Evolution\_terms.tex}, its $f(0)\log(1-\frac{k^+}{\vee})$ is the $\log(1-\frac{k^+}{\vee})$ line of the row, and all three together give
\begin{equation*}
f(0)\Big[\log\frac{\vee}{\Lambda}+\log\Big(1-\frac{k^+}{\vee}\Big)-\log\Big(\frac{\vee}{k^+}-1\Big)\Big]+\text{finite}=f(0)\log\frac{k^+}{\Lambda}+\text{finite}.
\end{equation*}
The phase $e^{-ik\cdot(y'-z)p^+/k^+}$ cuts the $p^+$ integral off at $p^+\sim k^+/|k\cdot(y'-z)|$, so this term really has $\log\frac{k^+}{\Lambda}$, not $\log\frac{\vee}{\Lambda}$; the $\log\frac{\vee}{k^+}$ of the $\left[\frac{1}{p^+}\right]_+$ term removes the difference. By contrast, Group I Row V has $\left[\frac{k^+}{p^+}\right]_+$ on $0<p^+<k^+$: no $\vee$, so nothing is taken out there.
""",
)


def bfkl_sections():
    pp = json.load(open(os.path.join(HERE, 'pplus_check.json')))['part1 kappa=0.7']
    rows = []
    for r in pp['rows']:
        e = int(round(np.log10(r['vee_over_k'])))
        rows.append(r"$10^{%d}$ & $%.4f%+.4f\,i$ & $%.5f%+.5f\,i$\\" % (e, r['I'][0], r['I'][1], r['I_plus_log'][0], r['I_plus_log'][1]))
    lim = pp['limit_minus_gamma_minus_log_ikappa']
    notes = dict(B_NOTES)
    notes['NOTE_G3C'] = notes['NOTE_G3C'].replace('@@PPTAB@@', '\n'.join(rows)).replace('@@PPLIM@@', r'%.5f%+.5f\,i' % (lim[0], lim[1]))
    for k, v in notes.items():
        setattr(MC, k, v)
    MC.RESULT, MC.GROUP4, MC.TOTAL = '', '', ''
    tmp = os.path.join(HERE, '_bfkl_tmp.tex')
    MC.OUT = tmp
    MC.assemble()
    t = open(tmp).read()
    os.remove(tmp)
    t = accept(t)
    t = t[t.index(r'\section{Group I (with c.c.)}'):t.index(r'\end{document}')]
    t = t.replace('\nRemove:\n', '\n').replace(r'\allowdisplaybreaks', '')
    t = t.replace(r"\textbf{Rows I, II.} The blue terms are the $\log\frac{\vee}{k^+}$ inside",
                  r"\textbf{Rows I, II.} The first term is from $\log\big[\frac{\vee}{k^++\Lambda}\big]$, the second from the $\log\frac{\vee}{k^+}$ inside")
    t = t.replace(r"\textbf{Row IV.} Both terms come from", r"\textbf{Row IV.} Both terms come from")
    assert 'blue' not in t and 'your term' not in t.lower(), [m.start() for m in re.finditer('blue|your term', t)]
    for g in ('Group I (with c.c.)', 'Group II', r'Group III\, (with c.c.)'):
        t = t.replace(r'\section{%s}' % g, r'\section{BFKL terms: %s}' % g)
    return t


SUM = r"""
\section{Sum of the BFKL terms}\label{sec:sum}
With $\mathcal N\equiv\frac{1}{(2\pi)^{3}}\frac{g^{4}}{16\pi^{5}}\frac{1}{k^{+}}$, $\mathcal J$ = the bracket of Sec.~3 of \texttt{Evolution\_terms.tex} (JIMWLK$\otimes$LO) and $\mathcal B$ = the bracket of its Sec.~4 (BFKL$\otimes$LO, your Eq.~(LOBFKL)):
\begin{enumerate}
\item The three-$\rho$ terms of Secs.~2--4 cancel after symmetrization (checked: \texttt{LogVeeKplus\_BFKL.pdf}, all 39 structures).
\item The two-$\rho$ terms of Secs.~2--4 add up to
\begin{equation}
\sum_{\rm Secs.~2-4}\frac{d^3N}{d^3k}=\mathcal N\,\log\Big[\frac{\vee}{k^{+}}\Big]\int_{x,x',z,w,w'}e^{-ik(w'-w)}\Big[\mathcal B-\mathcal J\Big].
\end{equation}
\item The reason for $-\mathcal J$: every row has the soft gluon integrated over $\Lambda<p^+<\vee$. With $T$ the integrand for $p^+\ll k^+$ and $P$ the integrand for $p^+\gg k^+$,
\begin{equation*}
\int_\Lambda^\vee\frac{dp^+}{p^+}\,G\Big(\frac{p^+}{k^+}\Big)=T\log\frac{k^+}{\Lambda}+P\log\frac{\vee}{k^+}+\dots=T\log\frac{\vee}{\Lambda}+(P-T)\log\frac{\vee}{k^+}+\dots
\end{equation*}
$T\log\frac{\vee}{\Lambda}$ is in \texttt{Evolution\_terms.tex}, and $\sum T=\mathcal J$ (Sec.~1 + Sec.~2 = Sec.~3, Sec.~1 above); $\sum P=\mathcal B$. So the $\log\frac{\vee}{k^+}$ terms are $\sum(P-T)=\mathcal B-\mathcal J$.
\item Both files together:
\begin{equation}
\boxed{\mathcal N\int_{x,x',z,w,w'}e^{-ik(w'-w)}\Big[\log\Big[\frac{\vee}{\Lambda}\Big]\mathcal J+\log\Big[\frac{\vee}{k^{+}}\Big]\big(\mathcal B-\mathcal J\big)\Big]=\mathcal N\int_{x,x',z,w,w'}e^{-ik(w'-w)}\Big[\log\Big[\frac{k^{+}}{\Lambda}\Big]\mathcal J+\log\Big[\frac{\vee}{k^{+}}\Big]\mathcal B\Big].}
\end{equation}
JIMWLK comes with $\log\frac{k^+}{\Lambda}$ and BFKL with $\log\frac{\vee}{k^+}$. In Sec.~4 of \texttt{Evolution\_terms.tex} the bracket is right, but its logarithm is $\log\frac{\vee}{k^+}$; in Sec.~3 the logarithm becomes $\log\frac{k^+}{\Lambda}$ once the $-\mathcal J\log\frac{\vee}{k^+}$ above is added to it. At mid rapidity, $k^+=\sqrt{\vee\Lambda}$, both logarithms are $\frac12\log\frac{\vee}{\Lambda}$.
\end{enumerate}
"""


# ---------------------------------------------------------------- Secs. 6-8: remaining terms
REM_INTRO = r"""
\section{Remaining terms: how they are obtained}
Every row of \texttt{Finite\_terms\_1.tex} is, exactly (any $\vee$, any $\Lambda$),
\begin{equation*}
\text{row}=\big(\text{its BFKL terms, Secs.~2--4}\big)+\big(\text{its remaining terms, Secs.~7--9}\big).
\end{equation*}
The remaining terms below are the rows of \texttt{Finite\_terms\_1.tex} in the same form, with these changes (blue: changed or new; red: removed):
\begin{itemize}
\item every logarithm that contains $\vee$ is replaced by what is left of it after $\log\frac{\vee}{k^+}$ is taken out ($A=(x-z)^2$, $B=(x-w)^2$, and the same with $x'$, $w'$):
\begin{align*}
\log\Big[\frac{(\vee-k^{+})A+k^{+}B}{\Lambda A+k^{+}B}\Big]&=\log\frac{\vee}{k^{+}}+\log\Big[\frac{\big(1-\frac{k^{+}}{\vee}\big)A+\frac{k^{+}}{\vee}B}{\frac{\Lambda}{k^{+}}A+B}\Big],\\
\log\Big[\frac{\vee}{k^{+}+\Lambda}\Big]&=\log\frac{\vee}{k^{+}}-\log\Big[1+\frac{\Lambda}{k^{+}}\Big],\qquad
\log\Big[\frac{\vee}{k^{+}-\Lambda}\Big]=\log\frac{\vee}{k^{+}}-\log\Big[1-\frac{\Lambda}{k^{+}}\Big],\\
\log\Big[\frac{k^{+}-\Lambda}{\vee}\Big]&=-\log\frac{\vee}{k^{+}}+\log\Big[1-\frac{\Lambda}{k^{+}}\Big],\qquad
\log\Big[\frac{\vee}{k^{+}}-1\Big]=\log\frac{\vee}{k^{+}}+\log\Big[1-\frac{k^{+}}{\vee}\Big],\\
\log\Big|\frac{\vee A-k^{+}B}{(\Lambda+k^{+})A-k^{+}B}\Big|&=\log\frac{\vee}{k^{+}}+\log\Big|\frac{A-\frac{k^{+}}{\vee}B}{\big(1+\frac{\Lambda}{k^{+}}\big)A-B}\Big|,\\
\log\Big|\frac{k^{+}B-\vee A}{k^{+}B-\Lambda A}\Big|&=\log\frac{\vee}{k^{+}}+\log\Big|\frac{\frac{k^{+}}{\vee}B-A}{B-\frac{\Lambda}{k^{+}}A}\Big|,\\
\log\Big|\frac{Bk^{+}-A\vee}{Bk^{+}-A(k^{+}+\Lambda)}\Big|&=\log\frac{\vee}{k^{+}}+\log\Big|\frac{B\frac{k^{+}}{\vee}-A}{B-A\big(1+\frac{\Lambda}{k^{+}}\big)}\Big|;
\end{align*}
\item Group II Row IV and Group III Rows III, V, VI: the $\xi$ and $p^+$ integrals that contain $\log\frac{\vee}{k^+}$ are done (notes under these rows).
\end{itemize}
All identities are checked numerically to $10^{-14}$ (\texttt{split\_logs.json}). In the remaining terms $\vee$ enters only through $\frac{k^+}{\vee}$, $\log(1-\frac{k^+}{\vee})$, $\int_{k^+/\vee}d\xi$, $\frac{1}{\vee}$ and the bounded oscillating functions $E_1\big(ik\cdot(y'-z)\frac{\vee-k^+}{k^+}\big)$, $\int_0^{\vee-k^+}dp^+$ of Group III: they have a finite limit for $\vee\to\infty$. The one exception is written out: $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ in Group I Row II, which vanishes like $\frac{k^+}{\vee}\log\frac{\vee}{k^+}$. $\Lambda\to0$ can be taken in all of them.
"""

G3_NEW = [r"\nc +\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,@@F@@}{@@C@@}\, \frac{1}{k^{+}}\,\int_{@@M@@}\, e^{-ik\cdot (y'-w)}\,"
          r"\Big[\gamma_{E}+\log\big(ik\cdot(y'-z)\big)+E_{1}\Big(ik\cdot(y'-z)\,\frac{\vee-k^{+}}{k^{+}}\Big)\Big]",
          r"\nc \times %s" % T2]
G3_NOTE = r"""
\new{The $\left[\frac{1}{p^+}\right]_+$ term (red) is integrated: $-\int_0^{\vee-k^+}\frac{dp^+}{2\pi}\left[\frac{1}{p^+}\right]_+e^{-ik\cdot(y'-z)\frac{p^+}{k^+}}=\frac{1}{2\pi}\Big[\log\frac{\vee}{k^+}+\log\big(1-\frac{k^+}{\vee}\big)+\gamma_E+\log\big(ik\cdot(y'-z)\big)+E_1\big(ik\cdot(y'-z)\frac{\vee-k^+}{k^+}\big)\Big]$. Its $\log\frac{\vee}{k^+}$ is in the BFKL terms; its $\log(1-\frac{k^+}{\vee})$ cancels the red last line; the rest is the blue term.}
"""
JB = (r"e^{ik\cdot(y'-z)\frac{(x-w)^{2}}{(x-z)^{2}}}\Big[E_{1}\Big(ik\cdot(y'-z)\,\frac{(x-w)^{2}}{(x-z)^{2}}\Big)"
      r"-E_{1}\Big(ik\cdot(y'-z)\Big[\frac{(x-w)^{2}}{(x-z)^{2}}+\frac{\vee-k^{+}}{k^{+}}\Big]\Big)\Big]")
G3_III2_NEW = [
    r"\nc +\, \frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\int_{w,x',y',z,x}\, e^{-ik\cdot (y'-w)}\,\bigg[\frac{(x-w)\cdot (x'-y')\,(z-w)\cdot (y'-z)}{(z-w)^{2}\,(x-w)^{2}\,(y'-z)^{2}\,(x'-y')^{2}}+\frac{(x-z)\cdot (y'-z)\,(x-w)\cdot (x'-y')}{2\,(x-z)^{2}\,(x-w)^{2}\,(y'-z)^{2}\,(x'-y')^{2}}\bigg]",
    r"\nc \times\Big[\gamma_{E}+\log\big(ik\cdot(y'-z)\big)+E_{1}\Big(ik\cdot(y'-z)\,\frac{\vee-k^{+}}{k^{+}}\Big)",
    r"\nc \qquad+" + JB + r"\Big]",
]
G3_III2_NOTE = r"""
\new{The $\left[\frac{1}{p^+}\right]_+$ term (red) is integrated, with $f(p^+)=\frac{e^{-ik\cdot(y'-z)p^+/k^+}}{p^+(x-z)^2+k^+(x-w)^2}$: $\int_0^{\vee-k^+}dp^+\frac{f(p^+)-f(0)}{p^+}=-\frac{1}{k^+(x-w)^2}\Big[\log\frac{\vee}{k^+}+\log\big(1-\frac{k^+}{\vee}\big)+\gamma_E+\log\big(ik\cdot(y'-z)\big)+E_1\big(ik\cdot(y'-z)\frac{\vee-k^+}{k^+}\big)+\mathcal J\Big]$, $\mathcal J$ = the last bracket of the blue term. Its $\log\frac{\vee}{k^+}$ is in the BFKL terms; its $\log(1-\frac{k^+}{\vee})$ cancels the red last line.}
"""
G2IV_NEW_INT = (r"\new{+\int_{\frac{k^{+}}{\vee}}^{1}d\xi\,\frac{\frac{1}{\xi}\,e^{-ik\cdot (y'-y)\,\frac{1}{\xi}}-e^{-ik\cdot (y'-y)}}{1-\xi}}")
G2IV_H4 = [
    r"\nc -\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]}\,\int_{x,x',z,w,w'}\,e^{-ik\cdot (w'-w)}\,\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}\,(w'-z)^{2}}\,\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}\,(x'-z)^{2}}",
    r"\nc \left[U^{e'c}(z) - U^{e'c}(x')\right] \times\, \left[U^{ec}(z) - U^{ce}(x)\right]\, \rho^{e'}(x') \rho^{e}(x)",
]
G2IV_NOTE = r"""
\new{\textbf{Row IV.} Your subtraction $\frac{1}{\xi^2}\frac{\xi}{(1-\xi)_+}e^{-ik\cdot(y'-y)\frac1\xi}=\frac{e^{-ik\cdot(y'-y)\frac1\xi}-e^{-ik\cdot(y'-y)}}{\xi(1-\xi)}$ (red) is written as $\frac{\frac1\xi e^{-ik\cdot(y'-y)\frac1\xi}-e^{-ik\cdot(y'-y)}}{1-\xi}-\frac{e^{-ik\cdot(y'-y)}}{\xi}$; the first part is the blue $\xi$ integral, the second gives $-e^{-ik\cdot(y'-y)}\log\frac{\vee}{k^+}$, which is in the BFKL terms. The first two lines are now the standard $\frac{\alpha_s}{2\pi}\int_{k^+/\vee}^1\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi$. At its lower end $\xi\to\frac{k^+}{\vee}$ the $\xi$ integral still produces $\log\frac{\vee}{k^+}$ times the integral of the last two (blue) lines, from the $\frac{1-\xi}{\xi}$ term at small $\xi$; this is the second BFKL term of this row (Sec.~3). The last two lines subtract it, so the sum is finite for $\vee\to\infty$.}
"""


def remaining_sections():
    s = open(SRC).read()
    body = s[s.index(r'\section{Group I (with c.c.)}'):s.index(r'\end{document}')]
    parts = re.split(r'(\\begin\{align\}.*?\\end\{align\})', body, flags=re.S)
    idx = [i for i, p in enumerate(parts) if p.startswith(r'\begin{align}')]
    assert len(idx) == 28
    al = {n: parts[i] for n, i in enumerate(idx)}
    new = {}

    def sub_logs(a):
        for orig, (c, R) in SL.SPLIT.items():
            a = a.replace(orig, r'\new{%s}' % R)
        for orig, (bf, R) in SL.CONTEXT.items():
            a = a.replace(orig, r'\new{%s}' % R)
        return a

    for n in al:
        new[n] = sub_logs(al[n])
    # Group I Row II 2a: all BFKL
    new[1] = r"\new{$d^3N_{2a}$: no remaining term (it is proportional to $\log\frac{\vee}{k^+}$; all of it is in Sec.~2).}" + "\n"

    def lines(a):
        body = a[len(r'\begin{align}'):a.rindex(r'\end{align}')]
        return body.split(r'\notag \\')

    def join(L):
        return r'\begin{align}' + r'\notag \\'.join(L).rstrip() + '\n' + r'\end{align}'

    def red(line):
        k = line.index('&')
        return line[:k + 1] + r' \rc ' + line[k + 1:]

    PP = r"-\left[\frac{1}{p^{+}}\right]_{+}\,\delta_{jm}\,\delta_{ik'}"
    for n, (coef, f, meas) in {24: (r"16\pi^{5}", "f^{c'd'a}", "w,x',y',z,y,x"), 26: (r"8\pi^{5}", "f^{c'ca}", "w,y',z,y,x,x'"),
                               27: (r"8\pi^{5}", "f^{c'd'a}", "w,y',z,y,x,x'")}.items():
        L = lines(al[n])
        assert L[1].count(PP) == 1 and r'\log{\left[1-\frac{k^{+}}{\vee}\right]}' in L[3]
        L[1] = L[1].replace(PP, r'{\color{red}%s}' % PP)
        col = L[1][L[1].index(r'\right]\,', L[1].index(PP)) + len(r'\right]\,'):]
        L[3], L[4] = red(L[3]), red(L[4])
        L += [' & ' + G3_NEW[0].replace('@@F@@', f).replace('@@C@@', coef).replace('@@M@@', meas), ' & ' + G3_NEW[1],
              ' & ' + NC_ + col.strip(), ' & ' + NC_ + L[2].split('&', 1)[1].strip()]
        new[n] = join(L) + G3_NOTE
    # Group III Row III part 2
    L = lines(al[25])
    PP2 = r"-\left[\frac{1}{p^{+}}\right]_{+}\,\bigg[\frac{(x-w)^{k'}\,(z-w)^{m}}{(z-w)^{2}}+\frac{(x-z)^{m}\,(x-w)^{k'}}{2\,(x-z)^{2}}\bigg]"
    assert L[3].count(PP2) == 1 and r'\log{\left[1-\frac{k^{+}}{\vee}\right]}' in L[6]
    L[3] = L[3].replace(PP2, r'{\color{red}%s}' % PP2)
    L[6], L[7], L[8] = red(L[6]), red(L[7]), red(L[8])
    col1 = L[5].split('&', 1)[1].strip()
    L += [' & ' + x for x in G3_III2_NEW] + [' & ' + NC_ + col1]
    new[25] = join(L) + G3_III2_NOTE
    # Group II Row IV
    L = lines(al[21])
    sub = r"+\frac{\xi}{(1-\xi)_{+}}\right]"
    assert L[1].count(sub) == 1
    L[1] = L[1].replace(sub, r"{\color{red}+\frac{\xi}{(1-\xi)_{+}}}\right]" + G2IV_NEW_INT)
    assert L[1].count(G2IV_NEW_INT + r"\Bigg]\times") == 1
    L[1] = L[1].replace(G2IV_NEW_INT + r"\Bigg]\times", G2IV_NEW_INT + r"\Bigg]\notag \\ & \times")
    L += [' & ' + G2IV_H4[0], ' & ' + G2IV_H4[1]]
    new[21] = join(L) + G2IV_NOTE

    for n in (16, 20):
        L = lines(new[n])
        for k, l in enumerate(L):
            for pat in (r"\times\,\new{\log", r"\times\, \new{\log"):
                j = l.find(pat)
                if j > 0 and l[:j].count(r'\left') == l[:j].count(r'\right') and len(l[:j].strip()) > 40:
                    L[k] = l[:j] + r"\notag \\ & " + l[j:]
                    break
        new[n] = join(L)
    out = []
    n = 0
    for i, p in enumerate(parts):
        if i in idx:
            out.append(new[n])
            n += 1
        else:
            out.append(p)
    t = ''.join(out)
    for g in ('Group I (with c.c.)', 'Group II', r'Group III\, (with c.c.)'):
        t = t.replace(r'\section{%s}' % g, r'\section{Remaining terms: %s}' % g)
    t = t.replace(r"\begin{align}" + "\n" + r"R^{i} =", r"\begin{align}" + "\n" + r"R^{i} =")
    # notes
    t = t.replace(r"\subsection{Row III}", r"\new{Part 1: the remaining $\log\big[1-\frac{k^+}{\vee}\big]\log\big[\frac{\vee}{k^+}\big]$ vanishes for $\vee\to\infty$; it is not a $\log\frac{\vee}{k^+}$ term. Part 2b: no change.}" + "\n" + r"\subsection{Row III}", 1)
    return t


NC_ = r"\nc "


def coverage(t):
    """every vee left in the remaining terms must be one of the allowed forms"""
    allowed = [r"\frac{k^{+}}{\vee}", r"\int_{\frac{k^{+}}{\vee}}", r"\frac{\vee-k^{+}}{k^{+}}", r"\int_{0}^{\vee-k^{+}}",
               r"\left[\frac{1}{\vee} - \frac{1}{k^{+} + \Lambda}\right]", r"\log{\left[\frac{\vee}{k^{+}}\right]}"]
    s = ''.join(re.findall(r'\\begin\{align\}.*?\\end\{align\}', t, flags=re.S))   # equations only, not the notes
    t = s
    for a in allowed:
        s = s.replace(a, '@')
    bad = [m.start() for m in re.finditer(r'\\vee', s)]
    left = [s[max(0, i - 80):i + 40] for i in bad]
    lv = t.count(r"\log{\left[\frac{\vee}{k^{+}}\right]}")
    return left, lv


def assemble():
    src = open(SRC).read()
    pre = src[:src.index(r'\begin{document}')]
    pre = pre.replace(r"\title{Finite terms}", r"\title{Finite terms: BFKL terms and remaining terms}")
    pre = pre.replace(r"\usepackage{comment}", r"\usepackage{comment}" + "\n" + r"\newcommand{\new}[1]{{\color{blue}#1}}" + "\n"
                      + r"\newcommand{\old}[1]{{\color{red}#1}}" + "\n" + r"\newcommand{\nc}{\color{blue}}" + "\n"
                      + r"\newcommand{\rc}{\color{red}}" + "\n" + r"\allowdisplaybreaks")
    rem = remaining_sections()
    left, lv = coverage(rem)
    assert not left, left
    assert lv == 2, lv          # log(1-k+/vee) log(vee/k+) of Group I Row II, and the H4 subtraction of Group II Row IV
    doc = (pre + r"\begin{document}" + "\n" + r"\maketitle" + "\n" + r"\sloppy" + "\n"
           + r"""\paragraph{Contents.} Sec.~1: \texttt{Evolution\_terms.tex}, Sec.~1 + Sec.~2 vs.\ Sec.~3. Secs.~2--4: the terms of \texttt{Finite\_terms\_1.tex} that enter BFKL (its $\log\frac{\vee}{k^+}$ terms), row by row. Sec.~5: their sum. Secs.~6--9: the remaining terms, row by row, in the form of \texttt{Finite\_terms\_1.tex}. For every row: row = BFKL terms + remaining terms, exactly.
""" + evol_section() + bfkl_sections() + SUM + REM_INTRO + rem + r"\end{document}" + "\n")
    open(OUT, 'w').write(doc)
    print('written', OUT)


if __name__ == '__main__':
    assemble()
