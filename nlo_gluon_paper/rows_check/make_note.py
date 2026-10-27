r"""Writes Group_check.tex from check_Group_{I,II,III}.json and measures.json."""
import json
import os
import check_all as CA
import eqnum as EN

HERE = os.path.dirname(os.path.abspath(__file__))
J = {g: json.load(open(os.path.join(HERE, 'check_%s.json' % g))) for g in ('Group_I', 'Group_II', 'Group_III')}


def tex_num(x):
    if x == 0:
        return '0'
    m, e = ('%.1e' % x).split('e')
    return r'%s\times10^{%d}' % (m, int(e))


def eqs(g, ranges):
    return ', '.join(EN.eq(g, a, b) for a, b in ranges)


def table():
    rows = []
    for g, gname in (('Group_I', 'I'), ('Group_II', 'II'), ('Group_III', 'III')):
        for name, kind, e1, e2, opt in CA.STEPS[g]:
            res = J[g][name]
            if kind == 'plus':
                r6 = max(r['rel_L1e-06'] for r in res)
                r7 = max(r['rel_L1e-07'] for r in res)
                val = r'$%s$ / $%s$' % (tex_num(r6), tex_num(r7))
            elif kind == 'int':
                k = [x for x in res[0] if x.startswith('rel')][0]
                val = '$%s$' % tex_num(max(r[k] for r in res))
            else:
                val = '$%s$' % tex_num(max(r['rel'] for r in res))
            mark = r' \textbf{(M2)}' if 'after the z integral' in name else ''
            disp = (name.replace('Eq. (111)', 'Eq. (1.16)').replace('(q+ = k+ - p+)', '($q^+=k^+-p^+$)')
                    .replace('P_gg', r'$\mathcal P_{gg}$').replace('f f U U = N_c U', '$ffUU=N_cU$')
                    .replace('p+ ', '$p^+$ ').replace('xi', r'$\xi$').replace("w'", "$w'$").replace('->', r'$\to$'))
            kname = {'int': 'integrals', 'plus': r'plus, $\Lambda\to0$', 'pt': r'fixed $p^+$'}[kind]
            rows.append(r'%s & %s%s & %s $\to$ %s & %s & %s\\' % (gname, disp, mark,
                                                             eqs(g, e1), eqs(g, e2), kname, val))
    return '\n'.join(rows)


UV = J['Group_I']['Row V 2b: UV term'][0]
U2, U3 = UV['UV_coefficient_L0.01'], UV['UV_coefficient_L0.001']

DOC = r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.55in]{geometry}
\usepackage{amsmath,amssymb,mathtools,booktabs,longtable,array}
\usepackage[dvipsnames]{xcolor}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\newcommand{\new}[1]{{\color{blue}#1}}
\newcommand{\old}[1]{{\color{red}#1}}
\allowdisplaybreaks
\title{Check of your \texttt{Group\_I.tex}, \texttt{Group\_II.tex}, \texttt{Group\_III.tex}}
\date{}
\begin{document}
\maketitle
\noindent Equation numbers are those of your files as they compile. Corrected copies (\new{blue} = new, \old{red} = remove; your own colours are shown in teal, olive, violet, brown there):\\ \texttt{Group\_I\_corrected.pdf}, \texttt{Group\_II\_corrected.pdf}, \texttt{Group\_III\_corrected.pdf}.

\section{Answer}
\paragraph{What was checked.} Every step that turns one of your equations into the next: the $p^+$, $q^+$ and $\xi$ integrals, the $w$ and $w'$ integrals over the $\delta^{(2)}$, the reductions of the Wilson lines, the splits (Group~I Row~V into 2a, 2b, and 2b into UV and finite; Group~II Row~I part~2 into 2a--2f), and every plus prescription. Each step is checked numerically with your equations as written (read directly from your \LaTeX), pointwise in all transverse coordinates, for random points, random adjoint SU(3) Wilson lines ($U^{bc}=U^{cb}$) and random charges, at two configurations; remaining $p^+$, $q^+$, $\xi$ integrals are done numerically (principal value where the integrand has a pole). The amplitude-level first lines ($\mathbb A$, $\mathbb B$, $\mathfrak C$) and the momentum integrals of Group~I Row~II part~2 are unchanged since they were checked against \texttt{Calculation\_main.tex} (\texttt{Group\_I\_v2\_check.pdf}, \texttt{Groups\_II\_III\_check.pdf}); compared with the files of \texttt{LL\_terms.pdf}, three things changed (Group~I Eq.~(1.32), a stray line in Group~II, colours in Group~III).

\paragraph{Mistakes.}
\begin{enumerate}
\item[\textbf{M1}] \textbf{Group I, Row V, 2b UV, Eqs.~(1.26) and (1.32) --- this one changes the result.} The bracket must be
\begin{equation*}
\Big[\old{\log\frac{k^+}{\Lambda}-\log\frac{\vee}{k^+}}\;\new{2\log\frac{k^+}{\Lambda}}+2\log\Big(1-\frac{\Lambda}{k^+}\Big)-\frac{11}{6}\Big].
\end{equation*}
There is no $\log\frac{\vee}{k^+}$ in Row V ($0<p^+<k^+$), and the coefficient of $\log\frac{k^+}{\Lambda}$ is 2. Your exact Eq.~(1.25) is right; the factor 2 is lost in the step ``Setting $\Lambda\to0$'' (Sec.~\ref{sec:M1}). As written, the poles $\times\log\frac{k^+}{\Lambda}$ of Group~I Row~II, Group~I Row~V 2b UV and Group~II Row~IV do not cancel, and a pole $\times\log\frac{\vee}{k^+}$ is left over that nothing cancels.
\item[\textbf{M2}] \textbf{Group II, Row IV, Eq.~(1.16):} in the second term the phase must be $e^{-ik\cdot(y'-y)\frac{p^++k^+}{k^+}}$, not $e^{-ik\cdot(y'-y)\frac{p^+}{(p^++k^+)}}$. It is a typo: Eq.~(1.14) before and Eq.~(1.17) after have the right phase. Numerically Eq.~(1.16)$\to$(1.17) fails by $@@M2A@@$ as written and agrees to $@@M2B@@$ with the phase corrected.
\item[\textbf{M3}] \textbf{Group III, Row III, first (unnumbered) line:} the prefactor is $\new{-}\frac{4}{(2\pi)^3}$. When the overall sign of Rows III, V, VI was corrected (they enter with $-\mathfrak C'^\dagger$, M1 of \texttt{Groups\_II\_III\_check.pdf}), the first lines of Rows V and VI were changed to $-\frac{2}{(2\pi)^3}$, that of Row III was not. The explicit forms below it, Eqs.~(1.5)--(1.10), are right.
\item[\textbf{M4}] \textbf{Group III, Row III, second part, Eqs.~(1.8)--(1.10) (C5):} the measure is $\int_{w,w',x',y',z,x}$ in Eq.~(1.8) and $\int_{w,x',y',z,x}$ after the $w'$ integral (Eqs.~(1.9), (1.10), both lines). $y$ does not appear in this part, $x'$ does, and $w'$ is already integrated.
\item[\textbf{M5}] \textbf{Group I, Row II, Eq.~(1.6):} $\int_{w,w'}\to\int_{w,w',x'}$ ($x'$ appears and is not integrated).
\item[\textbf{M6}] \textbf{Typos that do not change anything:} Group~II Eq.~(1.13): $\frac{(k^++p^+)}{k\{+\}}\to\frac{(k^++p^+)}{k^+}$ in the phase. Group~I Eq.~(1.9): the opening \texttt{\textbackslash Bigg\textbackslash\{} is never closed; Eqs.~(1.19), (1.20), (1.22), (1.24): a closing \texttt{\textbackslash Bigg\textbackslash\}} closes nothing. Group~I Eqs.~(1.16), (1.17), (1.30): $\int_{x,z,y,x',y'w'}\to\int_{x,z,y,x',y',w'}$; Eqs.~(1.1), (1.2): $\int_{w,w',x,x',}$ has a stray comma. Group~II Row~I: ``We write it in two parts'' twice.
\end{enumerate}

\paragraph{Still open from before (not errors in these files).}
\begin{itemize}
\item The logarithms of Group~I Row~III are principal values, $\log|\cdot|$: their arguments change sign at $p^+=k^+\frac{(x-w)^2}{(x-z)^2}$. With $\log|\cdot|$ every one of $d^3N_1$, $d^3N_2$, $d^3N_{3a}$, $d^3N_{3b}$ is the exact $p^+$ integral (Sec.~\ref{sec:table}).
\item The edge kernel of Group~I Row~III ($\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$ in $d^3N_{3a}$, the last line of $d^3N_{3b}$) is exact as a $p^+$ integral, but not integrable at $z\to w$ once $\Lambda\to0$: for the JIMWLK term it is replaced by $\mathcal D$ (F3), and its finite part is still open (C4) (\texttt{LL\_terms.pdf}, \texttt{Finite\_terms\_rows.pdf}).
\end{itemize}

\paragraph{The $\log\frac{\vee}{k^+}$ (BFKL) terms of Group III and of Group II Row IV are not missing.} They are inside line 1 of Group~III Rows~III (both parts), V, VI, Eqs.~(1.7), (1.10), (1.14), (1.17), and of Group~II Row~IV, Eq.~(1.21): in the $p^+$ (resp.\ $\xi$) integral, which you keep unevaluated, and they appear only after the $y'$ integral (at large $p^+$ only $y'\to z$ contributes; \texttt{Why\_log\_vee.pdf}). So your equations are exact and correct as they stand, and nothing has to be added to them. To make the BFKL term visible, the corrected copies now write each of these equations, in blue and without changing it, as
\begin{equation*}
\text{your equation}=\underbrace{\text{(B)}}_{\propto\log\frac{\vee}{k^+},\ \text{BFKL}}+\underbrace{\big[\text{your equation}-\text{(B)}\big]}_{\text{finite for }\vee\to\infty},
\end{equation*}
where $-$(B) is written as a $p^+$ integral, $-\text{(B)}=\frac{C}{2\pi}\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}\,F^z$, since $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$ exactly; with (B) as in \texttt{BFKL\_rows.pdf} (and $H_4$, your \texttt{BFKL.tex} line, for Group~II Row~IV). The visible $\log(\frac{\vee}{k^+}-1)$ of line~2 is not a BFKL term: it cancels the subtraction of $\left[\frac{1}{p^+}\right]_+$ (resp.\ $\frac{\xi}{(1-\xi)_+}$) of line~1. (Your \texttt{BFKL.tex} calls these Group~III Rows III, IV, V; in \texttt{Group\_III.tex} they are Rows III, V, VI, and Row~IV has only four $\rho$'s.)

\paragraph{Everything else is correct.} Group~I: Row~II part~1 (1.1)$\to$(1.2); Row~III, all four parts, exactly in $\Lambda$; Row~IV; Row~V part~1, part~2a, 2b finite, including the $q^+=k^+-p^+$ rewriting and every plus prescription. Group~II: Row~I (part~1, and part~2 = 2a+\dots+2f exactly in $\Lambda$); Row~IV ($w,w'$ integrals with their Jacobian, Wilson lines, $p^+\to\xi$, $\mathcal P_{gg}$, plus prescription in $\xi$; the $z$-integral formula is the one checked in \texttt{RowIV\_dperp.pdf}). Group~III: Rows~I, II, III (both parts), V (including the Wilson-line simplification), VI. Numbers in Sec.~\ref{sec:table}.

\section{M1: Group I Row V, 2b UV}\label{sec:M1}
\paragraph{Your exact form is right.} The UV term multiplies the $p^+$ integral of the splitting function over your range $\Lambda<p^+<k^+-\Lambda$; with $\xi=\frac{p^+}{k^+}$, $\epsilon=\frac{\Lambda}{k^+}$,
\begin{equation*}
\int_\epsilon^{1-\epsilon}d\xi\Big[\xi(1-\xi)+\frac{\xi}{1-\xi}+\frac{1-\xi}{\xi}\Big]=2\log\frac{1-\epsilon}{\epsilon}-\frac{11}{6}+4\epsilon-\epsilon^2+\frac{2}{3}\epsilon^3 ,
\end{equation*}
which is your Eq.~(1.25), since $2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^+}=2\log\frac{k^+}{\Lambda}$ exactly. There is no $\vee$ in it.
\paragraph{Where the 2 is lost.} In Eq.~(1.26) (``Setting $\Lambda\to0$''), $2\log\frac{\vee}{\Lambda}$ is written as $\log\frac{\vee}{k^+}+\log\frac{k^+}{\Lambda}$; it is $2\log\frac{\vee}{k^+}+2\log\frac{k^+}{\Lambda}$. With the $-2\log\frac{\vee}{k^+}$ next to it:
\begin{align*}
\frac{d^{3}N_{2b\,\rm UV\, NLO}}{d^{3}k} &= \frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}\left[\new{2\log\Big(\frac{k^{+}}{\Lambda}\Big)}+2\log\Big(1-\frac{\Lambda}{k^{+}}\Big)-\frac{11}{6}\right]\int_{y,z,x',w'}e^{-ik\cdot(w'-z)}\left[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y-z)^{2}\mu^{2}_{\overline{MS}}}{4}\right]\\
&\quad\times\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}\,(y-z)^{2}}\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[U^{af}(y)-U^{af}(z)\Big]\rho^{b'}(x')\,\rho^{f}(y).
\end{align*}
\paragraph{Your own 2b finite says the same.} 2b $=$ 2b UV $+$ 2b finite adds and subtracts one counterterm. In 2b finite, Eqs.~(1.27)--(1.29), it carries $\big[-\frac{11}{6}+2\log\frac{k^+-\Lambda}{\Lambda}\big]$, i.e.\ $2\log\frac{k^+}{\Lambda}$ again. Numerically, Eq.~(1.23) $=$ Eq.~(1.27) $+$ counterterm, pointwise in $p^+$, to $@@CT@@$; and with $\vee=20$, $k^+=1$:
\begin{center}
\begin{tabular}{lrr}\toprule
 & $\Lambda=10^{-2}$ & $\Lambda=10^{-3}$\\\midrule
$p^+$ integral (numerical) & $@@P2@@$ & $@@P3@@$\\
your Eq.~(1.25) & $@@O2@@$ & $@@O3@@$\\
$2\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}$ & $@@C2@@$ & $@@C3@@$\\
your Eqs.~(1.26), (1.32) & $@@Y2@@$ & $@@Y3@@$\\
your (1.32) minus the correct one & $@@D2@@$ & $@@D3@@$\\
$-\log\frac{\vee}{\Lambda}$ & $@@L2@@$ & $@@L3@@$\\\bottomrule
\end{tabular}
\end{center}
(The correct one differs from the exact $p^+$ integral by the $O(\Lambda)$ terms $4\epsilon-\epsilon^2+\frac23\epsilon^3$ only.)
\paragraph{Consequences.} The poles $\times\log\frac{k^+}{\Lambda}$ of Group~I Row~II part~1, Group~I Row~V 2b UV and Group~II Row~IV are in the ratio $1:-2:1$ and cancel (\texttt{LL\_terms.pdf}, \texttt{pole\_check.json}); with your coefficient they are $1:-1:1$ and do not. The new $-\log\frac{\vee}{k^+}$ adds a pole $\times\log\frac{\vee}{k^+}$ that no other row has. (In \texttt{BFKL\_rows.pdf} this line was red: it is to be removed, not moved into Group~I.)

\section{M2--M4 in your notation}
\paragraph{M2, Group II Eq.~(1.16), second term.}
\begin{equation*}
-\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{4\pi^{3}}\,\frac{N_{c}}{k^{+}}\,\int_{\Lambda}^{\vee-k^{+}}\frac{dp^{+}}{2\pi}\int_{x,y,y^{\prime},x'}\,e^{-ik\,\left(y'-y\right)\,\old{\frac{p^{+}}{(p^{+}+k^{+})}}\new{\frac{p^{+}+k^{+}}{k^{+}}}}\,\frac{(x-y)^{k}\,(x'-y')^{k'}\,(y'-y)^{m}\,(y'-y)^{m'}}{(x-y)^{2}\,(x'-y')^{2}\,(y'-y)^{2}}\bigg[\cdots\bigg]\times\cdots
\end{equation*}
\paragraph{M3, Group III Row III, first line.} $\displaystyle\frac{d^{3}N_{\rm NLO}}{d^{3}k} = \new{-}\,\frac{4}{(2\pi)^{3}} \int_{\Lambda}^{\vee -k^{+}} \frac{dp^{+}}{2\pi}\int_{z,w,w^{\prime}}\cdots$
\paragraph{M4, Group III Row III part 2.} Eq.~(1.8): $\int_{w,w',\new{x'},y',z,\old{y},x}$; Eqs.~(1.9), (1.10) (both lines): $\int_{w,\old{w'},\new{x'},y',z,\old{y},x}$.

\section{All steps, numerically}\label{sec:table}
Relative difference between the two sides of each step (largest of two configurations; $\vee=20$, $k^+=1$, $\Lambda=10^{-3}$ unless stated). ``fixed $p^+$'': both sides have the same longitudinal integral and the integrands are compared at $p^+=0.3,\,1.7,\,6$ (Group~I Row~V: $0.13,\,0.5,\,0.81$); a $\delta^{(2)}$ is done by moving the point it fixes, with its Jacobian. ``plus, $\Lambda\to0$'': the step uses your plus prescription; the right side starts at $p^+=0$ (or $\xi=1$) and has $\Lambda\to0$, so the difference must be $O(\Lambda)$; it drops by 10 from $\Lambda=10^{-6}$ to $10^{-7}$.
{\footnotesize
\begin{longtable}{@{}l>{\raggedright\arraybackslash}p{5.6cm}>{\raggedright\arraybackslash}p{4.0cm}lr@{}}
\toprule
Gr. & step & equations & kind & rel.\ difference (plus: $\Lambda=10^{-6}$ / $10^{-7}$)\\\midrule
\endhead
@@TABLE@@
\bottomrule
\end{longtable}}
\noindent Also: Group~I Row~V, Eq.~(1.23) $=$ Eq.~(1.27) $+$ counterterm at fixed $p^+$: $@@CT@@$. Group~I Row~III is evaluated with $\Lambda=10^{-2}$ (exact in $\Lambda$, with the principal value at $p^+=k^+\frac{(x-w)^2}{(x-z)^2}$, which lies inside the range for one configuration and outside for the other).

\paragraph{Measures.} A scan of every equation (points that appear vs.\ points integrated) finds only M4 and M5 (\texttt{measures.json}).

\section{How it was checked}
\texttt{rows\_check/pe.py} reads each equation from your file and evaluates it: transverse vectors and indices, Kronecker deltas, adjoint Wilson lines and $f^{abc}$ as tensors, classical charges, phases (with their $p^+$, $q^+$ or $1/\xi$ factors), $\log\to\log|\cdot|$, $\mathcal A$; repeated transverse indices are contracted term by term. Your plus prescriptions are used as you define them: $\left[\frac{1}{p^+}\right]_+f=\frac{1}{p^+}\big(f(p^+)-f(0)\big)$, the same for $\left[\frac{k^+}{p^+}\right]_+$ and $\left[\frac{1}{q^+}\right]_+$, and $\frac{d\xi}{\xi^2}\frac{\xi}{(1-\xi)_+}g=\frac{d\xi}{\xi(1-\xi)}\big(g(\xi)-g(1)\big)$. $d_\perp=2-\epsilon$; steps without a pole are evaluated at $\epsilon=0$, because the numerical index sums are two-dimensional. Scripts: \texttt{check\_all.py} (all steps), \texttt{measures.py}, \texttt{make\_corrected.py}, \texttt{make\_note.py}.
\end{document}
'''

if __name__ == '__main__':
    m2 = J['Group_II']['Row IV: after the z integral -> simplifying']
    m2b = J['Group_II']['Row IV: the same, phase of the 2nd term of Eq. (111) corrected']
    rep = {
        '@@M2A@@': '%.2f' % min(r['rel'] for r in m2), '@@M2B@@': tex_num(max(r['rel'] for r in m2b)),
        '@@CT@@': tex_num(max(x['two_b_equals_finite_plus_CT'] for x in J['Group_I']['Row V 2b: UV term'])),
        '@@P2@@': '%.4f' % U2['p_integral'], '@@P3@@': '%.4f' % U3['p_integral'],
        '@@O2@@': '%.4f' % U2['your_Eq_1_31'], '@@O3@@': '%.4f' % U3['your_Eq_1_31'],
        '@@C2@@': '%.4f' % U2['with_2log'], '@@C3@@': '%.4f' % U3['with_2log'],
        '@@Y2@@': '%.4f' % U2['your_final_bracket'], '@@Y3@@': '%.4f' % U3['your_final_bracket'],
        '@@D2@@': '%.4f' % U2['yours_minus_correct'], '@@D3@@': '%.4f' % U3['yours_minus_correct'],
        '@@L2@@': '%.4f' % U2['minus_log_vee_over_Lambda'], '@@L3@@': '%.4f' % U3['minus_log_vee_over_Lambda'],
        '@@TABLE@@': table(),
    }
    doc = DOC
    for k, v in rep.items():
        doc = doc.replace(k, v)
    assert '@@' not in doc
    open(os.path.join(HERE, 'Group_check.tex'), 'w').write(doc)
    print('written Group_check.tex')
