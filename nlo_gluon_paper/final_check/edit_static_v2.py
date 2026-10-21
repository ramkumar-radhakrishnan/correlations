"""One-off edit of final_static.tex: answers on the log(vee/Lambda) convention, the Lambda/vee dependence and the
running coupling (arXiv:2308.15545).  Every replacement is asserted to match once."""
p = 'final_static.tex'
t = open(p).read()


def rep(a, b):
    global t
    assert t.count(a) == 1, (t.count(a), a[:80])
    t = t.replace(a, b)


rep(r"""\title{Final check of the evolution terms:\\ three and four $\rho$'s, two $\rho$'s, JIMWLK, BFKL and DGLAP}""",
    r"""\title{Final notes on the evolution terms:\\ JIMWLK and BFKL with $\log\frac{\vee}{\Lambda}$, the scales $\Lambda$ and $\vee$, running coupling and DGLAP}""")

rep(r"""\section{Summary}
\paragraph{Results.}""", r"""\section{Summary}
\paragraph{Answers to your questions.}
\begin{enumerate}
\item \textbf{Your convention} $\log\frac{\Lambda}{k^+}=-\log\frac{\vee}{\Lambda}+\log\frac{\vee}{k^+}$ is exact, so the rewriting is always allowed (Sec.~\ref{sec:why}).
\begin{itemize}
\item \textbf{JIMWLK matches.} With this rule the coefficient of $\log\frac{\vee}{\Lambda}$ is JIMWLK$\otimes$LO, exactly. This needs the two-$\rho$ lists of \texttt{4rho.tex} (M1) and your JIMWLK formula corrected (M2).
\item \textbf{BFKL matches, but not as the coefficient of $\log\frac{\vee}{\Lambda}$ with this rule.} BFKL$\otimes$LO is exactly the coefficient of $\log\vee$. With your rule it sits in the left-over term $\log\frac{\vee}{k^+}(\mathcal B-\mathcal J)$, which is not finite: at mid rapidity it is as large as the evolution term.
\item \textbf{What gives both with $\log\frac{\vee}{\Lambda}$:} at mid rapidity, $k^+=\sqrt{\vee\Lambda}$, so $\log\frac{\vee}{k^+}=\log\frac{k^+}{\Lambda}=\frac12\log\frac{\vee}{\Lambda}$. Then, with nothing left over,
\begin{equation*}
\frac{d^3N_{\rm LL}}{d^3k}=\frac12\log\frac{\vee}{\Lambda}\,\big[{\rm JIMWLK}\otimes{\rm LO}+{\rm BFKL}\otimes{\rm LO}\big].
\end{equation*}
Away from mid rapidity the left-over is $\frac12\log\frac{k^{+2}}{\vee\Lambda}(\mathcal J-\mathcal B)$, with $\log\frac{k^{+2}}{\vee\Lambda}=2(y-y_{\rm mid})$.
\end{itemize}
\item \textbf{Scales.} Every term that depends on $\Lambda$ or $\vee$ is accounted for (Sec.~\ref{sec:scales}):
\begin{itemize}
\item the logarithms go into evolution (two $\rho$'s) or cancel (three and four $\rho$'s);
\item the $\frac1\Lambda$ terms cancel;
\item the poles multiplying $\log\frac{\vee}{\Lambda}$ cancel;
\item everything else is $O(\Lambda/k^+)$ or $O(k^+/\vee)$ and vanishes.
\end{itemize}
\item \textbf{Running coupling} (Sec.~\ref{sec:rc}, with arXiv:2308.15545):
\begin{itemize}
\item Your Group I Row II is exactly the NLO JIMWLK kernel $K_{JSJ}$ of that paper, including the constant $(\frac{67}{9}-\frac{\pi^2}{3})N_c$. It is twice the running-coupling logarithm.
\item Your Group I Row V is their dressed gluon: the same splitting weight and the same double-plus constant $-\frac{11}{6}=-\frac{\beta_0}{2N_c}$. It dresses both Wilson lines of the measured gluon (amplitude and conjugate), so it removes twice the running-coupling logarithm.
\item The two cancel exactly. With their scale choice, LO runs as $g(c_0^2/X^2)\,g(c_0^2/X'^2)$, with $X=x-w$, $X'=x'-w'$ the emitter--gluon distances. The rest is DGLAP-like and goes with the fragmentation function.
\end{itemize}
\end{enumerate}

\paragraph{Results.}""")

rep(r"""\item DGLAP for the fragmentation function: Eq.~\eqref{eq:DGLAP}.
\end{itemize}
""", r"""\item DGLAP for the fragmentation function: Eq.~\eqref{eq:DGLAP}.
\end{itemize}
At mid rapidity, $k^+=\sqrt{\vee\Lambda}$:
\begin{equation}
\boxed{\frac{d^3N_{\rm LL}}{d^3k}\bigg|_{\rm mid}=\frac12\,\NN\,\big(\mathcal J+\mathcal B\big)}
\label{eq:finalmid}
\end{equation}
""")

old_why = t[t.index(r"\section{Why $\log\frac{\vee}{\Lambda}$ carries JIMWLK and not BFKL}"):t.index(r"\section{The targets}")]
new_why = r"""\section{Your convention: BFKL and JIMWLK with $\log\frac{\vee}{\Lambda}$}\label{sec:why}
\subsection{Where the logarithms come from}
Every row is an exact $p^+$ integral of the form $\int dp^+\,F(p^+)$.
\begin{itemize}
\item $\Lambda$ enters only through the limit where the unobserved gluon is softer than the measured one ($p^+\to\Lambda$). There $F\to\mathfrak T/p^+$.
\item $\vee$ enters only through the upper limit, where the unobserved gluon is harder than the measured one ($p^+\to\vee$). There $F\to\mathfrak P/p^+$.
\end{itemize}
So, exactly and with unique coefficients,
\begin{equation}
\int dp^+F(p^+)=\mathfrak T\,\log\frac{k^+}{\Lambda}+\mathfrak P\,\log\frac{\vee}{k^+}+R ,
\label{eq:TP}
\end{equation}
where $R$ contains neither logarithm and stays finite for $\Lambda\to0$, $\vee\to\infty$. Summed over all rows, $\sum\mathfrak T=\mathcal J={\rm JIMWLK}\otimes{\rm LO}$ (Sec.~\ref{sec:two}) and $\sum\mathfrak P=\mathcal B={\rm BFKL}\otimes{\rm LO}$ (Sec.~\ref{sec:bfkl}), both exactly. Therefore
\begin{equation}
\frac{d^3N_{\rm LL}}{d^3k}=\mathcal N\Big[\log\frac{k^+}{\Lambda}\;\mathcal J+\log\frac{\vee}{k^+}\;\mathcal B\Big] .
\label{eq:LLexact}
\end{equation}
In rapidity, $\delta Y=\log\frac{\vee}{\Lambda}=y_P-y_T$, $\log\frac{\vee}{k^+}=y_P-y$ and $\log\frac{k^+}{\Lambda}=y-y_T$. The target evolves (JIMWLK) between $y_T$ and the gluon; the projectile evolves (BFKL) between the gluon and $y_P$. This is the standard structure of single inclusive gluon production.

\subsection{Writing everything with $\log\frac{\vee}{\Lambda}$}
Each of the following is an exact rewriting of Eq.~\eqref{eq:LLexact} (\texttt{sympy}):
\begin{align}
\text{(A) your rule, }\log\tfrac{k^+}{\Lambda}=\log\tfrac{\vee}{\Lambda}-\log\tfrac{\vee}{k^+}:&\quad\frac{d^3N_{\rm LL}}{d^3k}=\mathcal N\Big[\log\frac{\vee}{\Lambda}\;\mathcal J+\log\frac{\vee}{k^+}\,(\mathcal B-\mathcal J)\Big],\label{eq:A}\\
\text{(B) the other rule, }\log\tfrac{\vee}{k^+}=\log\tfrac{\vee}{\Lambda}-\log\tfrac{k^+}{\Lambda}:&\quad\frac{d^3N_{\rm LL}}{d^3k}=\mathcal N\Big[\log\frac{\vee}{\Lambda}\;\mathcal B+\log\frac{k^+}{\Lambda}\,(\mathcal J-\mathcal B)\Big],\label{eq:B}\\
\text{(C) both:}&\quad\frac{d^3N_{\rm LL}}{d^3k}=\mathcal N\Big[\frac12\log\frac{\vee}{\Lambda}\,(\mathcal J+\mathcal B)+\frac12\log\frac{k^{+2}}{\vee\Lambda}\,(\mathcal J-\mathcal B)\Big].\label{eq:C}
\end{align}
\begin{table}[H]\centering\small
\begin{tabular}{@{}lllll@{}}\toprule
rule & coefficient of $\log\frac{\vee}{\Lambda}$ & left-over & its log in rapidity & left-over small when\\\midrule
(A) yours & $\mathcal J$: JIMWLK matches & $\log\frac{\vee}{k^+}(\mathcal B-\mathcal J)$ & $y_P-y$ & gluon near the projectile\\
(B) & $\mathcal B$: BFKL matches & $\log\frac{k^+}{\Lambda}(\mathcal J-\mathcal B)$ & $y-y_T$ & gluon near the target\\
(C) & $\frac12(\mathcal J+\mathcal B)$: both & $\frac12\log\frac{k^{+2}}{\vee\Lambda}(\mathcal J-\mathcal B)$ & $2(y-y_{\rm mid})$ & gluon at mid rapidity (then $=0$)\\\bottomrule
\end{tabular}
\end{table}
\paragraph{Does it make sense?}
\begin{itemize}
\item The rule is exact, so the rewriting is always allowed. The only question is whether the left-over may be called a finite contribution. It may only if its logarithm is not large.
\item With your rule (A) the left-over is $\log\frac{\vee}{k^+}(\mathcal B-\mathcal J)$. At mid rapidity $\log\frac{\vee}{k^+}=\frac12\log\frac{\vee}{\Lambda}$, so it is $\frac12\log\frac{\vee}{\Lambda}(\mathcal B-\mathcal J)$: as large as the evolution term. Putting it into the finite part would drop the BFKL evolution of the projectile.
\item For a gluon at mid rapidity, (C) is the consistent choice. There $k^+=\sqrt{\vee\Lambda}$, the left-over vanishes identically, and Eq.~\eqref{eq:finalmid} holds:
\begin{equation}
\frac{d^3N_{\rm LL}}{d^3k}\bigg|_{\rm mid}=\frac12\,\mathcal N\log\frac{\vee}{\Lambda}\,\big({\rm JIMWLK}\otimes{\rm LO}+{\rm BFKL}\otimes{\rm LO}\big).
\end{equation}
In your normalization this is $\frac12\big[\sum_{i=1}^5\frac{dN^{(i)}}{d^2wd^2w'dy}+\frac{dN}{d^2wd^2w'dy}\big|_{\rm BFKL}\big]$ with $\delta Y=\log\frac{\vee}{\Lambda}$ in Eqs.~\eqref{eq:JIMc} and \eqref{eq:BFKL}. Both evolutions multiply $\log\frac{\vee}{\Lambda}$; each acts over half of the rapidity interval.
\item Away from mid rapidity, Eq.~\eqref{eq:LLexact}, or (C) with its left-over, must be kept.
\end{itemize}
\paragraph{So:}
\begin{itemize}
\item \textbf{JIMWLK matches}, exactly, as the coefficient of $\log\frac{\vee}{\Lambda}$ with your rule (Sec.~\ref{sec:two}).
\item \textbf{BFKL matches}, exactly, as the coefficient of $\log\vee$ (Sec.~\ref{sec:bfkl}). Equivalently it is the coefficient of $\log\frac{\vee}{\Lambda}$ with rule (B), or of $\frac12\log\frac{\vee}{\Lambda}$ at mid rapidity, rule (C).
\end{itemize}
The colours show this term by term. Every {\color{BFK}orange} kernel (the measured gluon emitted by the gluon at $z$) is absent from the $\log\frac{\vee}{\Lambda}$ terms, and every {\color{JIM}green} kernel is absent from the $\log\vee$ terms.

%==========================================================================================
"""
t = t.replace(old_why, new_why)

scales = r"""%==========================================================================================
\section{Every term that depends on $\Lambda$ or $\vee$}\label{sec:scales}
In your files $\Lambda$ and $\vee$ enter only through the limits of the $p^+$ integrals (and of the $\xi$ integral of Group II Row IV). After integration every row has the form of Eq.~\eqref{eq:TP}. Below are all the types of $\Lambda$- or $\vee$-dependence that occur in \texttt{Group\_I.tex}, \texttt{Group\_II.tex} and \texttt{Group\_III.tex}, and where each one goes.
\begin{table}[H]\centering\small
\begin{tabular}{@{}p{0.36\textwidth}p{0.25\textwidth}p{0.33\textwidth}@{}}\toprule
type & where & fate\\\midrule
$\log\frac{\vee}{\Lambda}$, two $\rho$'s & all groups & evolution: $\mathcal J$ (rule A) or $\frac12(\mathcal J+\mathcal B)$ (rule C), Sec.~\ref{sec:two}\\
$\log\frac{\vee}{\Lambda}$, three $\rho$'s & G1 II--V, G2 I, G3 I--III, V, VI & cancel after symmetrization, Sec.~\ref{sec:three}\\
$\log\frac{\vee}{\Lambda}$, four $\rho$'s & G1 I, G2 II, III, G3 IV & cancel after symmetrization, Sec.~\ref{sec:three}\\
$\log\frac{\vee}{k^+}$, two $\rho$'s: explicit ($\log\frac{\vee}{k^+}$, $\log\frac{\vee}{k^+\pm\Lambda}$, $\log\frac{k^+-\Lambda}{\vee}$, $\log(\frac{\vee}{k^+}-1)$) and hidden (Eqs.~\eqref{eq:asI}, \eqref{eq:asIII}) & all groups & $\mathcal B-\mathcal J$, Eq.~\eqref{eq:A}; at mid rapidity part of $\frac12\log\frac{\vee}{\Lambda}(\mathcal J+\mathcal B)$\\
$\log\frac{\vee}{k^+}$, three $\rho$'s / four $\rho$'s & & cancel / none (\texttt{LogVeeKplus\_BFKL.pdf})\\
$\frac1\Lambda$, $\frac{1}{\vee-k^+}$, $\frac{1}{k^+-\Lambda}$ & G1 III ($N_1$, $N_2$, $N_{3a}$) & cancel exactly, Eq.~\eqref{eq:green}\\
poles $\times\log\frac{\vee}{\Lambda}$ & G1 II, G1 V 2b UV, G2 IV & cancel, Sec.~\ref{sec:dglap}; the finite rest is part of $\mathcal J$\\
power corrections (list below) & all groups & vanish for $\Lambda/k^+\to0$, $k^+/\vee\to0$\\\bottomrule
\end{tabular}
\end{table}
\paragraph{The logarithms with $\vee$ and $\Lambda$ inside.} With $A=(x-z)^2$, $B=(x-w)^2$ (or the primed ones), exactly:
\begin{align}
\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}&=\log\frac{\vee}{k^+}+\log\frac{A}{B}+\log\Big(1-\frac{k^+}{\vee}+\frac{k^+B}{\vee A}\Big)-\log\Big(1+\frac{\Lambda A}{k^+B}\Big),\label{eq:asI}\\
\log\Big|\frac{k^+B-\vee A}{k^+B-(k^++\Lambda)A}\Big|&=\log\frac{\vee}{k^+}+\log\Big|\frac{A}{A-B}\Big|+O\Big(\frac{k^+}{\vee},\frac{\Lambda}{k^+}\Big).\label{eq:asIII}
\end{align}
The first occurs in Groups II and III and Group I Row IV; the second in Group I Row III ($N_1$, $N_2$). Both limits were checked with \texttt{sympy}. The other logarithms have no large part:
\begin{itemize}
\item $\log\big|\frac{k^+B-(k^+-\Lambda)A}{k^+B-\Lambda A}\big|=\log\big|\frac{B-A}{B}\big|+O(\frac{\Lambda}{k^+})$ ($N_1$, $N_{3a}$);
\item $\log\frac{[k^+(x-w)-(k^+-\Lambda)(x-z)]^2}{[k^+(x-w)-\Lambda(x-z)]^2}=\log\frac{(z-w)^2}{(x-w)^2}+O(\frac{\Lambda}{k^+})$ ($N_{3b}$);
\item $\mathcal A$ ($N_{3a}$, $N_{3b}$): $\Lambda$ only at $O(\frac{\Lambda}{k^+})$;
\item $\log\frac{\Lambda(x-z)^2+k^+(x-w)^2}{k^+(x'-w')^2+\Lambda(x'-z)^2}=\log\frac{(x-w)^2}{(x'-w')^2}+O(\frac{\Lambda}{k^+})$ (Group II Row I);
\item $R^i$ (Group I Row II part 2): $\vee$ only through $\frac{k^+}{\vee}$, which can be set to zero.
\end{itemize}
\paragraph{The $\frac1\Lambda$ terms.} With $D=(x-z)^2-(x-w)^2$ and the common factor $C=\frac{ig^4f^{abd}}{8\pi^5}\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2}\times$ colour, the four terms of $N_1$, $N_2$ (overall minus), $N_{3a}$ give
\begin{equation}
\frac{C}{D}\Big\{\Big[\frac1\Lambda-\frac{1}{\vee-k^+}\Big]+\Big[\frac1\Lambda-\frac{1}{k^+-\Lambda}\Big]+\Big[\frac{1}{\vee-k^+}-\frac1\Lambda\Big]-\Big[\frac1\Lambda-\frac{1}{k^+-\Lambda}\Big]\Big\}=0 .
\label{eq:green}
\end{equation}
\paragraph{Power corrections} (all vanish for $\Lambda\to0$, $\vee\to\infty$):
\begin{itemize}
\item $\frac{4\Lambda}{k^+}-\frac{\Lambda^2}{k^{+2}}+\frac{2\Lambda^3}{3k^{+3}}$ (Row V 2b), and the powers of $\delta=\Lambda/k^+$ in the exact Row II bracket;
\item $\log(1\pm\frac{\Lambda}{k^+})$;
\item $\log(1-\frac{k^+}{\vee})$, $\log^2(1-\frac{k^+}{\vee})$, and $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}=O\big(\frac{k^+}{\vee}\log\frac{\vee}{k^+}\big)$;
\item the lower limits $0$ instead of $\Lambda$ of the plus prescriptions: $O(\Lambda/k^+)$;
\item the remainder $-N_c\log(1-\frac{k^+}{\vee})$ in the $\delta(1-\xi)$ tally of Sec.~\ref{sec:dglap}: $O(k^+/\vee)$.
\end{itemize}
\paragraph{Result.}
\begin{equation}
\frac{d^3N}{d^3k}={\rm LO}\otimes D+\mathcal N\Big[\log\frac{k^+}{\Lambda}\,\mathcal J+\log\frac{\vee}{k^+}\,\mathcal B\Big]+\text{finite}+O\Big(\frac{\Lambda}{k^+},\frac{k^+}{\vee}\Big),
\end{equation}
with no other dependence on $\Lambda$ or $\vee$. At mid rapidity the bracket is $\frac12\log\frac{\vee}{\Lambda}(\mathcal J+\mathcal B)$. The scale $\mu$ cancels in the $\log\frac{\vee}{\Lambda}$ sector (Sec.~\ref{sec:dglap}) and between the $\beta_0$ terms (Sec.~\ref{sec:rc}); what is left of it is the scale of $D$.

The $\log\frac{\vee}{k^+}$ rows of the table rest on the earlier check of the same rows (\texttt{LogVeeKplus\_BFKL.pdf}). The fixes found there (P1, F3, F5, F6) are now in your files, but that check was not repeated row by row on these files.

"""
rep(r"""%==========================================================================================
\section{DGLAP}\label{sec:dglap}""", scales + r"""%==========================================================================================
\section{DGLAP}\label{sec:dglap}""")

rc = r"""%==========================================================================================
\section{Running coupling}\label{sec:rc}
Reference: A.~Kovner, M.~Lublinsky, V.~Skokov, Z.~Zhao, \emph{Not all that is $\beta_0$ is $\beta$-function: the DGLAP resummation and the running coupling in NLO JIMWLK}, arXiv:2308.15545 (KLSZ below).

\subsection{The $\beta_0$ terms of your cross section}
Write
\begin{itemize}
\item $X=x-w$, $X'=x'-w'$: the distances between the emitting charge and the measured gluon, in the amplitude and in the conjugate amplitude;
\item $E_X=\frac2\epsilon+\ell_X$, with $\ell_X=2\gamma_E+\log\frac{X^2\mu^2_{\overline{MS}}}{4}=\log\frac{X^2\mu^2_{\overline{MS}}}{c_0^2}$ and $c_0=2e^{-\gamma_E}$;
\item ${\rm LO}=\frac{1}{(2\pi)^3}\frac{g^2}{2\pi^2k^+}\int e^{-ik\cdot(w'-w)}\mathcal K(x-w)\cdot\mathcal K(x'-w')\,\mathcal O$, with $\mathcal O=[U(w')-U(x')]^{ab'}\rho^{b'}(x')[U(w)-U(x)]^{ab}\rho^b(x)$.
\end{itemize}
Relative to LO (inside the integrals), with the c.c.\ of Group I and $\beta_0=\frac{11}{3}N_c$ (no quarks):
\begin{itemize}
\item Group I Row II, Eq.~(1.2), the part without $\log\frac{\vee}{\Lambda}$, $\log\frac{\vee}{k^+}$: $\frac{\alpha_s}{2\pi}\big[\beta_0\,\frac12(E_X+E_{X'})+N_c(\frac{67}{9}-\frac{\pi^2}{3})\big]$;
\item Group I Row V 2b UV, Eq.~(1.28), after $y\to x$, $z\to w$: $-\frac{\alpha_s\beta_0}{4\pi}(E_X+E_{X'})$;
\item charge renormalization of the LO coupling, $\alpha_0=\alpha_s(\mu)\big[1-\frac{\alpha_s\beta_0}{4\pi}\frac2\epsilon\big]$: $-\frac{\alpha_s\beta_0}{4\pi}\frac{2}{\epsilon}$;
\item Group II Row IV: no $\beta_0$.
\end{itemize}

\subsection{Comparison with KLSZ}
\paragraph{Group I Row II is $K_{JSJ}$.} KLSZ Eq.~(7) gives, relative to the LO kernel $\frac{\alpha_s}{2\pi^2}\frac{X\cdot Y}{X^2Y^2}$,
\begin{equation*}
\frac{K_{JSJ}}{K_{\rm LO}}=\frac{\alpha_s}{2\pi}\Big[\beta_0\,\frac12\big(\ln X^2\mu^2+\ln Y^2\mu^2\big)+N_c\Big(\frac{67}{9}-\frac{\pi^2}{3}\Big)\Big],\qquad \mu^2=\frac{\mu^2_{\overline{MS}}}{4}e^{2\gamma_E}.
\end{equation*}
So $\ln X^2\mu^2=\ell_X$, and your Row II is identical, with $Y\to X'$: in the cross section the second emission is in the conjugate amplitude. The prefactor ratio was checked with \texttt{sympy}. This is twice the running-coupling logarithm (KLSZ, Sec.~II.C).

\paragraph{Group I Row V is the dressed gluon} (KLSZ Eqs.~(25)--(31)):
\begin{itemize}
\item \textbf{Colour.} After your identity, part 2a $[-f^{adc}f^{fbe}U^{db}(x)U^{ce}(z)+N_cU^{af}(x)]$ is KLSZ's $D(z_1,z_2)-N_cS$ (real minus virtual): UV finite. Part 2b, $N_c[U^{af}(y)-U^{af}(x)]$, carries the virtual $N_cS$ part with the LO colour.
\item \textbf{Longitudinal weight.} With $\xi=p^+/k^+$, $\Phi(\xi)=\xi(1-\xi)+\frac{\xi}{1-\xi}+\frac{1-\xi}{\xi}=\frac{\xi^2+(1-\xi)^2+\xi^2(1-\xi)^2}{\xi(1-\xi)}$. This is KLSZ's integrand (difference 0, \texttt{sympy}).
\item \textbf{Double plus prescription.} Exactly, $\int_{\lambda}^{1-\lambda}d\xi\,\Phi(\xi)=2\log\frac{1-\lambda}{\lambda}-\frac{11}{6}+4\lambda-\lambda^2+\frac23\lambda^3$, with $\lambda=\Lambda/k^+$. The constant is KLSZ Eq.~(28): $\int_0^1\frac{d\xi}{\xi_+(1-\xi)_+}[\cdots]=-\frac{11}{6}=-\frac{\beta_0}{2N_c}$ (\texttt{sympy}). KLSZ subtract both poles $\xi=0,1$ because they belong to the energy evolution. In your file the two endpoint logarithms are the $\log\frac{\vee}{\Lambda}$ terms of Row V, which go into JIMWLK$\otimes$LO (Sec.~\ref{sec:two}): the same separation.
\item \textbf{Transverse.} The counterterm kernel $K(z,z,x)-K(y,z,x)=\frac{(y-x)\cdot(y-z)}{(y-x)^2(x-z)^2}$ is $\frac{1}{(x-z)^2}$ for $|x-z|\ll|y-z|$ and is cut off at $|x-z|\sim|y-z|$. The size of the gluon pair runs from $1/\mu$ to the emitter distance: KLSZ's resolution $\tilde Q=1/|X|$ (their Eqs.~(35), (36)). Hence $E(y-z)$, i.e.\ KLSZ Eq.~(31) with $Q=1/|X|$.
\item \textbf{Result.} Row V 2b UV $=-\frac{\alpha_s\beta_0}{4\pi}(E_X+E_{X'})\,$LO: the dressing of the measured gluon, once on $U(w)$ (amplitude) and once on $U(w')$ (conjugate amplitude).
\end{itemize}

\paragraph{Counting.} In units $u=\frac{\alpha_s\beta_0}{8\pi}(\ell_X+\ell_{X'})\,$LO, which is what $\alpha_s(\mu)\to\sqrt{\alpha_s(c_0^2/X^2)\,\alpha_s(c_0^2/X'^2)}$ requires (\texttt{sympy}, one-loop running):
\begin{table}[H]\centering\small
\begin{tabular}{@{}lll@{}}\toprule
 & NLO JIMWLK (KLSZ) & your cross section\\\midrule
one-gluon NLO kernel & $K_{JSJ}$: $+2u$ & Group I Row II: $+2u$\\
dressing of the gluon & $-1u$ (one adjoint $S(z)$ in the real term) & Group I Row V 2b: $-2u$ (two Wilson lines $U(w)$, $U(w')$)\\
net & $+1u$: running coupling at $X$, $Y$ & $0$ (exact; numerically $2\times10^{-16}$, exact in $\Lambda$)\\\bottomrule
\end{tabular}
\end{table}
\noindent So in the cross section the $\beta_0$ logarithms of the emitter--gluon distance cancel at fixed order, and so do their poles and $\mu$. The coupling of LO must still be renormalized. Its counterterm leaves $-\frac{\alpha_s\beta_0}{4\pi}\frac2\epsilon\,$LO, which is collinear: it is the $\frac{\beta_0}{2}\delta(1-\xi)$ of the fragmentation counterterm (Sec.~\ref{sec:dglap}).

\subsection{The running coupling in your cross section}
At $O(\alpha_s^2)$, exactly (\texttt{sympy}),
\begin{equation}
\alpha_s(\mu)\,{\rm LO}+\big[\text{Row II}+\text{Row V}\big]_{\beta_0}={\rm LO}\big[\alpha_s\to g(c_0^2/X^2)\,g(c_0^2/X'^2)/4\pi\big]\;-\;\frac{\alpha_s\beta_0}{8\pi}\big(\ell_X+\ell_{X'}\big)\,{\rm LO}+O(\alpha_s^3),
\label{eq:rc}
\end{equation}
with $g^2(Q^2)=4\pi\alpha_s(Q^2)$ and one-loop running, $\alpha_s(c_0^2/X^2)=\alpha_s(\mu)\big[1+\frac{\alpha_s\beta_0}{4\pi}\ell_X\big]+O(\alpha_s^3)$. In words:
\begin{itemize}
\item \textbf{The first term is LO with the running coupling at the emitter--gluon distances.} This is KLSZ's choice, $\lambda=\frac12$ in their Eq.~(41): one $g$ at the emission vertex in the amplitude (distance $X$) and one in the conjugate amplitude (distance $X'$),
\begin{equation*}
\frac{d^3N}{d^3k}\bigg|_{\rm LO,\,rc}=\frac{1}{(2\pi)^3}\frac{1}{2\pi^2k^+}\int_{x,x',w,w'}e^{-ik\cdot(w'-w)}\,g\Big(\frac{c_0^2}{(x-w)^2}\Big)\,g\Big(\frac{c_0^2}{(x'-w')^2}\Big)\,\mathcal K(x-w)\cdot\mathcal K(x'-w')\,\mathcal O .
\end{equation*}
\item \textbf{The second term is DGLAP-like} (KLSZ Sec.~III): the dressing of the measured gluon with resolution $1/|X|$ ($1/|X'|$). It carries all the $\mu$ dependence, $-\frac{\alpha_s\beta_0}{4\pi}\log\mu^2\,$LO, which cancels against the $\frac{\beta_0}{2}\delta(1-\xi)$ term of the evolution of the fragmentation function, Eq.~\eqref{eq:DGLAP}. With the fragmentation scale $\mu_F$ it becomes $-\frac{\alpha_s\beta_0}{8\pi}\big[\log\frac{X^2\mu_F^2}{c_0^2}+\log\frac{X'^2\mu_F^2}{c_0^2}\big]\,$LO. This is a DGLAP-like logarithm between the resolution $1/|X|$ and $\mu_F\sim|\mathbf k|$; KLSZ resum such logarithms with their Eq.~(43).
\item \textbf{If $\alpha_s(\mu)$ is kept in LO,} the natural scale is the fragmentation scale, $\mu=\mu_F\simeq c_0/|y'-y|\sim|\mathbf k|$ (Sec.~\ref{sec:dglap}): then the $\beta_0$ part of the $\mu$ dependence cancels within LO$\otimes D$.
\item \textbf{Constant.} Row II also contains $\frac{\alpha_s}{2\pi}N_c(\frac{67}{9}-\frac{\pi^2}{3})\,$LO $=\frac{\alpha_s}{\pi}K\,$LO, with $K=N_c(\frac{67}{18}-\frac{\pi^2}{6})$: the same constant as in $K_{JSJ}$. It is not a logarithm and does not set the scale.
\item \textbf{Quarks} are not in your calculation, so $\beta_0=\frac{11}{3}N_c$. With quarks, $\beta_0=\frac{11N_c-2N_f}{3}$ and the $q\bar q$ dressing of KLSZ Eq.~(13) is added.
\end{itemize}

"""
rep(r"""%==========================================================================================
\section{Three and four $\rho$'s}\label{sec:three}""", rc + r"""%==========================================================================================
\section{Three and four $\rho$'s}\label{sec:three}""")

rep(r"""three-$\rho$ list: identical to v3; cancellation (\texttt{three\_cancel.py}) & yes; residual $@@TH_RES@@$\\\bottomrule""",
    r"""three-$\rho$ list: identical to v3; cancellation (\texttt{three\_cancel.py}) & yes; residual $@@TH_RES@@$\\
rules (A), (B), (C) equal Eq.~\eqref{eq:LLexact}; (C) at mid rapidity (\texttt{checks\_rc.py}) & exact\\
Row V weight $=$ KLSZ integrand; double-plus constant; exact endpoint form & 0; $-\frac{11}{6}$; $<10^{-100}$\\
Row II $=K_{JSJ}$ (prefactor ratio); $\sqrt{\alpha_s(X)\alpha_s(X')}$ at $O(\alpha_s^2)$; Row II $+$ Row V & 1; $u$; 0\\
limits of Eqs.~\eqref{eq:asI}, \eqref{eq:asIII} & 0, 0\\\bottomrule""")
open(p, 'w').write(t)
print('edited')
