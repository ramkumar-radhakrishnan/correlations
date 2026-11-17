import re
SRC="/home/user/correlations/nlo_review/final/Limiting_Cases.tex"
s=open(SRC).read()
pre,body=s.split(r"\begin{document}")
body=body.split(r"\maketitle")[1].split(r"\end{document}")[0]
pre=re.sub(r"\\title\{[^\n]*\n", r"\\title{Stage 2: check of the evolution terms, limiting finite terms, DGLAP and running coupling}"+"\n", pre, count=1)
def rep(t,old,new,n=1):
    assert t.count(old)==n,(t.count(old),old[:70]); return t.replace(old,new)
# demote sectioning of the embedded body
body=body.replace(r"\subsubsection*{",r"\paragraph{").replace(r"\subsubsection{",r"\paragraph{")
body=body.replace(r"\subsection{",r"\subsubsection{").replace(r"\section{",r"\subsection{")
# ---- notes inserted in the body ----
noteII=r"""
\rnote{From $N_{1}$ only the constant $-\frac{67}{9}+\frac{\pi^{2}}{3}$ survives: the two $\log(1-\frac{k^{+}}{\vee})$ terms written below vanish for $\vee\to\infty$. The pole bracket of $N_{1}$ times $\frac{11}{3}$ is not listed because it cancels exactly against the pole bracket times $-\frac{11}{6}$ of Row V, $2b\,{\rm UV}$ (after $x\to w$, $y\to x$): $\frac{1}{32}\cdot\frac{11}{3}=\frac{1}{16}\cdot\frac{11}{6}$ with opposite signs. $N_{2a}$ is entirely JIMWLK. In $R^{i}$ the limit gives $-2q^{i}\log\frac{p^{2}}{q^{2}}+\frac{(2p^{2}-2q\cdot p+q^{2})q^{i}+q^{2}p^{i}}{(q-p)^{2}}\log\frac{q^{2}}{(q-p)^{2}}$.}
"""
k0=body.index(r"\subsubsection{Row II}", body.index(r"\subsection{Group I (with c.c.)}"))+len(r"\subsubsection{Row II}")
body=body[:k0]+noteII+body[k0:]
noteIV=r"""
\chg{\paragraph{Constant left over from the $\xi\to0$ subtraction.} In \texttt{Group\_II.tex} the BFKL piece of this row is $\big[\log\frac{\vee}{k^{+}}-\big(1-\frac{k^{+}}{\vee}\big)\big]\,C$, and the second term of the DGLAP piece is exactly $-\big[\log\frac{\vee}{k^{+}}-\big(1-\frac{k^{+}}{\vee}\big)\big]C$, because $\int_{k^{+}/\vee}^{1}d\xi\,\frac{1-\xi}{\xi}=\log\frac{\vee}{k^{+}}-1+\frac{k^{+}}{\vee}$. If the DGLAP term is written with $\big[\frac{1-\xi}{\xi}\big]_{+}g(\xi)\equiv\frac{1-\xi}{\xi}\,[g(\xi)-g(0)]$ and only $\log\frac{\vee}{k^{+}}\,C$ is put into BFKL, the term}
\begin{align}
\chg{\frac{d^{3}N}{d^{3}k}=-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\int_{x,x',z,w,w'}e^{-ik\cdot(w'-w)}\,K(w,w',z)\,K(x,x',z)\left[U^{ce'}(z)-U^{ce'}(x')\right]\left[U^{ce}(z)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x)}
\end{align}
\chg{has to be added to the finite terms. It is in none of \texttt{Other\_terms}, \texttt{Limiting\_Cases}, \texttt{Final\_terms}. Alternatively define $\big[\frac{1-\xi}{\xi}\big]_{+}g(\xi)\equiv\frac{(1-\xi)\,g(\xi)-g(0)}{\xi}$, as you did for $\xi\to1$; then the subtraction is exactly $\log\frac{\vee}{k^{+}}\,g(0)$ and this term is absent. One of the two has to be chosen and stated.}
"""
i=body.index(r"\subsubsection{Row IV}", body.index(r"\subsection{Group II}"))
e=body.index(r"\end{align}",i)+len(r"\end{align}")
body=body[:e]+"\n"+noteIV+body[e:]
ruleIII=r"""
\chg{\paragraph{Rule for $\vee\to\infty$ in Rows III, V, VI.} These rows are left as integrals over $p^{+}$. In \texttt{BFKL.tex} each of them contributes three terms: $f(0)\log\frac{\vee}{k^{+}}$ from the blue logarithm, $-f(0)\log\frac{\vee}{k^{+}}$ from the large-$p^{+}$ tail of the plus-prescription integral, and the term from $-\frac{1}{k^{+}}\delta_{im}\delta_{jk'}$ (and $-\frac{p^{+}}{k^{+2}}[\ldots]$). The first two cancel, which is right: $\int_{\Lambda}^{\vee-k^{+}}\frac{dp^{+}}{p^{+}}f(p^{+})$ has no $\log\vee$ because $f$ falls off at large $p^{+}$. But then the plus-prescription integral with the upper limit $\vee-k^{+}$, as written below, still contains that $-f(0)\log\frac{\vee}{k^{+}}$, and the limit cannot be taken. The finite term is}
\begin{align}
\chg{\int_{0}^{\vee-k^{+}}dp^{+}\,\Big[\frac{1}{p^{+}}\Big]_{+}f(p^{+})+f(0)\log\frac{\vee}{k^{+}}\;\xrightarrow{\vee\to\infty}\;\int_{0}^{\infty}\frac{dp^{+}}{p^{+}}\Big[f(p^{+})-\theta(k^{+}-p^{+})\,f(0)\Big],}
\end{align}
\chg{i.e.\ in the formulas below read every upper limit $\vee-k^{+}$ as $\infty$ and every $\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}+\frac{1}{p^{+}+k^{+}}$ as $\frac{1}{p^{+}}$ with the subtraction $\theta(k^{+}-p^{+})f(0)/p^{+}$. The other $p^{+}$ integrals of these rows converge at large $p^{+}$ once the BFKL subtraction (the terms with $\mathcal{K}(x'-z)[U(z)\ldots]$) is made, so there $\vee-k^{+}\to\infty$ directly. The limit $\Lambda\to0$ of the same terms is the soft-sector problem of Sec.~\ref{sec:cannot}.}
"""
body=rep(body,r"\subsection{Group III\, (with c.c.)}",r"\subsection{Group III\, (with c.c.)}"+ruleIII)
MAIN1=r"""
\maketitle
\section{Check of the BFKL and JIMWLK terms}
The BFKL and JIMWLK formulas you sent with this request are character for character those of \texttt{Final\_terms.tex}, and the \texttt{Stage1\_split.tex} you sent back is the one I produced. The question is whether the blue terms of \texttt{Stage1\_split} (plus the four-$\rho$ terms) add up to these two formulas.

\textbf{Result: yes, with two qualifications.}
\begin{itemize}
\item The two formulas are the correct evolution of the leading order (BFKL equation applied to LO; $-H_{\rm JIMWLK}$ acting on the LO trace, all eight kernels).
\item $\log\frac{k^{+}}{\Lambda}$: every one of the 24 two-$\rho$ lists of \texttt{Evolution\_JIMWLK\_\_terms.tex} follows from the blue terms of the corresponding row (8 directly, 16 by symmetrising three and four $\rho$'s); their sum is the JIMWLK formula; the symmetrised three-$\rho$ terms cancel pointwise in all 12 classes.
\item $\log\frac{\vee}{k^{+}}$: for every row the bracket of \texttt{BFKL.tex} is the $\log\frac{\vee}{k^{+}}$ coefficient of the Group file and the two-$\rho$ list follows from it; the 18 lists add up to the BFKL formula; the symmetrised three-$\rho$ terms cancel pointwise.
\item \textbf{Qualification 1, Group I Row II.} For the structure
$X\propto\mathcal{K}^{i}(x'-w')[\mathcal{K}^{i}(x-w)-\mathcal{K}^{i}(y-w)]K(x,y,z)\,[U(x')-U(w')]\rho(x')\,U^{ah}(w)f^{hbc}\{\rho^{b}(x),\rho^{c}(y)\}$
the blue term of Row II is $+\log\frac{k^{+}}{\Lambda}X$, while the evolution files use $-\log\frac{\vee}{k^{+}}X$ instead. The cancellations above hold with the evolution files' version. The difference, $\log\frac{\vee}{\Lambda}X$, is not accounted for by Row II as written; its source has to be identified.
\item \textbf{Qualification 2, four $\rho$'s.} Eq.~(4rho) reduces to the three- and two-$\rho$ terms used in the files only after four Wilson lines in its second group are transposed ($U^{c'a}(w')\to U^{ac'}(w')$, $U^{bc'}(z)\to U^{c'b}(z)$, $U^{bd'}(z)\to U^{d'b}(z)$, $U^{bc}(z)\to U^{cb}(z)$). For $\log\frac{\vee}{k^{+}}$ the second two-$\rho$ block of \texttt{BFKL.tex} (line 69) was not reproduced.
\end{itemize}

\section{Finite terms}
\subsection{Definition and rules for the limit}
The finite terms are the black terms of \texttt{Stage1\_split}: everything that is not $\log\frac{k^{+}}{\Lambda}$ or $\log\frac{\vee}{k^{+}}$ times a coefficient, excluding the DGLAP and running-coupling terms of Group II Row IV, which are in Sec.~\ref{sec:dglap}. The limits are taken as follows.
\begin{center}
\begin{tabular}{p{0.52\textwidth}p{0.40\textwidth}}
\hline
Black piece in \texttt{Stage1\_split} & Limit $\Lambda\to0$, $\vee\to\infty$\\
\hline
$\log(1-\frac{k^{+}}{\vee})$, $\log(1+\frac{\Lambda}{k^{+}})$, $\log(1-\frac{\Lambda}{k^{+}})$ times a kernel & $0$\\
$\log\frac{(1-\frac{k^{+}}{\vee})A+\frac{k^{+}}{\vee}B}{B+\frac{\Lambda}{k^{+}}A}$, $A=(x-z)^{2}$, $B=(x-w)^{2}$ & $\log\frac{A}{B}$\\
$\log\frac{A-\frac{k^{+}}{\vee}B}{(1+\frac{\Lambda}{k^{+}})A-B}$, $\;\log\frac{\frac{k^{+}}{\vee}B-A}{B-\frac{\Lambda}{k^{+}}A}$, $\;\log\frac{\frac{k^{+}}{\vee}B-A}{B-(1+\frac{\Lambda}{k^{+}})A}$ & $\log\frac{A}{A-B}$, $\;\log\frac{A}{B}$, $\;\log\frac{A}{B-A}$ (principal value)\\
logs without large logarithm, arctan terms of Group I Row III & $\Lambda\to0$ inside, except the pair of Sec.~\ref{sec:cannot}\\
$[\frac{1}{\vee}-\frac{1}{k^{+}+\Lambda}]$ (Group II Row I, $2f$) & $-\frac{1}{k^{+}}$\\
plus-prescription integrals over $p^{+}\in(0,k^{+})$ (Group I Row V) & unchanged\\
integrals over $p^{+}\in(0,\vee-k^{+})$ (Group III Rows III, V, VI) & rule at the start of Sec.~\ref{sec:lim}, Group III\\
\hline
\end{tabular}
\end{center}

\subsection{Where the limit cannot be taken term by term}\label{sec:cannot}
\begin{enumerate}
\item \textbf{Group I Row III}, the term of $N_{3a}$ with $\frac{1}{(z-w)^{2}}$ and the last term of $N_{3b}$: only their sum with $\Lambda$ kept is integrable. The limit of the sum contains $\ln^{2}\frac{k^{+}}{\Lambda}$ and $\ln\frac{k^{+}}{\Lambda}$ pieces (remark in Sec.~\ref{sec:lim}).
\item \textbf{Group III Row V}, the $-\frac{1}{p^{+}}$ term: evaluated in closed form with $\Lambda$ kept, it is a pure squared logarithm, $\propto\ln^{2}\!\big[\frac{2e^{-\gamma_{E}}k^{+}}{\Lambda|k||x-w'|}\big]$. After adding the complex conjugate it cancels against item 1 completely, including the $\ln^{2}$ of the transverse scales. The two rows together leave a single term with the function $\mathcal{D}(k,b)=\int d^{2}v\,e^{-ik\cdot v}\frac{1}{v^{2}}\ln\frac{(b+v)^{2}}{b^{2}}$, which is listed in Group I Row III (red); nothing is listed in Group III Row V.
\item \textbf{Group III Rows III and VI and Group I Row V}, the plus-prescription terms that contain $U(z)$: the subtraction term and the plus-prescription integral are not separately convergent at large $z$. They are listed as they stand; they are defined only when the rows with the same soft gluon are added before the $z$ integration. This pairing is not verified.
\end{enumerate}

\subsection{Limiting finite terms, row by row}\label{sec:lim}
Changes with respect to your \texttt{Limiting\_Cases.tex} are in red, notes in teal.
"""
DGLAP=r"""
\section{DGLAP and running-coupling contributions}\label{sec:dglap}
Both come from Group II Row IV (\texttt{Group\_II.tex}, last formulas). With $\xi=\frac{k^{+}}{p^{+}+k^{+}}$ that row is, without approximation,
\[
\text{Row IV}=\text{BFKL piece}+\text{JIMWLK piece}+\text{DGLAP}+\text{RC}+\text{Finite}.
\]
\subsection{DGLAP}
\begin{align}
\frac{d^{3}N_{\rm DGLAP}}{d^{3}k}= &-\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{g^{2}}{8\pi^{2}}\frac{1}{k^{+}}\int_{0}^{1}\frac{d\xi}{\xi^{2}}\,\widetilde{\mathcal P}_{gg}(\xi)\,\int_{x,y,x',y'}e^{-ik\cdot(y'-y)\frac{1}{\xi}}\Big[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}}{4\,\xi^{2}}\Big]\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\notag\\
&\times\left[U^{ce'}(y') - U^{ce'}(x')\right]\left[U^{ce}(y) - U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x),\\
\widetilde{\mathcal{P}}_{gg}(\xi) =&\; 2N_{c}\left[\xi(1-\xi) + \frac{\xi}{(1-\xi)_{+}} + \Big[\frac{1-\xi}{\xi}\Big]_{+}\right]+\frac{\beta_{0}}{2}\,\delta(1-\xi),\qquad \beta_{0}=\frac{11N_{c}}{3}.\notag
\end{align}
\begin{itemize}
\item The coefficient of the splitting function agrees with \texttt{Group\_II.tex}. The collinear pole $\frac{2}{\epsilon}$ is kept here; \texttt{Final\_terms.tex} drops it, and it has to be said where it goes.
\item $\frac{\xi}{(1-\xi)_{+}}g\equiv\frac{\xi g(\xi)-g(1)}{1-\xi}$. Its subtraction term is exactly $\big[\log\frac{k^{+}}{\Lambda}+\log(1-\frac{k^{+}}{\vee})\big]g(1)$, which is the JIMWLK piece of this row; nothing is left over.
\item $\big[\frac{1-\xi}{\xi}\big]_{+}$ is not defined in \texttt{Final\_terms.tex}. With $\frac{1-\xi}{\xi}[g(\xi)-g(0)]$ the constant term of Sec.~\ref{sec:lim}, Group II Row IV, has to be added to the finite terms; with $\frac{(1-\xi)g(\xi)-g(0)}{\xi}$ it does not. Here $g(0)$ is the $\xi\to0$ limit after $y'-y=\xi u$, which gives the second term of the DGLAP piece of \texttt{Group\_II.tex}, proportional to $\frac{4\pi}{k^{2}}\int_{x,x',z}K(x,x',z)[U(z)-U(x')][U(z)-U(x)]$.
\item The lower limit of $\xi$ is $\frac{k^{+}}{\vee}\to0$.
\end{itemize}
\subsection{Running coupling}
\begin{align}
\frac{d^{3}N_{\rm RC}}{d^{3}k}= & +\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{g^{2}\beta_{0}}{16\pi^{2}}\frac{1}{k^{+}}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\Big[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}}{4}\Big]\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\notag\\
&\times\left[U^{ce'}(y') - U^{ce'}(x')\right]\left[U^{ce}(y) - U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x),
\end{align}
so that, with $\frac{g^{2}\beta_{0}}{16\pi^{2}}=\frac{\alpha_{s}\beta_{0}}{4\pi}$ and the pole removed,
\begin{align}
\frac{d^{3}N_{\rm LO}}{d^{3}k}+\frac{d^{3}N_{\rm RC}}{d^{3}k}&=\frac{1}{(2\pi)^{3}}\,\frac{g^{2}(\mu)}{2\pi^{2}}\,\frac{1}{k^{+}}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\left[1+\frac{\alpha_{s}(\mu)\,\beta_{0}}{4\pi}\log\frac{(y'-y)^{2}\,\mu^{2}_{\overline{MS}}\,e^{2\gamma_{E}}}{4}\right]\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}\,(x'-y')^{2}}\notag\\
&\times\left[U^{ce'}(y')-U^{ce'}(x')\right]\left[U^{ce}(y)-U^{ce}(x)\right]\rho^{e'}(x')\,\rho^{e}(x).
\end{align}
\textbf{What you should know about the $\beta_{0}$ terms.}
\begin{itemize}
\item The $\frac{\beta_{0}}{2}\delta(1-\xi)$ term of $\widetilde{\mathcal{P}}_{gg}$ is exactly minus the RC term (checked numerically, relative difference $3\times10^{-16}$). So in Group II Row IV the running-coupling term is added and subtracted: DGLAP $+$ RC contains no net $\beta_{0}$.
\item The only other $\beta_{0}$-type terms are virtual: the pole bracket times $\frac{11}{3}$ in Group I Row II and times $-\frac{11}{6}$ in Group I Row V ($2b\,{\rm UV}$). They cancel each other exactly.
\item Hence, as the files stand, the complete NLO result has no net $\beta_{0}\log\mu^{2}$ term: the $\mu$ dependence of ``LO $+$ RC'' is compensated by the $\delta(1-\xi)$ term of DGLAP. Whether this is intended (a choice of how to organise the collinear term) or whether a genuine running-coupling term is expected is a physics decision for you; the algebra of the files is consistent either way.
\end{itemize}
"""
doc=pre+r"\begin{document}"+MAIN1+body+DGLAP+"\n\\end{document}\n"
open("/home/user/correlations/nlo_review/stage2/Stage2_finite_terms.tex","w").write(doc)
print(len(doc))
