import json,sys; sys.path.insert(0,'.')
from lib_split import *
res=json.load(open("split_blocks.json"))
pre=open(UP+FILES["I"]).read().split(r"\begin{document}")[0]
pre=pre.replace(r"\title{", r"\title{Stage 1: BFKL and JIMWLK check, and what is left in each row\\[2pt]\large ").replace(r"\maketitle","")
import re
pre=re.sub(r"\\title\{[^\n]*\n", r"\\title{Stage 1: check of the BFKL and JIMWLK terms, and what is left in each row}"+"\n", pre, count=1)
MAC=r"""\newcommand{\chg}[1]{{\color{red}#1}}
\newcommand{\rnote}[1]{{\color{teal}[\textit{#1}]}}
"""
INTRO=r"""
\maketitle
\section{Result of the check}
\begin{enumerate}
\item \textbf{BFKL term.} The colour-singlet form is, term by term and with the same normalisation, the BFKL equation of your main file (Eq.~(v-bfkl) of \texttt{Calculation\_main.tex}) applied to the leading-order cross section, i.e.\ Eq.~(LOBFKL) with $\delta Y=\log\frac{\vee}{k^{+}}$ and $\alpha_{s}=g^{2}/4\pi$. The explicit form reduces to the singlet form under $\rho^{b}(x')\rho^{c}(x)\to\frac{\delta^{bc}}{N_{c}^{2}-1}\langle\rho^{f}(x')\rho^{f}(x)\rangle$ (all eight blocks checked). Eq.~(v-bfkl) itself is the standard BFKL equation: in momentum space its real part is $\frac{\alpha_{s}N_{c}}{\pi^{2}}\frac{q^{2}}{p^{2}(q-p)^{2}}$.
\item \textbf{JIMWLK term.} I acted with $H_{\rm JIMWLK}$ of Eq.~(HLO) on ${\rm Tr}\{(U^{\dagger}(x)-U^{\dagger}(w))(U(y)-U(w'))\}$ numerically (random adjoint Wilson lines, $SU(2)$ and $SU(3)$), using only the commutators $[J^{a}_{R}(z),U(x)]=\delta(x-z)U(x)T^{a}$. For each of the eight kernels $K(w,w',z)$, $K(x,w',z)$, $K(x,y,z)$, $K(y,w,z)$, $K(w',w',z)$, $K(w,w,z)$, $K(x,x,z)$, $K(y,y,z)$ the coefficient of $-H\,\mathcal{O}$ equals the coefficient in your singlet form and in your explicit form to $10^{-8}$; the kernels $K(w,x,z)$ and $K(w',y,z)$ have coefficient zero, as in your result. Sign and normalisation agree with $\delta Y=\log\frac{k^{+}}{\Lambda}$.
\item Only a notation point: after the JIMWLK formula the kernel is defined as $\mathcal{K}(x,y,z)$ but used as $K(x,y,z)$.
\end{enumerate}
So both expressions are the correct evolution of the leading-order cross section. What this check does \emph{not} cover is stated in Sec.~\ref{sec:open}.

\section{How to read the rest of this file}
For every row I took the last formula of \texttt{Group\_I.tex}, \texttt{Group\_II.tex}, \texttt{Group\_III.tex} and rewrote each logarithm \emph{exactly} as
(large logarithms, in {\color{blue}blue}) $+$ (remainder, in black). Nothing is expanded and no limit is taken: $\Lambda$ and $\vee$ are kept. The rewrites are
\begin{align}
\log\Big(\frac{\vee}{k^{+}}-1\Big)+\log\frac{k^{+}}{\Lambda}&={\color{blue}\log\frac{k^{+}}{\Lambda}+\log\frac{\vee}{k^{+}}}+\log\Big(1-\frac{k^{+}}{\vee}\Big),&
\log\frac{\vee}{k^{+}+\Lambda}&={\color{blue}\log\frac{\vee}{k^{+}}}-\log\Big(1+\frac{\Lambda}{k^{+}}\Big),\notag\\
\log\frac{k^{+}-\Lambda}{\Lambda}&={\color{blue}\log\frac{k^{+}}{\Lambda}}+\log\Big(1-\frac{\Lambda}{k^{+}}\Big),&
\log\frac{(\vee-k^{+})A+k^{+}B}{\Lambda A+k^{+}B}&={\color{blue}\log\frac{\vee}{k^{+}}}+\log\frac{\big(1-\frac{k^{+}}{\vee}\big)A+\frac{k^{+}}{\vee}B}{B+\frac{\Lambda}{k^{+}}A},\notag\\
\log\frac{\vee A-k^{+}B}{(\Lambda+k^{+})A-k^{+}B}&={\color{blue}\log\frac{\vee}{k^{+}}}+\log\frac{A-\frac{k^{+}}{\vee}B}{\big(1+\frac{\Lambda}{k^{+}}\big)A-B},&
\log\frac{k^{+}B-\vee A}{k^{+}B-\Lambda A}&={\color{blue}\log\frac{\vee}{k^{+}}}+\log\frac{\frac{k^{+}}{\vee}B-A}{B-\frac{\Lambda}{k^{+}}A},
\end{align}
with $A=(x-z)^{2}$, $B=(x-w)^{2}$ (or primed), and similarly for the other orderings. Logarithms that have no large logarithm at fixed transverse coordinates are left untouched, for example $\log\frac{[k^{+}(x-w)-(k^{+}-\Lambda)(x-z)]^{2}}{[k^{+}(x-w)-\Lambda(x-z)]^{2}}$.

\textbf{The blue terms are what goes into the BFKL and JIMWLK terms (together with the $4\rho$ terms of \texttt{4rho.tex}); everything in black is what is left.} Terms that carry only a blue factor are kept in place so that the bookkeeping can be followed. In rows written with the plus prescription, the black part is the plus-prescription integral and the blue part is the subtraction term.
"""
OPEN=r"""
\section{What is not covered, and points to watch}\label{sec:open}
\begin{enumerate}
\item \textbf{Sum of the blue terms.} I checked that your BFKL and JIMWLK expressions are the correct evolution of the leading order. I have not yet checked that the blue terms of all rows, plus the $4\rho$ terms, add up to them, i.e.\ the cancellation of the $4\rho$ and $3\rho$ logarithms and the $2\rho$ terms produced by reordering (\texttt{4rho.tex}). That is the next check.
\item \textbf{$p^{+}$ integrations.} Verified here by hand for Group I Row IV, Group II Row I (first part), Group III Rows I and II (they share one integral). Group I Row III was verified numerically before. Group II Row I, second part ($2a$--$2f$), is not rechecked yet.
\item \textbf{The limit $\Lambda\to0$ is not uniform in the soft sector.} A blue $\log\frac{k^{+}}{\Lambda}$ may only be split off if its coefficient is integrable. This fails in two places:
 \begin{itemize}
 \item Group I Row III: the term of $N_{3a}$ with kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}(z-w)^{2}}$ is not integrable at $z=w$; it must be kept together with the last term of $N_{3b}$ (black, $\Lambda$ kept). Their sum is $-\frac{ig^{4}f^{abd}}{16\pi^{5}k^{+}}\frac{1}{(z-w)^{2}}\log\frac{(k^{+}-\Lambda)^{2}[k^{+}(z-w)+\Lambda(x-z)]^{2}}{\Lambda^{2}[k^{+}(x-w)-\Lambda(x-z)]^{2}}$ times the same structure, and its $z$ integral contains $\ln^{2}\frac{k^{+}}{\Lambda}$.
 \item Group III Rows III (first part), V, VI: in the term $-\big[\frac{1}{p^{+}}\big]_{+}\delta_{jm}\delta_{ik'}$ the subtraction term (blue) has the kernel $\frac{(y'-z)\cdot(x-z)}{(y'-z)^{2}(x-z)^{2}}$, which is $1/z^{2}$ at large $z$. In Row V nothing else depends on $z$, so the blue coefficient and the plus-prescription integral are separately divergent; the unsplit term is $-\frac{ig^{4}f^{c'ca}}{8\pi^{4}k^{+}}\ln^{2}\!\big[\frac{2e^{-\gamma_{E}}k^{+}}{\Lambda|k||x-y'|}\big]$ times the colour structure (for $\vee\gg k^{+}$).
 \end{itemize}
 The $\ln^{2}$ of Group I Row III and the $\ln^{2}$ of the complex conjugate of Group III Row V cancel exactly, so there is no $\ln^{2}\frac{k^{+}}{\Lambda}$ in the sum; but these two rows cannot be reduced to ``blue $+$ finite'' separately.
\item \textbf{Lower limit.} \texttt{Group\_III.tex} confirms that the $p^{+}$ integrals of Group III start at $\Lambda$, which I had assumed.
\item \textbf{Group II Row I, first part.} Two factors $\log\frac{\vee}{\Lambda+k^{+}}$ are not coloured blue in \texttt{Group\_II.tex}; they do contain $\log\frac{\vee}{k^{+}}$ and are split here.
\item \chg{\textbf{Correction (added after comparing with \texttt{Other\_terms.tex}).}} \chg{In the rows that are left unintegrated in $p^{+}$ (Group III Rows III, V, VI) the black part still contains a $\log\frac{\vee}{k^{+}}$: it sits in the terms $-\frac{1}{k^{+}}\delta_{im}\delta_{jk'}$ (first part of Row III, Rows V, VI) and $-\frac{p^{+}}{k^{+2}}[\ldots]$ (second part of Row III), whose $p^{+}$ integrand falls only like $1/p^{+}$. Your \texttt{BFKL.tex} and \texttt{Other\_terms.tex} take it out by subtracting the same term with $\mathcal{K}(x'-y')[U(y')\ldots]\to\mathcal{K}(x'-z)[U(z)\ldots]$. I checked this: after $y'\to w'$ the subtracted term is exactly $\log\frac{\vee}{k^{+}}$ times the BFKL entry for the $-\frac{1}{k^{+}}$ terms, and for the $-\frac{p^{+}}{k^{+2}}$ term it is $\log\frac{\vee}{k^{+}}$ times the BFKL entry plus the remainder with $\frac{(x-w)^{2}}{(x-z)^{2}(p^{+}+k^{+})[p^{+}(x-z)^{2}+k^{+}(x-w)^{2}]}$ that you have in \texttt{Other\_terms.tex}.}
\item \textbf{Typos in the Group files (no effect on the results).} Group III Row II, first line: $\bar{\mathbb{A}}^{\dagger(1)a}_{i}(k^{+},-U(w')^{c'b}(z)$ should read $\bar{\mathbb{A}}^{\dagger(1)a}_{i}(k^{+},-w')\,U^{c'b}(z)$. Group I Row II, $N_{1}$: the product $\log\big[1-\frac{k^{+}}{\vee}\big]\log\frac{\vee}{k^{+}}$ vanishes for $\vee\to\infty$ and is not an evolution logarithm; it is left black.
\end{enumerate}
"""
NOTES={
 "I|Row II": r"\rnote{$N_{2a}$ is entirely blue. $N_{2b}$ (with $R^{i}$) is entirely left. In $N_{1}$ the blue piece multiplies the pole bracket $-\frac{2}{\epsilon}-2\gamma-\log\frac{(x-w)^{2}\mu^{2}}{4}$.}",
 "I|Row III": r"\rnote{See Sec.~\ref{sec:open}, item 3, for the term of $N_{3a}$ with $1/(z-w)^{2}$ and the last term of $N_{3b}$. The two red terms of $N_{2}$ cancel each other identically.}",
 "I|Row V, first part": r"\rnote{Black: the two plus-prescription integrals. Blue: their subtraction terms.}",
 "I|Row V, second part": r"\rnote{In $2b\,{\rm UV}$ the blue $2\log\frac{k^{+}}{\Lambda}$ multiplies the pole bracket; the $-\frac{11}{6}$ is left.}",
 "II|Row IV": r"\rnote{DGLAP, RC and Finite are entirely left. From the JIMWLK piece $\log\big(1-\frac{k^{+}}{\vee}\big)$ is left, from the BFKL piece $-\big(1-\frac{k^{+}}{\vee}\big)$ is left.}",
 "III|Row V": r"\rnote{See Sec.~\ref{sec:open}, item 3.}",
 "III|Row III, first part": r"\rnote{See Sec.~\ref{sec:open}, item 3.}",
 "III|Row VI": r"\rnote{See Sec.~\ref{sec:open}, item 3.}",
}
body=[]
cur=None
names={"I":"Group I (with c.c.)","II":"Group II","III":"Group III (with c.c.)"}
only4={"I":r"Row I: only $4\rho$.","II":r"Rows II and III: only $4\rho$.","III":r"Row IV: only $4\rho$."}
for g,title,spans in BLOCKS:
    if g!=cur:
        body.append(r"\section{"+names[g]+"}\n"+only4[g]+"\n"); cur=g
    body.append(r"\subsection{"+title+"}\n")
    if f"{g}|{title}" in NOTES: body.append(NOTES[f"{g}|{title}"]+"\n")
    for blk in res[f"{g}|{title}"]: body.append(blk+"\n")
doc=pre+MAC+r"\begin{document}"+INTRO+"\n".join(body)+OPEN+"\n\\end{document}\n"
open("/home/user/correlations/nlo_review/stage1/Stage1_split.tex","w").write(doc)
print(len(doc))
