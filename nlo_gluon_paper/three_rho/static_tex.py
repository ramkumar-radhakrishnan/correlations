PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.55in]{geometry}
\usepackage{amsmath,amssymb,amsfonts,mathtools}
\usepackage[dvipsnames]{xcolor}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue]{hyperref}
\allowdisplaybreaks
\numberwithin{equation}{section}
\setlength{\jot}{3pt}
\newcommand{\KK}{\mathcal K}
\newcommand{\NN}{\mathcal N}
\def\mf#1{\mathbb{#1}}
\setcounter{tocdepth}{2}

\title{Two-$\rho$ terms from the three-$\rho$ terms proportional to $\log(\vee/\Lambda)$}
\date{}

\begin{document}
\maketitle
\tableofcontents
"""

INTRO = r"""
\section{What is done here}\label{sec:intro}

The input consists of all three-$\rho$ terms of the NLO cross section that are proportional to $\log(\vee/\Lambda)$:
\begin{enumerate}
\item the remainders left after reordering the four-$\rho$ terms, i.e.\ the three-$\rho$ terms with two and with three Wilson lines of the notes (the remainder with four Wilson lines vanishes, as shown there);
\item the genuine three-$\rho$ terms of Groups I, II and III (file \texttt{Evolution\_three\_rho.tex}). Groups I and III come with their complex conjugate, Group II does not.
\end{enumerate}
In every term the three charges appear as an ordered product of operators. We write each product as its fully symmetrized part plus commutators, Eqs.~\eqref{eq:plain}--\eqref{eq:antiA} below. Each commutator removes one $\rho$ and produces a $\delta$-function. Those terms are the two-$\rho$ terms we are after. They are derived term by term in Secs.~\ref{sec:fourrho}--\ref{sec:G3} and collected in Sec.~\ref{sec:summary}.

Nothing else is changed: the kernels are those of the notes, evaluated at the points set equal by the $\delta$-functions. The two exceptions are a repeated colour index in Row II of Group I and the transformation of that row's kernel to coordinate space (Sec.~\ref{sec:GI.II} and App.~\ref{app:F}). The algebra was generated and checked by computer (App.~\ref{app:checks}).
"""

CONVENTIONS = r"""
\section{Conventions}\label{sec:conv}

\paragraph{Common prefactor.} Every term carries $\log(\vee/\Lambda)$ and a factor $g^4/(16\pi^5)=\frac{g^4}{8\pi^3}\cdot\frac{1}{2\pi^2}$. We define
\begin{equation}
\NN\equiv\frac{1}{(2\pi)^3}\,\frac{g^4}{16\pi^5}\,\frac{1}{k^+}\,\log\frac{\vee}{\Lambda},
\label{eq:NN}
\end{equation}
so that $\frac{1}{(2\pi)^3}\frac{ig^4}{16\pi^5}\frac{1}{k^+}\log\frac\vee\Lambda=i\NN$, $\frac{1}{(2\pi)^3}\frac{ig^4}{8\pi^5}\frac{1}{k^+}\log\frac\vee\Lambda=2i\NN$ and $\frac{1}{(2\pi)^3}\frac{ig^4}{32\pi^5}\frac{1}{k^+}\log\frac\vee\Lambda=\frac i2\NN$. All transverse positions except the two measured-gluon positions in the phase are integrated; we write $\int$ for all these integrals.

\paragraph{Kernels.} With $\KK(X)\equiv X/X^2$ (a two-vector),
\begin{equation}
\frac{(a-b)\cdot(c-d)}{(a-b)^2(c-d)^2}=\KK(a-b)\cdot\KK(c-d),\qquad \KK(X)\cdot\KK(X)=\frac{1}{X^2},\qquad \KK(-X)=-\KK(X).
\end{equation}
Each source has a kernel that depends on the three charge positions, e.g.\ $K_1(x',x,y)$. A $\delta$-function sets two of these equal, e.g.\ $K_1(x',x',y)$ means $K_1$ at $x=x'$.

\paragraph{Charges.} The valence charges obey
\begin{equation}
\big[\rho^a(x),\rho^b(y)\big]=if^{abc}\rho^c(x)\,\delta^{(2)}(x-y).
\label{eq:comm}
\end{equation}
We always integrate the $\delta$-function over the position of the \emph{right} operator. In $[\rho^a(x),\rho^b(y)]$ the point $y$ is replaced by $x$ everywhere: in the kernel, in the Wilson lines and in the remaining charge.

\paragraph{Colour.} $(T^a)_{bc}=-if^{abc}$, $f^{acd}f^{bcd}=N_c\delta^{ab}$. For adjoint Wilson lines we use
\begin{align}
&\text{(F)}\quad f^{acd}f^{bcd}=N_c\,\delta^{ab},\\
&\text{(U1)}\quad U^{ab}(x)\,U^{ac}(x)=\delta^{bc},\label{eq:U1}\\
&\text{(U2)}\quad f^{abc}\,U^{ad}(x)\,U^{be}(x)=f^{deh}\,U^{ch}(x).\label{eq:U2}
\end{align}
(U2) follows from $U^{aa'}U^{bb'}U^{cc'}f^{a'b'c'}=f^{abc}$ after multiplying by $U^{cc''}$ and using (U1) (written here for the transposed matrix). As written, (U1) and (U2) hold for any adjoint Wilson line: the contracted indices sit in the same slot of both $U$'s. Your notes use the convention
\begin{equation}
U^{\dagger\,bc}=U^{bc}=U^{cb},
\label{eq:Usym}
\end{equation}
and we follow it. With Eq.~\eqref{eq:Usym}, (U1) and (U2) can also be used when the contracted indices sit in different slots, e.g.\ $U^{ab}(x)U^{ca}(x)=\delta^{bc}$. Such steps are marked (U1$^\star$) and (U2$^\star$). They are exact only for a symmetric $U$; for a general adjoint matrix $U^{ab}U^{ca}=(U^2)^{cb}$. Finally, $f^{abc}$ contracted with an object symmetric in $a,b$ vanishes.
"""

IDENTITY = r"""
\section{The reordering identity}\label{sec:identity}

Let $A,B,C$ be three charges, e.g.\ $A=\rho^d(x')$, $B=\rho^a(x)$, $C=\rho^c(y)$. Their fully symmetrized product is
\begin{equation}
S(ABC)\equiv\frac16\big(ABC+ACB+BAC+BCA+CAB+CBA\big).
\end{equation}
\paragraph{Step 1: every ordering in terms of $ABC$.} Moving operators one step at a time,
\begin{align}
ACB&=ABC-A[B,C],\notag\\
BAC&=ABC-[A,B]C,\notag\\
BCA&=BAC-B[A,C]=ABC-[A,B]C-B[A,C],\notag\\
CAB&=ACB-[A,C]B=ABC-A[B,C]-[A,C]B,\notag\\
CBA&=BCA-[B,C]A=ABC-[A,B]C-B[A,C]-[B,C]A.
\end{align}
Adding the six orderings,
\begin{equation}
6\,S(ABC)=6\,ABC-3[A,B]C-2A[B,C]-2B[A,C]-[A,C]B-[B,C]A .
\label{eq:sixS}
\end{equation}
\paragraph{Step 2: symmetrize the products of two.} Each commutator is linear in $\rho$ (Eq.~\eqref{eq:comm}). For a commutator $X$ and a charge $Y$, $XY=\frac12\{X,Y\}+\frac12[X,Y]$ and $YX=\frac12\{X,Y\}-\frac12[X,Y]$. Here $[X,Y]$ is a double commutator, i.e.\ a single $\rho$. Inserting this in Eq.~\eqref{eq:sixS}, the anticommutators come with weights $3$ for $\{[A,B],C\}$, $2+1=3$ for $\{[A,C],B\}$ and $2+1=3$ for $\{[B,C],A\}$. Hence
\begin{equation}
\boxed{\;ABC=S(ABC)+\tfrac14\Big(\big\{[A,B],C\big\}+\big\{[A,C],B\big\}+\big\{[B,C],A\big\}\Big)+R_1\;}
\label{eq:plain}
\end{equation}
with the one-$\rho$ remainder
\begin{equation}
R_1=\tfrac1{12}\Big(3\big[[A,B],C\big]+2\big[A,[B,C]\big]+2\big[B,[A,C]\big]+\big[[A,C],B\big]+\big[[B,C],A\big]\Big).
\end{equation}
\paragraph{Step 3: drop the one-$\rho$ terms.} Each term of $R_1$ is a double commutator, $\propto f f\,\rho\,\delta\delta$. It is linear in $\rho$ and multiplies c-numbers (kernels and Wilson lines). Its expectation value vanishes for a colour-neutral projectile, $\langle\rho^a(x)\rangle=0$. We drop it.
\paragraph{Anticommutator forms.} Several rows of the notes contain $\{\rho,\rho\}$. Applying Eq.~\eqref{eq:plain} to $ABC$ and $ACB$ and adding, the terms $\{[B,C],A\}$ and $\{[C,B],A\}$ cancel, and $S(ACB)=S(ABC)$:
\begin{equation}
\boxed{\;A\{B,C\}=2S(ABC)+\tfrac12\big\{[A,B],C\big\}+\tfrac12\big\{[A,C],B\big\}+(\text{one }\rho)\;}
\label{eq:Aanti}
\end{equation}
When the single charge stands on the right, use $\{B,C\}A=A\{B,C\}-[A,\{B,C\}]$ and $[A,\{B,C\}]=\{[A,B],C\}+\{B,[A,C]\}$:
\begin{equation}
\boxed{\;\{B,C\}A=2S(ABC)-\tfrac12\big\{[A,B],C\big\}-\tfrac12\big\{[A,C],B\big\}+(\text{one }\rho)\;}
\label{eq:antiA}
\end{equation}
\paragraph{One piece.} For $X=\rho^x(p_X)$, $Y=\rho^y(p_Y)$, $Z=\rho^z(p_Z)$, Eq.~\eqref{eq:comm} gives
\begin{equation}
\big\{[X,Y],Z\big\}=i f^{xyh}\big\{\rho^h(p_X),\rho^z(p_Z)\big\}\,\delta^{(2)}(p_X-p_Y),
\end{equation}
and $p_Y\to p_X$ everywhere. The two-$\rho$ terms are written with the anticommutator $\{\rho,\rho\}$. A different ordering of two charges differs from it by one $\rho$, which we drop.

\paragraph{The symmetric three-$\rho$ part.} For every term, the part $S(\dots)$ is the same expression with $\rho\rho\rho$ replaced by $S(\rho\rho\rho)$ (and by $2S(\rho\rho\rho)$ for the anticommutator forms). These parts behave as products of commuting charges. We do not list them again, and we do not check their cancellation in this note.

\paragraph{Complex conjugate.} All prefactors below are $\pm i$, $\pm i/2$ or $\pm2i$ times $\NN$. The commutator supplies one more $i$, so every two-$\rho$ term has a \emph{real} coefficient. The kernels, $f$, $U$ and $\delta$ are real, and $\{\rho,\rho\}$ is Hermitian. The complex conjugate of a two-$\rho$ term therefore only turns $e^{-ik\cdot(P'-P)}$ into $e^{ik\cdot(P'-P)}$, where $P'$ and $P$ are the measured-gluon positions in the conjugate amplitude and the amplitude. Renaming $P\leftrightarrow P'$ gives back the original phase:
\begin{equation}
\text{c.c.}=\big(\text{same term with }P\leftrightarrow P'\big).
\end{equation}
This holds for Groups I and III, with $(P',P)=(w',w)$, except $(w',z)$ and $(w',x)$ in Row V of Group I and $(y',w)$ in Rows III, V and VI of Group III. Group II and the four-$\rho$ remainders get no c.c.
"""

APPF = r"""
\appendix
\section{The kernel of Row II of Group I in coordinate space}\label{app:F}

We need
\begin{equation}
F^i(w;x,y)\equiv\int d^2\mathbf k\,d^2\mathbf p\;e^{i\mathbf k\cdot(w-y)}e^{i\mathbf p\cdot(y-x)}\,\frac{B^i(\mathbf k,\mathbf p)}{\mathbf k^2\mathbf p^2},\qquad
B^i=-\frac{2\mathbf k^2\mathbf p^i+(\mathbf p^2-\mathbf k\cdot\mathbf p)\mathbf k^i}{(\mathbf k-\mathbf p)^2}-2\mathbf k^i .
\end{equation}
($B^i$ is the bracket of the notes: $-(\mathbf k^2\mathbf p^i-\mathbf k\cdot\mathbf p\,\mathbf k^i)-(\mathbf p^2\mathbf k^i+\mathbf k^2\mathbf p^i)=-[2\mathbf k^2\mathbf p^i+(\mathbf p^2-\mathbf k\cdot\mathbf p)\mathbf k^i]$.)

\paragraph{Step 1: partial fractions.} With $\mathbf q\equiv\mathbf k-\mathbf p$ and $\mathbf p^2-\mathbf k\cdot\mathbf p=-\mathbf p\cdot\mathbf q$,
\begin{equation}
\frac{B^i}{\mathbf k^2\mathbf p^2}=-\frac{2\mathbf p^i}{\mathbf q^2\mathbf p^2}+\frac{(\mathbf p\cdot\mathbf q)\,\mathbf k^i}{\mathbf q^2\mathbf k^2\mathbf p^2}-\frac{2\mathbf k^i}{\mathbf k^2\mathbf p^2}.
\end{equation}
\paragraph{Step 2: Fourier transforms.} We use $\int d^2\mathbf k\,e^{i\mathbf k\cdot X}\mathbf k^i/\mathbf k^2=2\pi i\,\KK^i(X)$; $\int d^2\mathbf q\,e^{i\mathbf q\cdot X}/\mathbf q^2=-\pi\log X^2+C_0$, where the constant $C_0$ is infrared divergent; and $\mathbf k^i/\mathbf k^2=\frac{i}{2\pi}\int d^2u\,e^{-i\mathbf k\cdot u}\KK^i(u)$. In the first term we change variables $(\mathbf k,\mathbf p)\to(\mathbf q,\mathbf p)$, so the phase is $e^{i\mathbf q\cdot(w-y)}e^{i\mathbf p\cdot(w-x)}$:
\begin{equation}
-2\cdot2\pi i\,\KK^i(w-x)\big(-\pi\log(w-y)^2+C_0\big)=-4\pi^2 i\,\KK^i(x-w)\log(w-y)^2+4\pi iC_0\KK^i(x-w).
\end{equation}
The third term factorizes directly:
\begin{equation}
-2\cdot2\pi i\,\KK^i(w-y)\big(-\pi\log(x-y)^2+C_0\big)=-4\pi^2 i\,\KK^i(y-w)\log(x-y)^2+4\pi iC_0\KK^i(y-w).
\end{equation}
In the second term we write $\mathbf k^i/\mathbf k^2$ as an integral over $u$, change variables to $(\mathbf q,\mathbf p)$, and do the two Fourier transforms. With $z'=w-u$,
\begin{equation}
\frac{i}{2\pi}(2\pi i)^2\int_{z'}\KK^i(w-z')\,\KK(z'-y)\cdot\KK(z'-x)=-2\pi i\int_{z'}K(x,y,z')\,\KK^i(w-z').
\end{equation}
\paragraph{Step 3: keep the odd part.} The colour structures of Row II are odd under $x\leftrightarrow y$, so only the odd part of $F^i$ contributes. The $z'$ term and the $C_0$ terms are even and drop out. With $F^i=2\pi i\,\mathfrak F^i$,
\begin{equation}
\mathfrak F^i(w;x,y)=-\pi\Big[\KK^i(x-w)\log\frac{(w-y)^2}{(x-y)^2}-\KK^i(y-w)\log\frac{(w-x)^2}{(x-y)^2}\Big].
\end{equation}
Using $\int_{z'}K(a,b,z')=-\pi\log(a-b)^2+\text{const}$, this is
\begin{equation}
\mathfrak F^i(w;x,y)=\KK^i(x-w)\int_{z'}\big[K(w,y,z')-K(x,y,z')\big]-\KK^i(y-w)\int_{z'}\big[K(w,x,z')-K(x,y,z')\big].
\end{equation}
The coordinate-space expression for $\mathbb A^{(3)}$ in your \texttt{main.tex} (the commented version), after $\log\frac{\Lambda}{k^+}=-\log\frac{\vee}{\Lambda}+\log\frac{\vee}{k^+}$, has the same structure: $\KK^i(x-w)\int_{z'}[K(w,y,z')-K(x,y,z')]-(x\leftrightarrow y)$.
"""
