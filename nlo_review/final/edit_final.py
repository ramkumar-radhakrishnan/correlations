import sys
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
OUT="/home/user/correlations/nlo_review/final/"
MAC=r"""\newcommand{\chg}[1]{{\color{red}#1}}
\newcommand{\rnote}[1]{{\color{teal}[\textit{#1}]}}
"""
def rep(s,old,new,n=1):
    c=s.count(old)
    assert c==n, (c,n,old[:90])
    return s.replace(old,new)

# ---------- common new pieces ----------
UU=r"\left[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\, \rho^{b'}(x')\right]"
# new last term of N_3b (inside the curly bracket of Other_terms / Limiting_Cases)
def new_last(pref):
    return (pref+r"""\chg{-\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\int_{x,y,z,x'}\frac{1}{(z-w)^{2}}\,\log{\left[\frac{(z-w)^{2}}{(x-w)^{2}}\right]}\,\frac{(x'-w')^{i}}{(x'-w')^{2}}\left[U^{ab'}(x')-U^{ab'}(w')\right]\rho^{b'}(x')}\notag \\ & \chg{\qquad\times U^{de}(x)\,\rho^{e}(x)\left\{\frac{(y-z)^{i}}{(y-z)^{2}}\left[U^{bc}(z)-U^{bc}(y)\right]-\frac{(x-w)^{2}}{(x-w)^{2}+(z-w)^{2}}\,\frac{(y-w)^{i}}{(y-w)^{2}}\left[U^{bc}(w)-U^{bc}(y)\right]\right\}\rho^{c}(y)}\notag \\ & \chg{-\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\frac{\pi^{3}}{6}\,\int_{x,y,x'}\frac{(x'-w')\cdot(y-w)}{(x'-w')^{2}\,(y-w)^{2}}\left[U^{ab'}(x')-U^{ab'}(w')\right]\rho^{b'}(x')}\notag \\ & \chg{\qquad\times U^{de}(x)\,\rho^{e}(x)\left[U^{bc}(w)-U^{bc}(y)\right]\rho^{c}(y)}""")

REMARK_III=r"""
{\color{red}
\paragraph{Remark on the last term of $N_{3b}$ (changed).}
In \texttt{Group\_I.tex} two terms multiply the kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}(z-w)^{2}}$: the second term of $N_{3a}$, with $\log\frac{k^{+}-\Lambda}{\Lambda}$, and the last term of $N_{3b}$. Each of them is $\propto 1/(z-w)^{2}$ and is not integrable at $z=w$ once $\Lambda\to0$ is taken; only their sum is. With $\Lambda$ kept the sum is exactly
\begin{align}
\mathcal{P}_{\Lambda}=-\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{z}\frac{1}{(z-w)^{2}}\,\log\frac{(k^{+}-\Lambda)^{2}\left[k^{+}(z-w)+\Lambda\,(x-z)\right]^{2}}{\Lambda^{2}\left[k^{+}(x-w)-\Lambda\,(x-z)\right]^{2}}\;\mathcal{F}^{i}(z),\qquad
\mathcal{F}^{i}(z)\equiv\frac{(y-z)^{i}}{(y-z)^{2}}\left[U^{bc}(z)-U^{bc}(y)\right]\rho^{c}(y),
\end{align}
(times $\frac{(x'-w')^{i}}{(x'-w')^{2}}[U^{ab'}(x')-U^{ab'}(w')]\rho^{b'}(x')\,U^{de}(x)\rho^{e}(x)$). The limit $\Lambda\to0$ is not uniform in $z$ (the logarithm changes behaviour at $|z-w|\sim\Lambda|x-w|/k^{+}$). Doing the $z$ integral first gives, up to $O(\Lambda)$ (checked numerically, agreement $10^{-5}$ at $\Lambda/k^{+}=10^{-6}$),
\begin{align}
\int_{z}\frac{1}{(z-w)^{2}}\log\frac{\left[k^{+}(z-w)+\Lambda(x-z)\right]^{2}}{\Lambda^{2}(x-w)^{2}}\,\mathcal{F}^{i}(z)
&=2\pi\,\ln^{2}\!\frac{k^{+}}{\Lambda}\;\mathcal{F}^{i}(w)
+2\ln\frac{k^{+}}{\Lambda}\int_{z}\frac{1}{(z-w)^{2}}\left[\mathcal{F}^{i}(z)-\frac{(x-w)^{2}}{(x-w)^{2}+(z-w)^{2}}\,\mathcal{F}^{i}(w)\right]\notag\\
&+\int_{z}\frac{1}{(z-w)^{2}}\ln\frac{(z-w)^{2}}{(x-w)^{2}}\left[\mathcal{F}^{i}(z)-\frac{(x-w)^{2}}{(x-w)^{2}+(z-w)^{2}}\,\mathcal{F}^{i}(w)\right]+\frac{\pi^{3}}{6}\,\mathcal{F}^{i}(w).
\end{align}
Under $\int_{w}e^{ik\cdot w}$ the single-log coefficient is exactly $\int_{z}\frac{1}{(z-w)^{2}}\big[\mathcal{F}^{i}(z)-\frac{(x-w)^{2}}{(x-w)^{2}+(z-w)^{2}}\mathcal{F}^{i}(w)\big]\to 2\pi\ln\frac{2e^{-\gamma_{E}}}{|k|\,|x-w|}\,\mathcal{F}^{i}(w)$, so the $\Lambda$ dependent part of this pair is $-\frac{ig^{4}f^{abd}}{8\pi^{4}k^{+}}\big[\ln^{2}\frac{k^{+}}{\Lambda}+2\ln\frac{k^{+}}{\Lambda}\ln\frac{2e^{-\gamma_{E}}}{|k||x-w|}\big]\mathcal{F}^{i}(w)$. It cancels exactly against the complex conjugate of the $-1/p^{+}$ term of Group III, Row V (see the remark there); the eikonal $\log\frac{k^{+}}{\Lambda}$ entries of these two rows in the evolution list ($\beta$ of Group I Row III, $\beta'$ of Group III Row V) are the images of these terms and cancel in the same way. What is left in this row are the last two terms of the expansion, which are the last three lines of $N_{3b}$ above. The previous form of this term, $\log\frac{(z-w)^{2}}{(x-w)^{2}}\,\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot(w-z)}{(x'-w')^{2}(y-w)^{2}(x-z)^{2}(w-z)^{2}}$, is not equal to it, and its $z$ integral diverges like $\ln^{2}$ at large $z$.
}
"""

def rowV_closed():
    return (r""" \notag \\
& \chg{-\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'ca}}{8\pi^{4}}\, \frac{1}{k^{+}}\,\int_{w,w',y,x,x'}\, e^{-ik\cdot (w'-w)}\,\ln^{2}\!\left[\frac{2\,e^{-\gamma_{E}}}{|k|\,|x-w'|}\right]\,\frac{(x'-w')\cdot (y-w)}{(x'-w')^{2}\,(y-w)^{2}}}\notag \\ & \chg{\qquad\times\left[U^{c'e'}(w')-U^{c'e'}(x')\right]\rho^{e'}(x')\,U^{ce}(x)\left[U^{ad}(y)-U^{ad}(w)\right]\rho^{e}(x)\rho^{d}(y)}""")

REMARK_V=r"""
{\color{red}
\paragraph{Remark on Row V (changed).}
The term $-\Big(\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}+\frac{1}{p^{+}+k^{+}}\Big)\delta_{jm}\delta_{ik'}$ has been taken out of the $p^{+}$ integral and replaced by the last two lines. Reason: before the plus prescription this term is $-\frac{1}{p^{+}}\,\delta_{jm}\delta_{ik'}$, and nothing else in this row depends on $z$, so the $z$ integral can be done exactly,
\begin{align}
J(\kappa,b)\equiv\int_{z}e^{-i\kappa\cdot(y'-z)}\,\frac{(y'-z)\cdot(x-z)}{(y'-z)^{2}(x-z)^{2}},\qquad \kappa=\frac{p^{+}}{k^{+}}\,k,\quad b=x-y',\qquad J\xrightarrow{\kappa\to0}2\pi\ln\frac{2e^{-\gamma_{E}}}{|\kappa|\,|b|}.
\end{align}
$J(0,b)$ is infinite (the integrand is $1/z^{2}$ at large $z$), so the subtraction term of the plus prescription, and with it the $\log$ term $f(0)\log\frac{\vee-k^{+}}{\Lambda}$ that was moved to the evolution list and the $\log(1-k^{+}/\vee)$ term, are not defined separately. Keeping $\Lambda$ as the lower limit and integrating over $p^{+}$ after $z$ gives (checked numerically)
\begin{align}
\int_{\Lambda}^{\infty}\frac{dp^{+}}{p^{+}}\,J\Big(\frac{p^{+}}{k^{+}}k,\,b\Big)=\pi\,\ln^{2}\!\left[\frac{2\,e^{-\gamma_{E}}\,k^{+}}{\Lambda\,|k|\,|b|}\right]+O(\Lambda),
\end{align}
with no further constant. For a finite upper limit $\vee-k^{+}$ one subtracts $\int_{\vee-k^{+}}^{\infty}\frac{dp^{+}}{p^{+}}J$, which vanishes for $(\vee-k^{+})|k||b|/k^{+}\gg1$; there is no $\log\vee$. Hence the whole term is
\begin{align}
&-\frac{1}{(2\pi)^{3}}\frac{ig^{4}f^{c'ca}}{8\pi^{4}}\frac{1}{k^{+}}\int_{w,w',y,x,x'} e^{-ik\cdot(w'-w)}\Big[\ln^{2}\frac{k^{+}}{\Lambda}+2\ln\frac{k^{+}}{\Lambda}\,\ln\frac{2e^{-\gamma_{E}}}{|k||x-w'|}+\ln^{2}\frac{2e^{-\gamma_{E}}}{|k||x-w'|}\Big]\frac{(x'-w')\cdot(y-w)}{(x'-w')^{2}(y-w)^{2}}\notag\\
&\qquad\times\big[U^{c'e'}(w')-U^{c'e'}(x')\big]\rho^{e'}(x')\,U^{ce}(x)\big[U^{ad}(y)-U^{ad}(w)\big]\rho^{e}(x)\rho^{d}(y).
\end{align}
Its complex conjugate (relabel $w\leftrightarrow w'$, $x'\leftrightarrow y$) has the same colour and operator structure as the pair $\mathcal{P}_{\Lambda}$ of Group I, Row III, with the opposite sign: the $\ln^{2}\frac{k^{+}}{\Lambda}$ and $\ln\frac{k^{+}}{\Lambda}$ pieces cancel exactly between the two rows, and only the $\ln^{2}\frac{2e^{-\gamma_{E}}}{|k||x-w'|}$ piece is kept here. (The lower limit $\Lambda$ is the one of \texttt{Group\_III.tex}.)
}
"""
NOTE_III_VI=r"""
\rnote{Same structure as Row V: the term with $\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}+\frac{1}{p^{+}+k^{+}}$ is $-1/p^{+}$ times $\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)$, and its subtraction term has a $1/z^{2}$ tail. The pieces of the colour factor that do not contain $U(z)$ could be integrated exactly as in Row V; taken together for Rows III and VI they reduce to a commutator, $-\frac{1}{2}\{A^{a},X^{d'}\}+A^{a}X^{d'}=\frac{1}{2}[A^{a},X^{d'}]$ with $A^{a}=U^{ad}(y)\rho^{d}(y)$, $X^{d'}=U^{d'e}(x)\rho^{e}(x)$, i.e.\ a two-$\rho$ term. The pieces with $U(z)$ cannot be integrated in closed form. For them the subtraction term and the plus-prescription integral are not separately convergent at large $z$; they are meaningful only when the rows with the same soft gluon are added before the $z$ integration (I expect the partner to be Group I Row V, first part, but have not verified this). This row is left as you wrote it.}
"""

# ================= Other_terms.tex =================
s=open(UP+"f10fa593-Other_terms.tex").read()
s=rep(s,r"\begin{document}",MAC+r"\begin{document}")
old=(r" & -\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\int_{x,y,z,x'}\log{\left[\frac{(z-w)^{2}}{(x-w)^{2}}\right]}\,\frac{(x'-w')\cdot(y-w)\, (x-z)\cdot (w-z)}{(x'-w')^{2}(y-w)^{2}(x-z)^{2}(w-z)^{2}}\notag \\ & "
     +UU+r"\times \left[U^{de}(x)\,U^{bc}(w)\, \rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]\Bigg\}"+"\n"+r"\end{align}"+"\n"+r"where: \begin{align}")
new=new_last(" & ")+r"\Bigg\}"+"\n"+r"\end{align}"+"\n"+REMARK_III+r"where: \begin{align}"
s=rep(s,old,new)
# Group III Row V
colV=r"\left[U^{c'e'}(y')\,\rho^{e'}(x')-U^{c'e'}(x')\,\rho^{e'}(x')\right]\times \notag \\ &  \left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right]"
brk=r"\left[\frac{\delta_{k'm}\,\delta_{ij}}{(p^{+}+k^{+})}-\left(\left[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\right]_{+}+\frac{1}{(p^{+}+k^{+})}\right)\,\delta_{jm}\,\delta_{ik'}\right]\, "
s=rep(s,brk+colV, r"\chg{\left[\frac{\delta_{k'm}\,\delta_{ij}}{(p^{+}+k^{+})}\right]}\, "+colV)
oldlog=(r" \notag \\"+"\n"+r"& -\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'ca}}{8\pi^{5}}\, \frac{1}{k^{+}}\,\log{\left(1-\frac{k^{+}}{\vee}\right)} \int_{w,y',z,y,x,x'}\, e^{-ik\cdot (y'-w)}\,\frac{(y'-z)\cdot (x-z)\,(x'-y')\cdot (y-w)}{(y'-z)^{2}\,(x'-y')^{2}\,(y-w)^{2}\,(x-z)^{2}}\notag \\ & \left[U^{c'e'}(y')\,\rho^{e'}(x')-U^{c'e'}(x')\,\rho^{e'}(x')\right]\times  \left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right]"+"\n"+r"\end{align}"+"\n"+r"\subsection{Row VI}")
s=rep(s,oldlog, rowV_closed()+"\n"+r"\end{align}"+"\n"+REMARK_V+r"\subsection{Row VI}")
# notes on Rows III and VI
s=rep(s,r"\subsection{Row IV}"+"\n"+r"\textsc{NA}"+"\n"+r"\subsection{Row V}", NOTE_III_VI+r"\subsection{Row IV}"+"\n"+r"\textsc{NA}"+"\n"+r"\subsection{Row V}")
i=s.rindex(r"\end{document}") if r"\end{document}" in s else len(s)
s=s[:i]+NOTE_III_VI+s[i:]
open(OUT+"Other_terms.tex","w").write(s)

# ================= Limiting_Cases.tex =================
s=open(UP+"56ca25e0-Limiting_Cases.tex").read()
s=rep(s,r"\begin{document}",MAC+r"\begin{document}")
old=(r" & -\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\int_{x,y,z,x'}\log{\left[\frac{(z-w)^{2}}{(x-w)^{2}}\right]}\,\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot (w-z)}{(x'-w')^{2}\,(y-w)^{2}\,(z-w)^{2}\,(x-z)^{2}}\notag \\ & "
     +UU+r"\times \left[U^{de}(x)\,U^{bc}(z)\, \rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]\Bigg\}"+"\n"+r"\end{align}"+"\n"+r"where: \begin{align}")
new=new_last(" & ")+r"\Bigg\}"+"\n"+r"\end{align}"+"\n"+REMARK_III+r"where: \begin{align}"
s=rep(s,old,new)
# summed form
old=(r"&-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{w,w'}e^{-ik\cdot(w'-w)}\int_{x,y,z,x'}\log\left(\frac{(z-w)^{2}}{(x-w)^{2}}\right)\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot(w-z)}{(x'-w')^{2}(y-w)^{2}(x-z)^{2}(w-z)^{2}}\notag\\"+"\n"
     r"&\times\left[U^{ab'}(x')\rho^{b'}(x')-U^{ab'}(w')\rho^{b'}(x')\right]\left[U^{de}(x)U^{bc}(w)\rho^{e}(x)\rho^{c}(y)-U^{de}(x)U^{bc}(y)\rho^{e}(x)\rho^{c}(y)\right]"+"\n"+r"\end{align}")
new=(r"&\chg{-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{w,w'}e^{-ik\cdot(w'-w)}\int_{x,y,z,x'}\frac{1}{(z-w)^{2}}\log\left(\frac{(z-w)^{2}}{(x-w)^{2}}\right)\frac{(x'-w')^{i}}{(x'-w')^{2}}\left[U^{ab'}(x')-U^{ab'}(w')\right]\rho^{b'}(x')}\notag\\"+"\n"
     r"&\chg{\qquad\times U^{de}(x)\rho^{e}(x)\left\{\frac{(y-z)^{i}}{(y-z)^{2}}\left[U^{bc}(z)-U^{bc}(y)\right]-\frac{(x-w)^{2}}{(x-w)^{2}+(z-w)^{2}}\,\frac{(y-w)^{i}}{(y-w)^{2}}\left[U^{bc}(w)-U^{bc}(y)\right]\right\}\rho^{c}(y)}\notag \\"+"\n"
     r"&\chg{-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\frac{\pi^{3}}{6}\int_{w,w'}e^{-ik\cdot(w'-w)}\int_{x,y,x'}\frac{(x'-w')\cdot(y-w)}{(x'-w')^{2}(y-w)^{2}}\left[U^{ab'}(x')-U^{ab'}(w')\right]\rho^{b'}(x')}\notag\\"+"\n"
     r"&\chg{\qquad\times U^{de}(x)\rho^{e}(x)\left[U^{bc}(w)-U^{bc}(y)\right]\rho^{c}(y)}"+"\n"+r"\end{align}"+"\n"
     r"\rnote{The last block is changed as explained in the remark above.}")
s=rep(s,old,new)
# Group III Row V
colV2=r"\left[U^{c'e'}(y')\,\rho^{e'}(x')-U^{c'e'}(x')\,\rho^{e'}(x')\right]\times \notag \\ &  \left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right]"
s=rep(s,brk+colV2, r"\chg{\left[\frac{\delta_{k'm}\,\delta_{ij}}{(p^{+}+k^{+})}\right]}\, "+colV2)
oldend=(r"\left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right] "+"\n"+r"\end{align}"+"\n"+r"\subsection{Row VI}")
s=rep(s,oldend, r"\left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right] "+rowV_closed()+"\n"+r"\end{align}"+"\n"+REMARK_V+r"\subsection{Row VI}")
s=rep(s,r"\subsection{Row IV}"+"\n"+r"\textsc{NA}"+"\n"+r"\subsection{Row V}", NOTE_III_VI+r"\subsection{Row IV}"+"\n"+r"\textsc{NA}"+"\n"+r"\subsection{Row V}")
i=s.rindex(r"\end{document}") if r"\end{document}" in s else len(s)
s=s[:i]+NOTE_III_VI+s[i:]
s=rep(s,r"R^{i} =& -2k^{i}", r"R^{i} =& \chg{-2q^{i}}")
s=rep(s,r"\Bigg]\Bigg\} \notag \\ & \Bigg[U^{ac'}(w')", r"\Bigg] \notag \\ & \Bigg[U^{ac'}(w')")
open(OUT+"Limiting_Cases.tex","w").write(s)

# ================= Final_terms.tex =================
s=open(UP+"690ba62f-Final_terms.tex").read()
s=rep(s,r"\begin{document}",MAC+r"\begin{document}")
old=(r"&-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}\,f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{x,y,z,w,x',w'}e^{-ik\cdot(w'-w)}\log\frac{(z-w)^{2}}{(x-w)^{2}}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)\,\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\notag\\"+"\n"
     r"&\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\,U^{de}(x)\Big[U^{bc}(w)-U^{bc}(y)\Big]\rho^{b'}(x')\rho^{e}(x)\rho^{c}(y)"+"\n"+r"\end{align}")
new=(r"&\chg{-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}\,f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{x,y,z,w,x',w'}e^{-ik\cdot(w'-w)}\,\frac{1}{(z-w)^{2}}\log\frac{(z-w)^{2}}{(x-w)^{2}}\;\mathcal{K}^{i}(x'-w')\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\,U^{de}(x)}\notag\\"+"\n"
     r"&\chg{\qquad\times\bigg\{\mathcal{K}^{i}(y-z)\Big[U^{bc}(z)-U^{bc}(y)\Big]-\frac{(x-w)^{2}}{(x-w)^{2}+(z-w)^{2}}\,\mathcal{K}^{i}(y-w)\Big[U^{bc}(w)-U^{bc}(y)\Big]\bigg\}\rho^{b'}(x')\rho^{e}(x)\rho^{c}(y)}\notag\\"+"\n"
     r"&\chg{-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}\,f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\frac{\pi^{3}}{6}\int_{x,y,w,x',w'}e^{-ik\cdot(w'-w)}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)}\notag\\"+"\n"
     r"&\chg{\qquad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\,U^{de}(x)\Big[U^{bc}(w)-U^{bc}(y)\Big]\rho^{b'}(x')\rho^{e}(x)\rho^{c}(y)}"+"\n"+r"\end{align}"+"\n"
     +REMARK_III.replace(r"the last three lines of $N_{3b}$ above",r"the last three lines above").replace(r"Remark on the last term of $N_{3b}$ (changed).",r"Remark on the last block (changed)."))
s=rep(s,old,new)
# Group III Row V
old=(r"\bigg[\frac{1}{p^{+}+k^{+}}\,\mathcal{K}(y'-z)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(y-w)\cdot\mathcal{K}(x-z)\notag \\ & -\bigg(\left[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\right]_{+}+\frac{1}{p^{+}+k^{+}}\bigg)\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(y-w)\cdot\mathcal{K}(x'-y')\bigg]\notag \\ & \times\Big[U^{c'e'}(y')-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ce}(x)\Big[U^{ad}(y)-U^{ad}(w)\Big]\rho^{e}(x)\rho^{d}(y)\notag\\")
new=(r"\chg{\frac{1}{p^{+}+k^{+}}}\,\mathcal{K}(y'-z)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(y-w)\cdot\mathcal{K}(x-z)\notag \\ & \times\Big[U^{c'e'}(y')-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ce}(x)\Big[U^{ad}(y)-U^{ad}(w)\Big]\rho^{e}(x)\rho^{d}(y)\notag\\"+"\n"
     r"&\chg{-\,\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}\,f^{c'ca}}{8\pi^{4}}\,\frac{1}{k^{+}}\int_{w,w',y,x,x'}e^{-ik\cdot(w'-w)}\,\ln^{2}\!\left[\frac{2\,e^{-\gamma_{E}}}{|k|\,|x-w'|}\right]\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)}\notag\\"+"\n"
     r"&\chg{\qquad\times\Big[U^{c'e'}(w')-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ce}(x)\Big[U^{ad}(y)-U^{ad}(w)\Big]\rho^{e}(x)\rho^{d}(y)}\notag\\")
s=rep(s,old,new)
old=(r"\times\Big[U^{c'e'}(z)-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ce}(x)\Big[U^{ad}(y)-U^{ad}(w)\Big]\rho^{e}(x)\rho^{d}(y)"+"\n"+r"\end{align}"+"\n"+r"\subsection{Row VI}")
s=rep(s,old, r"\times\Big[U^{c'e'}(z)-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ce}(x)\Big[U^{ad}(y)-U^{ad}(w)\Big]\rho^{e}(x)\rho^{d}(y)"+"\n"+r"\end{align}"+"\n"+REMARK_V.replace(r"\Big)\delta_{jm}\delta_{ik'}$ has been",r"\Big)\,\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(y-w)\cdot\mathcal{K}(x'-y')$ has been").replace(r"replaced by the last two lines",r"replaced by the second term (two lines)").replace(r"this term is $-\frac{1}{p^{+}}\,\delta_{jm}\delta_{ik'}$",r"this term is $-\frac{1}{p^{+}}\,\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(y-w)\cdot\mathcal{K}(x'-y')$")+r"\subsection{Row VI}")
s=rep(s,r"\subsection{Row IV}"+"\n"+r"\textsc{NA}"+"\n"+r"\subsection{Row V}"+"\n"+r"\begin{align}"+"\n"+r"\frac{d^{3}N_{\rm NLO}}{d^{3}k} & = \frac{1}{(2\pi)^{3}}\,\frac{ig^{4}\,f^{c'ca}}{4\pi^{4}}",
      NOTE_III_VI+r"\subsection{Row IV}"+"\n"+r"\textsc{NA}"+"\n"+r"\subsection{Row V}"+"\n"+r"\begin{align}"+"\n"+r"\frac{d^{3}N_{\rm NLO}}{d^{3}k} & = \frac{1}{(2\pi)^{3}}\,\frac{ig^{4}\,f^{c'ca}}{4\pi^{4}}")
# note after Row VI (three-rho Group III): find the end of its align
k=s.index(r"\subsection{Row VI}", s.index(r"f^{c'ca}}{4\pi^{4}}"))
e=s.index(r"\end{align}",k)+len(r"\end{align}")
s=s[:e]+"\n"+NOTE_III_VI+s[e:]
# notation typo
s=rep(s,r"where $\mathcal{K}(x,y,z) = ", r"where $\chg{K(x,y,z)} = ")
open(OUT+"Final_terms.tex","w").write(s)
print("ok")
