import re
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
def rep(s,old,new,n=1):
    assert s.count(old)==n,(s.count(old),old[:80]); return s.replace(old,new)
MAC=r"""\newcommand{\chg}[1]{{\color{red}#1}}
\newcommand{\rnote}[1]{{\color{teal}[\textit{#1}]}}
"""
s=open(UP+"56ca25e0-Limiting_Cases.tex").read()
s=re.sub(r"\\title\{[^\n]*\n", r"\\title{Limiting finite terms: $\\Lambda\\to0$ and $\\vee\\to\\infty$}"+"\n", s, count=1)
s=rep(s,r"\begin{document}",MAC+r"\begin{document}")
INTRO=r"""
\noindent\textbf{Procedure.} Start from \texttt{Finite\_terms.tex} (the black terms of \texttt{Stage1\_split} with $\Lambda\to0$) and take $\vee\to\infty$ at fixed transverse coordinates: terms with $\log(1-\frac{k^{+}}{\vee})$ vanish, $\frac{k^{+}}{\vee}\to0$ inside logarithms, $\frac{1}{\vee}\to0$. The result coincides with your \texttt{Limiting\_Cases.tex} except at the places marked in {\color{red}red}; notes are in {\color{teal}teal}.
"""
s=rep(s,r"\maketitle",r"\maketitle"+INTRO)
# Group I Row II N1
s=rep(s,r"\Bigg\{ -\frac{67}{9}+\frac{\pi^{2}}{3}+\log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}+\log{\left[1-\frac{k^{+}}{\vee}\right]}\,\log{\left[\frac{\vee}{k^{+}}\right]}\Bigg\}",
        r"\Bigg\{ \chg{\frac{11}{3}\left[-\frac{2}{\epsilon}-2\gamma-\log{\left[\frac{(x-w)^{2}\,\mu^{2}_{\overline{MS}}}{4}\right]}\right]}-\frac{67}{9}+\frac{\pi^{2}}{3}\Bigg\}")
s=rep(s,r"R^{i} =& -2k^{i}", r"R^{i} =& \chg{-2q^{i}}")
k=s.index(r"\subsection{Row III}", s.index(r"\section{Group I (with c.c.)}"))
s=s[:k]+r"""\rnote{$N_{1}$: the two terms with $\log(1-\frac{k^{+}}{\vee})$ that were still written here vanish and are removed; the pole bracket times $\frac{11}{3}$ is added (it cancels exactly against the $-\frac{11}{6}$ term of Row V, $2b\,{\rm UV}$). In $R^{i}$ the first term had $k^{i}$ for $q^{i}$.}
"""+s[k:]
# Group I Row III: last term of N3b, in the list and in the sum
UU=r"\left[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\, \rho^{b'}(x')\right]"
s=rep(s,r"\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot (w-z)}{(x'-w')^{2}\,(y-w)^{2}\,(z-w)^{2}\,(x-z)^{2}}\notag \\ & "+UU,
        r"\chg{\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}\,(y-z)^{2}\,(z-w)^{2}}}\notag \\ & "+UU)
old=(r"\log\left(\frac{(z-w)^{2}}{(x-w)^{2}}\right)\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot(w-z)}{(x'-w')^{2}(y-w)^{2}(x-z)^{2}(w-z)^{2}}\notag\\"+"\n"
     r"&\times\left[U^{ab'}(x')\rho^{b'}(x')-U^{ab'}(w')\rho^{b'}(x')\right]\left[U^{de}(x)U^{bc}(w)\rho^{e}(x)\rho^{c}(y)-U^{de}(x)U^{bc}(y)\rho^{e}(x)\rho^{c}(y)\right]"+"\n"+r"\end{align}")
new=(r"\log\left(\frac{(z-w)^{2}}{(x-w)^{2}}\right)\chg{\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}(z-w)^{2}}}\notag\\"+"\n"
     r"&\times\left[U^{ab'}(x')\rho^{b'}(x')-U^{ab'}(w')\rho^{b'}(x')\right]\left[U^{de}(x)\chg{U^{bc}(z)}\rho^{e}(x)\rho^{c}(y)-U^{de}(x)U^{bc}(y)\rho^{e}(x)\rho^{c}(y)\right]"+"\n"+r"\end{align}"+"\n"
     r"\rnote{Last term of $N_{3b}$, here and in the list above: this is what $\Lambda\to0$ gives for the last term of $N_{3b}$ of \texttt{Stage1\_split}. Your file had the kernel $\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot(w-z)}{(x'-w')^{2}(y-w)^{2}(x-z)^{2}(w-z)^{2}}$ (with $U(z)$ in the list and $U(w)$ in the sum). The term behaves like $1/(z-w)^{2}$ at $z\to w$, where the colour factor does not vanish.}")
s=rep(s,old,new)
# Group I Row V: 2b UV
uv=r"""\begin{align}
\chg{\frac{d^{3}N_{2b\,\rm UV}}{d^{3}k}}&\chg{=-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}\,\frac{11}{6}\int_{x,y,x',w'}e^{-ik\cdot(w'-x)}\,\left[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y-x)^{2}\,\bar\mu^{2}_{\overline{MS}}}{4}\right]\frac{(x'-w')\cdot(y-x)}{(x'-w')^{2}\,(y-x)^{2}}}\notag \\ & \chg{\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[U^{af}(y)-U^{af}(x)\Big]\rho^{b'}(x')\,\rho^{f}(y)}
\end{align}
"""
m=re.search(r"\\begin\{align\}\s*\n\s*\\frac\{d\^\{3\}N_\{2b\\,\\rm rest\}\}",s); assert m
s=s[:m.start()]+uv+s[m.start():]
# stray brace in Group II Row I
s=rep(s,r"\Bigg]\Bigg\} \notag \\ & \Bigg[U^{ac'}(w')", r"\Bigg] \notag \\ & \Bigg[U^{ac'}(w')")
# Group II Row IV: all pieces
i=s.index(r"\subsection{Row IV}", s.index(r"\section{Group II}"))
e=s.index(r"\end{align}",i)+len(r"\end{align}")
rowIV=r"""
\rnote{Your file listed only the Finite piece, which is unchanged. The other black pieces of this row are added in red.}
\begin{align}
\chg{\frac{d^{3}N_{\rm DGLAP}}{d^{3}k}=}&\chg{-\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{g^{2}}{8\pi^{2}}\frac{1}{k^{+}}\int_{0}^{1}\frac{d\xi}{\xi^{2}}\,\widetilde{\mathcal P}_{gg}(\xi)\,\int_{x,y,x',y'}e^{-ik\cdot(y'-y)\frac{1}{\xi}}\Big[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}}{4\,\xi^{2}}\Big]\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}}\notag\\
&\chg{\times\left[U^{ce'}(y') - U^{ce'}(x')\right]\left[U^{ce}(y) - U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x),}\\
&\chg{\widetilde{\mathcal{P}}_{gg}(\xi) = 2N_{c}\Big[\xi(1-\xi) + \frac{\xi}{(1-\xi)_{+}} + \Big[\frac{1-\xi}{\xi}\Big]_{+}\Big]+\frac{\beta_{0}}{2}\,\delta(1-\xi),}\notag\\
&\chg{\frac{\xi}{(1-\xi)_{+}}\,g(\xi)\equiv\frac{\xi\,g(\xi)-g(1)}{1-\xi},\qquad \Big[\frac{1-\xi}{\xi}\Big]_{+}g(\xi)\equiv\frac{1-\xi}{\xi}\,\big[g(\xi)-g(0)\big],}\notag\\[4pt]
\chg{\frac{d^{3}N_{\rm RC}}{d^{3}k}=}&\chg{+\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{g^{2}\beta_{0}}{16\pi^{2}}\frac{1}{k^{+}}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\Big[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}}{4}\Big]\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}}\notag\\
&\chg{\times\left[U^{ce'}(y') - U^{ce'}(x')\right]\left[U^{ce}(y) - U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x),}\\[4pt]
\chg{\frac{d^{3}N_{\rm rest\ of\ BFKL\ piece}}{d^{3}k}=}&\chg{-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\int_{x,x',z,w,w'}e^{-ik\cdot(w'-w)}\,\frac{(w-z)\cdot (w'-z)\,(x-z)\cdot (x'-z)}{(w-z)^{2}(x-z)^{2}(w'-z)^{2}(x'-z)^{2}}}\notag\\
&\chg{\times\left[U^{ce'}(z) - U^{ce'}(x')\right]\left[U^{ce}(z) - U^{ce}(x)\right]\rho^{e'}(x') \rho^{e}(x).}
\end{align}
\rnote{DGLAP: in \texttt{Group\_II.tex} the DGLAP piece has a second term, $2N_{c}\frac{1-\xi}{\xi}\frac{4\pi}{k^{2}}\int K(x,x',z)[\ldots]$, integrated from $\xi=\frac{k^{+}}{\vee}$; for $\vee\to\infty$ it diverges and has to be combined with the first term, which gives the plus distribution $[\frac{1-\xi}{\xi}]_{+}$ defined above ($g(0)$ is the $\xi\to0$ limit after $y'-y=\xi u$). The BFKL piece of that file is $[\log\frac{\vee}{k^{+}}-(1-\frac{k^{+}}{\vee})]$ times the kernel; after $\log\frac{\vee}{k^{+}}$ is absorbed, $-(1-\frac{k^{+}}{\vee})\to-1$ is left, which is the last line. The rest of the JIMWLK piece, $\log(1-\frac{k^{+}}{\vee})$ times its kernel, vanishes. The $\frac{\beta_{0}}{2}\delta(1-\xi)$ term of DGLAP is exactly minus the RC term. \textbf{The last line, taken alone, diverges logarithmically at large $z$}: the $w,w'$ integrals give $\frac{4\pi^{2}}{k^{2}}$ and $\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}(x'-z)^{2}}\sim\frac{1}{z^{2}}$, so it converges only if $\int_{x}[U(z)-U(x)]\rho(x)$ vanishes at large $z$. It is the same $z$ integral that multiplies $\log\frac{\vee}{k^{+}}$ in the BFKL piece and that appears in the second DGLAP term of \texttt{Group\_II.tex}, and the logarithmic part of the DGLAP term has the same tail for every $\xi$.}
"""
s=s[:e]+rowIV+s[e:]
# Group III: vee -> infinity in the p+ integrals
n=s.count(r"\int_{0}^{\vee-k^{+}}")
s=s.replace(r"\int_{0}^{\vee-k^{+}}",r"\int_{0}^{\chg{\infty}}")
k=s.index(r"\subsection{Row I}", s.index(r"\section{Group III\, (with c.c.)}"))
s=s[:k]+r"""\rnote{Rows III, V, VI: the upper limit $\vee-k^{+}$ of the $p^{+}$ integrals is now $\infty$ ("""+str(n)+r""" places, red). This is possible for all terms except the plus-prescription one: there $\int_{0}^{\vee-k^{+}}\frac{dp^{+}}{p^{+}}[f(p^{+})-f(0)]$ grows like $-f(0)\log\frac{\vee}{k^{+}}$. \texttt{BFKL.tex} has already taken this $-f(0)\log\frac{\vee}{k^{+}}$ out (it cancels the blue term there), so in these rows the plus prescription has to be read as}
\begin{align}
\chg{\Big(\Big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\Big]_{+}+\frac{1}{p^{+}+k^{+}}\Big)f(p^{+})\;\equiv\;\frac{f(p^{+})-\theta(k^{+}-p^{+})\,f(0)}{p^{+}},\qquad \Big[\frac{1}{p^{+}}\Big]_{+}f(p^{+})\equiv\frac{f(p^{+})-\theta(k^{+}-p^{+})\,f(0)}{p^{+}}.}
\end{align}
"""+s[k:]
open("/home/user/correlations/nlo_review/stage2/Finite_terms_limit.tex","w").write(s)
print("ok; p+ upper limits changed:",n," remaining 'vee':",len(re.findall(r"\\vee",s.split(r"\section{Group I (with c.c.)}")[1])))
