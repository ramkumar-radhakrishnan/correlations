import re,sys
sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/stage1")
from lib_split import get
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
def rep(s,old,new,n=1):
    assert s.count(old)==n,(s.count(old),old[:80]); return s.replace(old,new)
MAC=r"""\newcommand{\chg}[1]{{\color{red}#1}}
\newcommand{\rnote}[1]{{\color{teal}[\textit{#1}]}}
"""
s=open(UP+"f10fa593-Other_terms.tex").read()
s=re.sub(r"\\title\{[^\n]*\n", r"\\title{Finite terms: Stage1\\_split with the BFKL and JIMWLK terms removed, $\\Lambda\\to0$}"+"\n", s, count=1)
s=rep(s,r"\begin{document}",MAC+r"\begin{document}")
INTRO=r"""
\noindent\textbf{Procedure.} Start from \texttt{Stage1\_split}. Remove every blue term (they are absorbed into BFKL and JIMWLK). In what is left take $\Lambda\to0$ at fixed transverse coordinates; $\vee$ is kept. Nothing else is done. The result coincides with your \texttt{Other\_terms.tex} except at the places marked in {\color{red}red}; notes are in {\color{teal}teal}. The power terms $\propto\frac{1}{\Lambda}$ of Group I Row III cancel identically between $N_{1}$, $N_{2}$ and $N_{3a}$ and are not listed.
"""
s=rep(s,r"\maketitle",r"\maketitle"+INTRO)
# (b) Group I Row II, N1: black pole terms
s=rep(s,r"\Bigg\{ -\frac{67}{9}+\frac{\pi^{2}}{3}+\log^{2}",
        r"\Bigg\{ \chg{\left[-\frac{2}{\epsilon}-2\gamma-\log{\left[\frac{(x-w)^{2}\,\mu^{2}_{\overline{MS}}}{4}\right]}\right]\left[\frac{11}{3}+\log{\left(1-\frac{k^{+}}{\vee}\right)}\right]}-\frac{67}{9}+\frac{\pi^{2}}{3}+\log^{2}")
k=s.index(r"\subsection{Row III}", s.index(r"\section{Group I (with c.c.)}"))
s=s[:k]+r"""\rnote{$N_{1}$: the pole bracket times $\frac{11}{3}+\log(1-\frac{k^{+}}{\vee})$ is black in \texttt{Stage1\_split} and is restored here. Its $\frac{11}{3}$ part cancels exactly against the $-\frac{11}{6}$ term of Row V, $2b\,{\rm UV}$, below. $N_{2a}$ is entirely JIMWLK and is not listed.}
"""+s[k:]
# (a) Group I Row III, last term of N3b: the Lambda->0 limit of the term of Group_I.tex
UU=r"\left[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\, \rho^{b'}(x')\right]"
old=(r"\log{\left[\frac{(z-w)^{2}}{(x-w)^{2}}\right]}\,\frac{(x'-w')\cdot(y-w)\, (x-z)\cdot (w-z)}{(x'-w')^{2}(y-w)^{2}(x-z)^{2}(w-z)^{2}}\notag \\ & "
     +UU+r"\times \left[U^{de}(x)\,U^{bc}(w)\, \rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]\Bigg\}"+"\n"+r"\end{align}")
new=(r"\log{\left[\frac{(z-w)^{2}}{(x-w)^{2}}\right]}\,\chg{\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}\,(y-z)^{2}\,(z-w)^{2}}}\notag \\ & "
     +UU+r"\times \left[U^{de}(x)\,\chg{U^{bc}(z)}\, \rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]\Bigg\}"+"\n"+r"\end{align}"+"\n"
     r"\rnote{Last term of $N_{3b}$: this is the $\Lambda\to0$ limit of the last term of $N_{3b}$ in \texttt{Stage1\_split}, $\log\frac{[k^{+}(x-w)-(k^{+}-\Lambda)(x-z)]^{2}}{[k^{+}(x-w)-\Lambda(x-z)]^{2}}\to\log\frac{(z-w)^{2}}{(x-w)^{2}}$, with the kernel and the Wilson line of that term. Your file had $\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot(w-z)}{(x'-w')^{2}(y-w)^{2}(x-z)^{2}(w-z)^{2}}$ and $U^{bc}(w)$ here, which is not what the limit gives. The blue term absorbed into JIMWLK from this row is, correspondingly, $-\frac{ig^{4}f^{abd}}{8\pi^{5}k^{+}}\log\frac{k^{+}}{\Lambda}\,\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}(z-w)^{2}}$ with $U^{bc}(z)$. Both terms behave like $1/(z-w)^{2}$ at $z\to w$, where the colour factor does not vanish.}")
s=rep(s,old,new)
# (c) Group I Row V: 2b UV, black part
uv=r"""\begin{align}
\chg{\frac{d^{3}N_{2b\,\rm UV}}{d^{3}k}}&\chg{=-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}\,\frac{11}{6}\int_{x,y,x',w'}e^{-ik\cdot(w'-x)}\,\left[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y-x)^{2}\,\bar\mu^{2}_{\overline{MS}}}{4}\right]\frac{(x'-w')\cdot(y-x)}{(x'-w')^{2}\,(y-x)^{2}}}\notag \\ & \chg{\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[U^{af}(y)-U^{af}(x)\Big]\rho^{b'}(x')\,\rho^{f}(y)}
\end{align}
"""
m=re.search(r"\\begin\{align\}\s*\n\s*\\frac\{d\^\{3\}N_\{2b\\,\\rm rest\}\}",s); assert m
s=s[:m.start()]+uv+s[m.start():]
# (d) Group II Row IV: all black pieces of the row
blk=get("II",250,267)
blk=rep(blk,r"\Big[\log\frac{k^{+}}{\Lambda}+\log\Big(1-\frac{k^{+}}{\vee}\Big)\Big]",r"\chg{\log\Big(1-\frac{k^{+}}{\vee}\Big)}",2)
blk=rep(blk,r"\Big[\log\frac{\vee}{k^{+}}-\Big(1-\frac{k^{+}}{\vee}\Big)\Big]",r"\chg{\Big[-\Big(1-\frac{k^{+}}{\vee}\Big)\Big]}")
blk=rep(blk,r"\Big[\frac{1}{k^{+}}\log\frac{\vee}{k^{+}}-\frac{1}{k^{+}}+\frac{1}{\vee}\Big]",r"\chg{\Big[-\frac{1}{k^{+}}+\frac{1}{\vee}\Big]}")
blk=blk.replace(r"\frac{d^{3}N_{\rm JIMWLK}}{d^{3}k}=",r"\frac{d^{3}N_{\rm rest\ of\ JIMWLK\ piece}}{d^{3}k}=").replace(r"\frac{d^{3}N_{\rm BFKL}}{d^{3}k}=",r"\frac{d^{3}N_{\rm rest\ of\ BFKL\ piece}}{d^{3}k}=")
i=s.index(r"\subsection{Row IV}", s.index(r"\section{Group II}"))
a=s.index(r"\begin{align}",i); e=s.index(r"\end{align}",a)+len(r"\end{align}")
noteIV=r"""\rnote{All black pieces of this row are listed: DGLAP, RC and Finite as in \texttt{Group\_II.tex}, and what is left of the JIMWLK and BFKL pieces once $\log\frac{k^{+}}{\Lambda}$ and $\log\frac{\vee}{k^{+}}$ are removed (red). Your file listed only Finite, with the lower limit of $\xi$ set to $0$. The second term of DGLAP is minus the full BFKL piece of \texttt{Group\_II.tex}; it is the subtraction that turns $\frac{1-\xi}{\xi}$ into a plus distribution.}
"""
s=s[:a]+noteIV+blk+s[e:]
# Group III note on vee
k=s.index(r"\subsection{Row I}", s.index(r"\section{Group III\, (with c.c.)}"))
s=s[:k]+r"""\rnote{Rows III, V, VI are left as integrals over $p^{+}$, as in your file, including your subtraction terms with $\mathcal{K}(x'-z)[U(z)\ldots]$ that remove the $\log\frac{\vee}{k^{+}}$ of the $-\frac{1}{k^{+}}$ and $-\frac{p^{+}}{k^{+2}}$ terms. One point about $\vee$: \texttt{BFKL.tex} also takes $-f(0)\log\frac{\vee}{k^{+}}$ out of the plus-prescription integral (it cancels the blue term there). To match that, the plus-prescription integral here should be read as $\int_{0}^{\vee-k^{+}}\frac{dp^{+}}{p^{+}}[f(p^{+})-\theta(k^{+}-p^{+})f(0)]$; as printed, with the subtraction over the whole range, it still contains $-f(0)\log\frac{\vee}{k^{+}}$.}
"""+s[k:]
open("/home/user/correlations/nlo_review/stage2/Finite_terms.tex","w").write(s)
print("ok")
