import re
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
OUT="/home/user/correlations/nlo_review/final/"
MAC=r"""\newcommand{\chg}[1]{{\color{red}#1}}
\newcommand{\rnote}[1]{{\color{teal}[\textit{#1}]}}
"""
def rep(s,old,new,n=1):
    assert s.count(old)==n,(s.count(old),old[:80]); return s.replace(old,new)
# ---------------- BFKL.tex ----------------
L=open(UP+"9ea88190-BFKL.tex").read().split("\n")
cnt=0
for i in range(40,44):   # lines 41-44
    def fix(m):
        global cnt; cnt+=1
        return r"\chg{U^{"+m.group(1)+"}("+m.group(2)+")}"
    L[i]=re.sub(r"U\^\{([^}]*)\}\(([a-z]'?)(?=\s+U\^|\\,)",fix,L[i])
print("BFKL parentheses fixed:",cnt)
s="\n".join(L)
s=rep(s,r"\begin{document}",MAC+r"\begin{document}")
note4=r"""
\rnote{Status of this section. The two-$\rho$ list that follows from $N_{a}+N_{b}$ by symmetric reordering (first list below) is confirmed. Not confirmed: that $N_{a}+N_{b}$ themselves are the three-$\rho$ terms obtained by reordering the four-$\rho$ expression of \texttt{4rho.tex}, and the second two-$\rho$ list of this section. In lines 41--44 five closing brackets were missing; they are added in red.}
"""
s=rep(s,r"\subsection{From four $\rho$'s}",r"\subsection{From four $\rho$'s}"+note4)
noteII=r"""
\rnote{Status of this entry. \texttt{Group\_I.tex}, Row II, has no $\log\frac{\vee}{k^{+}}$ at fixed transverse momenta: in the momentum-space form the coefficient is $-2k^{i}+2k^{i}+\frac{(2p^{2}-2k\cdot p+k^{2})k^{i}+k^{2}p^{i}-p^{2}k^{i}-k^{2}p^{i}}{(k-p)^{2}}-k^{i}=0$, and $N_{1}$, $N_{2a}$ have none either. The entry below is equal to $-\frac{ig^{4}f^{hbc}}{32\pi^{5}k^{+}}\log\frac{\vee}{k^{+}}\,\mathcal{K}^{i}(x'-w')\big[\mathcal{K}^{i}(x-w)-\mathcal{K}^{i}(y-w)\big]K(x,y,z)$ times the same colour factor (use $\int_{z}[K(x,y,z)-K(w,y,z)]=\pi\log\frac{(w-y)^{2}}{(x-y)^{2}}$). The two-$\rho$ list written after it does follow from it, and the cancellation of all symmetrised three-$\rho$ $\log\frac{\vee}{k^{+}}$ terms holds only if this entry and $N_{a}+N_{b}$ are both included. Correspondingly, \texttt{Evolution\_three\_rho.tex} omits, for this colour structure, the $K(x,y,z)$ part of the $\log\frac{k^{+}}{\Lambda}$ term of $N_{2a}$ of \texttt{Group\_I.tex}. Net effect: relative to \texttt{Group\_I.tex} the evolution files contain $-\log\frac{\vee}{\Lambda}$ times this structure. With that, all three-$\rho$ logarithms cancel and the BFKL and JIMWLK terms are the correct evolution of the leading order, so the bookkeeping closes; but the term is not derived from Row II as written in \texttt{Group\_I.tex}, and its origin has to be found (a candidate is the ordering of the $\rho^{b}\rho^{b}\rho^{a}$ term of $\mathbb{A}^{(3)}$, which carries $\log\frac{\vee}{\Lambda}$).}
"""
s=rep(s,r"\subsection{Row II}"+"\n"+r"\begin{align}"+"\n"+r"\frac{d^{3}N_{2a\,\rm NLO}}{d^{3}k} & =  -\frac{1}{(2\pi)^{3}}\,\frac{i\,g^{4}\,f^{hbc}}{32\pi^{5}}",
        r"\subsection{Row II}"+noteII+r"\begin{align}"+"\n"+r"\frac{d^{3}N_{2a\,\rm NLO}}{d^{3}k} & =  -\frac{1}{(2\pi)^{3}}\,\frac{i\,g^{4}\,f^{hbc}}{32\pi^{5}}")
open(OUT+"BFKL.tex","w").write(s)
# ---------------- Group_II.tex ----------------
s=open(UP+"24c35375-Group_II.tex").read()
s=rep(s,r"\begin{document}",MAC+r"\begin{document}")
s=rep(s,r"U^{ac}(w)\,^{c'b}(z)",r"U^{ac}(w)\,\chg{U}^{c'b}(z)")
open(OUT+"Group_II.tex","w").write(s)
# ---------------- Group_III.tex ----------------
s=open(UP+"96b5b029-Group_III.tex").read()
s=rep(s,r"\begin{document}",MAC+r"\begin{document}")
s=rep(s,r"(k^{+},- U(w^{\prime})^{c^{\prime}b}(z)",r"\chg{(k^{+},-w^{\prime})\,U^{c^{\prime}b}(z)}")
open(OUT+"Group_III.tex","w").write(s)
print("ok")
