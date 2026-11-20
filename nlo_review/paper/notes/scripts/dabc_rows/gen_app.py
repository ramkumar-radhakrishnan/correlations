# Generates the appendix block app:fin_sym (row-by-row operator decomposition of the cubic rows).
import sys
D = r"\delta^{(2)}"
def S2(a, p, b, q): return r"\{\rho^{%s}(%s)\rho^{%s}(%s)\}_{\rm sym}" % (a, p, b, q)
def S3(p1, p2, p3): return r"\{\rho^a(%s)\rho^b(%s)\rho^c(%s)\}_{\rm sym}" % (p1, p2, p3)
def dl(p, q): return D + "(%s-%s)" % (p, q)
def dec(C, typ, p, m):
    p1, p2, p3 = p
    if typ == 'ord':
        lhs = C + r"\,\rho^a(%s)\rho^b(%s)\rho^c(%s)" % p
        L = [r"&=" + C + r"\," + S3(*p) + r"\notag\\",
             r"&\quad+\frac i2\," + C + r"\Big[f^{abr}\," + dl(p1, p2) + r"\," + S2('r', m[(1, 2)], 'c', p3) + r"\notag\\",
             r"&\qquad+f^{bcr}\," + dl(p2, p3) + r"\," + S2('a', p1, 'r', m[(2, 3)]) + r"\notag\\",
             r"&\qquad+f^{acr}\," + dl(p1, p3) + r"\," + S2('r', m[(1, 3)], 'b', p2) + r"\Big]\notag\\",
             r"&\quad-\frac16\," + C + r"\big[f^{abr}f^{rcs}-f^{bcr}f^{ras}\big]\notag\\" + "\n" + r"&\qquad\times" + dl(p1, p2) + r"\," + dl(p1, p3) + r"\,\rho^s(%s)" % m['all']]
    elif typ == 'A23':
        lhs = C + r"\,\rho^a(%s)\{\rho^b(%s),\rho^c(%s)\}" % p
        L = [r"&=2\," + C + r"\," + S3(*p) + r"\notag\\",
             r"&\quad+i\," + C + r"\,\Big[f^{abr}\," + dl(p1, p2) + r"\," + S2('r', m[(1, 2)], 'c', p3) + r"\notag\\",
             r"&\qquad+f^{acr}\," + dl(p1, p3) + r"\," + S2('r', m[(1, 3)], 'b', p2) + r"\Big]\notag\\",
             r"&\quad-\frac16\," + C + r"\big[f^{abr}f^{rcs}+f^{acr}f^{rbs}\big]\notag\\" + "\n" + r"&\qquad\times" + dl(p1, p2) + r"\," + dl(p1, p3) + r"\,\rho^s(%s)" % m['all']]
    elif typ == 'A12':
        lhs = C + r"\,\{\rho^a(%s),\rho^b(%s)\}\rho^c(%s)" % p
        L = [r"&=2\," + C + r"\," + S3(*p) + r"\notag\\",
             r"&\quad+i\," + C + r"\Big[f^{bcr}\," + dl(p2, p3) + r"\," + S2('a', p1, 'r', m[(2, 3)]) + r"\notag\\",
             r"&\qquad+f^{acr}\," + dl(p1, p3) + r"\," + S2('r', m[(1, 3)], 'b', p2) + r"\Big]\notag\\",
             r"&\quad+\frac16\," + C + r"\big[f^{bcr}f^{ras}+f^{acr}f^{rbs}\big]\notag\\" + "\n" + r"&\qquad\times" + dl(p1, p2) + r"\," + dl(p1, p3) + r"\,\rho^s(%s)" % m['all']]
    return lhs + "\n" + "\n".join(L)

# (heading, row-equation label, operator sentence, tensor name, tensor definition lines, type, points, merged points, label)
ROWS = [
 (r"Row I.II, $N_{2b}$", "eq:app_IIb",
  r"The colour operator is $\mathcal L^a\big[f^{hbc}\,U^{ah}(w)\,\{\rho^b(x),\rho^c(y)\}-f^{abc}\,U^{bd}(x)\,U^{ce}(y)\,\{\rho^d(x),\rho^e(y)\}\big]=\mathcal C^{abc}_{2b}\,\rho^a(x')\{\rho^b(x),\rho^c(y)\}$, type $\rho\{\rho,\rho\}$, $(p_1,p_2,p_3)=(x',x,y)$, with",
  r"\mathcal C^{abc}_{2b}", [r"\big[U^{ga}(x')-U^{ga}(w')\big]\big[f^{hbc}\,U^{gh}(w)-f^{gde}\,U^{db}(x)\,U^{ec}(y)\big]"],
  'A23', ("x'", "x", "y"), {(1, 2): "x", (1, 3): "y", 'all': "x"}, "app_sym_IIb"),
 (r"Row I.III, structures $C_D$ ($u=z$) and $C_A$ ($u=w$)", "eq:GIIII",
  r"The colour operators are $C_D,C_A=f^{abd}\,\mathcal L^a\,U^{de}(x)\big[U^{bc}(u)-U^{bc}(y)\big]\rho^e(x)\rho^c(y)=\mathcal C^{abc}_{\rm I.III}(u)\,\rho^a(x')\rho^b(x)\rho^c(y)$, ordered, $(p_1,p_2,p_3)=(x',x,y)$, with",
  r"\mathcal C^{abc}_{\rm I.III}(u)", [r"\big[U^{ga}(x')-U^{ga}(w')\big]\,f^{ghd}\,U^{db}(x)\big[U^{hc}(u)-U^{hc}(y)\big]"],
  'ord', ("x'", "x", "y"), {(1, 2): "x", (2, 3): "x", (1, 3): "y", 'all': "x"}, "app_sym_IIII"),
 (r"Row I.V, first part", "eq:app_V1",
  r"The colour operator is $\mathcal L^a\big[f^{adc}\,U^{db}(x)\,U^{ce}(z)\,\{\rho^e(y),\rho^b(y')\}-f^{abc}\,U^{cd}(y)\,U^{be}(y')\,\{\rho^d(y),\rho^e(y')\}\big]=\mathcal C^{abc}_{\rm I.V}\,\rho^a(x')\{\rho^b(y),\rho^c(y')\}$, type $\rho\{\rho,\rho\}$, $(p_1,p_2,p_3)=(x',y,y')$, with",
  r"\mathcal C^{abc}_{\rm I.V}", [r"\big[U^{ga}(x')-U^{ga}(w')\big]\,f^{gde}\big[U^{dc}(x)\,U^{eb}(z)-U^{eb}(y)\,U^{dc}(y')\big]"],
  'A23', ("x'", "y", "y'"), {(1, 2): "y", (1, 3): "y'", 'all': "y"}, "app_sym_IV1"),
 (r"Row II.I, $N_{1a}$", "eq:app_III1a",
  r"The colour operator, the product of the two brackets of Eq.~\eqref{eq:app_III1a}, is $\mathcal C^{abc}_{1a}\,\{\rho^a(x'),\rho^b(y)\}\rho^c(x)$, type $\{\rho,\rho\}\rho$, $(p_1,p_2,p_3)=(x',y,x)$, with",
  r"\mathcal C^{abc}_{1a}", [r"U^{ga}(w')\big[f^{kbc}\,U^{gk}(w)-f^{gkl}\,U^{kb}(z)\,U^{lc}(x)\big]",
                            r"-U^{ga}(x')\,U^{mh}(z)\,U^{mb}(y)\big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\big]"],
  'A12', ("x'", "y", "x"), {(2, 3): "x", (1, 3): "x", 'all': "x"}, "app_sym_III1a"),
 (r"Row II.I, $N_{1b}$", "eq:app_III1b",
  r"The colour operator, the product of the two brackets of Eq.~\eqref{eq:app_III1b}, is $\mathcal C^{abc}_{1b}\,\rho^a(x)\{\rho^b(x'),\rho^c(y)\}$, type $\rho\{\rho,\rho\}$, $(p_1,p_2,p_3)=(x,x',y)$, with",
  r"\mathcal C^{abc}_{1b}", [r"U^{gb}(w)\big[f^{kca}\,U^{gk}(w')-f^{gkl}\,U^{kc}(z)\,U^{la}(x)\big]",
                            r"-U^{gb}(x')\,U^{mh}(z)\,U^{mc}(y)\big[f^{kha}\,U^{gk}(w')-f^{gkl}\,U^{kh}(z)\,U^{la}(x)\big]"],
  'A23', ("x", "x'", "y"), {(1, 2): "x", (1, 3): "x", 'all': "x"}, "app_sym_III1b"),
 (r"Row III.II", "eq:app_IIIII",
  r"The colour operator, the product of the two brackets of Eq.~\eqref{eq:app_IIIII}, is $\mathcal C^{abc}_{\rm III.II}\,\rho^a(y)\rho^b(x')\rho^c(x)$, ordered, $(p_1,p_2,p_3)=(y,x',x)$, with",
  r"\mathcal C^{abc}_{\rm III.II}", [r"U^{ma}(y)\,U^{gb}(x')\,U^{mh}(z)\big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\big]",
                                    r"-U^{gb}(x')\big[f^{kac}\,U^{gk}(w)-f^{gkl}\,U^{ka}(z)\,U^{lc}(x)\big]"],
  'ord', ("y", "x'", "x"), {(1, 2): "y", (2, 3): "x", (1, 3): "x", 'all': "x"}, "app_sym_IIIII"),
 (r"Row III.III, first part, structures with $\mathcal L(v)$, $v=y'$ or $z$", "eq:app_III31",
  r"The colour operators are $f^{c'd'a}\,\mathcal L^{c'd'b}(v)\,\mathcal R_1^{ab}=\mathcal C^{abc}_{\rm III.III}(v)\,\rho^a(x')\{\rho^b(y),\rho^c(x)\}$, type $\rho\{\rho,\rho\}$, $(p_1,p_2,p_3)=(x',y,x)$, with",
  r"\mathcal C^{abc}_{\rm III.III}(v)", [r"f^{mnk}\big[U^{ma}(v)-U^{ma}(x')\big]", r"\times\big[U^{nc}(z)\,U^{kb}(w)-U^{kb}(y)\,U^{nc}(x)\big]"],
  'A23', ("x'", "y", "x"), {(1, 2): "y", (1, 3): "x", 'all': "x"}, "app_sym_IIIIII"),
 (r"Row III.V, structures with $\mathcal L_{\rm V}(v)$, $v=y'$ or $z$", "eq:app_III5",
  r"The colour operators are $f^{c'ca}\,\mathcal L_{\rm V}^{c'}(v)\,\mathcal R_{\rm V}^{ac}=\mathcal C^{abc}_{\rm III.V}(v)\,\rho^a(x')\rho^b(x)\rho^c(y)$, ordered, $(p_1,p_2,p_3)=(x',x,y)$, with",
  r"\mathcal C^{abc}_{\rm III.V}(v)", [r"f^{mkl}\big[U^{ma}(v)-U^{ma}(x')\big]\,U^{kb}(x)\big[U^{lc}(y)-U^{lc}(w)\big]"],
  'ord', ("x'", "x", "y"), {(1, 2): "x", (2, 3): "x", (1, 3): "y", 'all': "x"}, "app_sym_IIIV"),
 (r"Row III.VI, structures with $\mathcal L(v)$, $v=y'$ or $z$", "eq:app_III5",
  r"The colour operators are $f^{c'd'a}\,\mathcal L^{c'd'b}(v)\,\mathcal R_{\rm VI}^{ab}=\mathcal C^{abc}_{\rm III.VI}(v)\,\rho^a(x')\rho^b(y)\rho^c(x)$, ordered, $(p_1,p_2,p_3)=(x',y,x)$, with",
  r"\mathcal C^{abc}_{\rm III.VI}(v)", [r"f^{mnk}\big[U^{ma}(v)-U^{ma}(x')\big]\,U^{kb}(y)\big[U^{nc}(x)-U^{nc}(z)\big]"],
  'ord', ("x'", "y", "x"), {(1, 2): "y", (2, 3): "x", (1, 3): "x", 'all': "x"}, "app_sym_IIIVI"),
]

INTRO = r"""
\subsection{Cubic rows: symmetrised three-\texorpdfstring{$\rho$}{rho} part and commutator terms}\label{app:fin_sym}

Each cubic colour operator of the rows above is a $c$-number tensor times three ordered charges. Before any average over the projectile, and without any assumption on it, the ordered product is rewritten exactly as the totally symmetrised product plus single commutators, which are two-$\rho$ operators, plus double commutators, which are one-$\rho$ operators. Here $\{\dots\}_{\rm sym}$ is the average over all orderings of the operators inside the braces, the commutator is that of Eq.~\eqref{eq:rho}, $[\rho^a(\v x),\rho^b(\v y)]=if^{abc}\rho^c(\v x)\delta^{(2)}(\v x-\v y)$, and the indices $r$, $s$ are summed. For a tensor $C^{abc}$ and points $p_1,p_2,p_3$, Eq.~\eqref{eq:symmetrise} gives, for the three types that occur,
\begin{align}\label{eq:app_symord}
C^{abc}\rho^a(p_1)\rho^b(p_2)\rho^c(p_3)&=C^{abc}\{\rho^a(p_1)\rho^b(p_2)\rho^c(p_3)\}_{\rm sym}\notag\\
&\quad+\frac i2\,C^{abc}\Big[f^{abr}\,\delta^{(2)}(p_1-p_2)\,\{\rho^r(p_1)\rho^c(p_3)\}_{\rm sym}\notag\\
&\qquad+f^{bcr}\,\delta^{(2)}(p_2-p_3)\,\{\rho^a(p_1)\rho^r(p_2)\}_{\rm sym}\notag\\
&\qquad+f^{acr}\,\delta^{(2)}(p_1-p_3)\,\{\rho^r(p_1)\rho^b(p_2)\}_{\rm sym}\Big]\notag\\
&\quad-\frac16\,C^{abc}\big[f^{abr}f^{rcs}-f^{bcr}f^{ras}\big]\notag\\
&\qquad\times\delta^{(2)}(p_1-p_2)\,\delta^{(2)}(p_1-p_3)\,\rho^s(p_1)\,,\\
\label{eq:app_syma23}
C^{abc}\rho^a(p_1)\{\rho^b(p_2),\rho^c(p_3)\}&=2\,C^{abc}\{\rho^a(p_1)\rho^b(p_2)\rho^c(p_3)\}_{\rm sym}\notag\\
&\quad+i\,C^{abc}\Big[f^{abr}\,\delta^{(2)}(p_1-p_2)\,\{\rho^r(p_1)\rho^c(p_3)\}_{\rm sym}\notag\\
&\qquad+f^{acr}\,\delta^{(2)}(p_1-p_3)\,\{\rho^r(p_1)\rho^b(p_2)\}_{\rm sym}\Big]\notag\\
&\quad-\frac16\,C^{abc}\big[f^{abr}f^{rcs}+f^{acr}f^{rbs}\big]\notag\\
&\qquad\times\delta^{(2)}(p_1-p_2)\,\delta^{(2)}(p_1-p_3)\,\rho^s(p_1)\,,\\
\label{eq:app_syma12}
C^{abc}\{\rho^a(p_1),\rho^b(p_2)\}\rho^c(p_3)&=2\,C^{abc}\{\rho^a(p_1)\rho^b(p_2)\rho^c(p_3)\}_{\rm sym}\notag\\
&\quad+i\,C^{abc}\Big[f^{bcr}\,\delta^{(2)}(p_2-p_3)\,\{\rho^a(p_1)\rho^r(p_2)\}_{\rm sym}\notag\\
&\qquad+f^{acr}\,\delta^{(2)}(p_1-p_3)\,\{\rho^r(p_1)\rho^b(p_2)\}_{\rm sym}\Big]\notag\\
&\quad+\frac16\,C^{abc}\big[f^{bcr}f^{ras}+f^{acr}f^{rbs}\big]\notag\\
&\qquad\times\delta^{(2)}(p_1-p_2)\,\delta^{(2)}(p_1-p_3)\,\rho^s(p_1)\,.
\end{align}
Under (A1) and (A2) of Section~\ref{sec:fin_singlet}, $\langle\{\rho^a\rho^b\rho^c\}_{\rm sym}\rangle=S^{abc}$, $\langle\{\rho^a\rho^b\}_{\rm sym}\rangle=\delta^{ab}n/(N_c^2-1)$ and $\langle\rho^s\rangle=0$, which gives Eq.~\eqref{eq:cubicavg}: the first line of each identity yields the $(3\rho)$ equations of Section~\ref{sec:fin_singlet}, the bracket with the $\delta$-functions of two points the $(2\rho)$ (contact) equations, and the last line does not contribute. Below, each cubic row is written in this form with its own tensor. The indices $a,b,c$ are those of the charges in the order written, all other indices are summed, $\mathcal L^a=[U^{ab'}(x')-U^{ab'}(w')]\rho^{b'}(x')$, and the Group~III structures are those of Eq.~\eqref{eq:app_G3col}; in Rows~III.III and~III.VI $U(z)U^T(z)=1$ has been used. At a $\delta$-function the Wilson lines of the tensor are taken at the coincident point. The kernel and prefactor that multiply each colour operator are those of the row equation named in the heading; they are not affected by the decomposition.
"""

def block():
    out = [INTRO]
    for (head, src, sent, C, Cdef, typ, pts, m, lab) in ROWS:
        out.append(r"\paragraph{%s, Eq.~\eqref{%s}.} %s" % (head, src, sent))
        dd = C + r"&\equiv " + Cdef[0]
        for extra in Cdef[1:]:
            dd += r"\notag\\" + "\n" + r"&\qquad" + extra
        out.append(r"\begin{align}" + "\n" + dd + r"\,,\notag\\" + "\n" + dec(C, typ, pts, m) + r"\,.\label{eq:%s}" % lab + "\n" + r"\end{align}")
    return "\n".join(out) + "\n"

if __name__ == '__main__':
    p = sys.argv[1]
    s = open(p).read()
    end = "\\par\\endgroup"
    k = s.find("\n\\subsection{Cubic rows: symmetrised")
    if k < 0:
        k = s.rstrip().rindex(end)
    s = s[:k].rstrip('\n') + "\n" + block() + end + "\n"
    open(p, 'w').write(s)
