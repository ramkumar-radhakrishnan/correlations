# Parts D (what remains, symmetrised, row by row) and E (averaged, colour-singlet form).
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from texlib import *
from singlet_data import CUBIC, QUAD
nl = r"\notag\\"+"\n"
PT = {'x': 'x', 'xp': "x'", 'y': 'y', 'yp': "y'", 'z': 'z', 'w': 'w', 'wp': "w'"}
def c(color, s): return r"{\color{%s}%s}" % (color, s) if color else s
def cnum(z):
    z = complex(z); re, im = round(z.real, 12), round(z.imag, 12)
    def fr(v):
        from fractions import Fraction
        F = Fraction(v).limit_denominator(64)
        if F.denominator == 1: return '' if abs(F.numerator) == 1 else str(abs(F.numerator))
        return r'\frac{%d}{%d}' % (abs(F.numerator), F.denominator)
    if im == 0: return ('-' if re < 0 else '+')+(fr(re) or '1')
    if re == 0: return ('-' if im < 0 else '+')+fr(im)+'i'
    raise ValueError(z)
# ---------------- colour tensors ----------------
TEN = {
 'I.II 2b': r"\Big[U^{ga}(x')-U^{ga}(w')\Big]\Big[f^{hbc}\,U^{gh}(w)-f^{gde}\,U^{db}(x)\,U^{ec}(y)\Big]",
 'I.III C_D (u=z)': r"\Big[U^{ga}(x')-U^{ga}(w')\Big]f^{ghd}\,U^{db}(x)\Big[U^{hc}(z)-U^{hc}(y)\Big]",
 'I.III C_A (u=w)': r"\Big[U^{ga}(x')-U^{ga}(w')\Big]f^{ghd}\,U^{db}(x)\Big[U^{hc}(w)-U^{hc}(y)\Big]",
 'I.V first': r"\Big[U^{ga}(x')-U^{ga}(w')\Big]f^{gde}\Big[U^{dc}(x)\,U^{eb}(z)-U^{eb}(y)\,U^{dc}(y')\Big]",
 'II.I 1a': r"U^{ga}(w')\Big[f^{kbc}\,U^{gk}(w)-f^{gkl}\,U^{kb}(z)\,U^{lc}(x)\Big]-U^{ga}(x')\,U^{mh}(z)\,U^{mb}(y)\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]",
 'II.I 1b': r"U^{gb}(w)\Big[f^{kca}\,U^{gk}(w')-f^{gkl}\,U^{kc}(z)\,U^{la}(x)\Big]-U^{gb}(x')\,U^{mh}(z)\,U^{mc}(y)\Big[f^{kha}\,U^{gk}(w')-f^{gkl}\,U^{kh}(z)\,U^{la}(x)\Big]",
 'III.II': r"U^{ma}(y)\,U^{gb}(x')\,U^{mh}(z)\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]-U^{gb}(x')\Big[f^{kac}\,U^{gk}(w)-f^{gkl}\,U^{ka}(z)\,U^{lc}(x)\Big]",
 "III.III L(y')": r"f^{mnk}\Big[U^{ma}(y')-U^{ma}(x')\Big]\Big[U^{nc}(z)\,U^{kb}(w)-U^{kb}(y)\,U^{nc}(x)\Big]",
 'III.III L(z)': r"f^{mnk}\Big[U^{ma}(z)-U^{ma}(x')\Big]\Big[U^{nc}(z)\,U^{kb}(w)-U^{kb}(y)\,U^{nc}(x)\Big]",
 "III.V L(y')": r"f^{mkl}\Big[U^{ma}(y')-U^{ma}(x')\Big]U^{kb}(x)\Big[U^{lc}(y)-U^{lc}(w)\Big]",
 'III.V L(z)': r"f^{mkl}\Big[U^{ma}(z)-U^{ma}(x')\Big]U^{kb}(x)\Big[U^{lc}(y)-U^{lc}(w)\Big]",
 "III.VI L(y')": r"f^{mnk}\Big[U^{ma}(y')-U^{ma}(x')\Big]U^{kb}(y)\Big[U^{nc}(x)-U^{nc}(z)\Big]",
 'III.VI L(z)': r"f^{mnk}\Big[U^{ma}(z)-U^{ma}(x')\Big]U^{kb}(y)\Big[U^{nc}(x)-U^{nc}(z)\Big]",
}
# ---------------- kernels in the convention pref * kernel (lists of lines; first line, continuation lines) ----------------
KZ = [c('cPr', l) for l in KZ_24]+[c('cPr', '+'+KZ_LE), c('pD', '+'+KZ_E)]
III_Y = [G3_REG, None, G3_INST]   # Group III: regular, (soft remainder), instantaneous; phase e^{-ik(z-y')p+/k+} in the measure
PSI_Y = G3_SOFT
III_Z = [G3_INSTZ]
MEAS_G3 = lambda pts: DPI+r"\int_{%s}" % pts+PH_WY+PH_G3
def g3(pc):
    return {"y": [III_Y[0], c(pc, PSI_Y), III_Y[2]], "z": III_Z}
KT1a = ["-"+KT('w', "w'")[0], KT('w', "w'")[1]]
KT1b = ["+"+KT("w'", 'w')[0], KT("w'", 'w')[1]]
ROWS = [  # group, row, title, pref, measure, cc, structures [(struct name, kernel lines)], extra 2rho (list of (tex, avg tex))
 ('I', 'II', r"-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{64\pi^{6}}\,\frac{1}{k^{+}}", r"\int_{w,w'}"+PH_WW+r"\int_{x,y,x'}", True,
  [('I.II 2b', 'I.II 2b', [R_I2b])]),
 ('I', 'III', r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}", r"\int_{w,w'}"+PH_WW+r"\int_{x,y,z,x'}", True,
  [('I.III', 'I.III C_D (u=z)', KZ), ('I.III', 'I.III C_A (u=w)', [c('pA', KW)])]),
 ('I', 'III*', PF_DN, r"\int_{w,w'}"+PH_WW+r"\int_{x,y,z,x'}", True,
  [('I.III', 'I.III C_D (u=z)', KDZ), ('I.III', 'I.III C_A (u=w)', KDW)]),
 ('I', 'V', PF_I5, DPK+r"\int_{x,z,y,y',x',w'}"+PH_ZW, True,
  [('I.V first', 'I.V first', [KV1_J, c('pB', KV1_E1), c('pB', KV1_E2)])]),
 ('II', 'I', r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}", r"\int_{x,y,z,w,x',w'}"+PH_WW, False,
  [('II.I 1a', 'II.I 1a', KT1a), ('II.I 1b', 'II.I 1b', KT1b)]),
 ('III', 'II', r"-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}", r"\int_{x,y,x',z,w,w'}"+PH_WW, True,
  [('III.II', 'III.II', list(KT('w', "w'")))]),
 ('III', 'III', r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}", MEAS_G3("w,x',y',z,y,x"), True,
  [('III.III first', "III.III L(y')", g3('pB')['y']), ('III.III first', 'III.III L(z)', g3('pB')['z'])]),
 ('III', 'V', r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{4\pi^{4}}\,\frac{1}{k^{+}}", MEAS_G3("w,y',z,y,x,x'"), True,
  [('III.V', "III.V L(y')", g3('pA')['y']), ('III.V', 'III.V L(z)', g3('pA')['z'])]),
 ('III', 'VI', r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{4\pi^{4}}\,\frac{1}{k^{+}}", MEAS_G3("w,y',z,y,x,x'"), True,
  [('III.VI', "III.VI L(y')", g3('pD')['y']), ('III.VI', 'III.VI L(z)', g3('pD')['z'])]),
]
CHI = {'I.II 2b': 'I.II', 'I.III C_D (u=z)': 'I.III,z', 'I.III C_A (u=w)': 'I.III,w', 'I.V first': 'I.V', 'II.I 1a': 'II.I,1a', 'II.I 1b': 'II.I,1b',
       'III.II': 'III.II', "III.III L(y')": "III.III,y'", 'III.III L(z)': 'III.III,z', "III.V L(y')": "III.V,y'", 'III.V L(z)': 'III.V,z',
       "III.VI L(y')": "III.VI,y'", 'III.VI L(z)': 'III.VI,z'}
def Cn(name): return r"C^{abc}_{{\rm %s}}" % CHI[name]
def rho3(typ, p):
    a, b, cc = [PT[q] for q in p]
    if typ == 'ord': return r"\rho^{a}(%s)\,\rho^{b}(%s)\,\rho^{c}(%s)" % (a, b, cc)
    if typ == 'r{rr}': return r"\rho^{a}(%s)\big\{\rho^{b}(%s),\rho^{c}(%s)\big\}" % (a, b, cc)
    return r"\big\{\rho^{a}(%s),\rho^{b}(%s)\big\}\rho^{c}(%s)" % (a, b, cc)
def sym3(p): return r"\big\{\rho^{a}(%s)\rho^{b}(%s)\rho^{c}(%s)\big\}_{\rm sym}" % tuple(PT[q] for q in p)
def d2(a, b): return r"\delta^{(2)}(%s-%s)" % (PT[a], PT[b])
def contact(typ, p, C):
    p1, p2, p3 = p; P1, P2, P3 = [PT[q] for q in p]
    t12 = r"f^{abr}\,%s\big\{\rho^{r}(%s)\rho^{c}(%s)\big\}_{\rm sym}" % (d2(p1, p2), P1, P3)
    t23 = r"f^{bcr}\,%s\big\{\rho^{a}(%s)\rho^{r}(%s)\big\}_{\rm sym}" % (d2(p2, p3), P1, P2)
    t13 = r"f^{acr}\,%s\big\{\rho^{r}(%s)\rho^{b}(%s)\big\}_{\rm sym}" % (d2(p1, p3), P1, P2)
    if typ == 'ord': return [r"\frac{i}{2}\,%s\Big[%s" % (C, t12), "+"+t23, "+"+t13+r"\Big]"]
    if typ == 'r{rr}': return [r"i\,%s\Big[%s" % (C, t12), "+"+t13+r"\Big]"]
    return [r"i\,%s\Big[%s" % (C, t23), "+"+t13+r"\Big]"]
def onerho(typ, p, C):
    p1, p2, p3 = p
    dd = d2(p1, p2)+r"\,"+d2(p1, p3)+r"\,\rho^{s}(%s)" % PT[p1]
    if typ == 'ord': return r"-\frac{1}{6}\,%s\big[f^{abr}f^{rcs}-f^{bcr}f^{ras}\big]%s" % (C, dd)
    if typ == 'r{rr}': return r"-\frac{1}{6}\,%s\big[f^{abr}f^{rcs}+f^{acr}f^{rbs}\big]%s" % (C, dd)
    return r"+\frac{1}{6}\,%s\big[f^{bcr}f^{ras}+f^{acr}f^{rbs}\big]%s" % (C, dd)
def trace(kind, args):
    a = [PT[q] for q in args]
    if kind == 'D': return r"{\rm Tr}\{U^{\dagger}(%s)U(%s)\}" % tuple(a)
    if kind == 'F': return r"{\rm Tr}\{T^{e}U^{\dagger}(%s)U(%s)T^{e}U^{\dagger}(%s)U(%s)\}" % (a[0], a[1], a[1], a[2])
    if kind == 'Q': return r"{\rm Tr}\{T^{e}U^{\dagger}(%s)U(%s)T^{e}U^{\dagger}(%s)U(%s)\}" % tuple(a)
    if kind == 'G': return r"{\rm Tr}\{T^{h}U^{\dagger}(%s)U(%s)D^{h}U^{\dagger}(%s)U(%s)\}" % (a[0], a[1], a[2], a[0])
    if kind == 'Qd': return r"{\rm Tr}\{T^{h}U^{\dagger}(%s)U(%s)D^{h}U^{\dagger}(%s)U(%s)\}" % tuple(a)
def bracket_terms(br, flipQd=False):
    out = []
    for cf, kind, args, ncf in br:
        cf2 = -cf if (flipQd and kind == 'Qd') else cf
        sg = '-' if cf2 < 0 else '+'
        mag = abs(cf2); m = '' if mag == 1 else str(int(mag)) if mag == int(mag) else str(mag)
        out.append(sg+m+(r"N_{c}" if ncf else '')+trace(kind, args))
    out[0] = out[0].lstrip('+')
    return out
def bracket(br, flipQd=False): return ''.join(bracket_terms(br, flipQd))
def chunks(terms, per=2):
    return [''.join(terms[i:i+per]) for i in range(0, len(terms), per)]
def kerlines(K, first_prefix=r"\Big\{"):
    return K
STRUCT_DATA = {}
for row, R in CUBIC.items():
    for st in R['structs']:
        STRUCT_DATA[st[0]] = (R, st)
def eqblock(lhs, pref, measure, terms, cc, label=None, close=''):
    """terms: list of (kernel lines, colour-operator tex lines) ; returns align"""
    L = [lhs+r"&="+pref+measure+(r"\Bigg[" if len(terms) > 1 else "")]
    first = True
    for K, O in terms:
        for i, k in enumerate(K):
            if k is None: continue
            pre = (r"&\quad\times\bigg\{" if first else r"&\quad+\bigg\{") if i == 0 else r"&\qquad"
            L.append(pre+k)
            first = False
        L[-1] += r"\bigg\}"
        for io, o in enumerate(O): L.append((r"&\qquad\times " if io == 0 else r"&\qquad\quad ")+o)
    if len(terms) > 1: L[-1] += r"\Bigg]"
    if cc: L[-1] += r"+{\rm c.c.}"
    L[-1] += close
    body = (r"\notag\\"+"\n").join(L)
    body += (r"\label{%s}" % label) if label else r"\notag"
    return r"\begin{align}"+"\n"+body+"\n"+r"\end{align}"+"\n"
# ---------------- extra (original) two-rho terms of rows, Part D (as in Part B) and averages ----------------
def author_eq(a, b, nph, old, new, label):
    blk, n = phase(L(a, b), nph)
    return blk.replace(old, new)
EXTRA_D = {}
EXTRA_E = {}
FBR_II1 = r"\notag\\ &\quad\times\frac{n(x',x)}{N_{c}^{2}-1}\Bigg[N_{c}{\rm Tr}\{U^{\dagger}(x)U(x')\}+N_{c}{\rm Tr}\{U^{\dagger}(w)U(w')\}\notag\\ &\qquad-{\rm Tr}\{T^{e}U^{\dagger}(x)U(z)T^{e}U^{\dagger}(z)U(w')\}-{\rm Tr}\{T^{e}U^{\dagger}(x')U(z)T^{e}U^{\dagger}(z)U(w)\}\Bigg]"
def II1_second(avg):
    out = []
    for a, b, nm in [(281, 282, '2a'), (285, 287, '2b'), (290, 294, '2c'), (297, 303, '2d'), (306, 307, '2e')]:
        blk, n = phase(L(a, b), 1)
        blk = blk.replace(r"N_{%s\,\rm NLO}" % nm, r"N^{\rm II.I}_{%s}" % nm)
        blk = blk.replace(r"+\frac{(x-w)\cdot\mathcal{K}(x-z)\,\big[(x^{\prime}-z)^{2}-(x^{\prime}-w^{\prime})^{2}\big]}{4\,(z-w^{\prime})^{2}}\notag", r"\notag\\ &\qquad+\frac{(x-w)\cdot\mathcal{K}(x-z)\,\big[(x^{\prime}-z)^{2}-(x^{\prime}-w^{\prime})^{2}\big]}{4\,(z-w^{\prime})^{2}}\notag", 1)
        blk = blk.replace(r"+\frac{(x^{\prime}-w^{\prime})\cdot(z-w^{\prime})\,\big[(x-z)^{2}-(x-w)^{2}\big]}{2\,(z-w)^{2}\,(z-w^{\prime})^{2}}\notag", r"\notag\\ &\qquad+\frac{(x^{\prime}-w^{\prime})\cdot(z-w^{\prime})\,\big[(x-z)^{2}-(x-w)^{2}\big]}{2\,(z-w)^{2}\,(z-w^{\prime})^{2}}\notag", 1)
        if avg:
            i = blk.rindex(r"\Bigg[f^{c^{\prime}bd^{\prime}}")
            h = blk[:i].rstrip()
            if h.endswith(r"\times"): h = h[:-len(r"\times")].rstrip()
            if h.endswith('&'): h = h[:-1].rstrip()
            for t_ in [r"\notag \\", r"\notag\\"]:
                if h.endswith(t_): h = h[:-len(t_)].rstrip()
            blk = h+FBR_II1
            blk = blk.replace(r"\frac{d^{3}N^{\rm II.I}_{%s}}{d^{3}k}" % nm, r"\frac{d^{3}N^{\rm II.I}_{%s}}{d^{3}k}\bigg|_{\rm singlet}" % nm)
        out.append(blk)
    return r"\begin{align}"+"\n"+("\\\\\n".join(out))+r"\notag"+"\n"+r"\end{align}"+"\n"
I2a_D = author_eq(246, 247, 1, r"\frac{d^{3}N_{1\,\rm NLO}}{d^{3}k} & =", r"\frac{d^{3}N^{\rm I.II}_{1}}{d^{3}k} & =", None)
partB = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'partB.tex')).read()
def grab_eq(label):
    j = partB.index(r"\label{%s}" % label); i = partB.rindex(r"\begin{align}", 0, j); k = partB.index(r"\end{align}", j)
    return partB[i:k+len(r"\end{align}")]+"\n"
LO_REPL = r"\frac{n(x',x)}{N_{c}^{2}-1}\,{\rm Tr}\Big\{\big(U^{\dagger}(y')-U^{\dagger}(x')\big)\big(U(y)-U(x)\big)\Big\}"
# ---------------- generate ----------------
D = []; Ep = []
D.append(r"""\section{What remains, with every three-$\rho$ product symmetrised: Group by Group, Row by Row}\label{sec:D}
\paragraph{What remains after the cancellations of Sec.~\ref{sec:C}.} (1) Group I Row IV and Group III Row I are removed (exact). (2) The $\frac{11}{3}$ terms of Group I Rows II and V are absent (exact). (3) All other terms of Sec.~\ref{sec:B} remain: Group I Rows II, III, V; Group II Rows I, IV; Group III Rows II, III, V, VI. Among them, the soft remainders of the pairs A--D (coloured as in Sec.~\ref{sec:C}) remain as terms of their rows; they are integrable only in the pair combinations Eqs.~\eqref{eq:C-A}--\eqref{eq:C-D}, i.e.\ each Group I member has to be integrated together with its Group III partner, under the renaming $y'\to w'$ (pair D: and $x\leftrightarrow y$). (4) The 24-term kernel of Group I Row III and the part $\ell-\int_{0}^{k^{+}}dp^{+}\big[\frac{1}{p^{+}}\big]_{+}e^{-ik\cdot(z-w)p^{+}/k^{+}}$ of the third term of $N_{3b}$ (blue) are integrable only with the $w$ integral done before the $z$ integral, Eq.~\eqref{eq:C-Mhat}. (5) Group II Row IV, Group II Row I (second part) and Group III Row III (second part) have to be integrated together over the position of the gluon pair (joint tail, Sec.~\ref{sec:C6}). (6) In the $i\epsilon$ prescription the additional term $\Delta N$ of Group I Row III, Eq.~\eqref{eq:B-DN}, remains, without its Coulomb logarithm (cancelled exactly against the complex conjugate, Sec.~\ref{sec:C7}); it is listed after Row III. Nothing is dropped in this section; only the order of the charges is rewritten.

\paragraph{Symmetrisation (exact, no average, no model).} With $\{\dots\}_{\rm sym}$ the average over all orderings of the operators inside the braces, $\{AB\}_{\rm sym}=\frac12(AB+BA)$, and $[\rho^{a}(x),\rho^{b}(y)]=if^{abc}\rho^{c}(x)\,\delta^{(2)}(x-y)$, every ordered product of three charges with a $c$-number tensor $C^{abc}$ is, for the three types that occur (paper draft Eqs.~(app\_symord)--(app\_syma12)),
\begin{align}
C^{abc}\rho^{a}(p_{1})\rho^{b}(p_{2})\rho^{c}(p_{3})&=\underbrace{C^{abc}\{\rho^{a}(p_{1})\rho^{b}(p_{2})\rho^{c}(p_{3})\}_{\rm sym}}_{\rm (i)}\notag\\
&\quad+\underbrace{\frac{i}{2}C^{abc}\Big[f^{abr}\delta^{(2)}(p_{1}-p_{2})\{\rho^{r}(p_{1})\rho^{c}(p_{3})\}_{\rm sym}+f^{bcr}\delta^{(2)}(p_{2}-p_{3})\{\rho^{a}(p_{1})\rho^{r}(p_{2})\}_{\rm sym}}_{\rm (ii)}\notag\\
&\quad\underbrace{+f^{acr}\delta^{(2)}(p_{1}-p_{3})\{\rho^{r}(p_{1})\rho^{b}(p_{2})\}_{\rm sym}\Big]}_{\rm (ii)}\notag\\
&\quad\underbrace{-\frac{1}{6}C^{abc}\big[f^{abr}f^{rcs}-f^{bcr}f^{ras}\big]\delta^{(2)}(p_{1}-p_{2})\,\delta^{(2)}(p_{1}-p_{3})\,\rho^{s}(p_{1})}_{\rm (iii)}\,,\label{eq:D-ord}\\
C^{abc}\rho^{a}(p_{1})\{\rho^{b}(p_{2}),\rho^{c}(p_{3})\}&=2C^{abc}\{\rho^{a}(p_{1})\rho^{b}(p_{2})\rho^{c}(p_{3})\}_{\rm sym}\notag\\
&\quad+i\,C^{abc}\Big[f^{abr}\delta^{(2)}(p_{1}-p_{2})\{\rho^{r}(p_{1})\rho^{c}(p_{3})\}_{\rm sym}+f^{acr}\delta^{(2)}(p_{1}-p_{3})\{\rho^{r}(p_{1})\rho^{b}(p_{2})\}_{\rm sym}\Big]\notag\\
&\quad-\frac{1}{6}C^{abc}\big[f^{abr}f^{rcs}+f^{acr}f^{rbs}\big]\delta^{(2)}(p_{1}-p_{2})\,\delta^{(2)}(p_{1}-p_{3})\,\rho^{s}(p_{1})\,,\label{eq:D-r12}\\
C^{abc}\{\rho^{a}(p_{1}),\rho^{b}(p_{2})\}\rho^{c}(p_{3})&=2C^{abc}\{\rho^{a}(p_{1})\rho^{b}(p_{2})\rho^{c}(p_{3})\}_{\rm sym}\notag\\
&\quad+i\,C^{abc}\Big[f^{bcr}\delta^{(2)}(p_{2}-p_{3})\{\rho^{a}(p_{1})\rho^{r}(p_{2})\}_{\rm sym}+f^{acr}\delta^{(2)}(p_{1}-p_{3})\{\rho^{r}(p_{1})\rho^{b}(p_{2})\}_{\rm sym}\Big]\notag\\
&\quad+\frac{1}{6}C^{abc}\big[f^{bcr}f^{ras}+f^{acr}f^{rbs}\big]\delta^{(2)}(p_{1}-p_{2})\,\delta^{(2)}(p_{1}-p_{3})\,\rho^{s}(p_{1})\,.\label{eq:D-r21}
\end{align}
(i) is the symmetrised three-$\rho$ term, (ii) the commutator (contact, two-$\rho$) terms, (iii) the double-commutator (one-$\rho$) terms. At a $\delta$-function the Wilson lines of the tensor are at the coincident point. For every row below: the colour tensor $C^{abc}$ (the indices $a,b,c$ belong to the charges in the order written, all other colour indices summed) is written out; equation (i) is the symmetrised three-$\rho$ term; equation (ii) contains the contact terms of the row's three-$\rho$ products together with the row's original two-$\rho$ terms; the one-$\rho$ terms (iii) have the prefactor, integrals, phase and kernels of (i), with $t\,C^{abc}\{\rho^{a}\rho^{b}\rho^{c}\}_{\rm sym}$ replaced by the last line of Eq.~\eqref{eq:D-ord}, \eqref{eq:D-r12} or \eqref{eq:D-r21} (they have zero average over a colour-invariant projectile). That each row's tensor, contracted with the charges in the row's order, is the colour operator of the row as written in Sec.~\ref{sec:B}, and that (i)+(ii)+(iii) equals it, was checked at the operator level (Sec.~\ref{sec:F}). The two-$\rho$ products of the original two-$\rho$ terms are kept ordered (if wanted: $\rho^{a}(p)\rho^{b}(q)=\{\rho^{a}(p)\rho^{b}(q)\}_{\rm sym}+\frac{i}{2}f^{abr}\delta^{(2)}(p-q)\rho^{r}(p)$).
""")
Ep.append(r"""\section{Average over the projectile: colour-singlet form, Group by Group, Row by Row}\label{sec:E}
\noindent\textbf{Author's choices (2026-11-20).} Global colour invariance of the projectile gives $\langle\rho^{a}(x)\rangle=0$ and
\begin{align}
\big\langle\{\rho^{a}(x)\rho^{b}(y)\}_{\rm sym}\big\rangle=\big\langle\rho^{a}(x)\rho^{b}(y)\big\rangle_{x\neq y}=\delta^{ab}\,\frac{n(x,y)}{N_{c}^{2}-1}\,,\qquad n(x,y)=n(y,x)\,;\notag
\end{align}
the symmetrised three-point function is of $d$-type,
\begin{align}
\big\langle\{\rho^{a}(x)\rho^{b}(y)\rho^{c}(z)\}_{\rm sym}\big\rangle=d^{abc}\,\frac{N_{c}}{(N_{c}^{2}-4)(N_{c}^{2}-1)}\,o(x,y,z)\,,\qquad o(x,y,z)\equiv\big\langle d^{efg}\{\rho^{e}(x)\rho^{f}(y)\rho^{g}(z)\}_{\rm sym}\big\rangle\,,\notag
\end{align}
with the $f^{abc}$ part $g_{f}=0$; no model for $o$ ($o$ is real and totally symmetric; it is written with its arguments in the order of the charges); for $SU(2)$, $d^{abc}=0$ and the $o$ terms are absent. Hence, from Eqs.~\eqref{eq:D-ord}--\eqref{eq:D-r21} (the one-$\rho$ terms average to zero), with $C^{abc}$ the row's tensor,
\begin{align}
\langle C^{abc}\rho^{a}(p_{1})\rho^{b}(p_{2})\rho^{c}(p_{3})\rangle&=\frac{N_{c}}{(N_{c}^{2}-4)(N_{c}^{2}-1)}\,C^{abc}d^{abc}\,o(p_{1},p_{2},p_{3})\notag\\
&\quad+\frac{i}{2(N_{c}^{2}-1)}\Big[\Phi_{12}\,\delta^{(2)}(p_{1}-p_{2})\,n(p_{1},p_{3})+\Phi_{23}\,\delta^{(2)}(p_{2}-p_{3})\,n(p_{1},p_{2})-\Phi_{13}\,\delta^{(2)}(p_{1}-p_{3})\,n(p_{1},p_{2})\Big],\notag
\end{align}
and for the anticommutator types the factor $2$ on the first term and $i[\Phi_{12}\delta_{12}n-\Phi_{13}\delta_{13}n]$, resp.\ $i[\Phi_{23}\delta_{23}n-\Phi_{13}\delta_{13}n]$; here $\Phi_{ij}=C^{abc}f^{abc}$ with the Wilson lines at $p_{i}$ and $p_{j}$ equal. Every colour factor is written below in colour-singlet form with traces of adjoint Wilson lines ($U^{\dagger}=U^{T}$, $(T^{e})_{bc}=-if^{ebc}$, $(D^{h})_{bc}=d^{hbc}$, ${\rm Tr}\,T^{a}T^{b}=N_{c}\delta^{ab}$):
\begin{align}
&{\rm Tr}\{U^{\dagger}(a)U(b)\}\,,\quad{\rm Tr}\{T^{e}U^{\dagger}(a)U(b)T^{e}U^{\dagger}(b)U(c)\}=f^{ijk}f^{lmn}U^{il}(a)U^{jm}(b)U^{kn}(c)\,,\quad{\rm Tr}\{T^{e}U^{\dagger}(a)U(b)T^{e}U^{\dagger}(c)U(d)\}\,,\notag\\
&{\rm Tr}\{T^{h}U^{\dagger}(a)U(b)D^{h}U^{\dagger}(c)U(a)\}=i\,f^{ijk}d^{lmn}U^{il}(a)U^{jm}(b)U^{kn}(c)\,,\quad{\rm Tr}\{T^{h}U^{\dagger}(a)U(b)D^{h}U^{\dagger}(c)U(e)\}\notag
\end{align}
(the traces with $D^{h}$ are purely imaginary, C-odd and vanish for $SU(2)$). In the fundamental representation ($U^{ab}=2{\rm tr}(t^{a}Vt^{b}V^{\dagger})$; fix\_dabc, exact):
\begin{align}
-i\,{\rm Tr}\{T^{h}U^{\dagger}(a)U(b)D^{h}U^{\dagger}(c)U(a)\}&={\rm Im}\big[{\rm tr}(V^{\dagger}_{a}V_{b})\,{\rm tr}(V^{\dagger}_{b}V_{c})\,{\rm tr}(V^{\dagger}_{c}V_{a})\big]+{\rm Im}\,{\rm tr}\big(V^{\dagger}_{a}V_{b}V^{\dagger}_{c}V_{a}V^{\dagger}_{b}V_{c}\big)\,,\notag\\
i\,{\rm Tr}\{T^{h}U^{\dagger}(a)U(b)D^{h}U^{\dagger}(c)U(e)\}&=-{\rm Im}\big[{\rm tr}A\,{\rm tr}B\,{\rm tr}(A^{\dagger}B^{\dagger})\big]-{\rm Im}\,{\rm tr}\big(ABA^{\dagger}B^{\dagger}\big)\,,\qquad A=V^{\dagger}_{a}V_{b}\,,\ B=V^{\dagger}_{c}V_{e}\,.\notag
\end{align}
Each row below has the prefactor, integrals, phase and kernels of the same row in Sec.~\ref{sec:D}; only the colour operator is averaged. The $(3\rho)$ equations are the part of the finite terms odd in $k$ (configuration by configuration; they change sign with the phase convention, Sec.~\ref{sec:conv}), the $(2\rho)$ equations the even part; for the additional term $\Delta N$ of the $i\epsilon$ prescription it is the reverse (its $(3\rho)$ part is even, its contact part odd in $k$). The singlet colour factors of $\Delta N$ are those of Group I Row III (same colour tensors), with the kernels of $\Delta N$. These are the row equations of the paper draft (Eqs.~(sIIa)--(oIIIVI)) in the notation of this file, with the colour factors written as traces.
""")
GROUP_TITLES = {'I': r'Group I (with c.c.)', 'II': r'Group II', 'III': r'Group III (with c.c.)'}
order = [('I', 'I'), ('I', 'II'), ('I', 'III'), ('I', 'III*'), ('I', 'IV'), ('I', 'V'), ('II', 'I'), ('II', 'II'), ('II', 'III'), ('II', 'IV'),
         ('III', 'I'), ('III', 'II'), ('III', 'III'), ('III', 'IV'), ('III', 'V'), ('III', 'VI')]
rowmap = {(g, r): (g, r, pref, meas, cc, sts) for g, r, pref, meas, cc, sts in ROWS}
NA = {('I', 'I'): 'four-$\\rho$ row, evolution terms only; \\textsc{NA}', ('I', 'IV'): 'cancels Group III Row I exactly (Sec.~\\ref{sec:C1}); nothing remains',
      ('II', 'II'): '\\textsc{NA}', ('II', 'III'): '\\textsc{NA}', ('III', 'I'): 'cancels Group I Row IV exactly (Sec.~\\ref{sec:C1}); nothing remains', ('III', 'IV'): '\\textsc{NA}'}
curg = None
for g, r in order:
    if g != curg:
        D.append(r"\subsection{%s}" % GROUP_TITLES[g]+"\n"); Ep.append(r"\subsection{%s}" % GROUP_TITLES[g]+"\n"); curg = g
    ttl = r"Row III, $i\epsilon$ prescription ($\eta=+1$): the additional term $\Delta N$ of Eq.~\eqref{eq:B-DN}" if r == 'III*' else "Row %s" % r
    D.append(r"\subsubsection*{%s}" % ttl+"\n"); Ep.append(r"\subsubsection*{%s}" % ttl+"\n")
    if (g, r) in NA:
        D.append(NA[(g, r)]+"\n"); Ep.append(NA[(g, r)]+"\n"); continue
    tag = '%s.%s' % (g, r.replace('*', ''))
    NN = r'\Delta N' if r == 'III*' else 'N'
    lt = 'dN' if r == 'III*' else ''
    odd3, odd2 = ('(even in $k$)', '(odd in $k$)') if r == 'III*' else ('(odd in $k$)', '(even in $k$)')
    if (g, r) in rowmap:
        _, _, pref, meas, cc, sts = rowmap[(g, r)]
        # tensors
        tl = []
        for rowkey, sname, K in sts:
            Rd, st = STRUCT_DATA[sname]
            T_ = TEN[sname]
            if len(T_) > 160 and r"\Big]-U" in T_:
                i_ = T_.index(r"\Big]-U")+len(r"\Big]")
                T_ = T_[:i_]+r"\notag\\"+"\n"+r"&\qquad"+T_[i_:]
            tl.append(Cn(sname)+r"&\equiv "+T_)
        D.append(r"Colour tensor%s (the row's colour operator is $%s$%s):" % ('s' if len(sts) > 1 else '', Cn(sts[0][1])+rho3(STRUCT_DATA[sts[0][1]][0]['typ'], STRUCT_DATA[sts[0][1]][0]['pts']), ' and the same with the other tensor' if len(sts) > 1 and sts[0][0] == sts[1][0] else (', and $'+Cn(sts[1][1])+rho3(STRUCT_DATA[sts[1][1]][0]['typ'], STRUCT_DATA[sts[1][1]][0]['pts'])+'$' if len(sts) > 1 else ''))+"\n"
                 + r"\begin{align}"+"\n"+(r"\notag\\"+"\n").join(tl)+r"\notag"+"\n"+r"\end{align}"+"\n")
        # (i)
        terms = []
        for rowkey, sname, K in sts:
            Rd, st = STRUCT_DATA[sname]
            tt = '2' if Rd['t'] == 2 else ''
            terms.append((K, [tt+Cn(sname)+sym3(Rd['pts'])]))
        D.append(r"(i) Symmetrised three-$\rho$ term:"+"\n"+eqblock(r"\frac{d^{3}%s^{\rm %s}}{d^{3}k}\bigg|^{\rm (i)}_{3\rho}" % (NN, tag), pref, meas, terms, cc, label='eq:D-%s%s-3' % (tag, lt)))
        # (ii)
        terms = []
        for rowkey, sname, K in sts:
            Rd, st = STRUCT_DATA[sname]
            terms.append((K, contact(Rd['typ'], Rd['pts'], Cn(sname))))
        D.append(r"(ii) Two-$\rho$ terms: the contact terms of the symmetrisation"+(" and the row's original two-$\\rho$ terms" if (g, r) in [('I', 'II'), ('I', 'V'), ('II', 'I'), ('III', 'III')] else '')+":\n"
                 + eqblock(r"\frac{d^{3}%s^{\rm %s}}{d^{3}k}\bigg|^{\rm (ii)}_{2\rho,\ {\rm contact}}" % (NN, tag), pref, meas, terms, cc, label='eq:D-%s%s-2' % (tag, lt)))
        # (3rho) averaged
        terms = []
        for rowkey, sname, K in sts:
            Rd, st = STRUCT_DATA[sname]
            lab, contacts, c3tex, c3, br3, oargs = st
            coef = c3/Rd['A']*(-1j)
            ch = chunks(bracket_terms(br3, flipQd=True))
            O = [r"\frac{N_{c}}{(N_{c}^{2}-4)(N_{c}^{2}-1)}\,"+(r"\Big(%s\Big)" % cnum(coef))+r"\,o(%s)\Big[" % ','.join(PT[q] for q in oargs)+ch[0]]+ch[1:]
            O[-1] += r"\Big]"
            terms.append((K, O))
        Ep.append(r"Symmetrised three-$\rho$ part %s:" % odd3+"\n"+eqblock(r"\frac{d^{3}%s^{\rm %s}}{d^{3}k}\bigg|^{(3\rho)}_{\rm singlet}" % (NN, tag), pref, meas, terms, cc, label='eq:E-%s%s-3' % (tag, lt)))
        terms = []
        for rowkey, sname, K in sts:
            Rd, st = STRUCT_DATA[sname]
            lab, contacts, c3tex, c3, br3, oargs = st
            O = []
            for k_, (pair, fr, to, ctex, cc_, nargs, br) in enumerate(contacts):
                gco = cc_/Rd['A']
                ch = chunks(bracket_terms(br))
                O.append((r"\frac{1}{N_{c}^{2}-1}\bigg[" if k_ == 0 else '')+r"\Big(%s\Big)" % cnum(gco)+d2(fr, to)+r"\,n(%s)\Big[" % ','.join(PT[q] for q in nargs)+ch[0])
                O.extend(ch[1:]); O[-1] += r"\Big]"
            O[-1] += r"\bigg]"
            terms.append((K, O))
        Ep.append(r"Contact (two-$\rho$) part of the three-$\rho$ terms %s:" % odd2+"\n"+eqblock(r"\frac{d^{3}%s^{\rm %s}}{d^{3}k}\bigg|^{(2\rho)}_{\rm singlet,\ contact}" % (NN, tag), pref, meas, terms, cc, label='eq:E-%s%s-2' % (tag, lt)))
    # original two-rho terms
    if (g, r) == ('I', 'II'):
        D.append(r"and the original two-$\rho$ term of the row (unchanged):"+"\n"+r"\begin{align}"+"\n"+I2a_D+r"\label{eq:D-I.II-2o}"+"\n"+r"\end{align}"+"\n")
        e = I2a_D.replace(r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[U^{ab}(w)-U^{ab}(x)\Big]\rho^{b'}(x')\,\rho^{b}(x)",
                          r"\frac{n(x',x)}{N_{c}^{2}-1}\Big[{\rm Tr}\{U^{\dagger}(x')U(w)\}-{\rm Tr}\{U^{\dagger}(x')U(x)\}-{\rm Tr}\{U^{\dagger}(w')U(w)\}+{\rm Tr}\{U^{\dagger}(w')U(x)\}\Big]+{\rm c.c.}")
        assert e != I2a_D
        e = e.replace(r"\frac{d^{3}N^{\rm I.II}_{1}}{d^{3}k} & =", r"\frac{d^{3}N^{\rm I.II}_{1}}{d^{3}k}\bigg|_{\rm singlet} & =")
        Ep.append(r"Original two-$\rho$ term:"+"\n"+r"\begin{align}"+"\n"+e+r"\label{eq:E-I.II-2o}"+"\n"+r"\end{align}"+"\n")
    if (g, r) == ('I', 'V'):
        eq = grab_eq('eq:B-I5two').replace(r"\label{eq:B-I5two}", r"\label{eq:D-I.V-2o}")
        eq = eq.replace(r"\frac{d^{3}N^{\rm I.V}_{2}}{d^{3}k} &", r"\frac{d^{3}N^{\rm I.V}_{2}}{d^{3}k} &")
        S1 = r"-"+PLK+r"\frac{"+EP('x-z')+r"}{(y-z)^{2}}\Big[-\frac{(y-z)\cdot\mathcal{K}(x'-w')}{(x-z)^{2}}+\frac{1}{2}\,(y-z)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-x)\Big]"
        S2 = r"-"+PLKQ+r"\frac{"+EP('x-z')+r"}{(y-x)^{2}}\Big[-\frac{(y-x)\cdot\mathcal{K}(x'-w')}{(x-z)^{2}}-\frac{1}{2}\,(y-x)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-z)\Big]"
        assert eq.count(r"&\qquad"+S1) == 1 and eq.count(r"&\qquad"+S2+r"\Bigg\}") == 1
        eq = eq.replace(r"&\qquad"+S1, r"&\qquad{\color{pC}"+S1+"}")
        eq = eq.replace(r"&\qquad"+S2+r"\Bigg\}", r"&\qquad{\color{pC}"+S2+r"}\Bigg\}")
        D.append(r"and the original two-$\rho$ terms of the row (second part, Eq.~\eqref{eq:B-I5two}; soft remainders of pair C in violet):"+"\n"+eq)
        e = eq.replace(r"\label{eq:D-I.V-2o}", r"\label{eq:E-I.V-2o}").replace(r"\frac{d^{3}N^{\rm I.V}_{2}}{d^{3}k} &", r"\frac{d^{3}N^{\rm I.V}_{2}}{d^{3}k}\bigg|_{\rm singlet} &")
        old1 = r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[-f^{adc}\,f^{fbe}\,U^{db}(x)\,U^{ce}(z)+N_{c}\,U^{af}(y)\Big]\rho^{b'}(x')\,\rho^{f}(y)"
        new1 = (r"\frac{n(x',y)}{N_{c}^{2}-1}\Big[N_{c}{\rm Tr}\{U^{\dagger}(x)U(x')\}-N_{c}{\rm Tr}\{U^{\dagger}(x)U(w')\}\notag\\ &\qquad\qquad-{\rm Tr}\{T^{e}U^{\dagger}(x)U(x')T^{e}U^{\dagger}(x')U(z)\}+{\rm Tr}\{T^{e}U^{\dagger}(x)U(z)T^{e}U^{\dagger}(z)U(w')\}"
                r"\notag\\"+"\n"+r"&\qquad\qquad+N_{c}\Big(-{\rm Tr}\{U^{\dagger}(x)U(x')\}+{\rm Tr}\{U^{\dagger}(x)U(w')\}+{\rm Tr}\{U^{\dagger}(y)U(x')\}-{\rm Tr}\{U^{\dagger}(y)U(w')\}\Big)\Big]")
        old2 = r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[U^{af}(y)-U^{af}(x)\Big]\rho^{b'}(x')\,\rho^{f}(y)"
        new2 = r"\frac{n(x',y)}{N_{c}^{2}-1}\Big[-{\rm Tr}\{U^{\dagger}(x)U(x')\}+{\rm Tr}\{U^{\dagger}(x)U(w')\}+{\rm Tr}\{U^{\dagger}(y)U(x')\}-{\rm Tr}\{U^{\dagger}(y)U(w')\}\Big]"
        assert old1 in e and old2 in e
        e = e.replace(old1, new1).replace(old2, new2).replace(r",\label{eq:E-I.V-2o}", r"+{\rm c.c.}\,,\label{eq:E-I.V-2o}")
        Ep.append(r"Original two-$\rho$ terms (second part):"+"\n"+e)
    if (g, r) == ('II', 'I'):
        D.append(r"and the original two-$\rho$ terms of the row (second part, unchanged; to be integrated together with Group II Row IV and Group III Row III, second part, Sec.~\ref{sec:C6}):"+"\n"+II1_second(False))
        Ep.append(r"Original two-$\rho$ terms (second part):"+"\n"+II1_second(True))
    if (g, r) == ('II', 'IV'):
        eq = grab_eq('eq:B-II4')
        D.append(r"All terms are two-$\rho$ terms (no symmetrisation needed): the DGLAP and RC terms of Secs.~3, 4 and"+"\n"+eq.replace(r"\label{eq:B-II4}", r"\label{eq:D-II.IV}").replace(r"\label{eq:B-II4rem}", r"\label{eq:D-II.IVrem}"))
        e = eq.replace(r"\label{eq:B-II4}", r"\label{eq:E-II.IV}").replace(r"\label{eq:B-II4rem}", r"\label{eq:E-II.IVrem}")
        old = r"\Big[U^{ce'}(y') - U^{ce'}(x')\Big]\Big[U^{ce}(y) - U^{ce}(x)\Big]\rho^{e'}(x')\rho^{e}(x)"
        assert old in e
        e = e.replace(old, LO_REPL)
        old = r"\Big[U^{ce'}(z)-U^{ce'}(x')\Big]\Big[U^{ce}(z)-U^{ce}(x)\Big]\rho^{e'}(x')\rho^{e}(x)"
        assert old in e
        e = e.replace(old, r"\frac{n(x',x)}{N_{c}^{2}-1}\,{\rm Tr}\Big\{\big(U^{\dagger}(z)-U^{\dagger}(x')\big)\big(U(z)-U(x)\big)\Big\}")
        e = e.replace(r"\frac{d^{3}N^{\rm II.IV}_{\rm fin}}{d^{3}k}= &", r"\frac{d^{3}N^{\rm II.IV}_{\rm fin}}{d^{3}k}\bigg|_{\rm singlet}= &")
        Ep.append(r"All terms are quadratic; the DGLAP, RC (Secs.~3, 4) and LO+RC terms in singlet form follow by the replacement $[U^{ce'}(y')-U^{ce'}(x')][U^{ce}(y)-U^{ce}(x)]\rho^{e'}(x')\rho^{e}(x)\to"+LO_REPL+r"$; the finite terms are"+"\n"+e)
    if (g, r) == ('III', 'III'):
        eq = grab_eq('eq:B-III3two').replace(r"\label{eq:B-III3two}", r"\label{eq:D-III.III-2o}")
        PSI_LINE = r"&-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{4\pi^{4}}\,\frac{1}{k^{+}}"+DPI+PL3+r"\int_{w,x',y',z,x}"+PH_WY+PH_G3+r"\,\frac{1}{(x-w)^{2}}"
        PSI_L2 = r"&\qquad\times\Big[(x-w)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(z-w)\cdot\mathcal{K}(y'-z)+\frac{1}{2}\,\mathcal{K}(x-z)\cdot\mathcal{K}(y'-z)\,(x-w)\cdot\mathcal{K}(x'-y')\Big]"
        assert PSI_LINE in eq and PSI_L2 in eq
        eq = eq.replace(PSI_LINE, r"&{\color{pC}"+PSI_LINE[1:]+"}").replace(PSI_L2, r"&\qquad{\color{pC}"+PSI_L2[len(r"&\qquad"):]+"}")
        D.append(r"and the original two-$\rho$ terms of the row (second part, Eq.~\eqref{eq:B-III3two}; soft remainder of pair C in violet):"+"\n"+eq)
        e = eq.replace(r"\label{eq:D-III.III-2o}", r"\label{eq:E-III.III-2o}").replace(r"\frac{d^{3}N^{\rm III.III}_{2}}{d^{3}k} &", r"\frac{d^{3}N^{\rm III.III}_{2}}{d^{3}k}\bigg|_{\rm singlet} &")
        LY = r"\Big[f^{c'd'a}\,U^{c'e'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]"
        LZ = r"\Big[f^{c'd'a}\,U^{c'e'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]"
        LZa = r"\left[f^{c'd'a}\,U^{c'e'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\, U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]"
        R2 = r"\Big[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\Big]"
        R2a = r"\left[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x) -f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\right]"
        BY = r"\frac{n(x',x)}{N_{c}^{2}-1}\Big[-N_{c}{\rm Tr}\{U^{\dagger}(x)U(x')\}+N_{c}{\rm Tr}\{U^{\dagger}(x)U(y')\}\notag\\ &\qquad\qquad+{\rm Tr}\{T^{e}U^{\dagger}(x')U(z)T^{e}U^{\dagger}(z)U(w)\}-{\rm Tr}\{T^{e}U^{\dagger}(y')U(z)T^{e}U^{\dagger}(z)U(w)\}\Big]"
        BZ = r"\frac{n(x',x)}{N_{c}^{2}-1}\Big[-N_{c}{\rm Tr}\{U^{\dagger}(x)U(x')\}+N_{c}{\rm Tr}\{U^{\dagger}(x)U(z)\}\notag\\ &\qquad\qquad-N_{c}{\rm Tr}\{U^{\dagger}(z)U(w)\}+{\rm Tr}\{T^{e}U^{\dagger}(x')U(z)T^{e}U^{\dagger}(z)U(w)\}\Big]"
        LYR2 = LY+r"\notag\\"+"\n"+r"&\qquad\times"+R2
        assert e.count(LYR2) == 2, e.count(LYR2)
        e = e.replace(LYR2, BY)
        assert e.count(LY) == 1
        e = e.replace(LY, BY)
        LZR2 = LZ+r"\bigg\}\notag\\"+"\n"+r"&\qquad\times"+R2
        assert e.count(LZR2) == 1
        e = e.replace(LZR2, BZ+r"\bigg\}")
        old = LZa+r"\notag\\ &\qquad\times\left[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x) -f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\right]"
        assert e.count(old) == 1, 'klog'
        e = e.replace(old, BZ)
        assert 'f^{cbd}' not in e, e[e.index('f^{cbd}')-300:e.index('f^{cbd}')+50]
        e = e.replace(r",\label{eq:E-III.III-2o}", r"+{\rm c.c.}\,,\label{eq:E-III.III-2o}")
        Ep.append(r"Original two-$\rho$ terms (second part):"+"\n"+e)
    if (g, r) in rowmap and (g, r) not in [('I', 'II'), ('I', 'V'), ('II', 'I'), ('III', 'III')]:
        pass
open('partD.tex', 'w').write(''.join(D))
open('partE.tex', 'w').write(''.join(Ep))
print('ok')
