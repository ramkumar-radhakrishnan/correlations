r"""Finite_terms_1_1_corrected.tex: your log(vee/k+) list with corrections, red = remove, blue = add.
Built from input/Finite_terms_1_1.tex by explicit edits (every blue term is transcribed in corrected_lvk.py and
checked there against the exact coefficient of log(vee/k+) of the same row of Finite_terms_1.tex), plus the
JIMWLK bracket in colour-singlet form from ../check_v4/Evolution_terms_corrected.tex (Eq. (3.4))."""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'input', 'Finite_terms_1_1.tex')
EC = os.path.join(HERE, '..', 'check_v4', 'Evolution_terms_corrected.tex')
OUT = os.path.join(HERE, 'Finite_terms_1_1_corrected.tex')

LV = r"\log{\left[\frac{\vee}{k^{+}}\right]}"
NC = r"\nc "          # start of a blue cell inside an align

# ------------------------------------------------------------------ colour factors and brackets of your file
COL_III = (r"\left[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\, \rho^{b'}(x')\right]\times\, \left[U^{de}(x)\,U^{bc}(z)\, "
           r"\rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]")
COL_IV = (r"\left[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\,\rho^{b'}(x')\right] \left[f^{cdc'}\,U^{ac}(w)\,U^{bd}(z)\,"
          r"U^{be}(y)\,\rho^{e}(y)\rho^{c'}(x)-f^{abc}\,U^{be}(y)\,U^{cd}(x)\,\rho^{e}(y)\rho^{d}(x)\right]")
COL_A = (r"\Bigg[U^{ac'}(w')\,\{\rho^{c'}(x'),\rho^{b}(y)\}-U^{bc'}(z)\,U^{ad'}(x')\,U^{c'e'}(y)\, \{\rho^{d'}(x'),"
         r"\rho^{e'}(y)\}\Bigg]\times \Bigg[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\, U^{bc}(z)\,U^{de}(x)\,\rho^{e}(x)\Bigg]")
COL_B = (r"\Bigg[f^{c'bd'}\,U^{ac'}(w')\,\rho^{d'}(x)-f^{ac'd'}\, U^{bc'}(z)\,U^{d'e'}(x)\,\rho^{e'}(x)\Bigg] \times "
         r"\Bigg[U^{ac}(w)\,\{\rho^{c}(x'),\rho^{b}(y)\}-U^{bc}(z)\,U^{ad}(x')\,U^{ce}(y)\, \{\rho^{d}(x'),\rho^{e}(y)\}\Bigg]")
O58_1 = (r"\Bigg[f^{c'bd'}\,f^{cbd}\, U^{ac'}(w')\,U^{ac}(w)\,\rho^{d'}(x')\rho^{d}(x)- f^{c'bd'}\,f^{acd}\,U^{ac'}(w')\,"
         r"U^{bc}(z)\,U^{de}(x)\,\rho^{d'}(x')\rho^{e}(x)")
O58_2 = (r"-f^{ac'd'}\,f^{cbd}\,U^{ac}(w)\,U^{bc'}(z)\,U^{e'd'}(x')\,\rho^{e'}(x')\,\rho^{d}(x) +N_{c}\,U^{e'd}(x')\,"
         r"U^{de}(x)\,\rho^{e'}(x')\,\rho^{e}(x)\Bigg]")
COL_G3I_1 = (r"\Bigg[U^{ab'}(x')\,U^{c'd'}(y)\,U^{c'b}(z)\,\rho^{b'}(x')\rho^{d'}(y)-U^{bd'}(z)\,U^{d'e'}(y)U^{ac'}(w')\,"
             r"\rho^{c'}(x')\rho^{e'}(y)\Bigg]\times")
COL_R = r"\Bigg[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\, U^{bc}(z)\,U^{de}(x)\,\rho^{e}(x)\Bigg]"
COL_G3II_1 = (r"\Bigg[U^{c'd'}(y)\,U^{ab'}(x')\,U^{bc'}(z)\,\rho^{d'}(y)\rho^{b'}(x')-U^{ab'}(x')\,\rho^{b}(y)"
              r"\rho^{b'}(x')\Bigg]\times")
COL_LO = r"\left[U^{e'c}(y') - U^{e'c}(x')\right] \times\, \left[U^{ec}(y) - U^{ce}(x)\right]\, \rho^{e'}(x') \rho^{e}(x)"
COL_LOz = r"\left[U^{e'c}(z) - U^{e'c}(x')\right] \times\, \left[U^{ec}(z) - U^{ce}(x)\right]\, \rho^{e'}(x') \rho^{e}(x)"
COL_C1 = (r"\left[U^{e'c'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]\times")
COL_C1b = (r"\left[U^{ac}(w)\, \{\rho^{c}(y),\rho^{b}(x)\}-U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\{\rho^{d}(y),\rho^{e}(x)\}\right]")
COL_V = r"\left[U^{e'c'}(y')\,\rho^{e'}(x')-U^{c'e'}(x')\,\rho^{e'}(x')\right]\times"
COL_Vb = r"\left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right]"
COL_VIb = r"\left[U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\rho^{d}(y)\,\rho^{e}(x)-U^{ad}(y)\,\rho^{d}(y)\,\rho^{b}(x)\right]"
COL_32 = (r"\left[f^{c'd'a}\,U^{e'c'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\, U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')"
          r"\right]\times")
COL_32b = r"\left[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x) -f^{acd}\,U^{bc}(z)\,U^{de}(x)\,\rho^{e}(x)\right]"
T2 = r"\frac{(y'-z)\cdot (x-z)\,(x'-y')\cdot (y-w)}{(y'-z)^{2}\,(x'-y')^{2}\,(y-w)^{2}\,(x-z)^{2}}"


def tij(w):
    """your bracket T^{ij} (Group I Row IV, Group II Row I, Group III Rows I, II) at the point w"""
    return (r"\Bigg[\frac{\delta^{ij}}{2\,(z-W)^{2}}+\frac{(x-z)^{j}\,(z-W)^{i}}{(z-W)^{2}\,(x-z)^{2}}"
            r"-\frac{(x-z)^{j}\,(x-W)^{i}}{(x-W)^{2}\,(x-z)^{2}}-\frac{(x-W)^{i}\,(z-W)^{j}}{(z-W)^{2}\,(x-W)^{2}}\Bigg]"
            ).replace('W', w)


def kiv_term(sign, pref, xp_w, w, colour_lines):
    """a full blue K_IV term: sign (1/(2pi)^3) int e^{-ik(w'-w)} pref (1/k+) L_v X'^i Y^j T^{ij}(w) x colour"""
    s = (r"%s\frac{1}{(2\pi)^{3}}\,\int_{x,y,z,w,x',w'}e^{-ik\cdot (w^\prime-w)}\,%s\,\frac{1}{k^{+}}\,%s\,"
         r"\frac{(x'-%s)^{i}\,(y-z)^{j}}{(x'-%s)^{2}\,(y-z)^{2}}\notag \\ & %s\times %s" % (sign, pref, LV, xp_w, xp_w, NC, tij(w)))
    for cl in colour_lines:
        s += r"\notag \\ & %s%s" % (NC, cl)
    return s


def blue_lines(s):
    """make every cell of an align body blue"""
    return s.replace('&', '& ' + NC)


def strip_marks(t):
    """\new{a} -> a, \old{a} -> '' (balanced braces)"""
    out, i = [], 0
    while i < len(t):
        m = re.compile(r"\\(new|old)\{").match(t, i)
        if not m:
            out.append(t[i]); i += 1
            continue
        j, depth = m.end(), 1
        while depth:
            depth += {'{': 1, '}': -1}.get(t[j], 0)
            j += 1
        if m.group(1) == 'new':
            out.append(strip_marks(t[m.end():j - 1]))
        i = j
    return ''.join(out)


def jimwlk_singlet():
    L = open(EC).read().split('\n')
    i = next(k for k, l in enumerate(L) if l.startswith('Color singlet form'))
    j = next(k for k in range(i, len(L)) if L[k].strip() == r'\end{align}')
    body = '\n'.join(L[i + 1:j + 1])
    body = strip_marks(body)
    body = body.replace(r"\frac{d^{3}N_{\mathrm{NLO}}}{d^{3}k}", r"\frac{d^{3}N_{\rm IV}}{d^{3}k}")
    body = body.replace(r"\log\Big[\frac{\vee}{\Lambda}\Big]", LV)
    assert LV in body
    return body


# ------------------------------------------------------------------ the edits
def build():
    s = open(SRC).read()
    parts = re.split(r'(\\begin\{align\}.*?\\end\{align\})', s, flags=re.S)
    aligns = [p for p in parts if p.startswith(r'\begin{align}')]
    assert len(aligns) == 16
    new = {}

    # 1. G1 Row II part 1: coefficient 2
    a = aligns[0]
    old = r"\frac{g^{4}\,N_{c}}{32\pi^{4}}\,\frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]}\,"
    assert a.count(old) == 1
    new[0] = a.replace(old, r"\frac{g^{4}\,N_{c}}{32\pi^{4}}\,\frac{1}{k^{+}}\,\old{1\times}\new{2\times}"
                            r"\log{\left[\frac{\vee}{k^{+}}\right]}\,")
    # 2. G1 Row II 2a: spurious int_z
    a = aligns[1]
    assert a.count(r"\int_{x,y,x',z}") == 2
    new[1] = a.replace(r"\int_{x,y,x',z}", r"\int_{x,y,x'\old{,z}}")

    # 3. G1 Row III N1: your term red, the blue term
    a = aligns[2]
    old = (r"-\,\frac{ig^{4}\,f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]} \,\int_{x,y,z,x'}\, "
           r"\frac{(x'-w')\cdot (y-z)}{(x'-w')^{2}\,(y-z)^{2}}\,\frac{\left[(x-z)^{2}+(x-w)^{2}\right]}{\left[(x-z)^{2}"
           r"-(x-w)^{2}\right]^{2}}")
    assert a.count(old) == 1, 'N1'
    new[2] = a.replace(old, r"\old{%s}\notag \\ & %s-\,\frac{ig^{4}\,f^{abd}}{8\pi^{5}}\,\frac{1}{k^{+}}\,%s\,"
                            r"\int_{x,y,z,x'}\frac{(x'-w')\cdot (x-w)\, (x-z)\cdot (y-z)}{(x'-w')^{2}\, (x-w)^{2}\, "
                            r"(x-z)^{2}\, (y-z)^{2}}" % (old, NC, LV))

    # 4. G1 Row III N2: six blue terms inside the braces
    a = aligns[3]
    tail = r"\Bigg\} " + "\n" + r"\end{align}"
    assert a.endswith(r"\Bigg\} " + "\n" + r"\end{align}") or a.rstrip().endswith(r"\end{align}")
    den = r"(x'-w')^{2}\, (y-z)^{2}\, (z-w)^{2}\,\left[(x-w)^{2}-(x-z)^{2}\right]"
    six = [
        r"\frac{(x'-w')\cdot (z-w)\, (y-z)\cdot (x-w)}{%s}" % den,
        r"-\frac{(x'-w')\cdot (x-w)\, (y-z)\cdot (z-w)\,(x-z)^{2}}{%s\,(x-w)^{2}}" % den,
        r"+\frac{(x'-w')\cdot (y-z)\, (z-w)\cdot (x-w)\,(x-z)^{2}}{%s^{2}}" % den.replace(r"\right]", r"\right]"),
        r"-\frac{(x'-w')\cdot (z-w)\, (y-z)\cdot (x-z)\,(x-w)^{2}}{%s\,(x-z)^{2}}" % den,
        r"+\frac{(x'-w')\cdot (x-z)\, (y-z)\cdot (z-w)}{%s}" % den,
        r"-\frac{(x'-w')\cdot (y-z)\, (z-w)\cdot (x-z)\,(x-w)^{2}}{%s^{2}}" % den,
    ]
    blk = r"\notag \\ & %s+\frac{ig^{4}\, f^{abd}}{8\pi^{5}}\, \frac{1}{k^{+}}\,%s\,\int_{x,y,z,x'}\Bigg[%s" % (NC, LV, six[0])
    for t in six[1:]:
        blk += r"\notag \\ & %s\qquad%s" % (NC, t)
    blk += r"\Bigg]\notag \\ & %s%s" % (NC, COL_III)
    i = a.rfind(r"\Bigg\}")
    new[3] = a[:i] + blk + a[i:]

    # 6. G1 Row IV: blue T^{ij} term
    a = aligns[5]
    i = a.rfind(r"\end{align}")
    new[5] = a[:i] + r"\notag \\ & %s" % NC + kiv_term('+', r"\frac{ig^{4}}{8\pi^{5}}", "w'", 'w', [COL_IV]) + "\n" + a[i:]

    # 11. G2 Row I part 1: blue terms for blocks A and B
    a = aligns[10]
    i = a.rfind(r"\end{align}")
    new[10] = (a[:i] + r"\notag \\ & %s" % NC + kiv_term('+', r"\Big(-\frac{ig^{4}}{16\pi^{5}}\Big)", "w'", 'w', [COL_A])
               + r"\notag \\ & %s" % NC + kiv_term('+', r"\frac{ig^{4}}{16\pi^{5}}", "w", "w'", [COL_B]) + "\n" + a[i:])

    # 12-14. G2 Row I 2d, 2e, 2f: red
    for k in (11, 12, 13):
        new[k] = "{\\color{red}" + aligns[k] + "}"

    # 15, 16. G3 Rows I, II
    a = aligns[14]
    i = a.rfind(r"\end{align}")
    new[14] = (a[:i] + r"\notag \\ & %s" % NC
               + kiv_term('+', r"\Big(-\frac{ig^{4}}{8\pi^{5}}\Big)", "w'", 'w', [COL_G3I_1, COL_R]) + "\n" + a[i:])
    a = aligns[15]
    i = a.rfind(r"\end{align}")
    new[15] = (a[:i] + r"\notag \\ & %s" % NC
               + kiv_term('+', r"\Big(-\frac{ig^{4}}{8\pi^{5}}\Big)", "w'", 'w', [COL_G3II_1, COL_R]) + "\n" + a[i:])

    # reassemble
    out, n = [], 0
    for p in parts:
        if p.startswith(r'\begin{align}'):
            out.append(new.get(n, p))
            n += 1
        else:
            out.append(p)
    s = ''.join(out)
    return s


# ------------------------------------------------------------------ blue blocks inserted as text
NOTE_II = r"""
\new{\textbf{Row II.} Part 1: in \texttt{Finite\_terms\_1.tex} the pole bracket multiplies $\log\frac{\vee}{k^+}+\log\big(\frac{\vee}{k^+}-1\big)=2\log\frac{\vee}{k^+}+\log\big(1-\frac{k^+}{\vee}\big)$, so the coefficient is $2$. Part 2a: correct (nothing depends on $z$, so $\int_z$ is dropped).}
"""

NOTE_III = r"""
\new{\textbf{Row III.} The $\log\frac{\vee}{k^+}$ of this row also sits in three logarithms of \texttt{Finite\_terms\_1.tex} that you did not count:}
\begin{align*}
&\new{\log\Big|\frac{\vee\,(x-z)^{2}-k^{+}(x-w)^{2}}{(\Lambda+k^{+})(x-z)^{2}-k^{+}(x-w)^{2}}\Big|,\qquad \log\Big|\frac{k^{+}(x-w)^{2}-\vee\,(x-z)^{2}}{k^{+}(x-w)^{2}-\Lambda\,(x-z)^{2}}\Big|\quad(d^3N_1),}\\
&\new{\log\Big|\frac{(x-w)^{2}k^{+}-(x-z)^{2}\vee}{(x-w)^{2}k^{+}-(x-z)^{2}(k^{+}+\Lambda)}\Big|\quad(d^3N_2);\qquad\text{each}=\log\frac{\vee}{k^+}+O(1)\text{ for }\vee\to\infty.}
\end{align*}
\new{In $d^3N_1$ the first one cancels your term exactly (same kernel, opposite sign); the second one gives the blue term. In $d^3N_2$ the third one gives the six blue terms. $d^3N_{3a}$ is correct; $d^3N_{3b}$ has no $\vee$.}
"""

NOTE_IV = r"""
\new{\textbf{Row IV.} The blue term is the $\log\frac{\vee}{k^+}$ inside the first logarithm of the row in \texttt{Finite\_terms\_1.tex}:}
\begin{equation*}
\new{\log\Big[\frac{(\vee-k^{+})(x-z)^{2}+k^{+}(x-w)^{2}}{\Lambda(x-z)^{2}+k^{+}(x-w)^{2}}\Big]=\log\frac{\vee}{k^{+}}+\log\Big[\Big(1-\frac{k^+}{\vee}\Big)\frac{(x-z)^2}{(x-w)^2}+\frac{k^+}{\vee}\Big]+O\Big(\frac{\Lambda}{k^+}\Big).}
\end{equation*}
\new{Together, your term and the blue term have the kernel (with $\delta^{ij}$ the $\frac{1}{(z-w)^2}$ parts cancel)}
\begin{equation*}
\new{\frac{ig^4}{8\pi^5}\frac{1}{k^+}\Bigg[\frac{(x'-w')\cdot(z-w)\,(y-z)\cdot(x-z)}{(x'-w')^2(z-w)^2(y-z)^2(x-z)^2}-\frac{(x'-w')\cdot(x-w)\,(y-z)\cdot(x-z)}{(x'-w')^2(x-w)^2(y-z)^2(x-z)^2}-\frac{(x'-w')\cdot(x-w)\,(y-z)\cdot(z-w)}{(x'-w')^2(x-w)^2(y-z)^2(z-w)^2}\Bigg],}
\end{equation*}
\new{which is integrable at $z\to w$. Your term alone, $\propto\frac{1}{(z-w)^2}$, is not: its $z$ integral diverges, so the three-$\rho$ terms could not cancel. The same holds for Group II Row I part 1 and Group III Rows I, II.}
"""

NOTE_V = r"""
\new{\textbf{Row V.} All parts correct.}
"""

NOTE_G2I1 = r"""
\new{\textbf{Row I, part 1.} The blue terms are the $\log\frac{\vee}{k^+}$ inside $\log\Big[\frac{(\vee-k^{+})(x-z)^{2}+k^{+}(x-w)^{2}}{\Lambda(x-z)^{2}+k^{+}(x-w)^{2}}\Big]$ and its partner with $w\to w'$ (as in Group I Row IV).}
"""

G2I2_INTRO = r"""
\new{\textbf{Row I, part 2.} Remove 2d, 2e, 2f: in \texttt{Finite\_terms\_1.tex} their explicit $\log\frac{\vee}{k^++\Lambda}$ cancels exactly against the $\log\frac{\vee}{k^+}$ inside $\log\Big[\frac{(\vee-k^{+})(x-z)^{2}+k^{+}(x-w)^{2}}{\Lambda(x-z)^{2}+k^{+}(x-w)^{2}}\Big]$ and $\log\Big[\frac{(\vee-k^{+})(x'-z)^{2}+k^{+}(x'-w')^{2}}{\Lambda(x'-z)^{2}+k^{+}(x'-w')^{2}}\Big]$ in the same rows. The $\log\frac{\vee}{k^+}$ of part 2 comes from 2b and 2c, through the same two logarithms (2a: they cancel between numerator and denominator):}
"""

DEN2b = r"(x-w)^{2}\,(x'-w')^{2}"
DEN2c = r"(x-z)^{2}\,(x'-z)^{2}"
G2I2_BLUE = (r"""{\color{blue}
\begin{align}
\frac{d^{3}N_{2b\,\rm NLO}^{\rm}}{d^{3}k}  & = -\frac{1}{(2\pi)^{3}}\, \frac{g^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}\,LV\,\int_{z,w,w^{\prime},x,x'}\, e^{-ik(w^{\prime}-w)}\, \frac{1}{DEN2b} \notag \\
& \Bigg[\frac{(x'-w')\cdot (x-w)\, (z-w)\cdot (z-w')}{(z-w)^{2}\, (z-w')^{2}}+\frac{(x'-w')\cdot (x-w)\, (x'-z)\cdot (z-w)}{2\, (z-w)^{2}\, (x'-z)^{2}}\notag \\
& +\frac{(x-w)\cdot (x'-w')\, (x-z)\cdot (z-w')}{2\, (x-z)^{2}\, (z-w')^{2}}+\frac{(x-z)\cdot (x'-z)\, (x'-w')\cdot (x-w)}{4\, (x-z)^{2}\, (x'-z)^{2}}\Bigg]\notag \\
& O58_1\notag \\
& O58_2
\end{align}
\begin{align}
\frac{d^{3}N_{2c\,\rm NLO}^{\rm}}{d^{3}k}  & = \frac{1}{(2\pi)^{3}}\, \frac{g^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}\,LV\,\int_{z,w,w^{\prime},x,x'}\, e^{-ik(w^{\prime}-w)}\, \frac{1}{DEN2c} \notag \\
& \Bigg[\frac{(x'-z)\cdot (x-z)\, (z-w)\cdot (z-w')}{(z-w)^{2}\, (z-w')^{2}}-\frac{(x'-z)\cdot (x-z)\, (x'-w')\cdot (z-w)}{2\, (z-w)^{2}\, (x'-w')^{2}}\notag \\
& -\frac{(x-w)\cdot (z-w')\, (x-z)\cdot (x'-z)}{2\, (x-w)^{2}\, (z-w')^{2}}+\frac{(x-z)\cdot (x'-z)\, (x'-w')\cdot (x-w)}{4\, (x-w)^{2}\, (x'-w')^{2}}\Bigg]\notag \\
& O58_1\notag \\
& O58_2
\end{align}}
""").replace('LV', LV).replace('DEN2b', DEN2b).replace('DEN2c', DEN2c).replace('O58_1', O58_1).replace('O58_2', O58_2)

G2I2_NOTE = r"""
\new{The first term of $d^3N_{2c}$ has the kernel $\mathcal K(z-w)\cdot\mathcal K(z-w')\,\frac{(x-z)\cdot(x'-z)}{(x-z)^2(x'-z)^2}$: the measured gluon is emitted by the soft gluon $z$ in the amplitude and in the conjugate amplitude. This is the kernel of the second line of your Eq.~(LOBFKL); no term of your list had it.}
"""

G2IV = (r"""\old{\textsc{NA}}

{\color{blue}
\begin{align}
\frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & = -\,\frac{g^{2}}{4\pi^{2}}\, \frac{N_{c}}{(2\pi)^{3}}\,\frac{g^{2}}{2\pi^{2}}\, \frac{1}{k^{+}}\,LV\,\int_{x,y,x',y'}\,e^{-ik\cdot (y'-y)}\,\frac{(x-y)\cdot (x'-y')}{(x-y)^{2}\,(x'-y')^{2}}\notag \\
& \left[-\frac{d_{\perp}}{\epsilon}-d_{\perp}\,\gamma_{E}-\frac{d_{\perp}}{2}\,\log{\left(\frac{(y'-y)^{2}\,\mu^{2}_{\overline{MS}}}{4}\right)}-1\right]\,COL_LO\notag \\
& +\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,LV\,\int_{x,x',z,w,w'}\,e^{-ik\cdot (w'-w)}\,\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}\,(w'-z)^{2}}\,\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}\,(x'-z)^{2}}\notag \\
& COL_LOz
\end{align}}
\new{\textbf{Row IV.} Both terms come from the lower end $\xi\to\frac{k^+}{\vee}$ of $\int_{k^+/\vee}^1d\xi$ in \texttt{Finite\_terms\_1.tex}:}
\begin{itemize}
\item \new{first term: your subtraction is $\frac{1}{\xi^2}\frac{\xi}{(1-\xi)_+}e^{-ik\cdot(y'-y)\frac1\xi}=\frac{e^{-ik\cdot(y'-y)\frac1\xi}-e^{-ik\cdot(y'-y)}}{\xi(1-\xi)}$. With $\frac{1}{\xi(1-\xi)}=\frac1\xi+\frac{1}{1-\xi}$, the part $-\frac{e^{-ik\cdot(y'-y)}}{\xi}$ gives $-e^{-ik\cdot(y'-y)}\int_{k^+/\vee}^1\frac{d\xi}{\xi}=-e^{-ik\cdot(y'-y)}\log\frac{\vee}{k^+}$; the rest has no $\log\frac{\vee}{k^+}$. Together with Group I Row II part 1 (coefficient 2) and Row V 2b UV, its $\frac1\epsilon$ cancels: no pole is left in this list.}
\item \new{second term: the term $\frac{1-\xi}{\xi}$ of $P_{gg}$ at small $\xi$. There $e^{-ik\cdot(y'-y)/\xi}$ probes momentum $\frac{k}{\xi}$, and $\int d^2(y'-y)\,e^{-iq\cdot(y'-y)}\log(y'-y)^2\,h(y'-y)\to-\frac{4\pi}{q^2}h(0)$ puts $y=y'=z$ and gives a factor $\xi^2$, so $\int_{k^+/\vee}\frac{d\xi}{\xi}=\log\frac{\vee}{k^+}$. The factor $\frac{4\pi^2}{k^2}$ is written as $\int_{w,w'}e^{-ik\cdot(w'-w)}\mathcal K(w-z)\cdot\mathcal K(w'-z)$. This term is BFKL: the measured gluon is emitted by the gluon at $z$ on both sides.}
\end{itemize}
""").replace('LV', LV).replace('COL_LOz', COL_LOz).replace('COL_LO', COL_LO)


def crow_block(coef, f, cola, colb, measure):
    t = r"""\old{\textsc{NA}}

{\color{blue}
\begin{align}
\frac{d^{3}N_{\rm NLO}}{d^{3}k} & = \frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,@@F@@}{@@COEF@@}\, \frac{1}{k^{+}}\,@@LV@@\,\int_{@@MEAS@@}\, e^{-ik\cdot (y'-w)}\,@@T2@@\notag \\
& @@COLA@@\notag \\
& @@COLB@@
\end{align}}
"""
    for k, v in (('@@F@@', f), ('@@COEF@@', coef), ('@@LV@@', LV), ('@@MEAS@@', measure), ('@@T2@@', T2),
                 ('@@COLA@@', cola), ('@@COLB@@', colb)):
        t = t.replace(k, v)
    return t


G3_III1 = crow_block(r"16\pi^{5}", "f^{c'd'a}", COL_C1, COL_C1b, "w,x',y',z,y,x")
G3_V = crow_block(r"8\pi^{5}", "f^{c'ca}", COL_V, COL_Vb, "w,y',z,y,x,x'")
G3_VI = crow_block(r"8\pi^{5}", "f^{c'd'a}", COL_C1, COL_VIb, "w,y',z,y,x,x'")
G3_III2 = (r"""\old{\textsc{NA}}

{\color{blue}
\begin{align}
\frac{d^{3}N_{2\,\rm NLO}}{d^{3}k} & = \frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}\,LV\,\int_{w,x',y',z,x}\, e^{-ik\cdot (y'-w)}\,\bigg[\frac{(x-w)\cdot (x'-y')\,(z-w)\cdot (y'-z)}{(z-w)^{2}\,(x-w)^{2}\,(y'-z)^{2}\,(x'-y')^{2}}\notag \\
& +\frac{(x-z)\cdot (y'-z)\,(x-w)\cdot (x'-y')}{2\,(x-z)^{2}\,(x-w)^{2}\,(y'-z)^{2}\,(x'-y')^{2}}\bigg]\, COL32\notag \\
& COL32b
\end{align}}
""").replace('LV', LV).replace('COL32b', COL_32b).replace('COL32', COL_32)

NOTE_G3C = r"""
\new{\textbf{Group III Rows III, V, VI.} The blue terms come from the $\left[\frac{1}{p^+}\right]_+\delta_{jm}\delta_{ik'}$ term of \texttt{Finite\_terms\_1.tex}:}
\begin{equation*}
\new{\int_0^{\vee-k^+}\!dp^+\left[\frac{1}{p^+}\right]_+e^{-ik\cdot(y'-z)\frac{p^+}{k^+}}=-\Big[\log\frac{\vee}{k^+}+\log\Big(1-\frac{k^+}{\vee}\Big)+\gamma_E+\log\big(ik\cdot(y'-z)\big)+E_1\Big(ik\cdot(y'-z)\frac{\vee-k^+}{k^+}\Big)\Big]}
\end{equation*}
\new{(the same way in part 2, where $\frac{1}{p^+(x-z)^2+k^+(x-w)^2}\to\frac{1}{k^+(x-w)^2}$ at $p^+=0$). The kernels are those of the $\log(1-\frac{k^+}{\vee})$ lines of the same rows of \texttt{Finite\_terms\_1.tex}, with the opposite sign. The other terms of these rows give no $\log\frac{\vee}{k^+}$ (the oscillating $E_1$, $e^{-ik\cdot(y'-z)(\vee-k^+)/k^+}$ terms vanish after $\int_{y'}$).}
"""

RESULT = r"""
\section*{\new{Result}}
\new{Red: remove. Blue: add. Every blue term is the coefficient of $\log\frac{\vee}{k^+}$ in a logarithm of the same row of \texttt{Finite\_terms\_1.tex} that contains $\vee$ and that this list had left out (the notes after each row say which one). Every row of this list was checked against $\partial/\partial\log\vee$ of the same row of \texttt{Finite\_terms\_1.tex} at $\frac{\vee}{k^+}=10^7\dots10^9$: agreement to @@WORST@@ or better (\texttt{corrected\_check.json}).}

\paragraph{\new{Changes.}}
\begin{itemize}
\item \new{Group I Row II part 1: coefficient $2$, not $1$.}
\item \new{Group I Row III: $d^3N_1$ replaced; six terms added to $d^3N_2$.}
\item \new{Group I Row IV, Group II Row I part 1, Group III Rows I, II: the $\mathbb T^{ij}$ term added.}
\item \new{Group II Row I part 2: 2d, 2e, 2f removed (they add up to zero); 2b, 2c added.}
\item \new{Group II Row IV, Group III Row III (both parts), Rows V, VI: marked NA, but they contain $\log\frac{\vee}{k^+}$; added.}
\item \new{Group I Row II 2a, Row V (all parts): correct.}
\end{itemize}

\paragraph{\new{What the corrected list gives.}}
\begin{enumerate}
\item \new{\textbf{Three $\rho$'s} cancel after symmetrization (checked: \texttt{LogVeeKplus\_BFKL.pdf}, all 39 structures; \texttt{BFKL\_terms\_check.pdf}, Sec.~5). Without the blue $\mathbb T^{ij}$ terms they cannot: the kernel $\frac{1}{(z-w)^2}$ of your Group I Row IV, Group II Row I part 1 and Group III Rows I, II is not integrable at $z\to w$ (note after Group I Row IV).}
\item \new{\textbf{Two $\rho$'s} add up to}
\begin{equation}
\new{\sum_{\rm Groups\ I-III}\frac{d^3N}{d^3k}=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\log\Big[\frac{\vee}{k^{+}}\Big]\int_{x,x',z,w,w'}e^{-ik(w'-w)}\Big[\mathcal B-\mathcal J\Big],}
\end{equation}
\new{with $\mathcal B$ = BFKL$\otimes$LO (your Eq.~(LOBFKL)) and $\mathcal J$ = JIMWLK$\otimes$LO, the brackets of Secs.~4 and 3 of \texttt{Evolution\_terms\_corrected.tex}.}
\item \new{\textbf{Why $-\mathcal J$.} In every row the soft gluon is integrated over $\Lambda<p^+<\vee$. With $T$ = the integrand for $p^+\ll k^+$ and $P$ = the integrand for $p^+\gg k^+$,}
\begin{equation*}
\new{\int_\Lambda^\vee\frac{dp^+}{p^+}\,G\Big(\frac{p^+}{k^+}\Big)=T\log\frac{k^+}{\Lambda}+P\log\frac{\vee}{k^+}+\dots=T\log\frac{\vee}{\Lambda}+(P-T)\log\frac{\vee}{k^+}+\dots}
\end{equation*}
\begin{itemize}
\item \new{$T\log\frac{\vee}{\Lambda}$ is in your $\log\frac{\vee}{\Lambda}$ list, and $\sum T=\mathcal J$;}
\item \new{$\sum P=\mathcal B$;}
\item \new{so this list is $\sum(P-T)=\mathcal B-\mathcal J$.}
\end{itemize}
\item \new{\textbf{Both lists together:}}
\begin{equation}
\new{\boxed{\log\Big[\frac{\vee}{\Lambda}\Big]\mathcal J+\log\Big[\frac{\vee}{k^+}\Big]\big(\mathcal B-\mathcal J\big)=\log\Big[\frac{k^+}{\Lambda}\Big]\mathcal J+\log\Big[\frac{\vee}{k^+}\Big]\mathcal B.}}
\label{eq:both}
\end{equation}
\new{BFKL comes only with $\log\frac{\vee}{k^+}$, and JIMWLK only with $\log\frac{k^+}{\Lambda}$, as I said before. The $-\mathcal J$ of this list is not a new claim: your JIMWLK terms are written with $\log\frac{\vee}{\Lambda}=\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}$, and the $-\mathcal J\log\frac{\vee}{k^+}$ here removes the second piece.}
\item \new{\textbf{Group IV} (Sec.~4) moves $\log\frac{\vee}{k^+}\,\mathcal J$ from your $\log\frac{\vee}{\Lambda}$ list to this one. Then this list is exactly your Eq.~(LOBFKL) with $\delta Y=\log\frac{\vee}{k^+}$ (Sec.~5).}
\item \new{\textbf{At mid rapidity}, $k^+=\sqrt{\vee\Lambda}$: $\log\frac{\vee}{k^+}=\log\frac{k^+}{\Lambda}=\frac12\log\frac{\vee}{\Lambda}$. So Eq.~\eqref{eq:both} becomes $\frac12\log\frac{\vee}{\Lambda}(\mathcal J+\mathcal B)$, the form in \texttt{Evolution\_terms\_corrected.tex}. It is the same equation.}
\end{enumerate}
"""

LOBFKL = r"""
\begin{align}
    \frac{dN}{d^2 w\, d^2 w'\, dy}= - &\frac{\alpha_{s}^2 N_c  \old{\delta Y}\new{\log\big[\frac{\vee}{k^+}\big]} }{8\pi^{6} (N_c^2-1)}
    \int_{x,y,z} \langle \rho^f(x) \rho^f(y)\rangle
    \notag \\ &
    \Bigg[
    {\cal K} (x-w) \cdot {\cal K} (y-w') K_D(x,y,z) {\rm Tr} \left \{  (U^\dagger(x)-U^\dagger(w)) (U(y)-U(w')) \right\}
    \notag \\ &
    -2\, {\cal K} (z-w) \cdot {\cal K} (z-w') K(x,y,z) {\rm Tr} \left \{  (U^\dagger(z)-U^\dagger(w)) (U(z)-U(w')) \right\}
    \notag \\ &
    +
    {\cal K} (z-w) \cdot {\cal K} (y-w') (2K(x,y,z) - K(x,x,z)) {\rm Tr} \left \{  (U^\dagger(z)-U^\dagger(w)) (U(y)-U(w')) \right\}
    \notag \\ &
    +
    {\cal K} (x-w) \cdot {\cal K} (z-w') (2K(x,y,z) - K(y,y,z)) {\rm Tr} \left \{  (U^\dagger(x)-U^\dagger(w)) (U(z)-U(w')) \right\}
    \Bigg]
\label{eq:LOBFKL}
\end{align}
"""

GROUP4 = r"""
\section{\new{Group IV: the $\log\frac{\vee}{k^+}$ part of your $\log\frac{\vee}{\Lambda}$ list}}
\new{Your $\log\frac{\vee}{\Lambda}$ terms add up to $\frac{1}{(2\pi)^3}\frac{g^4}{16\pi^5}\frac{1}{k^+}\log\big[\frac{\vee}{\Lambda}\big]\int_{x,x',z,w,w'}e^{-ik(w'-w)}\,\mathcal J$ (\texttt{Evolution\_terms\_corrected.tex}, Sec.~3). With $\log\frac{\vee}{\Lambda}=\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}$:}
\begin{itemize}
\item \new{$\log\frac{k^+}{\Lambda}\,\mathcal J$ stays in the evolution list: JIMWLK, the evolution of the target from $\Lambda$ up to $k^+$;}
\item \new{$\log\frac{\vee}{k^+}\,\mathcal J$ belongs to this list. In colour-singlet form (Eq.~(3.4) of \texttt{Evolution\_terms\_corrected.tex}, with $\log\frac{\vee}{\Lambda}\to\log\frac{\vee}{k^+}$; $K(a,b,z)=\mathcal K(a-z)\cdot\mathcal K(b-z)$):}
\end{itemize}
{\color{blue}
@@JIMWLK@@}
"""

TOTAL = r"""
\section{\new{Total of this list}}
\new{Groups I--III give $\log\frac{\vee}{k^+}(\mathcal B-\mathcal J)$, Group IV gives $\log\frac{\vee}{k^+}\,\mathcal J$. Their sum is $\log\frac{\vee}{k^+}\,\mathcal B$. In colour-singlet form this is your Eq.~(LOBFKL) with $\delta Y=\log\frac{\vee}{k^+}$ (same normalization as your App.~E):}
@@LOBFKL@@
\new{and the evolution list keeps}
\begin{equation}
\new{\frac{d^3N}{d^3k}\bigg|_{\rm JIMWLK}=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\log\Big[\frac{k^{+}}{\Lambda}\Big]\int_{x,x',z,w,w'}e^{-ik(w'-w)}\,\mathcal J .}
\end{equation}
\new{The physical picture: the soft gluon $p^+$ has $\Lambda<p^+<\vee$. For $p^+>k^+$ it is part of the projectile and can emit the measured gluon (BFKL, interval $\log\frac{\vee}{k^+}$); for $p^+<k^+$ it is emitted next to the target (JIMWLK, interval $\log\frac{k^+}{\Lambda}$). Each soft gluon is counted once.}
"""


def assemble():
    s = build()
    chk = json.load(open(os.path.join(HERE, 'corrected_check.json')))
    worst = max(v['vs_dlog_file'] for v in chk.values() if v['vs_dlog_file'] is not None)
    m, e = ('%.0e' % worst).split('e')
    worst_tex = r'$%s\times10^{%d}$' % (m, int(e))

    # preamble, title, result
    s = s.replace(r"\usepackage{comment}", r"\usepackage{comment}" + "\n"
                  + r"\newcommand{\old}[1]{{\color{red}#1}}" + "\n" + r"\newcommand{\new}[1]{{\color{blue}#1}}" + "\n"
                  + r"\newcommand{\nc}{\color{blue}}", 1)
    s = s.replace(r"\title{BFKL Terms}", r"\title{BFKL Terms \new{(corrected)}}", 1)
    s = s.replace(r"\maketitle", r"\maketitle" + "\n" + r"\sloppy" + "\n" + RESULT.replace('@@WORST@@', worst_tex), 1)

    def after(anchor, text, start_anchor=None, nth=1):
        nonlocal s
        pos = 0
        if start_anchor:
            pos = s.index(start_anchor)
        for _ in range(nth):
            pos = s.index(anchor, pos) + len(anchor)
        s = s[:pos] + text + s[pos:]

    def before(anchor, text, start_anchor=None):
        nonlocal s
        pos = s.index(start_anchor) if start_anchor else 0
        i = s.index(anchor, pos)
        s = s[:i] + text + s[i:]

    # Group I notes (inserted before the next subsection heading)
    g1 = s.index(r"\section{Group I (with c.c.)}")
    before(r"\subsection{Row III}", NOTE_II, r"\section{Group I (with c.c.)}")
    before(r"\subsection{Row IV}", NOTE_III, r"\section{Group I (with c.c.)}")
    before(r"\subsection{Row V}", NOTE_IV, r"\section{Group I (with c.c.)}")
    before(r"\section{Group II}", NOTE_V)
    # Group II
    before(r"\subsubsection{Second part}", NOTE_G2I1, r"\section{Group II}")
    after(r"\subsubsection{Second part}", G2I2_INTRO, r"\section{Group II}")
    # blue 2b, 2c right after the intro, then the red 2d-2f follow
    after(G2I2_INTRO, G2I2_BLUE + G2I2_NOTE + "\n" + r"\new{Remove:}" + "\n", r"\section{Group II}")
    # Group II Row IV
    g2iv = s.index(r"\subsection{Row IV}", s.index(r"\section{Group II}"))
    na = s.index(r"\textsc{NA}", g2iv)
    s = s[:na] + G2IV + s[na + len(r"\textsc{NA}"):]
    # Group III NA rows: Row III part 1, part 2, Row V, Row VI (Row IV stays NA)
    g3 = s.index(r"\section{Group III\, (with c.c.)}")
    reps = [(r"\subsubsection{First part}", G3_III1), (r"\subsubsection{Second part}", G3_III2),
            (r"\subsection{Row V}", G3_V), (r"\subsection{Row VI}", G3_VI)]
    for head, blk in reps:
        h = s.index(head, g3)
        na = s.index(r"\textsc{NA}", h)
        s = s[:na] + blk + s[na + len(r"\textsc{NA}"):]
    before(r"\end{document}", NOTE_G3C + "\n\\allowdisplaybreaks\n" + GROUP4.replace('@@JIMWLK@@', jimwlk_singlet())
           + TOTAL.replace('@@LOBFKL@@', LOBFKL))
    # Group III note on Rows I, II
    h = s.index(r" \subsection{Row III}", g3)
    s = s[:h] + "\n" + r"\new{\textbf{Rows I, II.} The blue terms are the $\log\frac{\vee}{k^+}$ inside $\log\Big[\frac{(\vee-k^{+})(x-z)^{2}+k^{+}(x-w)^{2}}{k^{+}(x-w)^{2}+\Lambda(x-z)^{2}}\Big]$, as in Group I Row IV.}" + "\n" + s[h:]
    open(OUT, 'w').write(s)
    print('written', OUT)


if __name__ == '__main__':
    assemble()
