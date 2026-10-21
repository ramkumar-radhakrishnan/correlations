"""Writes Group_I_v2_corrected.tex: input/Group_I_v2.tex (the revised Group I file) with every correction of
Group_I_v2_check.pdf (R1-R5) and Group_I_check.pdf (G3, G4) marked in place: new text in blue, removed text in red.
The colour marks of the input are remapped first so that blue marks only the corrections
(blue -> teal, green -> red (removed), [Done] red -> gray, orange unchanged).  Every edit is asserted to match the
expected number of times.  New displayed equations are unnumbered, so the equation numbers of the input are kept."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(HERE, 'input', 'Group_I_v2.tex')).read()


def rep(old, new, count=1):
    global s
    n = s.count(old)
    assert n == count, (n, count, old[:100])
    s = s.replace(old, new)


R = lambda t: r'{\color{red}' + t + '}'
B = lambda t: r'{\color{blue}' + t + '}'


def red_term(green_snip, k=0):
    """Wrap the Row III term that contains the k-th green_snip (its sign, coefficient line and colour line) in red."""
    global s
    g = -1
    for _ in range(k + 1):
        g = s.index(green_snip, g + 1)
    i = s.rindex(r'\frac{ig^{4}', 0, g)
    j = len(s[:i].rstrip())
    if s[j - 1] in '+-':
        i = j - 1
    n1 = s.index(r'\notag', g)
    n2 = s.index(r'\notag', n1 + 1)
    line1, line2 = s[i:n1], s[n1:n2]
    sep = r'\notag \\ & '
    assert line2.startswith(sep), line2[:40]
    line1 = line1.replace(r'{\color{green}', '{')
    s = s[:i] + R(line1) + sep + R(line2[len(sep):]) + s[n2:]


# ---------------------------------------------------------------- colour marks of the input
rep(r'{\color{red} [Done\,$\checkmark$]}', r'{\color{gray} [Done\,$\checkmark$]}', count=5)
rep(r'{\color{blue}', r'{\color{teal}', count=22)

# ---------------------------------------------------------------- header, title
rep(r'\title{Group -I crossection derivation}', r'\title{Group -I ' + R('crossection') + B('cross-section') + ' derivation}')
rep(r'\maketitle', r'''\maketitle
\begingroup\color{blue}\noindent\textbf{Corrected copy of the revised Group I file, with all Group I terms.} New or changed text is in blue, removed text in red. So that blue marks only the corrections, the colour marks of your file were changed: your blue (the $\log\frac{\vee}{\Lambda}$ terms) is now teal, your green $\frac{1}{\Lambda}$ terms are red (they cancel, see the end of Row III), orange is unchanged, and the [Done] marks are gray. The reasons and the numerical checks are in \texttt{Group\_I\_v2\_check.pdf} (R1--R5) and \texttt{Group\_I\_check.pdf} (G1--G6). As in your file, $\Lambda\to0$ wherever possible. The equation numbers are those of your file. A summary of all terms and of the cancellations inside Group I is at the end.\par\endgroup''')

# ---------------------------------------------------------------- Row II, part 1 (R5)
rep(r"\int_{w, w',x,x',}", r"\int_{w, w',x,x'" + R(',') + '}', count=2)
rep(r'\log{\left[1-\frac{k^{+}}{\vee}\right]^{2}}',
    R(r'\log{\left[1-\frac{k^{+}}{\vee}\right]^{2}}') + r'\,' + B(r'\log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}'), count=2)
rep('Rewwriting, we get', R('Rewwriting') + B('Rewriting') + ', we get')

# ---------------------------------------------------------------- Row II, part 2a (R2 = G3)
old_2a = r'''-\frac{1}{(2\pi)^{3}}\,\frac{i\,g^{4}}{32\pi^{5}}\,\frac{1}{k^{+}}
\left[\log\left(\frac{\vee}{\Lambda}\right)-\log\left(\frac{\vee}{k^{+}}\right)\right]
\int_{w,w'} e^{-ik\cdot(w'-w)}\int_{x,y,x'}\frac{(x'-w')^{i}}{(x'-w')^{2}}\;\mathfrak{F}^{i}(w,x,y)\notag\\
&\times\Big[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\,\rho^{b'}(x')\Big]
\Big[U^{ah}(w)\,f^{hbc}\left\{\rho^{b}(x),\rho^{c}(y)\right\}-f^{abc}\,U^{bd}(x)\,U^{ce}(y)\left\{\rho^{d}(x),\rho^{e}(y)\right\}\Big]'''
new_2a = (r'''-\frac{1}{(2\pi)^{3}}\,\frac{i\,g^{4}}{32\pi^{5}}\,\frac{1}{k^{+}}
''' + R(r'\left[\log\left(\frac{\vee}{\Lambda}\right)-\log\left(\frac{\vee}{k^{+}}\right)\right]') + r'''
\int_{w,w'} e^{-ik\cdot(w'-w)}\int_{x,y,x'}\frac{(x'-w')^{i}}{(x'-w')^{2}}\;''' + R(r'\mathfrak{F}^{i}(w,x,y)') + r'''\notag\\
&\times\Big[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\,\rho^{b'}(x')\Big]\notag\\
&\times\Big['''
    + B(r'\Big[\log\Big(\frac{\vee}{\Lambda}\Big)\,\mathfrak{F}^{i}_{\beta}(w,x,y)-\log\Big(\frac{\vee}{k^{+}}\Big)\,\mathfrak{F}^{i}(w,x,y)\Big]')
    + r'''\,U^{ah}(w)\,f^{hbc}\left\{\rho^{b}(x),\rho^{c}(y)\right\}\notag\\
&\qquad-''' + B(r'\Big[\log\Big(\frac{\vee}{\Lambda}\Big)-\log\Big(\frac{\vee}{k^{+}}\Big)\Big]\,\mathfrak{F}^{i}(w,x,y)')
    + r'''\,f^{abc}\,U^{bd}(x)\,U^{ce}(y)\left\{\rho^{d}(x),\rho^{e}(y)\right\}\Big]''')
rep(old_2a, new_2a)
anchor = r'''\log\frac{(w-x)^{2}}{(x-y)^{2}}\right].
\end{align}'''
rep(anchor, anchor + r'''
\begingroup\color{blue}\noindent Here $\mathfrak{F}^{i}_{\beta}$ is the first line of $\mathfrak{F}^{i}$, and $\mathfrak{F}^{i}_{\alpha}=\mathfrak{F}^{i}-\mathfrak{F}^{i}_{\beta}$ its second line:
\begin{align*}
\mathfrak{F}^{i}_{\beta}(w,x,y)&=\frac{(x-w)^{i}}{(x-w)^{2}}\int_{z}\frac{(w-z)\cdot(y-z)}{(w-z)^{2}(y-z)^{2}}-\frac{(y-w)^{i}}{(y-w)^{2}}\int_{z}\frac{(w-z)\cdot(x-z)}{(w-z)^{2}(x-z)^{2}},\\
\mathfrak{F}^{i}_{\alpha}(w,x,y)&=-\left[\frac{(x-w)^{i}}{(x-w)^{2}}-\frac{(y-w)^{i}}{(y-w)^{2}}\right]\int_{z}\frac{(x-z)\cdot(y-z)}{(x-z)^{2}(y-z)^{2}}.
\end{align*}
\textbf{Why (G3, R2).} \texttt{Calculation\_main.tex} drops $-[\mathbb{N}^{(2)},\mathbb{A}^{(1)}_{2}]$ in the unitarity relation for $\mathbb{A}^{(3)\dagger}_{1}$. The two operators do not commute: $\big[\int_{x,y}K(x,y)\rho^{b}(x)\rho^{b}(y),\rho^{a}(u)\big]=if^{bac}\int_{x}K(x,u)\{\rho^{b}(x),\rho^{c}(u)\}$.
\begin{itemize}
\item The $-2k^{i}\log\frac{\vee}{\Lambda}$ term of the bracket in (1.3) is exactly this commutator; in coordinate space it is the $\mathfrak{F}_{\alpha}$ part of $\mathfrak{F}$.
\item With the virtual factor ordered in $p^{+}$ (gluons softer than $k^{+}$ to the left of $\mathbb{A}^{(1)}$, harder ones to the right), the $\log\frac{\vee}{\Lambda}$ part of the commutator drops out of the term with the charges $\rho$ not rotated, $U^{ah}(w)f^{hbc}\{\rho^{b}(x),\rho^{c}(y)\}$.
\item Only this three-$\rho$ term $\propto\log\frac{\vee}{\Lambda}$ changes. The symmetrized three-$\rho$ terms cancel only with this change (fix F4).
\item $\mathfrak{F}_{\alpha}$ and $\mathfrak{F}_{\beta}$ are separately logarithmically divergent at large $z$. The divergent parts are proportional to the total charge.
\end{itemize}\par\endgroup''')

# ---------------------------------------------------------------- Row II, part 2b (R5: int_{x'})
rep(r'\int_{x,y} e^{ik\cdot(w-y)}', r'\int_{x,y' + B(r",x'") + r'} e^{ik\cdot(w-y)}')

# ---------------------------------------------------------------- Row III: green 1/Lambda terms (they cancel), user's red pair
for g in (r'{\color{green}\left[\frac{1}{\Lambda}+\frac{1}{k^{+}-\vee}\right]}',
          r'{\color{green}\left[\frac{1}{(\vee-k^{+})}-\frac{1}{\Lambda}\right]}'):
    red_term(g)
gl = r'{\color{green}\left[\frac{1}{\Lambda}-\frac{1}{k^{+}-\Lambda}\right]}'
assert s.count(gl) == 2
red_term(gl, 1)                                   # N_3a
red_term(gl, 0)                                   # N_1, second green term
assert r'{\color{green}' not in s
rep(r'{\color{red}-\frac{ig^{4}\, f^{abd}}{8\pi^{5}}\, \frac{1}{k^{+}}\, {\color{teal}',
    r'{\color{red}-\frac{ig^{4}\, f^{abd}}{8\pi^{5}}\, \frac{1}{k^{+}}\, {')
rep(r'{\color{red}+ \frac{ig^{4}\, f^{abd}}{8\pi^{5}}\, \frac{1}{k^{+}}\,{\color{teal}',
    r'{\color{red}+ \frac{ig^{4}\, f^{abd}}{8\pi^{5}}\, \frac{1}{k^{+}}\,{')

# ---------------------------------------------------------------- Row III: principal values (R3 = G4)
logs = [(r'\left[\frac{\left[\vee\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}{\left[(\Lambda + k^{+})\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}\right]', 1),
        (r'\left[\frac{\left[(k^{+}-\Lambda)\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}{\left[\Lambda\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}\right]', 1),
        (r'\left[\frac{\left[k^{+}\,(x-w)^{2}-\vee\, (x-z)^{2}\right]}{\left[k^{+}\,(x-w)^{2}-\Lambda\, (x-z)^{2}\right]}\right]', 1),
        (r'\left[\frac{(x-w)^{2}\,k^{+}-(x-z)^{2}\,\vee}{(x-w)^{2}\,k^{+} - (x-z)^{2}\,(k^{+}+\Lambda)}\right]', 6),
        (r'\left[\frac{(x-w)^{2}\,k^{+}-(x-z)^{2}\,(k^{+}-\Lambda)}{(x-w)^{2}\,k^{+}-(x-z)^{2}\,\Lambda}\right]', 6)]
for arg, n in logs:
    rep(r'\log{' + arg + '}', B(r'\log{\left|' + arg[len(r'\left['):-len(r'\right]')] + r'\right|}'), count=n)

rep(r'\subsection{Row IV', r'''\begingroup\color{blue}\noindent\textbf{Row III: corrections.}
\begin{itemize}
\item \textbf{The $\frac{1}{\Lambda}$ terms cancel.} The four terms in red in $d^{3}N_{1}$, $d^{3}N_{2}$, $d^{3}N_{3a}$ (your green terms) multiply the same structure $C\equiv\frac{ig^{4}f^{abd}}{8\pi^{5}}\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}}\times$\,colour. With $D\equiv(x-z)^{2}-(x-w)^{2}$ and the overall signs of $d^{3}N_{2}$ ($-$) and of the denominators,
\begin{equation*}
\frac{C}{D}\Big\{\Big[\frac{1}{\Lambda}-\frac{1}{\vee-k^{+}}\Big]+\Big[\frac{1}{\Lambda}-\frac{1}{k^{+}-\Lambda}\Big]+\Big[\frac{1}{\vee-k^{+}}-\frac{1}{\Lambda}\Big]-\Big[\frac{1}{\Lambda}-\frac{1}{k^{+}-\Lambda}\Big]\Big\}=0 .
\end{equation*}
They are removed. Numerically the sum is zero to $6\times10^{-16}$ of one term.
\item \textbf{The two red terms of $d^{3}N_{2}$} that you had marked are equal and opposite; they cancel.
\item \textbf{Principal values (G4, R3).} The integrand of (1.10) has a pole at $p^{+}_{*}=k^{+}(x-w)^{2}/(x-z)^{2}$. When it lies inside the $p^{+}$ range the integral is a principal value, and the argument of the logarithm changes sign. Hence $\log[\cdots]\to\log|\cdots|$ in the 3 logarithms of $d^{3}N_{1}$, the 6 of $d^{3}N_{2}$ and the 6 of $d^{3}N_{3a}$. $\mathcal{A}$ and the logarithm of $d^{3}N_{3b}$ need no change.
\end{itemize}\par\endgroup
\subsection{Row IV''')

# ---------------------------------------------------------------- Row V, part 1: phases inside the integrals (R5)
ph_z = r'e^{-ik\cdot (w^\prime-z)}\,e^{-ik\cdot (z-x)\frac{p^{+}}{k^{+}}}'
ph_x = r'e^{-ik\cdot (w^\prime-x)}\,e^{-ik\cdot (x-z)\frac{q^{+}}{k^{+}}}'
xs = r"\int_{x,z,y,x',y'" + B(',') + "w'}"
for lim, n in ((r'\int_{\Lambda}^{k^{+}-\Lambda}', 1), (r'\int_{0}^{k^{+}}', 2)):
    rep(r"\int_{x,z,y,x',y'w'}\," + ph_z + r'\,\frac{ig^{4}}{16\pi^{4}}\, \frac{1}{k^{+\,2}}\,' + lim + r'\, \frac{dp^{+}}{2\pi}\,',
        xs + r'\,' + R(ph_z) + r'\,\frac{ig^{4}}{16\pi^{4}}\, \frac{1}{k^{+\,2}}\,' + lim + r'\, \frac{dp^{+}}{2\pi}\,' + B(ph_z) + r'\notag \\ & ', count=n)
rep(r"\int_{x,z,y,x',y'w'}\," + ph_x + r'\,\frac{ig^{4}}{16\pi^{4}}\, \frac{1}{k^{+}}\,\int_{0}^{k^{+}}\, \frac{dq^{+}}{2\pi}\,',
    xs + r'\,' + R(ph_x) + r'\,\frac{ig^{4}}{16\pi^{4}}\, \frac{1}{k^{+}}\,\int_{0}^{k^{+}}\, \frac{dq^{+}}{2\pi}\,' + B(ph_x) + r'\notag \\ & ', count=2)
rep(r"x',y'w'}", r"x',y'" + B(',') + "w'}", count=4)

# ---------------------------------------------------------------- Row V, part 2: phases inside the integrals, braces (R5)
ph2 = r"e^{-ik\cdot(w'-z)}\,e^{-ik\cdot (z-x)\,\frac{p^{+}}{k^{+}}}"
ph2q = r"e^{-ik\cdot(w'-x)}\,e^{-ik\cdot (x-z)\,\frac{q^{+}}{k^{+}}}"
rep(r"\int_{ w^\prime} " + ph2 + r" \Bigg\{\notag \\ & \frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\,\int_{\Lambda}^{k^{+}-\Lambda}\, \frac{dp^{+}}{2\pi} \, \int_{x,y,z,x'}\,",
    r"\int_{w^\prime}\, " + R(ph2) + r" \Bigg\{\notag \\ & \frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\,\int_{\Lambda}^{k^{+}-\Lambda}\, \frac{dp^{+}}{2\pi} \, \int_{x,y,z,x'}\,"
    + B(ph2) + r"\,", count=4)
rep(r"\int_{w^\prime} " + ph2 + r" \Bigg\{\notag \\ & \frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\, \int_{x,y,z,x'}\,",
    r"\int_{w^\prime}\, " + R(ph2) + r" \Bigg\{\notag \\ & \frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\, \int_{x,y,z,x'}\," + B(ph2) + r"\,")
rep(r"\int_{0}^{k^{+}}\, \frac{dp^{+}}{2\pi} \, " + ph2 + r" \Bigg\{\notag \\ & \frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\,\int_{x,y,z,x',w'}\,",
    r"\int_{0}^{k^{+}}\, \frac{dp^{+}}{2\pi} \, " + R(ph2) + r" \Bigg\{\notag \\ & \frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\,\int_{x,y,z,x',w'}\," + B(ph2) + r"\,", count=2)
rep(r"\int^{0}_{k^{+}}\, \frac{dq^{+}}{2\pi} \, " + ph2q + r"\, \frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\, \int_{x,y,z,x',w'}\,",
    r"\int^{0}_{k^{+}}\, \frac{dq^{+}}{2\pi} \, " + R(ph2q) + r"\, \frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\, \int_{x,y,z,x',w'}\," + B(ph2q)
    + r"\notag \\ & ", count=2)
Lor = r'{\color{orange}\left[\log{\left[\frac{\vee}{\Lambda}\right]}+\log{\left[\frac{k^{+}-\Lambda}{\vee}\right]}\right]}'
rep(r"\int_{w^\prime} e^{-ik\cdot(w'-z)}\, \frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\," + Lor + r" \int_{x,y,z,x'}\, ",
    r"\int_{w^\prime} " + R(r"e^{-ik\cdot(w'-z)}") + r"\, \frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\," + Lor + r" \int_{x,y,z,x'}\, "
    + B(r"e^{-ik\cdot(w'-z)}") + r"\notag \\ & ", count=2)
rep(r"\int_{w^\prime} e^{-ik\cdot(w'-x)}\, \frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\," + Lor + r"\, \int_{x,y,z,x'}\,",
    r"\int_{w^\prime} " + R(r"e^{-ik\cdot(w'-x)}") + r"\, \frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\," + Lor + r"\, \int_{x,y,z,x'}\,"
    + B(r"e^{-ik\cdot(w'-x)}") + r"\notag \\ & ", count=2)
rep(r"\Bigg]\,\Bigg\} \notag \\ & \Bigg[U^{ab'}(x')-U^{ab'}(w')\Bigg]",
    r"\Bigg]\," + R(r"\Bigg\}") + r" \notag \\ & \Bigg[U^{ab'}(x')-U^{ab'}(w')\Bigg]", count=7)
rep(r'''We write it in two parts:
\subsubsection{First part:''',
    B(r'''(R5) In Row V the phases $e^{-ik\cdot(w'-z)}e^{-ik\cdot(z-x)p^{+}/k^{+}}$ (and the $q^{+}$ ones) depend on $x$, $z$ and $p^{+}$; they were moved inside $\int dp^{+}$ and $\int_{x,y,z,\dots}$, and the unmatched brace $\Bigg\}$ was removed.''')
    + '\n\n' + r'''We write it in two parts:
\subsubsection{First part:''')

# ---------------------------------------------------------------- Row V, 2b UV: Lambda -> 0 (R4)
anchor = r'''\begin{align}
\frac{d^{3}N_{2b\,\rm finite\, NLO}^{\rm}}{d^{3}k} & ='''
rep(anchor, r'''\begingroup\color{blue}\noindent\textbf{Note (R4).} At the accuracy $\Lambda\to0$ of this equation, $2\log(1-\frac{\Lambda}{k^{+}})$ is $O(\Lambda/k^{+})$, the same order as the dropped $\frac{4\Lambda}{k^{+}}-\frac{\Lambda^{2}}{k^{+2}}+\frac{2\Lambda^{3}}{3k^{+3}}$; likewise $\log\frac{k^{+}-\Lambda}{\vee}=-\log\frac{\vee}{k^{+}}+O(\Lambda/k^{+})$ everywhere in Row V. Keeping these exact endpoint logarithms is harmless; dropping them gives the bracket $2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^{+}}-\frac{11}{6}$.\par\endgroup
''' + anchor)

# ---------------------------------------------------------------- Row V, 2b finite, 'Simplifying' (R1)
rep(r"\,\frac{k^{+}}{q^{+}}\times\notag \\ & \Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[N_{c}\,U^{af}(y)-N_{c}\,U^{af}(z)\Big]",
    r"\,\frac{k^{+}}{q^{+}}\times\notag \\ & \Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[N_{c}\,U^{af}(y)-N_{c}\,U^{af}(" + R('z') + B('x') + r")\Big]")
rep(r"{2(y-z)^{2}}\right]\notag \\ & \Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[N_{c}\,U^{af}(y)-N_{c}\,U^{af}(z)\Big]",
    r"{2(y-z)^{2}}\right]\notag \\ & \Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[N_{c}\,U^{af}(y)-N_{c}\,U^{af}(" + R('z') + B('x') + r")\Big]")
rep(r'''Now, we do plus prescription for''',
    B(r'''(R1) In the two $q^{+}$ blocks above (the $B_{2}$ term, $-\frac{p^{+}}{k^{+}-p^{+}}=-\frac{k^{+}}{q^{+}}+1$) the colour factor is $N_{c}[U(y)-U(x)]$: the $B_{2}$ term belongs to 2b. Only the counterterm has $U(z)$. With this change the equation equals the previous one exactly.''')
    + '\n\n' + r'''Now, we do plus prescription for''')

# ---------------------------------------------------------------- summary of all Group I terms
rep(r'\end{document}', r'''\begingroup\color{blue}
\section*{Summary: all Group I terms}
Equation numbers as above. $E(r)\equiv\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{r^{2}\mu^{2}_{\overline{MS}}}{4}$.
\begin{center}
\begin{tabular}{lll}
\hline
row & term & final form\\
\hline
Row I & four $\rho$'s & not treated here\\
Row II & part 1 (two $\rho$) & (1.2)\\
 & part 2a (three $\rho$) & (1.5), with the G3 change\\
 & part 2b (three $\rho$) & (1.8) with $R^{i}$ of (1.9), $\int_{x'}$ added\\
Row III & $d^{3}N_{1}$, $d^{3}N_{2}$, $d^{3}N_{3a}$, $d^{3}N_{3b}$ & (1.11)--(1.14) with $\mathcal{A}$; $\frac{1}{\Lambda}$ terms removed, $\log|\cdot|$\\
Row IV & & (1.16)\\
Row V & part 1 & (1.19) $=$ (1.32)\\
 & part 2a & (1.26) $=$ (1.33)\\
 & part 2b, UV & (1.28) $=$ (1.34)\\
 & part 2b, finite & (1.31) $=$ (1.35)\\
\hline
\end{tabular}
\end{center}
\paragraph{Cancellations inside Group I.}
\begin{enumerate}
\item \textbf{$\frac{1}{\Lambda}$ terms (Row III):} they cancel exactly (end of Row III).
\item \textbf{Running coupling (Row II part 1 against Row V 2b UV):} the $\beta_{0}$ terms are
\begin{align*}
\frac{d^{3}N_{\rm II}}{d^{3}k}\bigg|_{\beta_{0}}={}&\frac{1}{(2\pi)^{3}}\frac{g^{4}N_{c}}{32\pi^{4}}\frac{1}{k^{+}}\frac{11}{3}\int_{w,w',x,x'}e^{-ik\cdot(w'-w)}\frac{(x'-w')\cdot(x-w)}{(x'-w')^{2}(x-w)^{2}}\big[-E(x-w)\big]\\
&\times\big[U^{ab'}(x')-U^{ab'}(w')\big]\rho^{b'}(x')\big[U^{ab}(w)-U^{ab}(x)\big]\rho^{b}(x),\\
\frac{d^{3}N_{\rm V}}{d^{3}k}\bigg|_{\beta_{0}}={}&\frac{1}{(2\pi)^{3}}\frac{g^{4}N_{c}}{16\pi^{4}}\frac{1}{k^{+}}\Big(-\frac{11}{6}\Big)\int_{y,z,x',w'}e^{-ik\cdot(w'-z)}\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}}\,E(y-z)\\
&\times\big[U^{ab'}(x')-U^{ab'}(w')\big]\big[U^{af}(y)-U^{af}(z)\big]\rho^{b'}(x')\rho^{f}(y).
\end{align*}
Rename $y\to x$, $z\to w$, $f\to b$ in the second. Then
\begin{itemize}
\item the phase, the kernel and $E$ become those of the first;
\item $[U(x)-U(w)]=-[U(w)-U(x)]$;
\item $\frac{g^{4}N_{c}}{16\pi^{4}}\big(-\frac{11}{6}\big)=-\frac{g^{4}N_{c}}{32\pi^{4}}\frac{11}{3}$.
\end{itemize}
Hence
\begin{equation*}
\frac{d^{3}N_{\rm II}}{d^{3}k}\bigg|_{\beta_{0}}+\frac{d^{3}N_{\rm V}}{d^{3}k}\bigg|_{\beta_{0}}=0
\end{equation*}
per ordering and for the c.c. The pole, $\gamma_{E}$, $\mu$ and $\log(x-w)^{2}$ all cancel. With the exact brackets it holds exactly in $\Lambda$.
\item \textbf{What is left:}
\begin{itemize}
\item The $\log\frac{\vee}{\Lambda}$ terms (teal and orange) are the evolution (JIMWLK/BFKL) terms.
\item The remaining poles multiply only rapidity logarithms; they combine with Group II, Row IV.
\item The charge renormalization of LO gives the $\frac{\beta_{0}}{2}\delta(1-\xi)$ of the fragmentation counterterm.
\item Row II also contains the constant $-\frac{67}{9}+\frac{\pi^{2}}{3}$: per ordering this is $\frac{\alpha_{s}}{2\pi}K\times$LO with $K=N_{c}\big(\frac{67}{18}-\frac{\pi^{2}}{6}\big)$, the CMW constant.
\end{itemize}
\end{enumerate}
\endgroup
\end{document}''')

open(os.path.join(HERE, 'Group_I_v2_corrected.tex'), 'w').write(s)
print('wrote Group_I_v2_corrected.tex')
