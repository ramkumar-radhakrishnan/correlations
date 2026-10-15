"""Writes Calculation_main_corrected.tex: Calculation_main.tex with the corrections of Main_corrections.pdf marked in place
(red = as written / to be removed, blue = corrected or added).  Every edit must match exactly once.
Also writes Calculation_main_corrected_build.tex, a copy that compiles here: article class instead of revtex4-2
(not installed), and the tikzit diagrams replaced by a placeholder box."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, '..', 'input', 'Calculation_main.tex')).read()
end = src.index(r'\end{document}')
s, tail = src[:end], src[end:]
R = r'\textcolor{red}'
B = r'\textcolor{blue}'


def _pat(old):
    """regex for old in which every run of whitespace matches any (possibly empty) run of whitespace"""
    chunks = re.split(r'\s+', old.strip())
    lead = r'\s*' if old[:1].isspace() else ''
    return re.compile(lead + r'\s*'.join(re.escape(c) for c in chunks))


def rep(old, new, count=1, region=None):
    """replace old by new (whitespace-tolerant); region = (start_anchor, end_anchor) restricts the search"""
    global s
    pat = _pat(old)
    a, b = 0, len(s)
    if region is not None:
        a = s.index(region[0])
        b = s.index(region[1], a)
    part = s[a:b]
    n = len(pat.findall(part))
    if count is None:
        assert n >= 1, (n, old[:90])
        print('  %d x  %s' % (n, old[:70]))
    else:
        assert n == count, (n, count, old[:90])
    s = s[:a] + pat.sub(lambda m: new, part) + s[b:]


def after(anchor, text):
    """insert text right after the (unique, whitespace-tolerant) anchor"""
    global s
    ms = list(_pat(anchor).finditer(s))
    assert len(ms) == 1, (len(ms), anchor[:90])
    i = ms[0].end()
    s = s[:i] + text + s[i:]


def bnote(text):
    return '\n{\\color{blue}\\noindent\\textbf{Correction.} ' + text + '}\n'


P = r'\left[-\frac{2}{\epsilon}+\log{\left(\frac{k^{2}}{\mu^{2}_{\overline{MS}}}\right)}\right]'

# ------------------------------------------------------------------ header
after(r'\tableofcontents', '\n' + r'''{\color{blue}\noindent\textbf{Corrected copy.} Red: as written (to be removed); blue: corrected or added.
The corrections are explained in \texttt{Main\_corrections.pdf} (labels E1--E11, X1--X7, G1--G4, M1--M8 refer to it).
No limit $\Lambda\to0$ or $\vee\to\infty$ is taken; the only conditions are $\vee>p^+,k^+>\Lambda$. The two red terms in the first equations of Appendix F (the $\delta^2(u-v)$ term and the term it leads to) are red in the original and are not corrections.}
''')

# ------------------------------------------------------------------ X1: gtwoloop1
rep(r'Integrating over $\xi$ and setting $\Lambda$ to be $0$ wherever possible we get\\',
    r'Integrating over $\xi$ ' + R + r'{and setting $\Lambda$ to be $0$ wherever possible} we get\\')
after(r'''+2\ln^{2}\,\left(\frac{\Lambda}{k^{+}}\right)-\frac{67}{9}+\frac{2\pi^{2}}{3}\right]\, \ket{g^{a}_{i}(k^{+},\mathbf{k})}.
\end{align}''', bnote(r'''(X1) The exact $\xi$ integral, $\delta=\Lambda/k^+$, replaces the bracket by
\begin{align}
&\Big(\frac{11}{3}-8\delta+2\delta^2-\frac43\delta^3+4\ln\frac{\delta}{1-\delta}\Big)''' + P + r'''-\frac{67}{9}+\frac{2\pi^2}{3}+\frac{46}{3}\delta-\frac43\delta^2+\frac89\delta^3\notag\\
&+\Big(\frac{22}{3}-8\delta+2\delta^2-\frac43\delta^3\Big)\ln(1-\delta)-\Big(8\delta-2\delta^2+\frac43\delta^3\Big)\ln\delta-2\ln^2(1-\delta)+2\ln^2\delta-4\ln\delta\ln(1-\delta)-8\,{\rm Li}_2(\delta).\notag
\end{align}'''))

# ------------------------------------------------------------------ X2: instgg
rep(r'Integrating wrt $p^{+}$ and setting $\Lambda \rightarrow 0$ wherever possible we get',
    r'Integrating wrt $p^{+}$ ' + R + r'{and setting $\Lambda \rightarrow 0$ wherever possible} we get')
after(r'''\left(\frac{k^{+}}{\Lambda}-1-\ln\left(\frac{\Lambda}{k^{+}}\right)\right)\left|g_{i}^{b}(k^{+},\mathbf{k})\right\rangle
\end{align}''', bnote(r'''(X2) Exactly, $\frac{k^{+}}{\Lambda}-1-\ln\frac{\Lambda}{k^{+}}\to E\equiv\frac{k^{+}}{\Lambda}-\frac{k^{+}}{k^{+}-\Lambda}+\ln\frac{k^{+}-\Lambda}{\Lambda}$, here, in Eq.~\eqref{grho3}, in the one-$\rho$ term below, and in Eq.~\eqref{tworho2}.'''))

# ------------------------------------------------------------------ E8: part2
rep(r'+\delta^{ab}\delta^{ni}\delta^{(3)}(q-p)', r'+\delta^{ab}' + R + r'{\delta^{ni}}' + B + r'{\delta^{ji}}\delta^{(3)}(q-p)')

# ------------------------------------------------------------------ X3, X4: onerho2, tworho2
after(r'''+\frac{\pi^{2}}{3}+2\ln^{2}\left(\frac{\Lambda}{k^{+}}\right)\right)\left|g_{i}^{a}(k^{+},\mathbf{k})\right\rangle .
\end{align}''', bnote(r'''(X3) The $\xi$ integral is done with $\Lambda\to0$. Exactly ($\delta=\Lambda/k^+$) the bracket is
\begin{equation*}
3\ln\frac{\delta}{1-\delta}''' + P + r'''+\frac{\pi^{2}}{3}+2\ln^2\delta-2\ln^2(1-\delta)-2\ln\delta\ln(1-\delta)-4\,{\rm Li}_2(\delta).
\end{equation*}'''))
after(r'''+2\left(\frac{k^{+}}{\Lambda}-1-\ln\left(\frac{\Lambda}{k^{+}}\right)\right)(\mathbf{k}-\mathbf{p})^{2}\bigg)\mathbf{p}^{i}\bigg]\left|g_{i}^{a}(k^{+},\mathbf{k})\right\rangle.
\end{align}''', bnote(r'''(X4) Exactly, the anticommutator part of Eq.~\eqref{3grhorho} is
$\int\frac{dk^{+}}{\sqrt{k^{+}}}\int d^{2}\mathbf{k}\,d^2\mathbf{p}\,\frac{ig^{3}f^{abc}\{\rho^{b}(-\mathbf{p}),\rho^{c}(-\mathbf{k}+\mathbf{p})\}}{32\pi^{9/2}(\mathbf{k}-\mathbf{p})^{2}\mathbf{k}^{2}}V^i\ket{g^a_i(k^+,\mathbf{k})}$ with
\begin{equation*}
V^i=(\mathbf{k}-\mathbf{p})^i\Big[\frac{2\,\mathbf{k}\cdot(\mathbf{k}-\mathbf{p})}{\mathbf{p}^2}(\ell-{\cal L})-2\Big(\frac{k^+}{\Lambda}-\frac{k^+}{k^+-\Lambda}\Big)\Big]-\frac{2\,\mathbf{p}\cdot(\mathbf{k}-\mathbf{p})}{\Delta}{\cal L}\,\mathbf{k}^i+2(\ell+{\cal L})\mathbf{p}^i,
\end{equation*}
$\ell=\ln\frac{k^+-\Lambda}{\Lambda}$, $\Delta=(\mathbf{k}-\mathbf{p})^2-\mathbf{p}^2$, ${\cal L}=\ln\frac{k^+(\mathbf{k}-\mathbf{p})^2-\Lambda\Delta}{k^+\mathbf{p}^2+\Lambda\Delta}$. In the sum with Eq.~\eqref{grho3} (with $E$ exact), the part odd under $\mathbf p\to\mathbf k-\mathbf p$ (the only one that survives) is exactly the first term of Eq.~\eqref{combinedgrhorho} with $\ln\frac{\Lambda}{k^+}\to\ln\frac{\Lambda}{k^+-\Lambda}$.'''))

# ------------------------------------------------------------------ E9: stseve
rep(r'\frac{1}{\left(\frac{\mathbf{q}^{2}}{2q^{+}}\right)\left(\frac{\mathbf{r}^{2}}{2r^{+}}+\frac{\mathbf{p}^{2}}{2p^{+}}\right)\left(\frac{\mathbf{k}^{2}}{2k^{+}}\right)}',
    r'\frac{1}{' + R + r'{\left(\frac{\mathbf{q}^{2}}{2q^{+}}\right)\left(\frac{\mathbf{r}^{2}}{2r^{+}}+\frac{\mathbf{p}^{2}}{2p^{+}}\right)}' + B
    + r'{\left(\frac{\mathbf{p}^{2}}{2p^{+}}\right)\left(\frac{\mathbf{r}^{2}}{2r^{+}}+\frac{\mathbf{q}^{2}}{2q^{+}}\right)}\left(\frac{\mathbf{k}^{2}}{2k^{+}}\right)}')

# ------------------------------------------------------------------ X5: three rho
rep(r'''Setting $\Lambda = 0$, wherever possible ''',
    R + r'''{Setting $\Lambda = 0$, wherever possible} ''' + B + r'''{(nothing is set to zero: the next line is the previous one multiplied out)}''')
rep(r'''Using the argument when $\vee$ is very large or $\Lambda$ is very small, we get
\begin{align}\label{threerhoafter}''', R + r'''{Using the argument when $\vee$ is very large or $\Lambda$ is very small, we get} ''' + B
    + r'''{(X5) Not needed: keep $\frac{1}{\vee}-\frac{1}{\Lambda}$. It cancels exactly against the term $\frac{2k^{+}}{k^{2}}\big(\frac1\Lambda-\frac1\vee\big)$ of Eq.~\eqref{zetainteg}.}
\begin{align}\label{threerhoafter}''')
rep(r'''Using the argument when $\vee$ is very large or $\Lambda$ is very small, we get
\begin{align}
     \left|\Psi_{g\rho\rho\rho}^{7}\right\rangle & \approx''', R + r'''{Using the argument when $\vee$ is very large or $\Lambda$ is very small, we get} ''' + B
    + r'''{(X5) Not needed; Eq.~\eqref{3rhosmom} follows exactly from Eqs.~\eqref{zetainteg} and the exact line before Eq.~\eqref{threerhoafter}.}
\begin{align}
     \left|\Psi_{g\rho\rho\rho}^{7}\right\rangle & \approx''')

# ------------------------------------------------------------------ M6: normalisation
rep(r'\Bigg)\,-\,\frac{1}{8}\left\langle', R + r'{\Bigg)}' + B + r'{\Bigg]}\,' + R + r'{-\,\frac{1}{8}}' + B + r'{+\,\frac{3}{8}}\left\langle')

# ------------------------------------------------------------------ X6: grhoposition
after(r'''\log{\left[1-\frac{k^{+}}{\vee}\right]}\,\log{\left[\frac{\vee}{k^{+}}\right]}\bigg]\, \left|g_{i}^{a}(k^{+},\, \mathbf{z})\right\rangle
\end{align}''', bnote(r'''(X6) Written in the momentum $k^+$ of the outgoing gluon, the $\xi$ integral of the $\ket{\Psi^B_{g\rho}}$ part runs over $\frac{k^+}{\vee}<\xi<\frac{k^+}{k^++\Lambda}$ (not up to $1$ or $1-\frac{\Lambda}{k^+}$), and this part exists only for $k^+<\vee-\Lambda$. Exactly, for $\Lambda<k^+<\vee-\Lambda$ the two brackets become ($\delta=\Lambda/k^+$, $P\equiv-\frac2\epsilon-2\gamma-\ln\frac{(x-z)^2\bar\mu^2}{4}$)
\begin{align*}
&P\Big[\frac{11}{3}+\ln\frac{\Lambda^2}{k^{+2}-\Lambda^2}+\ln\Big(1-\frac{k^+}{\vee}\Big)-8\delta+2\delta^2-\frac43\delta^3\Big]+\frac{\pi^2}{3}-\frac{67}{9}+\ln^2\Big(1-\frac{k^+}{\vee}\Big)+\ln\Big(1-\frac{k^+}{\vee}\Big)\ln\frac{\vee}{k^+}+R_\delta,\\
&R_\delta=\frac{46}{3}\delta-\frac43\delta^2+\frac89\delta^3+\Big(\frac{22}{3}-8\delta+2\delta^2-\frac43\delta^3\Big)\ln(1-\delta)-\Big(8\delta-2\delta^2+\frac43\delta^3\Big)\ln\delta-2\ln\delta\ln(1-\delta)-4\,{\rm Li}_2(\delta)-\ln^2(1+\delta)+\ln(1+\delta)\ln\delta,
\end{align*}
(note $\frac{11}{3}+\ln\frac{\Lambda}{k^+}+\ln(\frac{\vee}{k^+}-1)-\ln\frac{\vee}{\Lambda}=\frac{11}{3}+2\ln\frac{\Lambda}{k^+}+\ln(1-\frac{k^+}{\vee})$). For $\vee-\Lambda<k^+<\vee$ only $\ket{\Psi^A_{g\rho}}$ contributes:
$P\big[\frac{11}{3}-8\delta+2\delta^2-\frac43\delta^3+\ln\frac{\delta}{1-\delta}-2\ln\frac{\vee}{\Lambda}\big]+\frac{\pi^2}{3}-\frac{67}{9}+\ln^2\delta-\ln^2\frac{\vee}{k^+}+R_\delta+\ln^2(1+\delta)-\ln(1+\delta)\ln\delta$.'''))

# ------------------------------------------------------------------ X7: combinedgrhorho
after(r'''Integrating wrt $\xi$ and applying the limits, we get the following''', ' ' + B + r'''{(X7: the exact limits are $\frac{k^+}{\vee}<\xi<\frac{k^+}{k^++\Lambda}$ for the last line, which exists only for $k^+<\vee-\Lambda$. In the result below $\ln\frac{\Lambda}{k^+}\to\ln\frac{\Lambda}{k^+-\Lambda}$ in the first term (X4), $-k^{i}\ln\frac{\vee}{k^+}\to-k^{i}\ln\frac{\vee}{k^++\Lambda}$, and add $-\frac{(2p^{2}-2\,k\cdot p+k^{2})k^{i}+k^{2}p^{i}}{(k-p)^{2}}\ln\big[1+\frac{\Lambda}{k^+}\frac{k^{2}}{(k-p)^{2}}\big]$.)}''')
rep(r'''Setting $\Lambda $ to be $0$ wherever possible, we get''', R + r'''{Setting $\Lambda $ to be $0$ wherever possible, we get} ''' + B
    + r'''{(X7) Not needed: keep $\log\frac{\Lambda k^{2}+k^{+}p^{2}}{\vee k^{2}+k^{+}p^{2}}$ as in the line above.}''')

# ------------------------------------------------------------------ E1, E10: onegrho12, onegrho22
rep(r'& = -\frac{ig^{2}f^{abc}}{8\pi^{3}}\,\int_{k^{+}+\Lambda}^{\vee} dp^{+}',
    r'& = ' + R + r'{-}' + B + r'{+}\frac{ig^{2}f^{abc}}{8\pi^{3}}\,\int_{k^{+}+\Lambda}^{\vee} dp^{+}')
rep(r'\bigg(\frac{k^{+}+p^{+}}{(k^{+}-p^{+})}+(k-p)^{l}-p^{l}-k^{l}\bigg)',
    r'\bigg(\frac{k^{+}+p^{+}}{(k^{+}-p^{+})}' + R + r'{+}(k-p)^{l}-p^{l}-k^{l}\bigg)')

# ------------------------------------------------------------------ E2: onegrhoX1
a = s.index(r'\label{Eq:onegrhoX1}')
b = s.index(r'\frac{ig^2 f^{abc}} {4 \pi^3}', a)
s = s[:b] + R + r'{+}' + B + r'{-}' + s[b:]

# ------------------------------------------------------------------ E5: unitarity
rep(r'''The fact that $A^{(1)}$ can at most be linear in $\rho$ at $g$ order we find ''',
    R + r'''{The fact that $A^{(1)}$ can at most be linear in $\rho$ at $g$ order we find} ''' + B
    + r'''{(E5) This is not correct: $\mf N^{(2)}\propto\int K(x,y)\rho^a(x)\rho^a(y)$ does not commute with $\rho^b(u)$, $[\int_{x,y}K\rho^a\rho^a,\rho^b(u)]=if^{abc}\int_xK(x,u)\{\rho^a(x),\rho^c(u)\}$. Keep $-\mf N^{\dagger(2)}\mf A^{(1)}_{2}+\mf A^{(1)}_{2}\mf N^{(2)}=-[\mf N^{(2)},\mf A^{(1)}_{2}]$ in the relation above; the two equations below and the simplified relation are to be removed.}''')
rep(r'''    \mf {N}^{(2)}\mf{A}^{(1)c}_{1k}(q) &= \mf{A}^{(1) c}_{1k}(q) \mf {N}^{(2)} \\
    \mf {N}^{(2)}\mf{A}^{(1)c}_{2k}(q) &= \mf{A}^{(1) c}_{2k}(q) \mf {N}^{(2)}''',
    r'''    \color{red}\mf {N}^{(2)}\mf{A}^{(1)c}_{1k}(q) &\color{red}= \mf{A}^{(1) c}_{1k}(q) \mf {N}^{(2)} \\
    \color{red}\mf {N}^{(2)}\mf{A}^{(1)c}_{2k}(q) &\color{red}= \mf{A}^{(1) c}_{2k}(q) \mf {N}^{(2)}''')
rep(r'''and thus
the relation for $\mf{A}^{(3)\, \dagger\, b}_{1j}(p)$ simplifies to:
\begin{align}
 \mf{A}^{(3)\, \dagger\, b}_{1j}(p) & = -\mf{A}^{(3)\, b}_{2j}(p) +\int''', R + r'''{and thus
the relation for $\mf{A}^{(3)\, \dagger\, b}_{1j}(p)$ simplifies to:}
\begin{align}
 \color{red}\mf{A}^{(3)\, \dagger\, b}_{1j}(p) & \color{red}= -\mf{A}^{(3)\, b}_{2j}(p) +\int''')

# ------------------------------------------------------------------ E11: A1coeff
rep(r'-\, \mf{A}^{a}_{2i}(k^{+},\mathbf{z}) =   \mf{A}^{\dagger\,a}_{1 i}(k^{+},\mathbf{z})',
    r'-\, \mf{A}^{a}_{2i}(k^{+},' + B + r'{-}\mathbf{z}) =   \mf{A}^{\dagger\,a}_{1 i}(k^{+},\mathbf{z})')

# ------------------------------------------------------------------ E4: A^(3) momentum space
old1 = (r'''\mf{A}^{(3)\, a}_{2i}(k^{+},k) & = \frac{g^{3}N_{c}\rho^{a}(-k)k^{i}}{8\sqrt{2k^{+}}\pi^{2}\,k^{2}}\,\Bigg[\left[-\frac{2}{\epsilon}+\log{\left(\frac{k^{2}}{\mu^{2}_{\overline{MS}}}\right)}\right]\times\, \left[\frac{11}{3}+\log{\left(\frac{\Lambda}{k^{+}}\right)}-2\,\log{\left(\frac{\vee}{\Lambda}\right)}\right] + \frac{\pi^{2}}{3} - \frac{67}{9} +\log^{2}\left(\frac{\Lambda}{k^{+}}\right)-\log^{2}\left(\frac{\vee}{k^{+}}\right)\Bigg]\notag \\ & +\frac{g^{3}N_{c}\rho^{a}(-k)k^{i}}{8\sqrt{2k^{+}}\pi^{2}\,k^{2}}\, \Bigg[\frac{\pi^{2}}{3} - \frac{67}{9} + \log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}+\log{\left[1-\frac{k^{+}}{\vee}\right]}\,\log{\left[\frac{\vee}{k^{+}}\right]}\Bigg]\notag \\ &''')
new1 = (r'''\mf{A}^{(3)\, a}_{2i}(k^{+},k) & = \text{\color{red}(remove)}\color{red}\frac{g^{3}N_{c}\rho^{a}(-k)k^{i}}{8\sqrt{2k^{+}}\pi^{2}\,k^{2}}\,\Bigg[\left[-\frac{2}{\epsilon}+\log{\left(\frac{k^{2}}{\mu^{2}_{\overline{MS}}}\right)}\right]\times\, \left[\frac{11}{3}+\log{\left(\frac{\Lambda}{k^{+}}\right)}-2\,\log{\left(\frac{\vee}{\Lambda}\right)}\right] + \frac{\pi^{2}}{3} - \frac{67}{9} +\log^{2}\left(\frac{\Lambda}{k^{+}}\right)-\log^{2}\left(\frac{\vee}{k^{+}}\right)\Bigg]\notag \\ & \text{\color{red}(remove)}\color{red}+\frac{g^{3}N_{c}\rho^{a}(-k)k^{i}}{8\sqrt{2k^{+}}\pi^{2}\,k^{2}}\, \Bigg[\frac{\pi^{2}}{3} - \frac{67}{9} + \log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}+\log{\left[1-\frac{k^{+}}{\vee}\right]}\,\log{\left[\frac{\vee}{k^{+}}\right]}\Bigg]\notag \\ & \text{\color{blue}(E4, new)}\color{blue}+\frac{g^{3}N_{c}\rho^{a}(-k)k^{i}}{8\sqrt{2k^{+}}\pi^{2}\,k^{2}}\,\Bigg['''
        + P + r'''\left[\frac{11}{3}+\log\frac{\Lambda}{k^{+}}+\log\Big(\frac{\vee}{k^{+}}-1\Big)-\log\frac{\vee}{\Lambda}\right]\notag \\ & \color{blue}\qquad +\frac{\pi^{2}}{3} - \frac{67}{9} + \log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}+\log{\left[1-\frac{k^{+}}{\vee}\right]}\,\log{\left[\frac{\vee}{k^{+}}\right]}\Bigg]\quad\text{\color{blue}(exact form: X6)}\notag \\ &''')
rep(old1, new1)
# coordinate A^(3): V -> vee (M3), first occurrence only (the rest of that line is commented out)
a = s.index(r'\mf{A}^{(3)\,a}_{2i}(k^{+},-z)  & =')
b = s.index(r'\log{\left(\frac{V}{\Lambda}\right)}', a)
s = s[:b] + r'\log{\left(\frac{' + R + '{V}' + B + r'{\vee}}{\Lambda}\right)}' + s[b + len(r'\log{\left(\frac{V}{\Lambda}\right)}'):]

# ------------------------------------------------------------------ E6: B2 ranges, (k+,p+) forms only
TH_OLD = r'\Theta\left(k^{+}-\Lambda-p^{+}\right)'
TH_NEW = R + r'{\Theta\left(k^{+}-\Lambda-p^{+}\right)}' + B + r'{\Theta(p^{+}-\Lambda)\Theta(k^{+}-\Lambda)\Theta(\vee-k^{+}-p^{+})}'
rep(TH_OLD, TH_NEW, count=4, region=(r'\mf{B}^{ba}_{2ji}(k^{+},k;p^{+},p) & =', r'\end{align}'))
rep(TH_OLD, TH_NEW, count=2, region=(r'\mf{B}^{bc}_{2ji}(k^{+}, -\v{w}; p^{+}, -\v{z}) & =', r'\end{align}'))

# ------------------------------------------------------------------ E3: B3
rep(r'\left[\frac{(k^{+}-p^{+})}{(k-p)^{2}}\,\Theta(p^{+}-k^{+}-\Lambda)+',
    r'\left[' + R + '{+}' + B + r'{-}\frac{(k^{+}-p^{+})}{(k-p)^{2}}\,\Theta(p^{+}-k^{+}-\Lambda)+')
rep(r'\left[\frac{1}{((x-z) -  (x-u))^2}\Theta\left(p^{+} - k^{+} - \Lambda \right)',
    r'\left[' + R + '{+}' + B + r'{-}\frac{1}{((x-z) -  (x-u))^2}\Theta\left(p^{+} - k^{+} - \Lambda \right)')
rep(r'\left[\frac{1}{((x-z) -  (x-u))^2}\,\Theta\left(\xi - 1 - \frac{\Lambda}{k^{+}}\right)',
    r'\left[' + R + '{+}' + B + r'{-}\frac{1}{((x-z) -  (x-u))^2}\,\Theta\left(\xi - 1 - \frac{\Lambda}{k^{+}}\right)')

# ------------------------------------------------------------------ G1-G3: amplitude placements (general U)
rep(r'U^{\dagger db}(v) \mf{A}^{\dagger\,(1)\, c}_{1i}', R + r'{U^{\dagger db}(v)}' + B + r'{U^{db}(v)} \mf{A}^{\dagger\,(1)\, c}_{1i}')
rep(r'U^{\dagger db}(v) \mf{A}^{(1)\, c}_{2i}', R + r'{U^{\dagger db}(v)}' + B + r'{U^{db}(v)} \mf{A}^{(1)\, c}_{2i}')
rep(r'U^{ bc}(v)\Bar{\mf{B}}^{\dagger ca}_{1ji}', R + r'{U^{bc}(v)}' + B + r'{U^{\dagger bc}(v)}\Bar{\mf{B}}^{\dagger ca}_{1ji}', count=2)
rep(r'-U^{bc}(v)\Bar{\mf{B}}^{ca}_{2ji}', '-' + R + r'{U^{bc}(v)}' + B + r'{U^{\dagger bc}(v)}\Bar{\mf{B}}^{ca}_{2ji}')
rep(r'w^{\prime}-v)U^{bd}(v)U^{ce}(z) \mf{A}^{\dagger\,(1)\,e}_{1k}', r'w^{\prime}-v)' + R + r'{U^{bd}(v)}' + B + r'{U^{\dagger bd}(v)}U^{ce}(z) \mf{A}^{\dagger\,(1)\,e}_{1k}')

# ------------------------------------------------------------------ G4: CS2 (general U)
CS = (r'\label{eq:CS2}', r'\end{align}')
rep(r'U^{ c^{\prime}a}(w^{\prime})', R + r'{U^{c^{\prime}a}(w^{\prime})}' + B + r'{U^{ac^{\prime}}(w^{\prime})}', count=3, region=CS)   # with and without the space
rep(r'U^{bc^{\prime}}(z)', R + r'{U^{bc^{\prime}}(z)}' + B + r'{U^{c^{\prime}b}(z)}', count=4, region=CS)
rep(r'U^{bd^{\prime}}(z)', R + r'{U^{bd^{\prime}}(z)}' + B + r'{U^{d^{\prime}b}(z)}', count=2, region=CS)
rep(r'U^{bc}(z)\Bar', R + r'{U^{bc}(z)}' + B + r'{U^{cb}(z)}\Bar', count=5, region=CS)
rep(r'U^{bd}(z)', R + r'{U^{bd}(z)}' + B + r'{U^{db}(z)}', count=2, region=CS)
after(r'''\begin{align}\label{eq:CS2}''', '')
s = s.replace(r'''\begin{align}\label{eq:CS2}''', bnote(r'''(G4) The red/blue Wilson-line placements below matter only for a general adjoint $U$ ($U^\dagger=U^T\neq U$); with the convention $U^{bc}=U^{cb}$ both forms are equal.''') + r'''\begin{align}\label{eq:CS2}''', 1)

# ------------------------------------------------------------------ E7: Appendix F
def fix_line(kernel, block_anchor, old_mat, new_mat):
    global s
    a = s.index(block_anchor)
    k = s.index(kernel, a)
    e = s.index(r'\right\}', k)
    line = s[k:e]
    assert line.count(old_mat) == 1, (kernel, line)
    s = s[:k] + line.replace(old_mat, R + '{' + old_mat + '}' + B + '{' + new_mat + '}') + s[e:]


M1 = r'( U^\dagger(x) - U^\dagger(t) )  ( U(y) - U(t) )'
M2 = r'( U^\dagger(y) - U^\dagger(t) )  ( U(x) - U(t) )'
mat = lambda u, v: r'( U^\dagger(%s) - U^\dagger(t) )  ( U(%s) - U(t) )' % (u, v)
A1 = r"&\frac{dN^{(1)}}{d^2 w d^2 w' d y}"
A2 = r"&\frac{dN^{(2)}}{d^2 w d^2 w' d y}"
A3 = r"&\frac{dN^{(3)}}{d^2 w d^2 w' d y}"
A4 = r"&\frac{dN^{(4)}}{d^2 w d^2 w' d y}"
fix_line(r"{\cal K}(y,x,t)", A1, M1, mat('y', 'x'))
fix_line(r"{\cal K}(w',x,t)", A1, M1, mat("w'", 'x'))
fix_line(r"{\cal K}(y,w,t)", A1, M1, mat('y', 'w'))
fix_line(r"{\cal K}(w',w,t)", A1, M1, mat("w'", 'w'))
fix_line(r"{\cal K}(x, w',t)", A2, M2, mat("w'", 'x'))
fix_line(r"{\cal K}(w,y,t)", A2, M2, mat('y', 'w'))
fix_line(r"{\cal K}(w, w',t)", A2, M2, mat("w'", 'w'))
fix_line(r"{\cal K}(y,y,t)", A3, M1, mat('y', 'y'))
fix_line(r"{\cal K}(w', w',t)", A3, M1, mat("w'", "w'"))
fix_line(r"{\cal K}(x,x,t)", A4, M2, mat('x', 'x'))
fix_line(r"{\cal K}(w, w,t)", A4, M2, mat('w', 'w'))
rep(r'''&\frac{dN^{(5)}}{d^2 w d^2 w' d y} =  -\frac{\alpha_s^2 \delta Y}''', r'''&\frac{dN^{(5)}}{d^2 w d^2 w' d y} =  ''' + R + '{-}' + B + r'''{+}\frac{\alpha_s^2 \delta Y}''')
a = s.index(r"&\frac{dN^{(5)}}")
part = s[a:]
n0 = part.count(r'\right\} \left\{')
assert n0 == 2, n0
part = part.replace(r'\right\} \left\{', r'\right\} ' + B + r'{{\rm Tr}}\left\{')
n1 = part.count('\\right\\} \n \\left\\{')
assert n1 == 2, n1
part = part.replace('\\right\\} \n \\left\\{', '\\right\\} \n ' + B + r'{{\rm Tr}}\left\{')
s = s[:a] + part
after(r'''The sum of the above gives the needed target evolution of the leading order single inclusive production. ''', '')
s = s.replace(r'''Note that ${\cal K}(x,y,z) = {\cal K}(y,x,z)$. ''', r'''Note that ${\cal K}(x,y,z) = {\cal K}(y,x,z)$. ''' + B
              + r'''{(E7) In $dN^{(1)}$--$dN^{(4)}$ the matrix must sit at the points $u,v$ of ${\cal K}(u,v,t)$ of each line: $(U^\dagger(u)-U^\dagger(t))(U(v)-U(t))$ in $dN^{(1)}$, its transpose $(U^\dagger(v)-U^\dagger(t))(U(u)-U(t))$ in $dN^{(2)}$, and $(U^\dagger(v)-U^\dagger(t))(U(v)-U(t))$ for ${\cal K}(v,v,t)$ in $dN^{(3)}$, $dN^{(4)}$.} ''', 1)

out = s + tail
open(os.path.join(HERE, 'Calculation_main_corrected.tex'), 'w').write(out)

# ------------------------------------------------------------------ build copy for this machine
b = out
b = b.replace(r'\documentclass[aps,prd,amsmath,amssymb,nofootinbib]{revtex4-2}',
              '\\documentclass[10pt]{article}\n\\usepackage[a3paper,margin=1.5cm]{geometry}\n\\usepackage{amsmath,amssymb}')
b = b.replace(r'\usepackage{tikzit}', '\\usepackage[T1]{fontenc}\n' + r'\newcommand{\tikzfig}[1]{\fbox{\small diagram \texttt{\detokenize{#1}} (not included)}}')
b = b.replace(r'\input{stiles.tikzstyles}', '')
b = b.replace(r'\usepackage[inline]{showlabels}', '')
b = b.replace(r'\begin{document}', '\\begin{document}\n\\maketitle\n\\allowdisplaybreaks', 1)
open(os.path.join(HERE, 'Calculation_main_corrected_build.tex'), 'w').write(b)
print('wrote Calculation_main_corrected.tex (%d chars) and the build copy' % len(out))
