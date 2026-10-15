"""Write Group_III_corrected.tex: the notes with the corrections in blue (removed text in red)."""
L = open('input/Group_III.tex').read().split('\n')
B = r"{\color{blue}%s}"


def rep(n, old, new, count=1):
    """replace in line n (1-based); assert presence"""
    s = L[n - 1]
    assert s.count(old) >= count, (n, old, s.count(old))
    L[n - 1] = s.replace(old, new, count)


def neg_first(n, old):
    rep(n, old, old.replace('-', B % '+', 1))


# ---- measures (Rows I, II)
rep(44, r"\int_{z,w,w^{\prime}}", r"\int_{%s z,w,w^{\prime}}" % (B % 'x,y,x^{\\prime},'))
rep(55, r"\int_{z,w,w^{\prime}}", r"\int_{%s z,w,w^{\prime}}" % (B % 'x,y,x^{\\prime},'))
# ---- Row III: overall sign of the row
rep(63, r"& =   \frac{4}{(2\pi)^{3}}", r"& = %s\,\frac{4}{(2\pi)^{3}}" % (B % '-'))
rep(68, r"& = -\, \frac{1}{(2\pi)^{3}}", r"& = %s\, \frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(72, r"& = -\,\frac{1}{(2\pi)^{3}}", r"& = %s\,\frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(76, r"& = -\,\frac{1}{(2\pi)^{3}}", r"& = %s\,\frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(76, r"\int_{0}^{\vee-k^{+}}", r"\int_{%s}^{\vee-k^{+}}" % (B % r"\Lambda"))
rep(76, r"& +\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{16\pi^{5}}", r"& %s\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{16\pi^{5}}" % (B % '-'))
rep(80, r"& = -\, \frac{1}{(2\pi)^{3}}", r"& = %s\, \frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(85, r"& = -\, \frac{1}{(2\pi)^{3}}", r"& = %s\, \frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(85, r"\int_{w,w',y',z,y,x}", r"\int_{w,%s,y',z,y,x}" % (B % "x'"))
rep(90, r"& = -\, \frac{1}{(2\pi)^{3}}", r"& = %s\, \frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(90, r"\int_{0}^{\vee-k^{+}}", r"\int_{%s}^{\vee-k^{+}}" % (B % r"\Lambda"))
rep(90, r"\int_{w,w',y',z,y,x}", r"\int_{w,%s,y',z,y,x}" % (B % "x'"))
rep(91, r"&  +\, \frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}", r"&  %s\, \frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}" % (B % '-'))
rep(91, r"\int_{w,y',z,y,x}", r"\int_{w,%s,y',z,y,x}" % (B % "x'"))
# ---- Row V
rep(98, r"& =   \frac{2}{(2\pi)^{3}}", r"& = %s\,\frac{2}{(2\pi)^{3}}" % (B % '-'))
rep(101, r"& = -\, \frac{1}{(2\pi)^{3}}", r"& = %s\, \frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(105, r"& = -\,\frac{1}{(2\pi)^{3}}", r"& = %s\,\frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(109, r"& = -\,\frac{1}{(2\pi)^{3}}", r"& = %s\,\frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(109, r"\int_{0}^{\vee-k^{+}}", r"\int_{%s}^{\vee-k^{+}}" % (B % r"\Lambda"))
rep(109, r"& +\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{8\pi^{5}}",
    r"& %s\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{8\pi^{5}}" % (B % '-'))
rep(113, r"& = -\,\frac{1}{(2\pi)^{3}}", r"& = %s\,\frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(113, r"\int_{0}^{\vee-k^{+}}", r"\int_{%s}^{\vee-k^{+}}" % (B % r"\Lambda"))
rep(113, r"& +\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'ca}}{8\pi^{5}}", r"& %s\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'ca}}{8\pi^{5}}" % (B % '-'))
# ---- Row VI
rep(117, r"& =   \frac{2}{(2\pi)^{3}}", r"& = %s\,\frac{2}{(2\pi)^{3}}" % (B % '-'))
rep(120, r"& = -\, \frac{1}{(2\pi)^{3}}", r"& = %s\, \frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(120, r"\int_{w,w',y',z,y,x}", r"\int_{w,w',y',z,y,x,%s}" % (B % "x'"))
rep(124, r"& = -\,\frac{1}{(2\pi)^{3}}", r"& = %s\,\frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(128, r"& = -\,\frac{1}{(2\pi)^{3}}", r"& = %s\,\frac{1}{(2\pi)^{3}}" % (B % '+'))
rep(128, r"\int_{0}^{\vee-k^{+}}", r"\int_{%s}^{\vee-k^{+}}" % (B % r"\Lambda"))
rep(128, r"&  +\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{8\pi^{5}}", r"&  %s\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{8\pi^{5}}" % (B % '-'))
# ---- plus prescription definition (after each "Using + prescription")
PLUS = (r"{\color{blue}Here $\int_{\Lambda}^{\vee-k^{+}}dp^{+}\big[\frac{1}{p^{+}}\big]_{+}g(p^{+})\equiv\int_{\Lambda}^{\vee-k^{+}}dp^{+}\,\frac{g(p^{+})-g(0)}{p^{+}}$, "
        r"so that $\int_{\Lambda}^{\vee-k^{+}}\frac{dp^{+}}{p^{+}}g(p^{+})=\int_{\Lambda}^{\vee-k^{+}}dp^{+}\big[\frac{1}{p^{+}}\big]_{+}g(p^{+})+g(0)\log\frac{\vee-k^{+}}{\Lambda}$ exactly. "
        r"The lower limit stays $\Lambda$ for all terms: extending it to $0$ drops $\int_0^\Lambda dp^+(\dots)$, which is not small unless $\Lambda\ll k^+$ and $\Lambda|\mathbf k||y'-z|\ll k^+$.}")
for n in (74, 88, 107, 126):
    assert 'prescription' in L[n - 1], (n, L[n - 1])
    L[n - 1] = L[n - 1] + ' ' + PLUS
# ---- header with the list of corrections
HEAD = r"""{\color{blue}\noindent\textbf{Corrections (blue; checked against Calculation\_main.tex, see Groups\_II\_III\_check.pdf):}
\begin{itemize}
\item Rows III, V, VI: the overall sign. In Eq.~(CS2) of Calculation\_main.tex these rows enter as $-4\,\mathfrak C'^\dagger\mathfrak B$, $-2\,\mathfrak C'^\dagger\mathbb X_1'$, $-2\,\mathfrak C'^\dagger\mathbb X_2'$, because the splitting amplitude after the target is $\mathfrak C=-\int_y[\bar{\mathbb C}U^\dagger U\mathbb A-U^\dagger\bar{\mathbb C}\bar{\mathbb A}]$ ($\bar{\mathbb C}^\dagger_2=-\bar{\mathbb C}_1$, $\mathbb A^{\dagger}_1=-\mathbb A_2$). Every line of these rows changes sign, including the $\log(\vee/\Lambda)$ terms.
\item Plus prescription: the lower limit of the $p^+$ integral stays $\Lambda$ (exact); with $0$ the result is only approximate.
\item Missing or extra integration variables in the measures.
\item Not changed here: index placements of the adjoint Wilson lines. They are correct with your convention $U^{bc}=U^{cb}$; for a general adjoint $U$ ($U^\dagger=U^T\neq U$) see Sec.~5 of the check note.
\end{itemize}}"""
i = L.index(r"\maketitle")
L.insert(i + 1, HEAD)
T = '\n'.join(L)
T = T.replace(r"\title{Group - III crossection derivation}", r"\title{Group - III crossection derivation {\color{blue}(corrected)}}")
open('Group_III_corrected.tex', 'w').write(T)
print('ok')
