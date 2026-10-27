r"""Writes Group_I_corrected.tex, Group_II_corrected.tex, Group_III_corrected.tex: your files (input/) with every
correction of Group_check.pdf marked in place.  Red = remove, blue = new.  Your own colour marks are remapped so that
they cannot be confused with these: blue -> teal, green -> olive, red -> violet, orange -> brown.  Every edit is asserted."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = r'{\color{blue}%s}'
OLD = r'{\color{red}%s}'


def remap(t):
    for a, b in (('blue', 'teal'), ('green', 'olive'), ('red', 'violet'), ('orange', 'brown')):
        t = t.replace(r'\color{%s}' % a, r'\color{%s}' % b)
    return t


def head(title, items):
    return (r'{\color{blue}\noindent\textbf{Corrections (blue = new, red = remove; see Group\_check.pdf):}' + '\n'
            r'\begin{itemize}' + '\n' + '\n'.join(r'\item ' + x for x in items) + '\n' + r'\end{itemize}}' + '\n')


class File:
    def __init__(self, g):
        self.g = g
        self.L = [remap(x) for x in open(os.path.join(HERE, 'input', g + '.tex')).read().split('\n')]

    def rep(self, n, old, new, count=1):
        s = self.L[n - 1]
        assert s.count(old) >= count, (self.g, n, old[:60], s.count(old))
        self.L[n - 1] = s.replace(old, new, count)

    def braces(self, a, b):
        """mark a \\Bigg\\} that closes nothing (red) and close a \\Bigg\\{ that is never closed (blue)"""
        t = '\n'.join(self.L[a - 1:b])
        out, opened, k, n = '', [], 0, 0
        while k < len(t):
            if t.startswith(r'\Bigg\{', k):
                opened.append(k)
                out += r'\Bigg\{'
                k += 7
                continue
            if t.startswith(r'\Bigg\}', k):
                if opened:
                    opened.pop()
                    out += r'\Bigg\}'
                else:
                    out += OLD % r'\Bigg\}'
                    n += 1
                k += 7
                continue
            out += t[k]
            k += 1
        if opened:
            j = out.rindex(r'\end{align}')
            out = out[:j].rstrip() + (NEW % (r'\Bigg\}' * len(opened))) + '\n' + out[j:]
            n += len(opened)
        assert n, (a, b)
        new = out.split('\n')
        assert len(new) == b - a + 1
        self.L[a - 1:b] = new

    def write(self, items):
        i = self.L.index(r'\maketitle')
        self.L.insert(i + 1, head(self.g, items))
        t = '\n'.join(self.L)
        t = t.replace(r'\usepackage{xcolor}', r'\usepackage[dvipsnames]{xcolor}')
        open(os.path.join(HERE, self.g + '_corrected.tex'), 'w').write(t)




# ------------------------------------------------------------------ the BFKL terms hidden in line 1 (shown explicitly)
LVK = r'\log{\left(\frac{\vee}{k^{+}}\right)}'


def hidden_B():
    """the terms (B) of Group III Rows III (both parts), V, VI, as in BFKL_rows.pdf (from ../bfkl_final/static.tex)"""
    S = open(os.path.join(HERE, '..', 'bfkl_final', 'static.tex')).read()
    i = S.index(r'\paragraph{The new terms}')
    blk = S[S.index(r'\begin{align}', i) + len(r'\begin{align}'):S.index(r'\end{align}', i)]
    parts = re.split(r'\\text\{Row ([IV()12]+):\}\\quad&', blk)[1:]
    out = {}
    for k in range(0, len(parts), 2):
        t = parts[k + 1].replace(r'\Lv', LVK).replace(r'\nc', '')
        out[parts[k]] = re.sub(r'[,.]\s*(\\\\)?\s*$', '', t.strip())
    return out


def as_integral(t):
    """-(B) with log(vee/k+) written as (1/2pi) int_0^{vee-k+} dp+ / (p+ + k+)  (exact)"""
    t = t.lstrip()
    assert t.startswith('-')
    t = '+' + t[1:]
    for a, b in ((r'{16\pi^{5}}', r'{8\pi^{4}}'), (r'{8\pi^{5}}', r'{4\pi^{4}}')):
        if a in t:
            t = t.replace(a, b)
            break
    assert t.count(LVK) == 1
    return t.replace(LVK, r'\int_{0}^{\vee-k^{+}}\frac{dp^{+}}{2\pi}\,\frac{1}{(p^{+}+k^{+})}\\ &\quad\times')


def blue_align(lhs, rhs):
    cells = re.split(r'(?<!\\)&', rhs)
    body = r'\color{blue}' + cells[0] + ''.join(r'& \color{blue} ' + c for c in cells[1:])
    return '\\begin{align*}\n{\\color{blue}' + lhs + '} &= ' + body + '\n\\end{align*}'


def bfkl_block(eq, where, B):
    txt = (r'{\color{blue}\noindent\textbf{The BFKL term of this row, written out.} It is inside line 1 of Eq.~' + eq + r', in the ' + where
           + r' term: at large $p^+$ the phase $e^{-ik\cdot\frac{p^+}{k^+}(y^\prime-z)}$ lets only $y^\prime\to z$ contribute to the $y^\prime$ integral, '
           r'which gives $\frac{k^+}{p^++k^+}$ times the value at $y^\prime=z$, and $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$ '
           r'(\texttt{Why\_log\_vee.pdf}). Your equation stays as it is; it equals (B) $+$ [your equation $-$ (B)], with}' + '\n')
    b1 = blue_align(r'\text{(B)}', B)
    b2 = blue_align(r'\text{Eq.~' + eq + r'} - \text{(B)}', r'\text{Eq.~' + eq + r'}\ ' + as_integral(B))
    tail = (r'{\color{blue}\noindent (B) $\propto\log\frac{\vee}{k^+}$ is the BFKL term (it is in \texttt{BFKL\_rows.pdf}); '
            r'Eq.~' + eq + r'$\,-\,$(B) has no $\log\frac{\vee}{k^+}$: it is finite for $\vee\to\infty$. '
            r'The $\log(\frac{\vee}{k^+}-1)$ of your line 2 is not a BFKL term: it cancels the subtraction in $\left[\frac{1}{p^+}\right]_+$ of line 1.}')
    return '\n' + txt + b1 + '\n' + b2 + '\n' + tail + '\n'


H4 = (r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]}\,\int_{x,x',z,w,w'}\,e^{-ik\cdot (w'-w)}\,"
      r"\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}\,(w'-z)^{2}}\,\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}\,(x'-z)^{2}}\notag \\ & \left[U^{e'c}(z) - U^{e'c}(x')\right] \times\, "
      r"\left[U^{ec}(z) - U^{ce}(x)\right]\, \rho^{e'}(x') \rho^{e}(x)")


def h4_block(eq):
    txt = (r'{\color{blue}\noindent\textbf{The BFKL term of this row, written out.} It is inside the $\xi$ integral of line 1 of Eq.~' + eq
           + r', in the $\frac{(1-\xi)}{\xi}$ term: for $\xi\to0$ the phase $e^{-ik\cdot(y^\prime-y)\frac{1}{\xi}}$ lets only $y^\prime\to y$ contribute, '
           r'the $\log\frac{(y^\prime-y)^2\mu^2}{4}$ of the bracket gives $\propto\frac{d\xi}{\xi}$, and $\int_{k^+/\vee}^1\frac{d\xi}{\xi}=\log\frac{\vee}{k^+}$ '
           r'(\texttt{Why\_log\_vee.pdf}). Your equation stays as it is; it equals $H_4$ $+$ [your equation $-$ $H_4$], with}' + '\n')
    b1 = blue_align(r'H_4', H4)
    b2 = blue_align(r'\text{Eq.~' + eq + r'} - H_4', r'\text{Eq.~' + eq + r'}\ -' + H4.replace(r'\log{\left[\frac{\vee}{k^{+}}\right]}', r'\int_{\frac{k^{+}}{\vee}}^{1}\frac{d\xi}{\xi}\\ &\quad\times'))
    tail = (r'{\color{blue}\noindent $H_4\propto\log\frac{\vee}{k^+}$ is the BFKL term (your \texttt{BFKL.tex}, Group~II Row~IV, line 2). '
            r'The $\log(\frac{\vee}{k^+}-1)$ of line 2 of Eq.~' + eq + r' is not a BFKL term: it cancels the subtraction in $\frac{\xi}{(1-\xi)_+}$ of line 1.}')
    return '\n' + txt + b1 + '\n' + b2 + '\n' + tail + '\n'


# ------------------------------------------------------------------ Group I
G1 = File('Group_I')
# Row II 2b: x' is integrated
G1.rep(75, r"\int_{w,w'} e^{-ik\cdot(w'-w)}", r"\int_{w,w'%s} e^{-ik\cdot(w'-w)}" % (NEW % ",x'"))
# Row II part 1: stray comma in the measure
G1.rep(46, r"\int_{w, w',x,x',}", r"\int_{w, w',x,x'%s}" % (OLD % ','))
G1.rep(50, r"\int_{w, w',x,x',}", r"\int_{w, w',x,x'%s}" % (OLD % ','))
# Row V part 1: missing comma in the measure y'w'
for n in (157, 161, 245):
    s = G1.L[n - 1]
    G1.L[n - 1] = s.replace(r"\int_{x,z,y,x',y'w'}", r"\int_{x,z,y,x',y'%sw'}" % (NEW % ','))
    assert G1.L[n - 1] != s, n
# \Bigg\{ never closed (Row III d^3N_1) or \Bigg\} that closes nothing (Row V part 2): display only
for a, b in ((95, 97), (168, 170), (172, 174), (181, 183), (188, 190)):
    G1.braces(a, b)
# Row V 2b UV, "Setting Lambda -> 0" and the final answer: E1
G1.rep(202, r"{\color{teal} \log{\left(\frac{\vee}{k^{+}}\right)} +\log{\left(\frac{k^{+}}{\Lambda}\right)}}-2\log\Big(\frac{\vee}{k^{+}}\Big)",
       (OLD % r"\log{\left(\frac{\vee}{k^{+}}\right)} +\log{\left(\frac{k^{+}}{\Lambda}\right)}-2\log\Big(\frac{\vee}{k^{+}}\Big)")
       + (NEW % r"2\log{\left(\frac{k^{+}}{\Lambda}\right)}"))
G1.rep(253, r"{\color{teal} \log{\left(\frac{k^{+}}{\Lambda}\right)}-\log{\left(\frac{\vee}{k^{+}}\right)}}",
       (OLD % r"\log{\left(\frac{k^{+}}{\Lambda}\right)}-\log{\left(\frac{\vee}{k^{+}}\right)}") + (NEW % r"2\log{\left(\frac{k^{+}}{\Lambda}\right)}"))
G1_ITEMS = [
    r"\textbf{Row V, 2b UV (E1), Eqs.\ after ``Setting $\Lambda\to0$'' and in the final answer:} the bracket is "
    r"$\big[2\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}\big]$. Your exact bracket just before it has "
    r"$2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^+}=2\log\frac{k^+}{\Lambda}$; there is no $\log\frac{\vee}{k^+}$ in Row V ($0<p^+<k^+$).",
    r"Row II, $d^3N_{2b}$: $x'$ is integrated, $\int_{w,w',x'}$; Row II part 1: a stray comma in $\int_{w,w',x,x'}$.",
    r"Display only: in Row III $d^3N_1$ the opening \texttt{\textbackslash Bigg\textbackslash\{} is never closed; in Row V part 2 (``Integrating w.r.t.\ $w$'', ``This can be written as'', $d^3N_{2a}$ and its plus-prescription form) a closing \texttt{\textbackslash Bigg\textbackslash\}} closes nothing.",
    r"Row V part 1: $\int_{x,z,y,x',y',w'}$ (comma).",
    r"Not marked (notation): the logarithms of Row III are principal values, $\log|\cdot|$; the arguments change sign at $p^+=k^+\frac{(x-w)^2}{(x-z)^2}$.",
]
G1.write(G1_ITEMS)

# ------------------------------------------------------------------ Group II
G2 = File('Group_II')
G2.rep(44, 'We write it in two parts:', OLD % 'We write it in two parts:')
G2.rep(98, r"\frac{(k^{+}+p^{+})}{k{+}}", r"\frac{(k^{+}+p^{+})}{%s}" % (OLD % 'k{+}') + '')
G2.L[97] = G2.L[97].replace(r"{%s}" % (OLD % 'k{+}'), r"{%s%s}" % (OLD % 'k{+}', NEW % 'k^{+}'))
assert NEW % 'k^{+}' in G2.L[97]
G2.rep(112, r"e^{-ik\,\left(y'-y\right)\,\frac{p^{+}}{(p^{+}+k^{+})}}",
       r"e^{-ik\,\left(y'-y\right)\,%s%s}" % (OLD % r"\frac{p^{+}}{(p^{+}+k^{+})}", NEW % r"\frac{p^{+}+k^{+}}{k^{+}}"))
G2_ITEMS = [
    r"Row IV, after the $z$ integral (``This gives us''): the phase of the second term is $e^{-ik\cdot(y'-y)\frac{p^++k^+}{k^+}}$, "
    r"as in the line before and in the next equation (only a typo: the next equation is right).",
    r"Row IV, ``Integrating $w$ and $w'$'': $\frac{(k^++p^+)}{k^+}$ in the phase (typo $k\{+\}$).",
    r"Row I: the sentence ``We write it in two parts'' appears twice.",
]
assert G2.L[132].strip() == r'\end{align}', G2.L[132]
G2.L[132] = G2.L[132] + h4_block('(1.21)')
G2_ITEMS.append(r"Shown in blue after Eq.~(1.21) (nothing changes): the BFKL term $\log\frac{\vee}{k^+}H_4$ that sits inside the $\xi$ integral of line 1.")
G2.write(G2_ITEMS)

# ------------------------------------------------------------------ Group III
G3 = File('Group_III')
G3.rep(63, r"& =   \frac{4}{(2\pi)^{3}}", r"& = %s\,\frac{4}{(2\pi)^{3}}" % (NEW % '-'))
G3.rep(80, r"\int_{w,w',y',z,y,x}", r"\int_{w,w',%sy',z,%sx}" % (NEW % "x',", OLD % 'y,'))
G3.rep(85, r"\int_{w,w',y',z,y,x}", r"\int_{w,%s%sy',z,%sx}" % (OLD % "w',", NEW % "x',", OLD % 'y,'))
G3.rep(90, r"\int_{w,w',y',z,y,x}", r"\int_{w,%s%sy',z,%sx}" % (OLD % "w',", NEW % "x',", OLD % 'y,'))
G3.rep(91, r"\int_{w,y',z,y,x}", r"\int_{w,%sy',z,%sx}" % (NEW % "x',", OLD % 'y,'))
G3_ITEMS = [
    r"Row III, first line ($\mathfrak A^\dagger\bar{\mathfrak C}^\dagger\mathfrak B$ form): the overall sign is $-$, as you already have in "
    r"Rows V and VI (Rows III, V, VI enter with $-\mathfrak C'^\dagger$). The explicit forms below it are right.",
    r"Row III, second part (C5): the measure is $\int_{w,w',x',y',z,x}$ before and $\int_{w,x',y',z,x}$ after the $w'$ integral "
    r"($y$ does not appear in this part; $x'$ does).",
]
HB = hidden_B()
for end, eq, key, where in ((77, '(1.7)', 'III(1)', r'$-\frac{1}{k^+}\,\delta_{im}\delta_{jk^\prime}$'),
                            (93, '(1.10)', 'III(2)', r'$-\frac{p^+}{k^{+2}}$'),
                            (114, '(1.14)', 'V', r'$-\frac{1}{k^+}\,\delta_{im}\delta_{jk^\prime}$'),
                            (129, '(1.17)', 'VI', r'$-\frac{1}{k^+}\,\delta_{im}\delta_{jk^\prime}$')):
    assert G3.L[end - 1].strip() == r'\end{align}', (end, G3.L[end - 1])
    G3.L[end - 1] = G3.L[end - 1] + bfkl_block(eq, where, HB[key])
G3_ITEMS.append(r"Shown in blue after Eqs.~(1.7), (1.10), (1.14), (1.17) (nothing changes): the BFKL term (B) $\propto\log\frac{\vee}{k^+}$ "
                r"of Rows III (both parts), V, VI, which sits inside line 1, and each row as (B) $+$ a part that is finite for $\vee\to\infty$.")
G3.write(G3_ITEMS)

if __name__ == '__main__':
    print('written Group_I_corrected.tex, Group_II_corrected.tex, Group_III_corrected.tex')
