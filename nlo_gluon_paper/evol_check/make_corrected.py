r"""Your Evolution_two_rho.tex, Evolution_three_rho.tex and Evolution_terms.tex (input/) with every correction applied.
  python3 make_corrected.py          -> corrected/<file>.tex        clean: corrections applied, your colours kept
  python3 make_corrected.py marked   -> corrected/<file>_marked.tex red = remove, blue = new; your own colours remapped
                                                                     blue -> teal, green -> olive, red -> violet,
                                                                     orange -> brown
Every edit is asserted (the text it replaces must be there, the given number of times).  The clean copies keep the
line numbers of your files, so cmp_two.py, cmp_three.py and cmp_list.py can be re-run on them (EVOL_DIR=corrected)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'corrected')
MARKED = 'marked' in sys.argv
R0, B0 = '@@RED@@', '@@BLUE@@'      # placeholders: replaced after your own colours are remapped


def old(s):
    return (r'{\color{%s}%s}' % (R0, s)) if MARKED else ''


def new(s):
    return (r'{\color{%s}%s}' % (B0, s)) if MARKED else s


def remap(t):
    for a, b in (('blue', 'teal'), ('green', 'olive'), ('red', 'violet'), ('orange', 'brown')):
        t = t.replace(r'\color{%s}' % a, r'\color{%s}' % b)
    return t


def balanced(m):
    d = 0
    for ch in m.replace('\\{', '').replace('\\}', ''):
        d += (ch == '{') - (ch == '}')
        if d < 0:
            return False
    return d == 0 and m.count('\\left') == m.count('\\right') and not m.endswith('\\')


def mark(a, b):
    """common prefix + old(differing middle of a) + new(differing middle of b) + common suffix; the middles are
    brace-balanced and do not cut a command name, else the whole of a and b is marked"""
    if not MARKED:
        return b
    p = 0
    while p < min(len(a), len(b)) and a[p] == b[p]:
        p += 1
    q = 0
    while q < min(len(a), len(b)) - p and a[len(a) - 1 - q] == b[len(b) - 1 - q]:
        q += 1
    while p > 0 or q > 0:
        ma, mb = a[p:len(a) - q], b[p:len(b) - q]
        cut = (p > 0 and (a[p - 1] in '\\^_' or (a[p - 1].isalpha() and ((ma[:1] or 'x').isalpha() or (mb[:1] or 'x').isalpha())))) or \
              (q > 0 and ((ma[-1:] or ' ').isalpha() or (mb[-1:] or ' ').isalpha()) and a[len(a) - q].isalpha())
        if balanced(ma) and balanced(mb) and not cut:
            return a[:p] + (old(ma) if ma else '') + (new(mb) if mb else '') + a[len(a) - q:]
        if p >= q and p > 0:
            p -= 1
        else:
            q -= 1
    return old(a) + new(b)


class File:
    def __init__(self, name):
        self.name = name
        self.t = open(os.path.join(HERE, 'input', name + '.tex')).read()
        self.n_lines = self.t.count('\n')
        self.items = []

    def edit(self, a, b, count=1, line=None, keep='', brk=False):
        """replace the text a by b (count times; or only on the given 1-based line); keep is not used (mark() finds
        the common prefix and suffix itself)"""
        if line is not None:
            L = self.t.split('\n')
            assert L[line - 1].count(a) >= 1, (self.name, line, a[:80])
            m = (old(a) + r'\notag \\ & ' + new(b)) if (brk and MARKED) else mark(a, b)
            L[line - 1] = L[line - 1].replace(a, m, 1)
            self.t = '\n'.join(L)
            return
        assert self.t.count(a) == count, (self.name, a[:80], self.t.count(a))
        self.t = self.t.replace(a, mark(a, b))

    def write(self):
        t = self.t
        if MARKED:
            t = remap(t).replace(R0, 'red').replace(B0, 'blue')
            head = (r'{\color{blue}\noindent\textbf{Corrections (blue = new, red = remove; see Evolution\_check.pdf):}'
                    '\n' r'\begin{itemize}' '\n' + '\n'.join(r'\item ' + x for x in self.items) + '\n'
                    r'\end{itemize}}' '\n')
            i = t.index(r'\maketitle') + len(r'\maketitle')
            t = t[:i] + '\n' + head + t[i:]
        else:
            assert t.count('\n') == self.n_lines, self.name     # same line numbers as your file
        t = t.replace(r'\usepackage{xcolor}', r'\usepackage[dvipsnames]{xcolor}')
        os.makedirs(OUT, exist_ok=True)
        open(os.path.join(OUT, self.name + ('_marked' if MARKED else '') + '.tex'), 'w').write(t)


# ------------------------------------------------------------------------------------------------ Evolution_two_rho
two = File('Evolution_two_rho')
two.edit(r"\int_{w, w',x,x',}", r"\int_{w, w',x,x'}", keep=r"\int_{w, w',x,x'")
two.edit(r"\int_{w,y',z,y,x}\, e^{-ik\cdot (y'-w)}", r"\int_{w,x',y',z,x}\, e^{-ik\cdot (y'-w)}")
two.items = [r"Group I, Row II: stray comma in the measure, $\int_{w,w',x,x',}\to\int_{w,w',x,x'}$.",
             r"Group III, Row III, second part: measure $\int_{w,y',z,y,x}\to\int_{w,x',y',z,x}$ (there is no $y$; "
             r"$x'$ is integrated).",
             r"Everything else is correct: every row equals the $\log(k^+/\Lambda)$ two-$\rho$ terms of your Group rows "
             r"exactly (Evolution\_check.pdf, Table~1)."]
two.write()

# ---------------------------------------------------------------------------------------------- Evolution_three_rho
thr = File('Evolution_three_rho')
# Group I Row II: the log form has no z left
thr.edit(r"\int_{x,y,x',z}\frac{(x'-w')^{i}}{(x'-w')^{2}}\,\Bigg[\frac{(x-w)^{i}}{(x-w)^{2}}\,\log",
         r"\int_{x,y,x'}\frac{(x'-w')^{i}}{(x'-w')^{2}}\,\Bigg[\frac{(x-w)^{i}}{(x-w)^{2}}\,\log", keep=r"\int_{x,y,x'")
# Group I Row III, d^3N_2 (first form): signs of the first two terms; the edge kernel -> D (F3)
thr.edit(r"\notag \\ &  \frac{ig^{4}\, f^{abd}}{8\pi^{5}}", r"\notag \\ &  -\frac{ig^{4}\, f^{abd}}{8\pi^{5}}",
         keep=r"\notag \\ &  ")
thr.edit(r"\notag \\ &   - \frac{ig^{4}\, f^{abd}}{8\pi^{5}}", r"\notag \\ &   + \frac{ig^{4}\, f^{abd}}{8\pi^{5}}",
         keep=r"\notag \\ &   ")
thr.edit(r"\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}(z-w)^{2}}",
         r"\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot(w-z)}{(x'-w')^{2}\,(y-w)^{2}\,(x-z)^{2}\,(w-z)^{2}}", line=59)
thr.edit(r"U^{bc}(z)", r"U^{bc}(w)", line=59)
# Group I Row III, d^3N_2 (rewritten form) = D
thr.edit(r"\frac{ig^{4}\, f^{abd}}{16\pi^{5}}", r"\frac{ig^{4}\, f^{abd}}{8\pi^{5}}", line=68)
thr.edit(r"\mathcal{K}(x'-w')\cdot \mathcal{K}(y-z)\, \frac{\left[(x-w)^{2}+(x-z)^{2}\right]}{\left[(x-w)^{2}-(x-z)^{2}\right]^{2}}",
         r"\mathcal{K}(x'-w')\cdot \mathcal{K}(y-w)\, \mathcal{K}(x-z)\cdot \mathcal{K}(w-z)", line=68, brk=True)
thr.edit(r"U^{bc}(z)", r"U^{bc}(w)", line=68)
# Group I Row V: comma
thr.edit(r"\int_{w,z,y,x',y'w'}", r"\int_{w,z,y,x',y',w'}", keep=r"\int_{w,z,y,x',y'")
# Group II Row I, rewritten form: K.K.K.K -> K.K K.K
thr.edit(r"\mathcal{K}(x'-w')\cdot \mathcal{K}(w-z)", r"\mathcal{K}(x'-w')\, \mathcal{K}(w-z)", count=4,
         keep=r"\mathcal{K}(x'-w')")
# consolidated list: four-rho part with log(k+/Lambda) only (this is the JIMWLK file)
thr.edit(r"\left[\log{\left[\frac{k^{+}}{\Lambda}\right]}+\log{\left[\frac{\vee}{k^{+}}\right]}\right]",
         r"\log{\left[\frac{k^{+}}{\Lambda}\right]}", count=2)
# consolidated list, Group I Row II d^3N_2: the closing bracket
L = thr.t.split('\n')
assert L[206].strip() == r'\end{align}' and r'\Bigg]' not in L[205][-40:], L[205][-60:]
L[205] = L[205].rstrip() + new(r'\Bigg]')
thr.t = '\n'.join(L)
# consolidated list, Group III Row III: '\-' prints nothing
thr.edit(r"&\-\mathcal{K}", r"& -\mathcal{K}", count=2, keep='&')
thr.items = [
    r"Group I, Row II: the second (logarithmic) term has no $z$ left: $\int_{x,y,x',z}\to\int_{x,y,x'}$.",
    r"Group I, Row III, $d^3N_2$: the first two terms have the wrong sign (in Group\_I.tex the bracket has the "
    r"prefactor $-\frac{1}{(2\pi)^3}$). With the right signs they cancel the fourth and fifth terms exactly (same "
    r"kernels). The remaining edge-kernel term is $\mathcal D$ (fix F3): $U^{bc}(z)\to U^{bc}(w)$ and kernel "
    r"$\mathcal K(x'-w')\cdot\mathcal K(y-w)\,\mathcal K(x-z)\cdot\mathcal K(w-z)$.",
    r"Group I, Row III, rewritten $d^3N_2$: as written it is the $[(x-w)^2+(x-z)^2]$ term, which cancels inside Row III "
    r"($+\log-\log$). It must be $\mathcal D$, which is what your list ``All three $\rho$ term with c.c.'' already uses.",
    r"Group I, Row V: comma in the measure, $\int_{w,z,y,x',y'w'}\to\int_{w,z,y,x',y',w'}$.",
    r"Group II, Row I, rewritten form: four terms read $\mathcal K\cdot\mathcal K\cdot\mathcal K\cdot\mathcal K$; it is "
    r"$\mathcal K(x-w)\cdot\mathcal K(x'-w')\,\mathcal K(w-z)\cdot\mathcal K(y-z)$.",
    r"``From four $\rho$'s'': only $\log(k^+/\Lambda)$ belongs in this (JIMWLK) file; the $\log(\vee/k^+)$ part of these "
    r"terms goes to the BFKL list, where it cancels too. Your cancellation section already uses $\log(k^+/\Lambda)$ only.",
    r"Genuine terms, Group I Row II, $d^3N_2$: the closing \texttt{\textbackslash Bigg]} is missing.",
    r"Genuine terms, Group III Row III: two terms start with \texttt{\textbackslash-}, which prints nothing; it is a "
    r"minus sign.",
    r"Everything else is correct (Evolution\_check.pdf, Tables~2--4)."]
thr.write()

# ---------------------------------------------------------------------------------------------------- Evolution_terms
trm = File('Evolution_terms')
trm.edit(r"\int_{x,,z,w,w'}", r"\int_{x,y,z,w,w'}", count=2, keep=r"\int_{x,")
trm.edit(r"where $K(y,y,z) = \frac{(y-z)\cdot (y-z)}{(y-z)^{4}}$",
         r"where $K(y,y,z) = \frac{(y-z)\cdot (y-z)}{(y-z)^{4}}$" +
         r", and in general $K(a,b,z)=\frac{(a-z)\cdot(b-z)}{(a-z)^{2}\,(b-z)^{2}}$",
         keep=r"where $K(y,y,z) = \frac{(y-z)\cdot (y-z)}{(y-z)^{4}}$")
i = trm.t.index(r'\begin{comment}')
j = trm.t.index(r'\end{comment}')
blk = trm.t[i:j]
LK = r"{\color{blue}\log{\left(\frac{k^{+}}{\Lambda}\right)}}"
assert blk.count(LK) == 3, blk.count(LK)
trm.t = trm.t[:i] + blk.replace(LK, r"{\color{blue}\log{\left(\frac{\vee}{k^{+}}\right)}}") + trm.t[j:]
trm.items = [
    r"Colour-singlet form: measure $\int_{x,,z,w,w'}\to\int_{x,y,z,w,w'}$.",
    r"Define $K(a,b,z)=\frac{(a-z)\cdot(b-z)}{(a-z)^{2}(b-z)^{2}}$ (only $K(y,y,z)$ is defined).",
    r"In the commented-out BFKL section the logarithm is $\log(\vee/k^+)$, not $\log(k^+/\Lambda)$ (no change in the "
    r"PDF).",
    r"Everything else is correct: Sec.~1 + Sec.~2 = Sec.~3 = JIMWLK$\otimes$LO exactly (Evolution\_check.pdf, Sec.~4)."]
trm.write()
print('written', 'marked' if MARKED else 'clean', 'copies to', OUT)
