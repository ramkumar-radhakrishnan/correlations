r"""Your BFKL.tex (input/) with the corrections applied.
  python3 make_corrected.py          -> corrected/BFKL.tex          clean (your colours kept)
  python3 make_corrected.py marked   -> corrected/BFKL_marked.tex   red = remove, blue = new (your blue -> teal)
Every edit is asserted."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'corrected')
MARKED = 'marked' in sys.argv
CHK = open(os.path.join(HERE, '..', 'full_check', 'Checks_JIMWLK_BFKL.tex')).read()


def old(s):
    return (r'{\color{red}%s}' % s) if MARKED else ''


def new(s):
    return (r'{\color{blue}%s}' % s) if MARKED else s


def g1ii_b():
    """the log(vee/k+) terms of Group I Row II (eq:B:G1II2a of Checks_JIMWLK_BFKL.tex); the log form has no z left"""
    i = CHK.index(r'\label{eq:B:G1II2a}')
    a = CHK.rindex(r'\begin{align}', 0, i)
    t = CHK[a + len(r'\begin{align}'):CHK.index(r'\end{align}', i)].replace(r'\label{eq:B:G1II2a}', '')
    t = t.replace(r'\nc', '').strip()
    k = t.index(r'\log\frac{(w-y)^{2}}{(x-y)^{2}}')
    j = t.rindex("\\int_{x,y,x',z}", 0, k)
    t = t[:j] + "\\int_{x,y,x'}" + t[j + len("\\int_{x,y,x',z}"):]
    return t


t = open(os.path.join(HERE, 'input', 'BFKL.tex')).read()
items = []

# C1: Group I Row II, first form
L = t.split('\n')
assert L[104].startswith(r'\subsection{Row II}') and L[105].strip() == r'\begin{align}' and L[111].strip() == r'\end{align}'
body_old = '\n'.join(L[106:111])
body_new = g1ii_b()
if MARKED:
    blk = r'\begin{align}' + '\n' + body_new + '\n' + r'\end{align}'
    rep = ['{\\color{red}\\noindent Remove (this is Group I Row III, $d^3N_1$):}',
           r'{\color{red}\begin{align}', body_old, r'\end{align}}',
           r'{\color{blue}\noindent The terms $\propto\log\frac{\vee}{k^+}$ of Group I Row II:', blk + '}']
    L[105:112] = rep
else:
    L[106:111] = body_new.split('\n')
t = '\n'.join(L)
items.append(r"Group I, Row II: the first form was a copy of Row III, $d^3N_1$. It is replaced by the $\log\frac{\vee}{k^+}$ "
             r"terms of Row II (Checks\_JIMWLK\_BFKL.pdf, eq.\ (B:G1II2a)); the logarithmic term has no $z$ left: "
             r"$\int_{x,y,x'}$. Your ``This can be written as'' list is the reordering of these terms (exact).")

# C2: Group II Row IV, third line: + sign
a = r"&   \frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}"
assert t.count(a) == 1
t = t.replace(a, r"&   " + new('+') + r"\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}")
items.append(r"Group II, Row IV, third line ($H_4$): the sign $+$ is missing at the start of the line (it was also missing "
             r"in my Checks\_JIMWLK\_BFKL.pdf; corrected there).")

# C3: singlet form measure
a = r"\int_{x,,z,w,w'}"
assert t.count(a) == 1
t = t.replace(a, r"\int_{x," + new('y') + r",z,w,w'}")
items.append(r"Colour-singlet form: measure $\int_{x,,z,w,w'}\to\int_{x,y,z,w,w'}$.")
items.append(r"Everything else is correct: the sum of all two-$\rho$ lists and both forms of the BFKL contribution equal "
             r"BFKL$\otimes$LO (see BFKL\_check.pdf).")

if MARKED:
    for a_, b_ in (('blue', 'teal'), ('green', 'olive'), ('orange', 'brown')):
        t = t.replace(r'\color{%s}' % a_, r'\color{%s}' % b_)
    # restore the correction colours
    t = t.replace(r'{\color{teal}\noindent', r'{\color{blue}\noindent')
    t = t.replace(r'\color{teal}+', r'\color{blue}+').replace(r'\color{teal}y', r'\color{blue}y')
    head = (r'{\color{blue}\noindent\textbf{Corrections (blue = new, red = remove; see BFKL\_check.pdf):}' '\n'
            r'\begin{itemize}' '\n' + '\n'.join(r'\item ' + x for x in items) + '\n' r'\end{itemize}}' '\n')
    i = t.index(r'\maketitle') + len(r'\maketitle')
    t = t[:i] + '\n' + head + t[i:]
t = t.replace(r'\usepackage{xcolor}', r'\usepackage[dvipsnames]{xcolor}')
os.makedirs(OUT, exist_ok=True)
open(os.path.join(OUT, 'BFKL_marked.tex' if MARKED else 'BFKL.tex'), 'w').write(t)
print('written', 'marked' if MARKED else 'clean')
