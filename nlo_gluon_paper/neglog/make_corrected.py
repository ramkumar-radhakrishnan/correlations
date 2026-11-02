r"""Your Other_terms.tex (input/) with every logarithm whose argument can change sign written as log|.|
(principal value of the p+ integral, as in your Fourier transforms (i3)-(i5)).
  python3 make_corrected.py          -> corrected/Other_terms.tex          (clean)
  python3 make_corrected.py marked   -> corrected/Other_terms_marked.tex   (red = removed, blue = new)
Every replacement is counted and asserted."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MARKED = 'marked' in sys.argv
t = open(os.path.join(HERE, 'input', 'Other_terms.tex')).read()

# (old text, inner fraction, expected count): the log argument changes sign where the p+ pole is inside the range
REPL = [
    (r"\log{\left[\frac{(x-z)^{2}-\frac{k^{+}}{\vee}\,(x-w)^{2}}{(x-z)^{2}-(x-w)^{2}}\right]}",
     r"\frac{(x-z)^{2}-\frac{k^{+}}{\vee}\,(x-w)^{2}}{(x-z)^{2}-(x-w)^{2}}", 1),
    (r"\log{\left[\frac{(x-w)^{2}-(x-z)^{2}}{(x-w)^{2}}\right]}",
     r"\frac{(x-w)^{2}-(x-z)^{2}}{(x-w)^{2}}", 7),
    (r"\log{\left[\frac{(x-z)^{2}-\frac{k^{+}}{\vee}\,(x-w)^{2}}{(x-w)^{2}}\right]}",
     r"\frac{(x-z)^{2}-\frac{k^{+}}{\vee}\,(x-w)^{2}}{(x-w)^{2}}", 1),
    (r"\log{\left[\frac{(x-w)^{2}\,\frac{k^{+}}{\vee}-(x-z)^{2}}{(x-w)^{2}-(x-z)^{2}}\right]}",
     r"\frac{(x-w)^{2}\,\frac{k^{+}}{\vee}-(x-z)^{2}}{(x-w)^{2}-(x-z)^{2}}", 6),
]
n_tot = 0
for old, frac, n in REPL:
    c = t.count(old)
    assert c == n, (old[:60], c, n)
    new = r"\log{\left|" + frac + r"\right|}"
    if MARKED:
        new = (r"\log{\color{red}\Bigg[}\!{\color{blue}\Bigg|}" + frac
               + r"{\color{blue}\Bigg|}\!{\color{red}\Bigg]}")
    t = t.replace(old, new)
    n_tot += c
if MARKED:
    note = (r"{\color{blue}\noindent\textbf{Logarithms with arguments that change sign} (red brackets = removed, blue bars = new): "
            r"%d logarithms of Group I Row III are now $\log|\cdot|$, the principal value of their $p^+$ integral "
            r"(see Negative\_logs.pdf).}" % n_tot) + '\n'
    i = t.index(r'\maketitle') + len(r'\maketitle')
    t = t[:i] + '\n' + note + t[i:]
os.makedirs(os.path.join(HERE, 'corrected'), exist_ok=True)
fn = 'Other_terms_marked.tex' if MARKED else 'Other_terms.tex'
open(os.path.join(HERE, 'corrected', fn), 'w').write(t)
print('written', fn, n_tot, 'logs')
