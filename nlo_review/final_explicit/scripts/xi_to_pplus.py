# Post-processing step (author's request 2026-11-22): in Group I Row V (second part), all three
# appearances (explicit list, symmetrised list, colour-singlet list) write xi = p+/k+ explicitly.
# Run after assemble.py (and colourfree_sec10.py):  python3 -I scripts/xi_to_pplus.py Final_terms_explicit.tex
import re, sys
f = sys.argv[1]; s = open(f).read()
s = re.sub(r'\\frac\{(e\^\{-ik\\cdot(?:\\ph)?\{\(x-z\)\}\\frac\{p\^\{\+\}\}\{k\^\{\+\}\}\})\}\{\\xi\\,\(y-x\)\^\{2\}\+\(1-\\xi\)\\,\(y-z\)\^\{2\}\}',
           lambda m: r'\frac{k^{+}\,' + m.group(1) + r'}{p^{+}\,(y-x)^{2}+(k^{+}-p^{+})\,(y-z)^{2}}', s)
s = s.replace(r'\xi(1-\xi)\,\frac{(y-x)', r'\frac{p^{+}(k^{+}-p^{+})}{k^{+\,2}}\,\frac{(y-x)')
s = s.replace(r'+(1-\xi)\,\mathcal{K}(x\'-w\')', r'+\frac{k^{+}-p^{+}}{k^{+}}\,\mathcal{K}(x\'-w\')')
s = s.replace(r'+\xi\,\mathcal{K}(x\'-w\')', r'+\frac{p^{+}}{k^{+}}\,\mathcal{K}(x\'-w\')')
s = s.replace(r'With $\xi\equiv p^{+}/k^{+}$: the', 'The')
s = s.replace(r'/[\xi(y-x)^{2}+(1-\xi)(y-z)^{2}]', r'\,k^{+}/[p^{+}(y-x)^{2}+(k^{+}-p^{+})(y-z)^{2}]')
s = s.replace(r'[\xi(1-\xi)-2]=-\frac{11}{6}', r'\Big[\frac{p^{+}(k^{+}-p^{+})}{k^{+\,2}}-2\Big]=-\frac{11}{6}')
open(f, 'w').write(s)
