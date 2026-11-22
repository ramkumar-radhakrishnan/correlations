# Post-processing step (author's request 2026-11-22): remove colours from Section 10
# ("What remains, with every three-rho product symmetrised") of Final_terms_explicit.tex.
# Run after assemble.py:  python3 -I scripts/colourfree_sec10.py Final_terms_explicit.tex
import re, sys
f = sys.argv[1]; s = open(f).read()
a = s.index(r'\section{What remains, with every three-$\rho$ product symmetrised')
b = s.index(r'\section{Average over the projectile: colour-singlet form')
sec = re.sub(r'\\color\{[A-Za-z]+\}\s?', '', s[a:b]).replace(r'\ph{', '{')
sec = sec.replace(r' (coloured as in Sec.~\ref{sec:C})', r' (Sec.~\ref{sec:C})').replace(r' (blue) ', ' ')
if 'This section is printed without colours' not in sec:
    i = sec.index('\n') + 1
    sec = sec[:i] + (r"\noindent\textit{This section is printed without colours. Terms that are integrable only in combination are identified in the paragraph below and in Sec.~\ref{sec:C}; the remnant term of Group~II Row~IV is not in the author's file (Sec.~\ref{sec:A}, A5); the phases are $e^{-ik\cdot(\text{unprimed}-\text{primed})}$ throughout.}\par" + "\n") + sec[i:]
open(f, 'w').write(s[:a] + sec + s[b:])
