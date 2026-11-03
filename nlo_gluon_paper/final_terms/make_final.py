r"""Builds Final_result.tex: front.tex (new: checks of Other_terms.tex, DGLAP/RC, leading power, final result,
exact-a dictionary, open points) + the verified simplified forms (../limiting/Limiting_terms.tex, Sec. 3) + the
Gaussian/colour-singlet forms (../limiting/Gaussian_singlet.tex, Secs. 2-5) + the AABKL comparison, with labels
prefixed so that nothing collides.   python3 make_final.py"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
LIM = os.path.join(HERE, '..', 'limiting')


def chunk(fn, start, stop):
    t = open(os.path.join(LIM, fn)).read()
    i = t.index(start)
    j = t.index(stop, i) if stop else t.index(r'\end{document}')
    return t[i:j]


def prefix(t, p):
    t = re.sub(r'\\label\{([^}]*)\}', lambda m: r'\label{%s:%s}' % (p, m.group(1)), t)
    t = re.sub(r'\\(eqref|ref)\{([^}]*)\}', lambda m: r'\%s{%s:%s}' % (m.group(1), p, m.group(2)), t)
    return t


lt = prefix(chunk('Limiting_terms.tex', r'\subsection{The c.c.}', r'\subsection{Final result}'), 'LT')
lt = lt.replace(r'\subsection{', r'\subsection{').replace('(LT Eq.', '(Eq.')
aab = prefix(chunk('Limiting_terms.tex', r'\section{Comparison with forward rapidity', r'\section{Files}'), 'LT')
aab = aab.replace(r'\section{', r'\section{')
gs = prefix(chunk('Gaussian_singlet.tex', r'\section{Rules}', r'\section{Final result in colour-singlet form'), 'GS')
rep = {r'(LT Eq.~(3.13), weights': r'(Eq.~\eqref{LT:eq:G3}, weights',
       r'(LT (3.6), first two lines)': r'(Eq.~\eqref{LT:eq:2ab}, first two lines)',
       r'With LT Eq.~(3.4) the': r'With Eq.~\eqref{LT:eq:V1cf} the',
       r'\section{Kernels (from LT, $d_\perp=2$)}': r'\section{Kernels ($d_\perp=2$, leading power)}',
       r'the kernels are those of LT.': r'the kernels are those of Sec.~\ref{sec:simp}.',
       r'(LT: unchanged)': '(unchanged)',
       r'(from LT, ': '(', }
for k, v in rep.items():
    gs = gs.replace(k, v)
assert 'LT (' not in gs and 'LT Eq' not in gs, [m.start() for m in re.finditer('LT', gs)]
gs = gs.replace(r'\section{', r'\section{')

front = open(os.path.join(HERE, 'front.tex')).read()
back = open(os.path.join(HERE, 'back.tex')).read()
doc = front.replace('%%SIMP%%', lt).replace('%%SINGLET%%', gs).replace('%%AABKL%%', aab) + back
open(os.path.join(HERE, 'Final_result.tex'), 'w').write(doc)
print('written', len(doc))
