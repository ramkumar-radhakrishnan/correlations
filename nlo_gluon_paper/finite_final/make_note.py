"""Builds Finite_terms_final.tex from static.tex:
 - the table of Sec. 2 from compare_file.json,
 - the C3 and C4 numbers from checks_xi.json and edge_check.json,
 - the running-coupling and DGLAP sections from ../dglap_rc/Finite_terms_DGLAP_RC.tex (checked again for your file in
   checks_rc.json),
 - the finite terms from the corrected copy (make_corrected.py): red terms dropped, notes shortened, Group I Row II
   part 1 reduced to V_vee, Row V 2b UV and Group II Row IV lines 1, 2 (DGLAP) taken out."""
import json
import os
import re

import make_corrected as MC

HERE = os.path.dirname(os.path.abspath(__file__))
S = open(os.path.join(HERE, 'static.tex')).read()


def brace_end(t, j):
    """index after the brace group that opens at t[j] == '{'"""
    d, k = 1, j + 1
    while d:
        d += {'{': 1, '}': -1}.get(t[k], 0)
        k += 1
    return k


def drop_cmd(t, cmd, keep):
    while cmd + '{' in t:
        i = t.index(cmd + '{')
        k = brace_end(t, i + len(cmd))
        t = t[:i] + (t[i + len(cmd) + 1:k - 1] if keep else '') + t[k:]
    return t


def once(t, old, new):
    assert t.count(old) == 1, (old[:70], t.count(old))
    return t.replace(old, new)


# ---------------------------------------------------------------- Sec. 2 table
CF = json.load(open(os.path.join(HERE, 'compare_file.json')))
CH = {'G2 IV': 'C1', 'G3 III(1)': 'C2, C3', 'G3 III(2)': 'C2, C3, C5', 'G3 V': 'C2, C3', 'G3 VI': 'C2, C3',
      'G1 III N3a': 'C4 (first term)', 'G1 III N3b': 'C4 (last term)'}
LAB = {'G1': 'Group I', 'G2': 'Group II', 'G3': 'Group III'}
rows = []
for name, res in CF.items():
    g, r = name.split(' ', 1)
    what = 'identical' if res == 'identical' else 'differs: ' + ('the $-\\log\\frac{\\vee}{k^+}H_4$ line is missing'
                                                                if name == 'G2 IV' else '$^{(k^+)}$ missing')
    rows.append('%s Row %s & %s & %s\\\\' % (LAB[g], r.replace('^i', '$^i$'), what, CH.get(name, '---')))
NID = sum(1 for v in CF.values() if v == 'identical')
S = S.replace('@@TABLE@@', '\n'.join(rows))
NCH = len(set(CH) | {n for n, v in CF.items() if v != 'identical'})
S = S.replace('@@NID@@', str(len(CF) - NCH)).replace('@@NBLK@@', str(len(CF)))

# ---------------------------------------------------------------- C3, C4, xi numbers
X = json.load(open(os.path.join(HERE, 'checks_xi.json')))
sci = lambda v: ('%.0e' % v).replace('e-0', 'e-').replace('e+0', 'e') if v else '0'
tex = lambda v: sci(v).replace('e', r'\times10^{') + '}' if v and 'e' in sci(v) else sci(v)
S = S.replace('@@ALONE@@', ', '.join('%.2f' % v for v in X['-(1/k+) term alone for vee = 1e2, 1e4, 1e6 (grows like -log)']))
S = S.replace('@@SUM@@', ', '.join('%.4f' % v for v in X['[1/(1-xi)]_+ term for vee = 1e2, 1e4, 1e6 (finite)']))
S = S.replace('@@D1@@', tex(X['[1/p+]_+^(k+) -> [1/xi]_+']))
S = S.replace('@@D2@@', tex(X['-(1/k+) F3 + (B) line -> -[1/(1-xi)]_+']))
S = S.replace('@@DB@@', ', '.join(tex(v) for v in X['xi -> 1: (1-xi) x term = (B) integrand, max rel. diff at 1-xi = 1e-4, 1e-6, 1e-8']))
S = S.replace('@@DV@@', tex(X['Row V (III(1), VI): k+ x p+ integrand = xi integrand, max rel. diff']))
S = S.replace('@@D32@@', tex(X['Row III(2): k+ x p+ integrand = xi integrand, max rel. diff']))
E = json.load(open(os.path.join(HERE, 'edge_check.json')))['z_to_w']['3']
S = S.replace('@@EDGE@@', '\n'.join(r'$10^{%d}$ & $%s$ & $%s$ & $%s$\\' % (round(__import__('math').log10(q['r'])),
                                    tex(q['r2_avg']), tex(q['r2_avg_over_log']), tex(q['without_last_N3b_r_avg']))
                                    for q in E))

# ---------------------------------------------------------------- running coupling, DGLAP
D = open(os.path.join(HERE, '..', 'dglap_rc', 'Finite_terms_DGLAP_RC.tex')).read()
D = D[D.index(r'\section{Running coupling}'):D.index(r'\section{The finite terms}')]
D = D.replace(r'\texttt{checks.json}', r'\texttt{checks\_rc.json}')
D = once(D, 'The poles of the remaining terms', 'The poles of your rows')
D = once(D, 'Group II Row IV, lines 1 and 2 of the remaining terms', 'Group II Row IV, lines 1 and 2 of your file')
D = once(D, r'together with the $\log\frac{\vee}{k^+}H_4$ subtraction:',
         r'together with the $\log\frac{\vee}{k^+}H_4$ subtraction (C1):')
S = S.replace('@@RCDGLAP@@', D)

# ---------------------------------------------------------------- the finite terms
T = '\n'.join(MC.L)
T = T[T.index(r'\section{Group I (with c.c.)}'):T.index(r'\end{document}')]
# notes of the corrected copy -> removed (short notes are added below)
while r'\new{\textbf{' in T:
    i = T.index(r'\new{\textbf{')
    T = T[:i] + T[brace_end(T, i + 4):]
# C4: drop the red terms
T, n = re.subn(r'& \\rc -\\frac\{ig\^\{4\}f\^\{abd\}\}\{8\\pi\^\{5\}\}.*?& \\rc .*?\\notag \\\\ ', '', T, count=1, flags=re.S)
assert n == 1
T, n = re.subn(r'\\notag \\\\ & \\rc -\\frac\{ig\^\{4\}f\^\{abd\}\}\{16\\pi\^\{5\}\}.*?\\color\{black\}', '', T, count=1,
               flags=re.S)
assert n == 1 and r'\rc' not in T
T = drop_cmd(T, r'\old', False)
# Group I Row II part 1: only V_vee is left
i = T.index(r'\Bigg[\left[-\frac{2}{\epsilon}')
j = T.index(r'+\log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}', i)
T = T[:i] + r'\Bigg[\new{\log^{2}' + T[j + len(r'+\log^{2}'):]
T = once(T, r'\log{\left[1-\frac{k^{+}}{\vee}\right]}\,\log{\left[\frac{\vee}{k^{+}}\right]}\Bigg]\times',
         r'\log{\left[1-\frac{k^{+}}{\vee}\right]}\,\log{\left[\frac{\vee}{k^{+}}\right]}}\Bigg]\times')
T = once(T, r'''\subsection{Row III}
\begin{align}
   \frac{d^{3}N_{1\,\rm NLO}^{\rm}}{d^{3}k}  & = \frac{1}{(2\pi)^{3}}  \int_{w, w^\prime} e^{-ik(w^\prime-w)}''',
         r'''\new{Part 1: only $V_\vee=-\frac{\alpha_s}{2\pi}N_c\log(1-\frac{k^+}{\vee})\log(\frac{\vee}{k^+}-1)\,{\rm LO}$ (with c.c.) is left; it vanishes for $k^+\ll\vee$. The $\frac{11}{3}$ cancels against Row V 2b UV, the constant is in Eq.~\eqref{eq:cmw}, the pole bracket times $\log(1-\frac{k^+}{\vee})$ is Q1, Eq.~\eqref{eq:q1}. Part 2b: no change.}
\subsection{Row III}
\begin{align}
   \frac{d^{3}N_{1\,\rm NLO}^{\rm}}{d^{3}k}  & = \frac{1}{(2\pi)^{3}}  \int_{w, w^\prime} e^{-ik(w^\prime-w)}''')
T = once(T, r'''\subsection{Row IV}
\begin{align}
    \frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & =  \frac{1}{(2\pi)^{3}}''', r'''\new{The edge-kernel terms are removed (C4). Their place is taken by the finite part of the $\mathcal D$ region (Q3, open).}
\subsection{Row IV}
\begin{align}
    \frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & =  \frac{1}{(2\pi)^{3}}''')
# Group I Row V 2b UV: removed
i = T.index(r'\frac{d^{3}N_{2b\,\rm UV\, NLO}^{\rm}}{d^{3}k}')
i = T.rindex(r'\begin{align}', 0, i)
j = T.index(r'\end{align}', i) + len(r'\end{align}')
T = T[:i] + r'''\new{$d^3N_{2b\,\rm UV}$: nothing is left. Its $-\frac{11}{6}$ cancels the $\frac{11}{3}$ of Row II part 1 (Sec.~\ref{sec:rc}); $2\log(1-\frac{\Lambda}{k^+})\to0$.}''' + T[j:]
# Group II Row IV: only line 3
i = T.index(r'\frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & = \frac{g^{2}}{4\pi^{2}}\, \frac{N_{c}}{(2\pi)^{3}}')
j = T.index(r'& +\,\frac{g^{2}}{4\pi^{2}}\, \frac{1}{(2\pi)^{3}}\,\frac{g^{2}}{2\pi^{2}}\, \frac{N_{c}}{k^{+}}\,\int_{\frac{k^{+}}{\vee}}', i)
T = T[:i] + r'\frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & = ' + T[j + len(r'& +\,'):]
T = once(T, r"\int_{x,y,x',y'}\,e^{-ik\cdot (y'-y)\frac{1}{\xi}}\,\Bigg[\frac{(x-y)\cdot (x'-y')}",
         r"\int_{x,y,x',y'}\,e^{-ik\cdot (y'-y)\frac{1}{\xi}}\notag \\ & \times\Bigg[\frac{(x-y)\cdot (x'-y')}")
i = T.index(r'\notag \\ & \nc -\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}')
j = T.index(r'\end{align}', i)
T = T[:i].rstrip() + '\n' + T[j:]
T = once(T, r'''\section{Group III\, (with c.c.)}''', r'''\new{Only line 3 (linear polarization) is left. Lines 1, 2 and $-\log\frac{\vee}{k^+}H_4$ (C1) are the DGLAP term, Eq.~\eqref{eq:dglap}.}
\section{Group III\, (with c.c.)}''')
# Group III notes (each note goes after its row)
NOTE3 = (r'\new{C2: $\left[\frac{1}{p^+}\right]_+^{(k^+)}$; C3: the last line.@@ In the variables of AABKL: '
         r'Sec.~\ref{sec:aabkl}.}')
for anchor, extra in ((r'\subsubsection{Second part}' + '\n' + r'\begin{align}' + '\n' + r' \frac{d^{3}N_{2\,\rm NLO}}', ''),
                      (r'\subsection{Row IV}' + '\n' + r'\textsc{NA}' + '\n' + r'\subsection{Row V}', ' C5: the measure.'),
                      (r'\subsection{Row VI}', '')):
    assert T.count(anchor) == 1, anchor
    T = T.replace(anchor, NOTE3.replace('@@', extra) + '\n' + anchor)
T = T.rstrip() + '\n' + NOTE3.replace('@@', '') + '\n'
T = T.replace(r'\section{Group', r'\section{Finite terms: Group')
S = S.replace('@@FINITE@@', T)
assert '@@' not in S

if __name__ == '__main__':
    open(os.path.join(HERE, 'Finite_terms_final.tex'), 'w').write(S)
    print('written Finite_terms_final.tex')
