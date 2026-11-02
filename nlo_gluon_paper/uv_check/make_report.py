r"""Writes UV_divergences.tex from the scan outputs (uvscan2_*.json, uvtarget_*.json, surfscan.json)."""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(pat):
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, pat))):
        d.update(json.load(open(f)))
    return d


SCAN = load('uvscan2_*.json')
TARG = load('uvtarget_*.json')
SURF = json.load(open(os.path.join(HERE, 'surfscan.json')))
PN = {'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w', 'xp': "x'", 'wp': "w'", 'yp': "y'"}


def pair_tex(p):
    a, b = p.split('-')
    return '$%s\\to %s$' % (PN[a], PN[b])


def short(key):
    k = key.replace('(with c.c.) ', '').replace('  ', ' ').replace(' / ', ', ')
    k = k.replace('ROW TOTAL', '\\textbf{row total}')
    k = k.replace(' #', ', eq.~').replace('Second part', '2nd part').replace('First part', '1st part')
    return k


def worst(res):
    """(min exponent, pair) over all pairs and all seeds / p+ values"""
    best = (99.0, '')
    for v in res.values():
        if not isinstance(v, dict):
            continue
        for p, x in v.items():
            if isinstance(x, dict) and 's' in x and x['s'] < best[0]:
                best = (x['s'], p)
    return best


def verdict(s):
    return '\\textcolor{red}{divergent}' if s <= -1.7 else 'finite'


rows = []
for key, res in SCAN.items():
    if 'skipped' in res:
        continue
    s, p = worst(res)
    divs = sorted({pp for v in res.values() for pp, x in v.items() if x['s'] <= -1.7})
    rows.append((key, '%.2f' % s, pair_tex(p) if p else '', ', '.join(pair_tex(q) for q in divs) or '---',
                 verdict(s), 'all pairs, 2 configurations'))
INT = json.load(open(os.path.join(HERE, 'intcheck.json'))) if os.path.exists(os.path.join(HERE, 'intcheck.json')) else {}
for key, res in TARG.items():
    if 'skipped' in res:
        continue
    if 'Row V / Second part #3' in key and INT:
        ss = min(v['sum']['s'] for v in INT.values())
        ls = [min(v['lines'][i]['s'] for v in INT.values()) for i in range(len(next(iter(INT.values()))['lines']))]
        rows.append((key, '%.2f' % ss, pair_tex('x-z'), '---', verdict(ss),
                     'p$^+$, q$^+$ integrals done; lines %s diverge at $x\\to z$' % ', '.join(str(i + 1) for i, l in enumerate(ls) if l <= -1.7)))
        continue
    if 'terms' in res:      # lines with different longitudinal variables: judge line by line
        src = res['terms']
        s, p = min((worst({'a': t}) for t in src), key=lambda z: z[0])
        divs = sorted({pp for t in src for pp, x in t.items() if x['s'] <= -1.7})
    else:
        s, p = worst({'a': res['sum']})
        divs = sorted({pp for pp, x in res['sum'].items() if x['s'] <= -1.7})
    tdiv = []
    for i, t in enumerate(res.get('terms', [])):
        d = [pp for pp, x in t.items() if x['s'] <= -1.7]
        if d:
            tdiv.append('line %d: %s' % (i + 1, ', '.join(pair_tex(q) for q in d)))
    if 'Group II' in key:   # no longitudinal integral: the sum is meaningful
        s, p = worst({'a': res['sum']})
        divs = sorted({pp for pp, x in res['sum'].items() if x['s'] <= -1.7})
    rows.append((key, '%.2f' % s, pair_tex(p), ', '.join(pair_tex(q) for q in divs) or '---', verdict(s),
                 ('selected pairs; ' + '; '.join(tdiv)) if tdiv else 'selected pairs, $p^+=0.37k^+$, line by line'))


def order(r):
    k = r[0]
    g = 1 if 'Group I (' in k else (2 if 'Group II /' in k else 3)
    return (g, k.replace('ROW TOTAL', 'zz'))


rows.sort(key=order)
merged = {}
for r in rows:
    if r[0] in merged:      # scanned twice: keep the all-pairs verdict, add the line-by-line information
        old = merged[r[0]]
        extra = r[5].replace('selected pairs; ', '')
        merged[r[0]] = old[:5] + (old[5] + '; ' + extra,)
    else:
        merged[r[0]] = r
table = []
for r in sorted(merged.values(), key=order):
    table.append('%s & %s & %s & %s & %s & {\\footnotesize %s}\\\\' % (short(r[0]), r[1], r[2], r[3], r[4], r[5]))

srows = []
NAMEC = {'D1': '$(x-z)^2=(x-w)^2$', 'D2': "$(x'-w')^2(x-z)^2=(x-w)^2(x'-z)^2$", 'D3': "$(x'-w')^2=(x'-z)^2$",
         'S': '$w$ on the segment $[x,z]$ ($\\mathcal A_0$)'}
for k, v in SURF.items():
    key, sname = k.split(' | ')
    cls = sorted({r['sum'][2] for r in v.values()})
    tcl = sorted({t[2] for r in v.values() for t in r['terms'].values()})
    ex = '; '.join('$e_E=%.2f$, $e_O=%.2f$' % (r['sum'][0], r['sum'][1]) for r in list(v.values())[:1])
    col = (lambda c: '\\textcolor{red}{%s}' % c if c == 'DIVERGENT' else c)
    srows.append('%s & %s & %s & %s & %s\\\\' % (short(key), NAMEC[sname], ', '.join(col(c) for c in tcl), ex,
                                                 ', '.join(col(c) for c in cls)))

T = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.5in]{geometry}
\usepackage{amsmath,amssymb,amsfonts,mathtools}
\usepackage[dvipsnames]{xcolor}
\usepackage{array,booktabs,longtable,float}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\setlength{\emergencystretch}{3em}
\newcommand{\KK}{\mathcal K}
\title{Transverse divergences of the terms in your \texttt{Other\_terms.tex}}
\date{}
\begin{document}
\maketitle
\noindent Your file without Group I Row II part 1 and Group II Row IV (both inside \verb|\begin{comment}|). Question: does any term diverge when transverse coordinates come together (UV), or anywhere else in the transverse plane? $\KK(X)=X/X^2$, $a=\frac{k^+}{\vee}$, $U^{bc}=U^{cb}$.

\section{Answer}
\begin{center}\fbox{\begin{minipage}{0.96\textwidth}
\textbf{Flagged.}
\begin{enumerate}
\item[\textbf{F1}] \textbf{Group I Row V, $d^3N_{2b\,\rm UV}$: UV divergent.} The $\frac{2}{\epsilon}$ is the dimensionally regularized UV divergence of 2b at $x\to z$ (the two gluons of the loop at one point). After it, $\log\frac{(y-z)^2\mu^2}{4}$ times $\frac{(y-z)^i}{(y-z)^2}[U(y)-U(z)]$ is integrable at $y\to z$. In this file its partner is gone: the $\frac{11}{3}$ term of Group I Row II part 1, which cancels its $\beta_0$ pole, is commented out.
\item[\textbf{F2}] \textbf{Group I Row III: each of the parts 2, 3a, 3b is divergent; only the sum of the four parts is finite.}
\begin{itemize}
\item At $z\to w$, each of parts 2, 3a, 3b grows like $r^{-2}$ (logarithmic divergence).
\item Across the curve $(x-z)^2=(x-w)^2$ (it passes through $z=w$), the separate lines go like $[(x-z)^2-(x-w)^2]^{-2}$, which is not integrable.
\item Across the segment where $w$ lies between $x$ and $z$, the $\mathcal A_0$ lines go like $1/{\rm distance}$, which is not integrable either.
\item The sum of parts 1, 2, 3a, 3b is regular on both curves, and its worst UV behaviour is $r^{-1}$ (at $x\to w$), which is integrable.
\end{itemize}
Row III must be kept as one expression.
\item[\textbf{F3}] \textbf{Group I Row V, $d^3N_{2b\,\rm finite}$: separate lines divergent at $x\to z$, the sum finite.} The $d_\perp$ line, the line without the plus prescription, and the $\frac{11}{6}$ line each go like $(x-z)^{-2}$. The $\frac{11}{6}$ line is the subtraction that removes this divergence.
\item[\textbf{F4}] \textbf{Group II Row I, $d^3N_{2f}$: each of its three lines is divergent, the sum finite.} Line 2, $\big[\frac1\vee-\frac1{k^+}\big]\frac{1}{(z-w)^2(z-w')^2}$, diverges logarithmically at $z\to w$ and at $z\to w'$. Line 1 diverges at $z\to w$ and line 3 at $z\to w'$, and they cancel it. Keep the three lines together.
\item[\textbf{F5}] \textbf{Group I Row II, $d^3N_{2b}$ (momentum form, $R^i$): UV finite, but infrared divergent} (not UV). Two regions are affected.
\begin{itemize}
\item At $p\to0$: $R^i\to-2k^i\log\frac{p^2}{k^2}$.
\item At $q=k-p\to0$: $R^i\to\frac{2k^2k^i}{q^2}\log\frac{k^2}{q^2}$.
\end{itemize}
With the $\frac{1}{p^2k^2}$, $\int d^2p$ diverges like $\log^2$ at both points. These are the regions $|x-y|\to\infty$. In addition, $\int d^2k$ uses the same letter $k$ as the measured gluon.
\end{enumerate}
\textbf{Not divergent} (UV and on curves): every other equation, i.e.\ Group I Rows III (sum) and IV, Row V first part and $d^3N_{2a}$, Group II Row I (every equation; separate lines of the second part need a principal value across the curve $(x'-w')^2(x-z)^2=(x-w)^2(x'-z)^2$), Group III Rows I, II, III, V, VI.
\end{minipage}}\end{center}

\section{Method}
\begin{itemize}
\item Your equations are evaluated pointwise as written (each equation split into its lines, \texttt{terms.py}). Wilson lines and charges are smooth functions of position:
\begin{itemize}
\item $U(p)$ is the adjoint of $V(p)=W(p)\,{\rm diag}(1,-1,-1)\,W(p)^\dagger$, with $W(p)=e^{i\theta_a(p)t^a}$; this is an SU(3) matrix and $U=U^T$.
\item $\theta_a(p)$ and $\rho^a(p)$ are sums of sines with wave numbers of order 1.
\end{itemize}
A colour factor that vanishes when two points coincide therefore vanishes linearly, as it does in the real integral.
\item \textbf{UV.} For every pair of points $P,Q$, set $P=Q+r\,\hat n(\varphi)$ and average the integrand over $\varphi$, giving $A(r)$ at $r=10^{-1.5}\dots10^{-4.5}$. The exponent $s$ in $A\sim r^s$ decides: $\int d^2r$ diverges iff $s\le-2$.
\begin{itemize}
\item Where a singular curve crosses the circle, the $\varphi$ average is a principal value.
\item Every equation was scanned, and every row with several equations as a whole.
\item Two field and point configurations were used. Terms with $p^+$ or $q^+$ integrals are evaluated at fixed $p^+$ ($0.37k^+$; the plus prescriptions subtract at $p^+=0$), on the pairs of points that can diverge (Sec.~\ref{sec:pc}).
\item In Group I Row V $d^3N_{2b\,\rm finite}$ the lines have different longitudinal variables ($p^+$, $q^+$, none), so they can cancel only after the integrals. There the $p^+$ and $q^+$ integrals are done (16-point Gauss--Legendre on $[0,k^+]$, with the orientation and $\frac{1}{2\pi}$ of each), \texttt{intcheck.py}.
\end{itemize}
\item \textbf{Singular curves.} A point is moved across the curve; from $g(\pm\delta)$ take $E=\frac{g_++g_-}{2}\sim\delta^{e_E}$ and $O=\frac{g_+-g_-}{2}\sim\delta^{e_O}$.
\begin{itemize}
\item Integrable: $e_E>-1$ and $e_O>-1$.
\item Principal value: $e_O=-1$.
\item Divergent: $e_E\le-1$ or $e_O\le-2$.
\end{itemize}
\end{itemize}

\section{UV: coincident points}
$s$ = exponent of $A(r)$ at the worst pair (all configurations and $p^+$ values); $s\le-2$ = divergent.
{\small
\begin{longtable}{p{5.2cm}cp{1.6cm}p{2.4cm}p{1.4cm}p{4.2cm}}
\toprule
equation & worst $s$ & pair & divergent pairs & verdict & scan\\
\midrule
""" + '\n'.join(table) + r"""
\bottomrule
\end{longtable}}
Exponents near $0$ or $-1$ are bounded or $1/r$ singularities, both integrable. Group I Row II $d^3N_{2b}$ (momentum integrals) and Row V $d^3N_{2b\,\rm UV}$ (explicit pole) are not scanned; see F5 and F1.

\section{Singular curves}
{\small
\begin{longtable}{p{4.3cm}p{3.6cm}p{3.4cm}p{3cm}p{2.2cm}}
\toprule
equation & curve & separate lines & sum: $e_E$, $e_O$ & sum\\
\midrule
""" + '\n'.join(srows) + r"""
\bottomrule
\end{longtable}}

\section{Analytic checks}\label{sec:pc}
\paragraph{Power counting.} In two dimensions $\int d^2r\,r^{-n}$ converges at $r\to0$ for $n<2$. When $m$ points collapse together the measure is $r^{2m-3}dr$. Each $\KK(u-v)$ contributes $r^{-1}$ when $u\to v$. A colour bracket like $U(z)-U(y)$ adds a factor $r$ when $y\to z$.
\begin{itemize}
\item \textbf{Group III Rows III, V, VI, line 1} ($p^+$ integral) has $\KK^m(y'-z)\,\KK^{k'}(x'-y')\,\KK^i(y-w)\,\KK^j(x-z)$. Each kernel joins a different pair, so every two-point limit is $r^{-1}$. A collapse of $m$ points brings at most $m-1$ kernels, i.e.\ $r^{-(m-1)}$ against $r^{2m-3}$. The same holds for line 2 ($\KK(x'-z)$ in place of $\KK(x'-y')$) and line 3. Hence these rows are finite.
\item \textbf{Group III Row III, second part.} $p^+(x-z)^2+k^+(x-w)^2$ vanishes only if $x,z,w$ all coincide. Its $\frac{1}{(z-w)^2}$ always multiplies $(z-w)$ or $(x-z)^2-(x-w)^2$, which is $O(r)$ at $z\to w$.
\item \textbf{Group I Row V, $d^3N_{2a}$.} At $x\to z$ the $d_\perp$ line goes like $\frac{(z-x)^i}{(x-z)^2}\frac{(y-x)^2-(y-z)^2}{(x-z)^2}\sim r^{-2}$, with a non-zero angular average. But its colour factor vanishes there: $f^{adc}U^{db}(z)U^{ce}(z)=f^{hbe}U^{ah}(z)$ and $f^{fbe}f^{hbe}=N_c\delta^{fh}$. So
\begin{equation*}
\big[-f^{adc}f^{fbe}U^{db}(x)U^{ce}(z)+N_cU^{af}(x)\big]_{x=z}=-N_cU^{af}(z)+N_cU^{af}(z)=0,
\end{equation*}
and the integrand is $O(r^{-1})$. The same holds for the $q^+$ line.
\item \textbf{Group I Row V, $d^3N_{2b\,\rm finite}$.} Its colour factor $N_c[U(y)-U(x)]$ does not vanish at $x\to z$.
\begin{itemize}
\item The $d_\perp$ line and the $q^+$ line without the plus prescription go like $(x-z)^{-2}$, with a non-zero average.
\item In the $[\frac{k^+}{q^+}]_+$ line the leading term does not depend on $q^+$ at $x=z$, so the plus prescription removes it.
\item The $\frac{11}{6}$ line, with $\frac{(y-x)\cdot(y-z)}{(y-x)^2(x-z)^2}\to\frac{1}{(x-z)^2}$, cancels the rest. This is the subtraction of \texttt{rowV/Part2b\_MSbar.pdf}: the sum is finite for every $p^+$.
\end{itemize}
\end{itemize}

\paragraph{Group I Row III, the $\mathcal A_0$ segment.} $\mathcal A_0=\int_0^1\frac{dt}{|(x-w)-t(x-z)|^2}$, which is $\simeq\frac{\pi}{|(x-w)\times(x-z)|}$ when $w$ lies on the segment from $x$ to $z$: a $1/{\rm distance}$ singularity across a line, not integrable. On the segment $x-w=t\,e$, $z-w=-(1-t)e$, with $e=x-z$. Let $Q\equiv\frac{[(x'-w')\cdot e][(y-z)\cdot e]}{e^2}$, and drop the common $\frac{1}{(x'-w')^2(y-z)^2}$. The four $\mathcal A_0$ lines (3 in part 3a, 1 in part 3b) become
\begin{equation*}
\underbrace{-\frac{t}{1-t}Q-\frac{1}{1-t}Q+\frac{t}{1-t}Q}_{\rm 3a}\;\underbrace{+\frac{1}{1-t}Q}_{\rm 3b}=0.
\end{equation*}
So their sum vanishes linearly on the segment and $\mathcal A_0\times$(sum) is bounded. Part 3a or part 3b alone is not: they must stay together.

\paragraph{Group I Row II, $R^i$.} Limits, checked numerically to $10^{-6}$:
\begin{itemize}
\item $p\to\infty$: $R^i\to k^i[4\log a-\log(1-a)]$, giving $\int\frac{d^2p}{p^2}e^{ip\cdot(y-x)}\sim-2\pi\log|x-y|$ (integrable).
\item $k\to\infty$: $R^i\to2k^i\log\frac{k^2}{p^2}$, giving $\frac{(w-y)^i}{(w-y)^2}\log|w-y|$ (integrable).
\item $x,y\to w$ together: the transform is homogeneous of degree $-1$, against $r^3dr$ (integrable).
\item $p\to0$ and $p\to k$: the $\log^2$ infrared divergences of F5.
\end{itemize}

\section*{Reproducing}
In \texttt{nlo\_gluon\_paper/uv\_check/}: \texttt{uvscan2.py} (all pairs, equations and rows), \texttt{uvtarget.py} (rows with $p^+$ integrals), \texttt{intcheck.py} (Row V $d^3N_{2b\,\rm finite}$ with the integrals done), \texttt{surfscan.py} (singular curves), \texttt{make\_report.py} (this file).
\end{document}
"""
open(os.path.join(HERE, 'UV_divergences.tex'), 'w').write(T)
print('written', len(table), 'rows,', len(srows), 'curve rows')
