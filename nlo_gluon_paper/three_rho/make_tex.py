"""Generate the LaTeX note 'Two-rho terms from the three-rho terms'."""
import pickle
import sympy as sp
from collections import OrderedDict
from engine import tex_colour, tex_coeff, Nc
from kernels import subst_kernel, tex_kernel, subst_label
import process

RULE_TEX = {'UU': r'\text{(U1)}', 'UU(sym)': r'\text{(U1$^\star$)}', 'fUU': r'\text{(U2)}',
            'fUU(sym)': r'\text{(U2$^\star$)}', 'ff': r'\text{(F)}', 'delta': None, 'trd': r'\delta^{aa}=N_c^2-1'}

PHASE_SWAP = {"w'": 'w'}

REMARKS = {
    'A1': r"""In the notes this term reads $-\frac{i}{(2\pi)^3}\frac{g^4}{8\pi^3}\frac{\log(\vee/\Lambda)}{2\pi^2}\frac1{k^+}\int(\dots)$, i.e.\ its prefactor is $-i\NN$. It is the complete remainder of reordering the four-$\rho$ terms, so no c.c.\ is added.""",
    'A2': r"""Same prefactor $-i\NN$ and no c.c.""",
    'B1': r"""In the notes the three-Wilson-line remainder has prefactor $+\frac{i}{(2\pi)^3}\frac{g^4}{8\pi^3}\frac{\log(\vee/\Lambda)}{2\pi^2}\frac1{k^+}=+i\NN$. Its first and fourth terms share the kernel $K_1$ and are treated here; no c.c.\ is added.""",
    'B2': r"""The second and third terms of the three-Wilson-line remainder share the kernel $K_3$; prefactor $+i\NN$, no c.c.""",
    'GI.II': r"""\emph{Input.} In the notes,
\begin{align*}
\frac{d^3N_2}{d^3k}={}&-\frac{1}{(2\pi)^3}\frac{g^4f^{abc}}{64\pi^6}\frac1{k^+}\int_{w,w'}e^{-ik\cdot(w'-w)}\int d^2\mathbf p\,d^2\mathbf k\int_{x,y}e^{i\mathbf k\cdot(w-y)}e^{i\mathbf p\cdot(y-x)}\frac{(x'-w')^i}{(x'-w')^2\,\mathbf p^2\mathbf k^2}\,B^i(\mathbf k,\mathbf p)\,\log\frac{\vee}{\Lambda}\\
&\times\Big[U^{ab'}(x')\rho^{b'}(x')-U^{ab'}(w')\rho^{b'}(x')\Big]\Big[U^{ab}(w)\{\rho^b(x),\rho^c(y)\}-U^{bd}(x)U^{ce}(y)\{\rho^d(x),\rho^e(y)\}\Big]+{\rm c.c.},\\
B^i(\mathbf k,\mathbf p)\equiv{}&-\frac{\mathbf k^2\mathbf p^i-(\mathbf k\cdot\mathbf p)\mathbf k^i}{(\mathbf k-\mathbf p)^2}-2\mathbf k^i-\frac{\mathbf p^2\mathbf k^i+\mathbf k^2\mathbf p^i}{(\mathbf k-\mathbf p)^2}.
\end{align*}
Two remarks before reordering.
\begin{enumerate}
\item \emph{Colour.} As written, the indices $a$ and $b$ each appear three times in the first term of the second bracket: in $f^{abc}$, $U^{ab'}(x')$, $U^{ab}(w)$ and $\rho^b(x)$. This term is $U^{ab}(w)\mathbb A^{(3)b}$ with $\mathbb A^{(3)b}\propto f^{bcd}\{\rho^c(x),\rho^d(y)\}$, so we read it as $U^{ah}(w)f^{hbc}\{\rho^b(x),\rho^c(y)\}$. The second term, $f^{abc}U^{bd}(x)U^{ce}(y)\{\rho^d(x),\rho^e(y)\}=\bar{\mathbb A}^{(3)a}$, is unchanged.
\item \emph{Kernel.} The Fourier transform is done in App.~\ref{app:F}. With $\int d^2\mathbf p\,d^2\mathbf k\,e^{i\mathbf k\cdot(w-y)}e^{i\mathbf p\cdot(y-x)}B^i/(\mathbf p^2\mathbf k^2)=2\pi i\,\mathfrak F^i(w;x,y)$, and since both colour structures are odd under $x\leftrightarrow y$, only the part of $\mathfrak F^i$ odd under $x\leftrightarrow y$ contributes. It is
\begin{equation*}
\mathfrak F^i(w;x,y)=\KK^i(x-w)\int_{z'}\big[K(w,y,z')-K(x,y,z')\big]-\KK^i(y-w)\int_{z'}\big[K(w,x,z')-K(x,y,z')\big],
\end{equation*}
with $K(a,b,z')=\KK(a-z')\cdot\KK(b-z')$; equivalently $\mathfrak F^i(w;x,y)=-\pi\big[\KK^i(x-w)\log\frac{(w-y)^2}{(x-y)^2}-\KK^i(y-w)\log\frac{(w-x)^2}{(x-y)^2}\big]$. The prefactor becomes $-\frac{g^4}{64\pi^6}\cdot2\pi i=-\frac{i}{2}\cdot\frac{g^4}{16\pi^5}$, i.e.\ $-\frac i2\NN$.
\end{enumerate}""",
    'GI.III': r"""\emph{Input.} Part 1 of the notes has three terms with the same colour structure. The first two (red) are equal and opposite and cancel. The third has prefactor $\frac{ig^4f^{abd}}{16\pi^5}=i\NN f^{abd}$ and kernel $K^{(1)}=\KK(x'-w')\cdot\KK(x-w)\,\KK(x-z)\cdot\KK(y-z)$. Part 2 has five terms, again with the same colour structure and prefactor $\pm\frac{ig^4f^{abd}}{8\pi^5}=\pm2i\NN f^{abd}$; its terms 1 and 5 are equal, and so are terms 2 and 4. Together,
\begin{gather*}
\text{Row III}=i\NN\,f^{abd}\int K_{\rm I.III}(x',x,y)\,[\dots][\dots]+{\rm c.c.},\qquad K_{\rm I.III}=K^{(1)}+2\,\KK(x'-w')\cdot\KK(y-z)\,\Phi(x),\\
\Phi(x)\equiv\frac{1}{(z-w)^2}\bigg[\frac{2\,(z-w)\cdot(x-z)\,(x-w)^2-2\,(z-w)\cdot(x-w)\,(x-z)^2}{\big[(x-w)^2-(x-z)^2\big]^2}-1\bigg].
\end{gather*}""",
    'GI.IV': r"""\emph{Input.} Prefactor $\frac{1}{(2\pi)^3}\frac{ig^4}{8\pi^5}\frac1{k^+}\log\frac\vee\Lambda=2i\NN$. The kernel is the sum of the two kernels of the notes, $K_{\mathfrak B}\equiv\KK(x-w)\cdot\KK(x'-w')\,\big[\KK(z-w)\cdot\KK(y-z)+\tfrac12\KK(x-z)\cdot\KK(y-z)\big]$.""",
    'GI.V1': r"""\emph{Input.} Prefactor $-\frac{1}{(2\pi)^3}\frac{ig^4}{32\pi^5}\frac1{k^+}\log\frac\vee\Lambda=-\frac i2\NN$. After the $w$ integration with the $\delta$-function of $\mathbb C_1$ the measured gluon in the amplitude sits at $z$: the phase is $e^{-ik\cdot(w'-z)}$. The charges sit at $x'$, $y$ and $y'$.""",
    'GI.V2': r"""\emph{Input.} Prefactor $-\frac i2\NN$; the measured gluon in the amplitude sits at $x$ (phase $e^{-ik\cdot(w'-x)}$). Same colour structure as the first term.""",
    'GII.I1': r"""\emph{Input.} The prefactor is $-\frac{1}{(2\pi)^3}\frac{ig^4}{16\pi^5}\frac1{k^+}\log\frac\vee\Lambda=-i\NN$ and the kernel is $K_{\mathfrak B}$, as in Row IV of Group I. The charges appear as $\{\rho(x'),\rho(y)\}\rho(x)$. Group II has no c.c.; the second term of the row is its conjugate.""",
    'GII.I2': r"""\emph{Input.} Prefactor $+i\NN$ and kernel
\begin{equation*}
\bar K_{\mathfrak B}=K_{\mathfrak B}\big|_{w\leftrightarrow w'}=\KK(x-w')\cdot\KK(x'-w)\,\big[\KK(z-w')\cdot\KK(y-z)+\tfrac12\KK(x-z)\cdot\KK(y-z)\big].
\end{equation*}
The charges appear as $\rho(x)\{\rho(x'),\rho(y)\}$.""",
    'GIII.I': r"""\emph{Input.} Prefactor $-\frac{1}{(2\pi)^3}\frac{ig^4}{8\pi^5}\frac1{k^+}\log\frac\vee\Lambda=-2i\NN$, kernel $K_{\mathfrak B}$. Order of the charges: $\rho(x')\rho(y)\rho(x)$.""",
    'GIII.II': r"""\emph{Input.} Prefactor $-2i\NN$, kernel $K_{\mathfrak B}$. Order of the charges: $\rho(y)\rho(x')\rho(x)$.""",
    'GIII.III': r"""\emph{Input.} Prefactor $-\frac{1}{(2\pi)^3}\frac{ig^4f^{c'd'a}}{16\pi^5}\frac1{k^+}\log\frac\vee\Lambda=-i\NN f^{c'd'a}$. After the $w'$ integration the measured gluon of the conjugate amplitude sits at $y'$: the phase is $e^{-ik\cdot(y'-w)}$, and the kernel is $K_{\mathfrak C}\equiv\KK(y'-z)\cdot\KK(x-z)\,\KK(x'-y')\cdot\KK(y-w)$.""",
    'GIII.V': r"""\emph{Input.} Prefactor $-\frac{1}{(2\pi)^3}\frac{ig^4f^{c'ca}}{8\pi^5}\frac1{k^+}\log\frac\vee\Lambda=-2i\NN f^{c'ca}$, phase $e^{-ik\cdot(y'-w)}$, kernel $K_{\mathfrak C}$. Order of the charges: $\rho(x')\rho(x)\rho(y)$.""",
    'GIII.VI': r"""\emph{Input.} Prefactor $-2i\NN f^{c'd'a}$, phase $e^{-ik\cdot(y'-w)}$, kernel $K_{\mathfrak C}$. Order of the charges: $\rho(x')\rho(y)\rho(x)$.""",
}

SECTION_OF = OrderedDict([
    ('fourrho', ('Three-$\\rho$ terms from reordering the four-$\\rho$ terms', ['A1', 'A2', 'B1', 'B2'])),
    ('G1', ('Group I (with c.c.)', ['GI.II', 'GI.III', 'GI.IV', 'GI.V1', 'GI.V2'])),
    ('G2', ('Group II (no c.c.)', ['GII.I1', 'GII.I2'])),
    ('G3', ('Group III (with c.c.)', ['GIII.I', 'GIII.II', 'GIII.III', 'GIII.V', 'GIII.VI'])),
])


def tex_word(word, form):
    r = [r'\rho^{%s}(%s)' % (i, p) for i, p in word]
    if form == 'plain':
        return r'\,'.join(r)
    if form == 'A{BC}':
        return r[0] + r'\,\big\{' + r[1] + ',' + r[2] + r'\big\}'
    return r'\big\{' + r[0] + ',' + r[1] + r'\big\}\,' + r[2]


def tex_pref(c):
    """prefactor such as -i, i/2, 2i as LaTeX (with sign)."""
    s, b = tex_coeff(c)
    b = b.replace(r'\,', '')
    if b == '':
        b = '1'
    return ('-' if s == '-' else '') + b


def piece_rule_tex(rules):
    out = [RULE_TEX[r] for r in rules if RULE_TEX.get(r)]
    return r'\quad' + r',\ '.join(out) if out else ''


def kargs_after(src, subst):
    old, new = subst
    return [subst_label(l, old, new) for l in src.kargs]


def phase_tex(src):
    P1, P0 = src.phase
    return r'e^{-ik\cdot(%s-%s)}' % (P1, P0)


def cc_tex(src):
    if not src.cc:
        return ''
    P1, P0 = src.phase
    return r'+\big(%s\leftrightarrow %s\big)' % (P1, P0)


def grouped_result(src, sres):
    groups = OrderedDict()
    for sub in sres['subs']:
        for pc in sub['pieces']:
            if pc['coeff_after'] == 0:
                continue
            args = tuple(kargs_after(src, pc['subst']))
            col = tex_colour(pc['facs_after'], rho_mode='anti')
            key = (args, col)
            groups[key] = groups.get(key, 0) + sp.expand(src.pref * pc['coeff_after'])
    return [(k, sp.simplify(v)) for k, v in groups.items() if sp.simplify(v) != 0]


def coeff_line(c, first):
    s, b = tex_coeff(c)
    if first and s == '+':
        s = ''
    return s + b


def gen_source(src, sres):
    L = []
    L.append(r'\subsection{%s}\label{sec:%s}' % (src.title, src.key))
    L.append(REMARKS.get(src.key, ''))
    # input
    L.append(r'\paragraph{The three-$\rho$ term.}')
    L.append(r'\begin{align}')
    head = r'\frac{d^3N}{d^3k}\bigg|_{\rm %s}={}&%s\,\NN\int %s\,%s(%s)\,\Big\{' % (
        src.key, tex_pref(src.pref), phase_tex(src), src.kname, ','.join(src.kargs))
    L.append(head + r'\notag\\')
    for n, (sign, colour, word) in enumerate(src.subterms):
        sg = '+' if sign > 0 else '-'
        line = r'&\quad %s\,%s\,%s' % (sg, tex_colour(colour, keep_order=True), tex_word(word, src.form))
        if n == len(src.subterms) - 1:
            line += r'\Big\}' + (r'+{\rm c.c.}' if src.cc else '')
            L.append(line + r'\label{eq:in-%s}' % src.key)
        else:
            L.append(line + r'\notag\\')
    L.append(r'\end{align}')
    L.append(r'with $%s(%s)=%s$.' % (src.kname, ','.join(src.kargs), tex_kernel(src.kexpl)))
    # identity used
    form_ref = {'plain': r'Eq.~\eqref{eq:plain}', 'A{BC}': r'Eq.~\eqref{eq:Aanti}', '{BC}A': r'Eq.~\eqref{eq:antiA}'}[src.form]
    w0 = src.subterms[0][2]
    if src.form == '{BC}A':
        abc = (w0[2][1], w0[0][1], w0[1][1])
    else:
        abc = (w0[0][1], w0[1][1], w0[2][1])
    L.append(r'\paragraph{Reordering.} In every sub-term (one product of the two brackets) the charges appear in the order shown, so we apply %s with $A$, $B$, $C$ the charges at $%s$, $%s$, $%s$. The sub-terms are numbered 1--%d in the order above. For each sub-term we list the commutators, and each resulting two-$\rho$ piece before and after the colour algebra.' % (form_ref, abc[0], abc[1], abc[2], len(src.subterms)))
    for n, sub in enumerate(sres['subs']):
        sign = sub['sign']
        L.append(r'\paragraph{Sub-term %d:} $%s\,%s\,%s$.' % (
            n + 1, '+' if sign > 0 else '-', tex_colour(sub['colour'], keep_order=True), tex_word(sub['word'], src.form)))
        L.append(r'\begin{align}')
        rows = []
        for pc in sub['pieces']:
            ncol = len(sub['colour'])
            fnew = pc['facs_before'][ncol]
            hh = fnew[3]
            (xi, xp), (yi, yp), (zi, zp) = pc['X'], pc['Y'], pc['Z']
            comm = r'[\rho^{%s}(%s),\rho^{%s}(%s)]=i f^{%s%s%s}\rho^{%s}(%s)\,\delta^{(2)}(%s-%s)' % (
                xi, xp, yi, yp, xi, yi, hh, hh, xp, xp, yp)
            args = kargs_after(src, pc['subst'])
            kn = r'%s(%s)' % (src.kname, ','.join(args))
            s0, b0 = tex_coeff(pc['coeff_before'])
            before = r'%s%s\,%s\,%s' % ('' if s0 == '+' else '-', b0, kn, tex_colour(pc['facs_before'], rho_mode='anti', keep_order=True))
            if pc['coeff_after'] == 0:
                after = '0'
            else:
                s1, b1 = tex_coeff(pc['coeff_after'])
                after = r'%s%s\,%s\,%s' % ('' if s1 == '+' else '-', b1, kn, tex_colour(pc['facs_after'], rho_mode='anti'))
            tag = pc['tag']
            fac = {'plain': r'\tfrac14', 'A{BC}': r'\tfrac12', '{BC}A': r'-\tfrac12'}[src.form]
            X = r'%s' % tag[0]
            Y = r'%s' % tag[1]
            other = ({'AB': 'C', 'AC': 'B', 'BC': 'A'}[tag])
            rows.append(r'&[%s,%s]:\ %s,\quad %s\to %s\notag\\' % (X, Y, comm, yp, xp))
            rows.append(r'&\qquad %s\{[%s,%s],%s\}=%s\notag\\' % (fac, X, Y, other, before))
            rule = piece_rule_tex(pc['rules'])
            rows.append(r'&\qquad\phantom{%s\{[%s,%s],%s\}}=%s%s' % (fac, X, Y, other, after, rule))
            rows.append(r'\label{eq:%s-%d-%s}\\' % (src.key, n + 1, tag))
        rows[-1] = rows[-1].rstrip('\\')
        L += rows
        L.append(r'\end{align}')
    # result
    res = grouped_result(src, sres)
    L.append(r'\paragraph{Two-$\rho$ terms of this source.} Multiplying by the prefactor $%s\,\NN$ and collecting equal terms (dummy colour indices renamed),' % tex_pref(src.pref))
    L.append(r'\begin{align}')
    L.append(r'\frac{d^3N^{(2\rho)}}{d^3k}\bigg|_{\rm %s}={}&\NN\int %s\,\Big\{\notag\\' % (src.key, phase_tex(src)))
    for j, ((args, col), c) in enumerate(res):
        ln = r'&\quad %s%s(%s)\,%s' % (coeff_line(c, False), src.kname, ','.join(args), col)
        if j == len(res) - 1:
            L.append(ln + r'\Big\}' + cc_tex(src) + r'\label{eq:res-%s}' % src.key)
        else:
            L.append(ln + r'\notag\\')
    L.append(r'\end{align}')
    # kernels
    pats = OrderedDict()
    for (args, col), c in res:
        pats[args] = True
    L.append(r'The kernel at the coincident points:')
    L.append(r'\begin{align*}')
    lines = []
    for args in pats:
        # find substitution that produced it
        old = [l for l in src.kargs if l not in args]
        if old:
            old = old[0]
            new = [a for a, l in zip(args, src.kargs) if l == old][0]
            ke = subst_kernel(src.kexpl, old, new)
        else:
            ke = src.kexpl
        lines.append(r'%s(%s)&=%s' % (src.kname, ','.join(args), tex_kernel(ke)))
    L.append(r'\\'.join(lines))
    L.append(r'\end{align*}')
    return '\n'.join(L), res


def explicit_kernel(src, args):
    old = [l for l in src.kargs if l not in args]
    if not old:
        return tex_kernel(src.kexpl)
    old = old[0]
    new = [a for a, l in zip(args, src.kargs) if l == old][0]
    return tex_kernel(subst_kernel(src.kexpl, old, new))


def main():
    res = process.process()
    body = []
    allres = {}
    srcmap = {r['src'].key: r for r in res}
    for key, (title, keys) in SECTION_OF.items():
        body.append(r'\section{%s}\label{sec:%s}' % (title, key))
        for k in keys:
            txt, rr = gen_source(srcmap[k]['src'], srcmap[k])
            body.append(txt)
            allres[k] = rr
    # summary
    S = [r'\section{Summary: all two-$\rho$ terms}\label{sec:summary}']
    S.append(r"""Below are all two-$\rho$ terms generated by the reordering, with the kernels written out. Every coefficient is real. For Groups I and III the complex conjugate is included as the exchange of the two measured-gluon positions in the phase, shown at the end of each line. The four-$\rho$ remainders and Group II have no c.c.\ (Sec.~\ref{sec:conv}). The symmetric three-$\rho$ remainders are not listed (Sec.~\ref{sec:identity}).""")
    for key, (title, keys) in SECTION_OF.items():
        S.append(r'\subsection*{%s}' % title)
        for k in keys:
            src = srcmap[k]['src']
            S.append(r'\begin{align}')
            S.append(r'\frac{d^3N^{(2\rho)}}{d^3k}\bigg|_{\rm %s}={}&\NN\int %s\,\Big\{\notag\\' % (k, phase_tex(src)))
            rr = allres[k]
            for j, ((args, col), c) in enumerate(rr):
                kt = explicit_kernel(src, args)
                S.append(r'&\quad %s%s\notag\\' % (coeff_line(c, False), kt))
                ln = r'&\qquad\times %s' % col
                if j == len(rr) - 1:
                    S.append(ln + r'\Big\}' + cc_tex(src) + r'\label{eq:sum-%s}' % k)
                else:
                    S.append(ln + r'\notag\\')
            S.append(r'\end{align}')
    pickle.dump(allres, open('allres.pkl', 'wb'))
    return '\n'.join(body), '\n'.join(S), res


if __name__ == '__main__':
    body, summ, res = main()
    open('_body.tex', 'w').write(body)
    open('_summary.tex', 'w').write(summ)
    print('written')
