import re, sys, os
B = '/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/fixedbuild'
S = B + '/src/'
R = '/home/user/correlations/nlo_review/'

def rd(f):
    return open(S + f).read()

def removed(lines, what, color='orange'):
    t = '\n'.join(lines)
    t = t.replace(r'\begin{align}', r'\begin{align*}').replace(r'\end{align}', r'\end{align*}')
    t = t.replace(r'{\color{red}', '{').replace(r'{\color{orange}', '{').replace(r'{\color{teal}', '{')
    t = re.sub(r'\\subsubsection\*?\{([^}]*)\}', r'\\par\\noindent\\textit{\1}\\par', t)
    head = r'{\color{%s}\sloppy\noindent[\textbf{Removed} (%s):]\par' % (color, what)
    return head + '\n' + t + '\n}\\par\n'

def g2notes():
    t = rd('g2_notes.tex').split('%G2R4\n')
    return t[0].replace('%G2R1\n', ''), t[1]

def final():
    L = open(R + 'final_checked/Final_terms_checked.tex').read().split('\n')
    ln = lambda a, b: L[a - 1:b]          # 1-based inclusive
    n1, n4 = g2notes()
    out = []
    out += ln(1, 35)
    out.append(rd('header_final.tex'))
    out += ln(36, 229)
    out.append(rd('dglap_add.tex'))
    out += ln(230, 234)
    out.append(rd('rc_add.tex'))
    out += ln(235, 239)
    out.append(rd('lorc_add.tex'))
    # finite terms intro: line 241 contains the old red note on vee-k+
    l241 = L[240]
    old = r"{\color{red}[Group III Rows III, V, VI keep $\vee-k^{+}$ as the upper limit of the $p^{+}$ integral: in this form the integrals do not depend on $\vee$ up to corrections of order $k^{+}/\vee$.]}"
    assert old in l241
    l241 = l241.replace(old, r"{\color{orange}[Removed: \sout{Group III Rows III, V, VI keep $\vee-k^{+}$ as the upper limit of the $p^{+}$ integral: in this form the integrals do not depend on $\vee$ up to corrections of order $k^{+}/\vee$.}]} " + rd('finite_intro_add.tex').strip())
    out += ln(240, 240) + [l241] + ln(242, 242)
    out.append(rd('notation.tex'))
    out += ln(243, 250)
    out.append(rd('g1r2_add.tex'))
    out += ln(251, 263)
    out.append(removed(ln(264, 313), r'the separate equations $N_{1}$, $N_{2}$, $N_{3a}$, $N_{3b}$ and the notes of the earlier check; replaced by Eq.~\eqref{eq:fx-GIIII} and the pairs A, D of the block $\mathcal{S}$'))
    out.append(rd('g1r3_new.tex'))
    out += ln(314, 317)
    out.append(removed(ln(318, 337), r'first and second part of Row V; replaced by Eqs.~\eqref{eq:fx-V1}, \eqref{eq:fx-V2} and the pairs B, C of the block $\mathcal{S}$'))
    out.append(rd('g1r5_new.tex'))
    out += ln(338, 374)
    out.append(n1)
    out += ln(375, 383)
    out.append(n4)
    out += ln(384, 390)
    out.append(removed(ln(391, 404), r'first and second part of Row III; replaced by Eqs.~\eqref{eq:fx-G3R31}, \eqref{eq:fx-G3R32} and the pairs B, C of the block $\mathcal{S}$'))
    out.append(rd('g3r3_new.tex'))
    out += ln(405, 407)
    out.append(removed(ln(408, 410), r'Row V; replaced by Eq.~\eqref{eq:fx-G3R5} and the pair A of the block $\mathcal{S}$'))
    out.append(rd('g3r5_new.tex'))
    out += ln(411, 411)
    out.append(removed(ln(412, 414), r'Row VI; replaced by Eq.~\eqref{eq:fx-G3R6} and the pair D of the block $\mathcal{S}$'))
    out.append(rd('g3r6_new.tex'))
    out.append(rd('soft_block.tex'))
    if os.path.exists(S + 'p5_final.tex'):
        out.append(rd('p5_final.tex'))
    out.append(rd('singlet.tex'))
    if os.path.exists(S + 'notes_final.tex'):
        out.append(rd('notes_final.tex'))
    assert L[414] == r'\end{document}', L[414]
    out += ln(415, len(L))
    txt = '\n'.join(out)
    txt = txt.replace(r'\usepackage{cancel}', '\\usepackage{cancel}\n\\usepackage[normalem]{ulem}', 1)
    return txt

def limiting():
    L = open(R + 'limits_checked/Limiting_cases_final_1_checked.tex').read().split('\n')
    ln = lambda a, b: L[a - 1:b]
    n1, n4 = g2notes()
    out = []
    out += ln(1, 36)
    out.append(rd('header_lim.tex'))
    out += ln(37, 37)
    out.append(rd('finite_intro_add.tex'))
    out.append(rd('notation.tex'))
    out.append(rd('lim_defs.tex'))
    out += ln(38, 45)
    out.append(rd('g1r2_add.tex'))
    out += ln(46, 58)
    out.append(removed(ln(59, 108), r'the separate equations $N_{1}$, $N_{2}$, $N_{3a}$, $N_{3b}$ and the notes of the earlier check; replaced by Eq.~\eqref{eq:fx-GIIII} and the pairs A, D of the block $\mathcal{S}$'))
    out.append(rd('g1r3_new.tex'))
    out += ln(109, 118)
    out.append(rd('lim_rowIV_note.tex'))
    out += ln(119, 119)
    out.append(removed(ln(120, 139), r'first and second part of Row V; replaced by Eqs.~\eqref{eq:fx-V1}, \eqref{eq:fx-V2} and the pairs B, C of the block $\mathcal{S}$'))
    out.append(rd('g1r5_new.tex'))
    out += ln(140, 176)
    out.append(n1)
    out += ln(177, 185)
    out.append(n4)
    out += ln(186, 197)
    out.append(removed(ln(198, 211), r'first and second part of Row III; replaced by Eqs.~\eqref{eq:fx-G3R31}, \eqref{eq:fx-G3R32} and the pairs B, C of the block $\mathcal{S}$'))
    out.append(rd('g3r3_new.tex'))
    out += ln(212, 214)
    out.append(removed(ln(215, 217), r'Row V; replaced by Eq.~\eqref{eq:fx-G3R5} and the pair A of the block $\mathcal{S}$'))
    out.append(rd('g3r5_new.tex'))
    out += ln(218, 218)
    out.append(removed(ln(219, 221), r'Row VI; replaced by Eq.~\eqref{eq:fx-G3R6} and the pair D of the block $\mathcal{S}$'))
    out.append(rd('g3r6_new.tex'))
    out.append(rd('soft_block.tex'))
    if os.path.exists(S + 'p5_lim.tex'):
        out.append(rd('p5_lim.tex'))
    out.append(rd('singlet.tex'))
    if os.path.exists(S + 'notes_lim.tex'):
        out.append(rd('notes_lim.tex'))
    assert L[221] == r'\end{document}', L[221]
    out += ln(222, len(L))
    txt = '\n'.join(out)
    txt = txt.replace(r'\usepackage{cancel}', '\\usepackage{cancel}\n\\usepackage[normalem]{ulem}', 1)
    for i, lab in enumerate(['eq:fx-fun', 'eq:fx-split', 'eq:fx-DGLAP']):
        assert txt.count('\\label{%s}' % lab) == 1, lab
        txt = txt.replace('\\label{%s}' % lab, '\\tag{N.%d}\\label{%s}' % (i + 1, lab))
    txt = txt.replace('(see the NOTE at the DGLAP term)', r'(see the NOTE at the DGLAP term in \texttt{Final\_terms\_fixed})')
    return txt


def div():
    L = open(R + 'limits_checked/Divergent_terms_checked.tex').read().split('\n')
    ln = lambda a, b: L[a - 1:b]
    ref = r'Sec.~\ref{sec:dv-table}'
    def fix(s):
        return s.replace(r'\texttt{Divergent\_terms\_fixed}', ref)
    out = []
    out += ln(1, 35)
    out.append(rd('div_header.tex'))
    out += ln(36, 42)
    out.append(rd('div_rowII.tex'))
    out += ln(43, 53)
    out.append(rd('div_rowV.tex'))
    out += ln(54, 60)
    out.append(rd('div_add.tex'))
    out += ln(61, 68)
    out.append(rd('div_zero.tex'))
    out += ln(69, 85)
    l86 = L[85]
    old = 'is not cancelled by any other term.'
    assert l86.count(old) == 1
    l86 = l86.replace(old, r'{\color{teal}is not cancelled by any other term.} {\color{red}[\textbf{Superseded} (fix P3, P4; fix\_poles\_dglap Steps 4--6, 9): this pole is collinear (IR), multiplies the real $\mathcal{P}_{gg+}$ and is absorbed, together with the virtual $\frac{\beta_{0}}{2}\delta(1-\xi)$ of the dressed-gluon normalisation, by the $\overline{\rm MS}$ fragmentation counterterm with the full $\mathcal{P}_{gg}$; see the next note and ' + ref + '.]}')
    out.append(l86)
    dg = fix(rd('dglap_add.tex'))
    o2 = r'is exactly the remnant ($-\frac{1}{k^{+}}\times$BFKL structure) of Group II Row IV (its second term)'
    assert o2 in dg
    dg = dg.replace(o2, r'is exactly the remnant ($-\frac{1}{k^{+}}\times$BFKL structure) of Group II Row IV (its second term, in the finite terms, \texttt{Final\_terms\_fixed}; note (ii) above)')
    out.append(dg)
    out += ln(87, 90)
    out.append(fix(rd('rc_add.tex')))
    out += ln(91, 111)
    l112 = L[111]
    o3 = r'so the order of the limits $\vee\to\infty$ and $\epsilon\to0$ does not arise;'
    assert l112.count(o3) == 1
    l112 = l112.replace(o3, r'{\color{teal}so the order of the limits $\vee\to\infty$ and $\epsilon\to0$ does not arise;} {\color{red}[\textbf{Superseded} (fix P3; fix\_poles\_dglap Step 6): the order of the limits does arise: with all counterterms the poles add up to $-\frac{\ell}{2}\frac{2}{\epsilon}$ (in units $\mathfrak{a}\times$LO), zero only if $\vee\to\infty$ is taken before $\epsilon\to0$; this order is assumed and must be stated, ' + ref + '.]}')
    o4 = r'so the statement holds only together with the DGLAP term with $\delta(1-\xi)$ and the same $\mu$.'
    assert l112.count(o4) == 1
    l112 = l112.replace(o4, o4 + r' {\color{red}[\textbf{Fix P3} (fix\_poles\_dglap Steps 3, 5, 7): the presupposition is fulfilled, but not by the RC pole: the coupling counterterm $-\frac{\alpha_{s}\beta_{0}}{4\pi}\frac{2}{\epsilon_{\rm UV}}\times$LO is matched by the UV half of the (scaleless) normalisation of the dressed produced gluon ($Z_{3}$, Eq.~\eqref{eq:dv-norm}); the RC pole is collinear and cancels the IR half of that normalisation. Only in the single-$\epsilon$ bookkeeping, where the normalisation is zero, does the coupling counterterm cancel the RC pole formally.]}')
    out.append(l112)
    out.append(rd('lorc_add.tex'))
    out.append(rd('div_table.tex'))
    out.append(rd('notes_div.tex'))
    assert L[112] == r'\end{document}', L[112]
    out += ln(113, len(L))
    txt = '\n'.join(out)
    txt = txt.replace(r'\usepackage{cancel}', '\\usepackage{cancel}\n\\usepackage[normalem]{ulem}', 1)
    return txt

os.makedirs(B + '/out', exist_ok=True)
which = sys.argv[1]
if which == 'final':
    open(B + '/out/Final_terms_fixed.tex', 'w').write(final())
elif which == 'div':
    open(B + '/out/Divergent_terms_fixed.tex', 'w').write(div())
elif which == 'lim':
    open(B + '/out/Limiting_cases_final_fixed.tex', 'w').write(limiting())
