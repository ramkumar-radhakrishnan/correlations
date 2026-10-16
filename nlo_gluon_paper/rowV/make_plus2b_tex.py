"""Writes PlusPrescription_2b.tex from plus2b_static.tex, input_2b_finite.tex, plus2b_checks.json, plus2b_ir.json."""
import json


def e(x, d=0):
    """x in scientific notation, as LaTeX math"""
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


def lab(k):
    for a, b in (('->', r'$\to$'), ('p+', '$p^+$'), ('k+', '$k^+$'), ('xi', r'$\xi$'), ('F4', '$F_4$'),
                 ('I0', '$I_0$'), ('I1', '$I_1$')):
        k = k.replace(a, b)
    return k.replace('$$', '')


def cplx(v):
    return r'$%.10f%+.10f\,i$' % tuple(v)


s = open('plus2b_static.tex').read()
inp = open('input_2b_finite.tex').read().replace(r'\begin{align}', r'\begin{align}\label{eq:input}', 1)
C = json.load(open('plus2b_checks.json'))
IR = json.load(open('plus2b_ir.json'))
rows = []
for k, v in C.items():
    if k.startswith('config'):
        rows.append(r'%s: endpoint coefficients vs.\ limits of the integrand at $p^+\to0$, $k^+$ & rel.\ $%s$, $%s$ '
                    r'(limit taken at $10^{-8}$)\\' % (k, e(v['endpoint0_closed_vs_limit']), e(v['endpoint1_closed_vs_limit'])))
        rows.append(r'%s: plus-prescribed integrand of Eq.~(\ref{eq:plus}) vs.\ its definition & rel.\ $%s$\\'
                    % (k, e(v['note_integrand_vs_definition'])))
        rows.append(r'%s: your $p^+$ integral vs.\ plus + endpoint & %s\newline vs.\ %s\newline (rel.\ $%s$)\\'
                    % (k, cplx(v['original']), cplx(v['plus_form']), e(v['rel_diff'])))
for k, v in C['UV r^2 angular average'].items():
    rows.append(r'UV: $r^2\langle\cdot\rangle$, %s, $r=10^{-2},10^{-3},10^{-4}$ & %s\\'
                % (lab(k), ', '.join('$%s$' % e(x, 2) for x in v)))
c = C['combined endpoint form vs I0+I1 (grid, symmetric weight)']
rows.append(r'combined form \eqref{eq:endc} vs.\ Eq.~\eqref{eq:end}, $x,z$ summed on a grid with a symmetric weight '
            r'& rel.\ $%s$\\' % e(c['rel_diff']))
irrows = [r'%s & %s\\' % (lab(k), ', '.join('$%s$' % e(x, 3) for x in v)) for k, v in IR.items()]
s = s.replace('@@INPUT@@', inp).replace('@@CHECKS@@', '\n'.join(rows)).replace('@@IR@@', '\n'.join(irrows))
open('PlusPrescription_2b.tex', 'w').write(s)
print('wrote PlusPrescription_2b.tex')
