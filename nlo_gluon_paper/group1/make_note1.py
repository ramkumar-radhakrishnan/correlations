"""Writes Group_I_check.tex from check_static.tex and checks1.json (plus ../groups23/main/main_checks.json for the
non-instantaneous entries of B_3)."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'checks1.json')))
M = json.load(open(os.path.join(HERE, '..', 'groups23', 'main', 'main_checks.json')))
s = open(os.path.join(HERE, 'check_static.tex')).read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


def rng(vals):
    return r'$%s$ to $%s$' % (e(min(vals)), e(max(vals)))


# B3 table
B = C['B3_vs_wavefunction']
lt = [r for r in B if r['pp_over_kp'] < 1]
gt = [r for r in B if r['pp_over_kp'] > 1]
mean = lambda rows, key: sum((r[key][0] if isinstance(r[key], list) else r[key]) for r in rows) / len(rows)
ni = M['B3 / wave function (8 pi^3 = 248.05021344)']
rows = [
    r'instantaneous (m5), new & $%+.6f$ & $%+.6f$\\' % (mean(lt, 'inst'), mean(gt, 'inst')),
    r'two-$\rho$, $\{\rho,\rho\}$ (m1, m1), new & $%+.6f$ & $%+.6f$\\' % (mean(lt, 'anti'), mean(gt, 'anti')),
    r'two-$\rho$, $[\rho,\rho]$ (m1, m1), new & $%+.6f$ & $%+.6f$\\' % (mean(lt, 'comm'), mean(gt, 'comm')),
    r'non-instantaneous (m2, m3), \texttt{groups23/main} & $%+.6f$ & $\mathbf{%+.6f}$ (E3)\\'
    % (ni['Theta(k+-Lambda-p+) entry']['mean_real'], ni['Theta(p+-k+-Lambda) entry']['mean_real']),
]
s = s.replace('@@B3TABLE@@', '\n'.join(rows))

s = s.replace('@@V1RANGE@@', rng([r['written_0_k'] for r in C['RowV_part1_plus']]))
s = s.replace('@@V2RANGE@@', rng([r['written_0_k'] for r in C['RowV_2a_plus']]))
s = s.replace('@@BR1@@', 'rel. $%s$' % e(C['RowII_part2_notes_minus_main_minus_k_logvk']))
s = s.replace('@@BR2@@', 'rel. $%s$' % e(C['RowII_part2_split_2a_2b_vs_main']))
s = s.replace('@@BR3@@', 'rel. $%s$' % e(C['RowII_part2_prefactor_colour']))
s = s.replace('@@BR4@@', 'rel. $%s$' % e(C['RowII_part2_exact_split']))
s = s.replace('@@DIMREG@@', '$%s$' % C['dimreg_int_K_x_w_z_over_pi_minus_bracket'])
s = s.replace('@@FOUR@@', 'rel. $%s$' % e(max(max(r['RowIII'], r['RowIV']) for r in C['four_rho'])))

III = C['RowIII_integrand']
w_all = max(r['written'] for r in III)
e3_lt = max(r['E3'] for r in III if r['pp'] < 1.3)
e3_gt = min(r['E3'] for r in III if r['pp'] > 1.3)
ints = C['RowIII_p_integrals']
pvmax = max(max(r[L] for L in ('L1', 'L2', 'L3', 'L4', 'L5')) for r in ints)
cl = [r[L + '_complex_log'] for r in ints for L in ('L1', 'L2', 'L3', 'L4', 'L5') if r[L + '_pole_inside']]
IV = C['RowIV']
V = C['RowV_integrands']
all_rows = [
    (r'$\mf B_3$ entries vs.\ light-cone wave functions (Table~\ref{tab:B3})', r'ratio $8\pi^3$, except $-8\pi^3$ for the non-inst.\ $p^+>k^+$ entry'),
    (r'$\frac1\pi\int_z[K(x,w,z)-K(w,w,z)]$ in $d_\perp=2-\epsilon$, $\overline{\rm MS}$, minus the Row II bracket', r'$%s$ (\texttt{sympy})' % C['dimreg_int_K_x_w_z_over_pi_minus_bracket']),
    (r'Row II part 1 (1.1) vs.\ first principles (any bracket; symmetric and general $U$)', r'rel.\ $%s$' % e(C['RowII_part1_vs_first_principles'])),
    (r'Row II rewriting (1.2)', r'$%s$' % e(C['RowII_rewriting'])),
    (r'Row II part 2: (1.3) as written minus \texttt{Calculation\_main.tex} bracket minus $k^i\log\frac{\vee}{k^+}$', r'rel.\ $%s$' % e(C['RowII_part2_notes_minus_main_minus_k_logvk'])),
    (r'Row II part 2: split 2a $+$ 2b vs.\ full bracket (random $\bk,\bp$, $k^+$, $\Lambda$, $\vee$)', r'rel.\ $%s$' % e(C['RowII_part2_split_2a_2b_vs_main'])),
    (r'Row II part 2: prefactor and colour vs.\ first principles', r'rel.\ $%s$' % e(C['RowII_part2_prefactor_colour'])),
    (r'Row II part 2: exact bracket (X4, X7) $=(E+2k+D)[\dots]+R_{\rm exact}$', r'rel.\ $%s$' % e(C['RowII_part2_exact_split'])),
    (r'Row III (1.10) vs.\ first principles with $\mf B_3$ as written (48 points)', r'rel.\ $%s$' % e(w_all)),
    (r'Row III (1.10) vs.\ first principles with E3', r'$p^+<k^+$: rel.\ $%s$; $p^+>k^+$: differ by $\geq%.0f\%%$' % (e(e3_lt), 100 * e3_gt)),
    (r'Row III: $N_1,N_2,N_{3a},N_{3b}$ vs.\ numerical PV integrals, logs as $\log|\cdot|$ (6 configurations)', r'rel.\ $%s$' % e(pvmax)),
    (r'Row III: the same with complex logarithms, pole inside the range', r'rel.\ %s' % rng(cl)),
    (r'Row III: green terms with G1', r'rel.\ $%s$' % e(max(r['rational_E3'] for r in ints))),
    (r'Row IV (1.15) vs.\ first principles (symmetric and general $U$)', r'rel.\ $%s$' % e(max(r['integrand'] for r in IV))),
    (r'Row IV (1.16) vs.\ numerical $\int_\Lambda^{\vee-k^+}$ of (1.15)', r'rel.\ $%s$' % e(max(r['integral'] for r in IV))),
    (r'Row V part 1 (1.17) vs.\ first principles', r'rel.\ $%s$' % e(max(r['part1'] for r in V))),
    (r'Row V part 2 (1.20) vs.\ first principles', r'rel.\ $%s$' % e(max(r['part2'] for r in V))),
    (r'Four-$\rho$ parts of Rows III and IV (Sec.~\ref{sec:four})', r'rel.\ $%s$' % e(max(max(r['RowIII'], r['RowIV']) for r in C['four_rho']))),
]
s = s.replace('@@ALLTABLE@@', '\n'.join(r'%s & %s\\' % a for a in all_rows))

trows = []
for r in ints:
    lab = r'%s $U$, $k^+=%.1f$, $\Lambda=%.2f$, $\vee=%.0f$, $p^+_*=%.3g$' % ('symm.' if r['sym'] else 'general', r['kp'], r['Lam'], r['vee'], r['pole'])
    parts = ', '.join(r'$L_%d%s$: $%s$' % (i, r'^\ast' if r['L%d_pole_inside' % i] else '', e(r['L%d' % i])) for i in range(1, 6))
    parts += r', G1: $%s$' % e(r['N2_E3'])
    trows.append(r'%s & %s\\' % (lab, parts))
s = s.replace('@@IIITABLE@@', '\n'.join(trows))

vr = []
for a, b in zip(C['RowV_part1_plus'], C['RowV_2a_plus']):
    vr.append(r'%d & %.3f & $%s$ & $%s$ & $%s$ \quad / \quad $%s$\\' % (a['seed'], a['Lam'] / 1.3, e(a['written_0_k']), e(a['limits_Lam']),
                                                                         e(b['written_0_k']), e(b['limits_Lam'])))
s = s.replace('@@VPLUSTABLE@@', '\n'.join(vr))
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'Group_I_check.tex'), 'w').write(s)
print('wrote Group_I_check.tex')
