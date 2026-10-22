"""Writes Finite_terms_final.tex from final_static.tex and checks.json."""
import json

C = json.load(open('checks.json'))
W = json.load(open('crows_w.json'))
s = open('final_static.tex').read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'$%s\times10^{%d}$' % (m, int(p))


v = C['verify_spread']
closed = ['G1-III', 'G1-IV', 'G2-I(1)', 'G2-I(2)', 'G3-I', 'G3-II']
plus = ['G1-V(1)', 'G3-III(1)', 'G3-III(2)', 'G3-V', 'G3-VI']
simp = C['simp_rel_at_1e6_1e9']
sub = {
    'G1III': e(v['G1-III']), 'G1IV': e(v['G1-IV']), 'G1V1': e(v['G1-V(1)']), 'G2I1': e(v['G2-I(1)']),
    'G2I2': e(v['G2-I(2)']), 'G3I': e(v['G3-I']), 'G3II': e(v['G3-II']),
    'G3III1': e(max(v['G3-III(1)'], v['G3-V'], v['G3-VI'])), 'G3III2': e(v['G3-III(2)']),
    'CLOSEDMAX': e(max(v[k] for k in closed)), 'PLUSMAX': e(max(v[k] for k in plus)),
    'SIMP': '%s / %s' % (e(max(x[0] for x in simp.values())), e(max(x[1] for x in simp.values()))),
    'LIMITS': r'$\le$ ' + e(max(C['Lv_coefficient_vs_P_minus_T'].values())),
    'TPLIM': r'$\le$ ' + e(max(max(r['A'], r['P']) for r in W)),
    'CW': '%s / %s' % (e(max(r['sub_1e3'] for r in W)), e(max(r['sub_1e5'] for r in W))),
    'IVI': e(C['G1IV_plus_G3I_integrand_max_rel']), 'COMB': e(C['G2I1A_plus_G3II_combined_max_rel']),
    'N2': r'%.1f--%.1f' % (min(C['N2_network_G1III_plus_G2I1G3II_rel']), max(C['N2_network_G1III_plus_G2I1G3II_rel'])),
    'LV2': r'$%s$, $%s$, $%s$, $%s$, $%s$' % (C['Lv_coef_2b'], C['Lv_coef_2c'], C['Lv_coef_2d'], C['Lv_coef_2e'], C['Lv_coef_2f']),
    'B0': r'$%s$, $%s$; $%s$' % (C['rowII_beta0_part'], C['rowV_beta0_part'], C['tally_beta0_coefficient_of_Ehat']),
    'TALLY': r'$%s$; $\frac{\alpha_s}{2\pi}N_c(\frac{67}{9}-\frac{\pi^2}{3})$' % C['tally_Lv_coefficient_of_Ehat'],
}
for k, val in sub.items():
    s = s.replace('@@%s@@' % k, val)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open('Finite_terms_final.tex', 'w').write(s)
print('wrote Finite_terms_final.tex')
