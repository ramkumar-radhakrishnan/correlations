"""Writes Group_I_v2_check.tex from v2_static.tex and check_v2.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'check_v2.json')))
s = open(os.path.join(HERE, 'v2_static.tex')).read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


R5 = C['RowV_2b_finite']
sub = {
    'R2': '$%s$' % e(C['RowII_part2_new_bracket_vs_main']),
    'R3': '$%s$' % e(C['RowIII_new_integrand_vs_first_principles_E3']),
    'N2': '$%s$' % e(C['RowIII_new_N2_vs_PV_integral']),
    'GREEN': '$%s$' % e(C['RowIII_green_terms_sum_over_one_term']),
    'SW': r'$%.2f$--$%.2f$' % (min(r['simplifying_written'] for r in R5), max(r['simplifying_written'] for r in R5)),
    'SF': '$%s$' % e(max(r['simplifying_fixed'] for r in R5)),
    'FIN3': '$%s$' % e(max(r['final_written'] for r in R5 if r['Lam_over_kp'] == 1e-3), 1),
    'FIN5': '$%s$' % e(max(r['final_written'] for r in R5 if r['Lam_over_kp'] == 1e-5), 1),
    'B0': '$%s$' % e(C['beta0_RowII_plus_RowV2b_per_ordering']),
    'B0X': '$%s$' % C['beta0_and_powers_exact_in_Lambda'],
}
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'Group_I_v2_check.tex'), 'w').write(s)
print('wrote Group_I_v2_check.tex')
