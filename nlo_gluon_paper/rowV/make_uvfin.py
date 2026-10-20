"""Writes Part2b_UV_finite.tex from uvfin_static.tex, part2b_dperp_checks.json and part2b_checks.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
N = json.load(open(os.path.join(HERE, 'part2b_dperp_checks.json')))
O = json.load(open(os.path.join(HERE, 'part2b_checks.json')))
s = open(os.path.join(HERE, 'uvfin_static.tex')).read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


sav = max(abs(r['r2_avg'] - r['expected']) for r in O['S_angular_average'])
dr = 'all 0' if all(O[k] == '0' for k in ('dimreg_theta_minus_target', 'dimreg_K_minus_target', 'dimreg_smooth_minus_target')) else 'FAIL'
r1 = max(r['r2_rest'] for r in N['rest_small_r'] if r['r'] == 1e-2)
r2 = max(r['r2_rest'] for r in N['rest_small_r'] if r['r'] == 1e-3)
sub = {
    'SD': '$%s$' % e(N['short_distance_dperp_ne_d_max_rel']),
    'SAV': r'max.\ dev.\ $%s$' % e(sav),
    'DR': dr,
    'PI': '$%s$' % e(N['p_integral_minus_closed_form_max']),
    'EXP': '$%s$' % N['dperp_over_d_expansion_O(eps^0)'],
    'CDR': '$%s$' % N['CDR_minus_previous_result_O(eps^0)'],
    'PREF': '$%s$' % N['prefactor_minus_g4Nc_16pi4'],
    'R1': '$%s$' % e(r1, 1),
    'R2': '$%s$' % e(r2, 1),
}
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'Part2b_UV_finite.tex'), 'w').write(s)
print('wrote Part2b_UV_finite.tex')
