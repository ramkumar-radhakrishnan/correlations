"""Writes RowIV_FF_userform.tex from user_ff_static.tex and user_ff_checks.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'user_ff_checks.json')))
s = open(os.path.join(HERE, 'user_ff_static.tex')).read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


up1 = [r['upper1_minus_direct'] for r in C['user_form_vs_exact'] if r['Lam_over_kp'] < 1e-3 and r['kp'] == 1.0][0]
sub = {
    'INT': '$%s$' % C['int_dxi_over_xi(1-xi)_minus_user_bracket'],
    'EXACT': '$%s$' % e(C['user_form_exact_max_rel']),
    'UP1': '$%s$' % e(up1, 1),
    'WCDR': '$%s$' % C['W_CDR_minus_(-2/eps-2gamma-log)'],
    'PREF': '$%s$' % C['g2Nc_4pi2_minus_alpha_2pi_2Nc'],
    'CONV': '$%s$' % e(C['user_plus_to_standard_plus_max_rel']),
    'SHIFT': r'$2\log\frac{\vee}{k^+}+2\log(1-\frac{k^+}{\vee})$' if C['conversion_shift'] == '2*logvk + 2*log(1 - xi0)' else C['conversion_shift'],
    'TALLY': r'$\frac{11}{6}-\log(1-\xi_0)$' if C['tally_std'] == '11/6 - log(1 - xi0)' else C['tally_std'],
    'MOM': r'$\frac{11N_c}{6}$' if C['delta_coefficient_from_sum_rule'] == '11*N_c/6' else C['delta_coefficient_from_sum_rule'],
    'BA': '$%s$' % e(C['B_equals_A_with_plus_and_cuts_max_rel']),
    'MU': '$%s$' % C['dlogmu2_of_alpha_times_minus2_over_eps_over_alpha'],
}
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'RowIV_FF_userform.tex'), 'w').write(s)
print('wrote RowIV_FF_userform.tex')
