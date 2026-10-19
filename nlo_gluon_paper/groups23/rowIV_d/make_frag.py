"""Writes RowIV_FF_DGLAP.tex from frag_static.tex and frag_checks.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'frag_checks.json')))
s = open(os.path.join(HERE, 'frag_static.tex')).read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


sub = {
    'USER': '$%s$' % C['user_xi_form_minus_ours'],
    'CONV': '$%s$' % C['LO_conversion_minus_alpha_Nc_over_2pi_xi3'],
    'COMPACT': '$%s$' % C['convolution_form_minus_compact'],
    'CDR': '$%s$' % C['CDR_weight_minus_(1/epsbar-ell)'],
    'PLUS': '$%s$' % e(C['plus_prescription_exact_limits_max_rel']),
    'RL': ', '.join('$%.5f$' % v for v in C['R_Lambda_over_delta']),
    'RLE': '$%.5f$' % C['R_Lambda_over_delta_expected'],
    'MOM': r'$-\frac{11N_c}{6}$' if C['int_xi_Preal_plus'] == '-11*N_c/6' else C['int_xi_Preal_plus'],
    'DELTA': r'$\frac{11N_c}{6}\delta(1-\xi)$' if C['delta_coefficient_from_sum_rule'] == '11*N_c/6' else C['delta_coefficient_from_sum_rule'],
    'JAC': '$%s$' % e(C['jacobian_dz_over_z_with_dxi_over_xi3_rel']),
    'MU': '$%s$' % C['dlogmu2_of_alpha_over_epsbar_over_alpha'],
    'CII': '$%+.0f$' % C['colour_RowII_over_LO'],
    'CV': '$%+.0f$' % C['colour_RowV2b_over_LO'],
    'PII': r'$\frac{N_c}{2}$' if C['RowII_prefactor_over_(alpha_s/2pi)'] == 'N_c/2' else C['RowII_prefactor_over_(alpha_s/2pi)'],
    'PV': '$N_c$' if C['RowV2b_prefactor_over_(alpha_s/2pi)'] == 'N_c' else C['RowV2b_prefactor_over_(alpha_s/2pi)'],
    'TALLY': '$%s$' % C['delta_tally_minus_(11/6+3log((1+d)/(1-d))-log(1-xi0))'],
    'SER': r'6\delta+2\delta^3' if C['leftover_series'] == '2*delta**3 + 6*delta' else C['leftover_series'],
}
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'RowIV_FF_DGLAP.tex'), 'w').write(s)
print('wrote RowIV_FF_DGLAP.tex')
