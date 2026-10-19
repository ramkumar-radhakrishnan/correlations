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


def c(v):
    return r'$%.5f%+.5f\,i$' % (v[0], v[1])


sub = {
    'USER': '$%s$' % C['user_xi_form_minus_ours'],
    'PREF': '$%s$' % C['prefactor_g4_16pi4_minus_alpha_2pi_g2_2pi2'],
    'BR': '$%s$' % C['braces_explicit_form_minus_RowIV_dperp'],
    'BRC': '$%s$' % C['braces_explicit_form_minus_compact'],
    'CDR': '$%s$' % C['CDR_weight_minus_(1/epsbar-ell)'],
    'PLUS': '$%s$' % e(C['plus_prescription_explicit_phase_max_rel']),
    'A': '%.1f' % C['R_Lambda_a'],
    'RL': ', '.join(c(v) for v in C['R_Lambda_over_delta']),
    'RLE': c(C['R_Lambda_over_delta_expected']) + r' $=2(1-ia)e^{-ia}$',
    'FULL': '$%s$' % e(C['full_RowIV_plus_form_vs_direct_max_rel']),
    'MOM': r'$-\frac{11N_c}{6}$' if C['int_xi_Preal_plus'] == '-11*N_c/6' else C['int_xi_Preal_plus'],
    'DELTA': r'$\frac{11N_c}{6}\delta(1-\xi)$' if C['delta_coefficient_from_sum_rule'] == '11*N_c/6' else C['delta_coefficient_from_sum_rule'],
    'COMP': '$%s$' % e(C['composition_explicit_phase_max_rel']),
    'MU': '$%s$' % C['dlogmu2_of_alpha_over_epsbar_over_alpha'],
    'CII': '$%+.0f$' % C['colour_RowII_over_O'],
    'CV': '$%+.0f$' % C['colour_RowV2b_over_O'],
    'PII': r'$\frac{N_c}{2}$' if C['RowII_prefactor_over_unit'] == 'N_c/2' else C['RowII_prefactor_over_unit'],
    'PV': '$N_c$' if C['RowV2b_prefactor_over_unit'] == 'N_c' else C['RowV2b_prefactor_over_unit'],
    'TALLY': '$%s$' % C['delta_tally_minus_(11/6+3log((1+d)/(1-d))-log(1-xi0))'],
    'SER': r'6\delta+2\delta^3' if C['leftover_series'] == '2*delta**3 + 6*delta' else C['leftover_series'],
}
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'RowIV_FF_DGLAP.tex'), 'w').write(s)
print('wrote RowIV_FF_DGLAP.tex')
