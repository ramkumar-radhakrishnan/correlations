"""RowIV_dperp.tex from rowIV_static.tex and check_d.json"""
import json
C = json.load(open('check_d.json'))
s = open('rowIV_static.tex').read()


def e(x):
    if x == 0:
        return '$0$'
    m, p = ('%.0e' % x).split('e')
    return r'$%s\times10^{%d}$' % (m, int(p))


rep = {
    '@@POL@@': 'max.\\ difference ' + e(C['polarization_sum_max_diff']),
    '@@J2@@': '$%.12f$' % C['d2_traceless_numeric'],
    '@@J3T@@': '$%.12f$' % C['d3_trace_numeric'],
    '@@J3TC@@': '$%.12f$' % C['d3_trace_closed'],
    '@@J33@@': '$%.1e$' % C['d3_J33_numeric'],
    '@@J3R@@': e(C['d3_max_rel']),
    '@@SYM@@': 'difference $%s$' % C['contraction_minus_result'],
    '@@RED@@': 'difference $%s$' % C['reduction_to_earlier_dperp_2_minus_eps'],
    '@@XI@@': 'differences $%s$, $%s$, $%s$' % (C['xi_Phi'], C['xi_c2'], C['xi_c1']),
    '@@END@@': e(C['endpoint_rel']),
}
for k, v in rep.items():
    s = s.replace(k, v)
assert '@@' not in s
open('RowIV_dperp.tex', 'w').write(s)
print('wrote RowIV_dperp.tex')
