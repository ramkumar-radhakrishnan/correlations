"""Writes RowIV_plus.tex from plus_static.tex and plus_check.json."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'plus_check.json')))
s = open(os.path.join(HERE, 'plus_static.tex')).read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


G0 = np.array(R['C0'])
asym = []
for r in R['asymptotics']:
    if r['q'] in (10.0, 20.0, 30.0, 45.0):
        asym.append(r'$%.1f$ & $%g$ & $%.5f$ & $%.5f$ & $%.4f$ & $%.4f$ & $%.4f$\\' % (
            r['ang'], r['q'], r['q2F_T'][0], r['pred_T'], r['q2F_T'][0] / r['pred_T'], r['q2F_E'][0], r['q2F_E'][0] / (4 * np.pi)))
X = R['xi_integrals']
xi = []
for r in X['rows']:
    if 'dT' in r:
        xi.append(r'$%g$ & $%.6f$ & $%.6f$ & $%.6f$ & $%.4f$ & $%.4f$ & $%.4f$\\' % (
            r['xi0'], r['T'][0], r['dT'][0], r['dT_pred'], r['TE'][0], r['dTE'][0], r['dTE_pred']))
    else:
        xi.append(r'$%g$ & $%.6f$ & -- & -- & $%.4f$ & -- & --\\' % (r['xi0'], r['T'][0], r['TE'][0]))
up = []
for r in R['upper_end']['rows']:
    up.append(r'$%g$ & $%s$ & $%s$ & $%.3f$\\' % (r['Lambda_over_kp'], e(r['strip'][0], 3), e(r['pred'], 3), r['strip'][0] / r['pred']))
sub = {
    'G0': r'\begin{pmatrix}%.4f&%.4f\\%.4f&%.4f\end{pmatrix}' % (G0[0, 0], G0[0, 1], G0[1, 0], G0[1, 1]),
    'CAL': e(R['calibration_max_abs_err_times_q2'], 1),
    'ASYM': '\n'.join(asym),
    'K': '%g' % X['k'],
    'XI': '\n'.join(xi),
    'T0': '%.5f' % X['T_extrapolated_to_0'],
    'F1': '%.5f' % R['upper_end']['F_at_xi1'][0],
    'UP': '\n'.join(up),
}
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'RowIV_plus.tex'), 'w').write(s)
print('wrote RowIV_plus.tex')
