"""Writes Part2b_UV_phase.tex from uvphase_static.tex and part2b_phase_checks.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'part2b_phase_checks.json')))
s = open(os.path.join(HERE, 'uvphase_static.tex')).read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


ab = max(max(r['A_rel'], r['B_rel']) for r in C['A_B'])
lim = r[-1] if False else C['A_B'][-1]
fin = C['finite_term_small_r']
sub = {
    'J': '$%s$' % C['J_minus_pi(2/eps-log(q2/mu2))'],
    'AB': '$%s$' % e(ab),
    'LIM': r'at $\lambda=10^{-4}$ the differences are $%s$ and $%s$' % (e(lim['m2A_minus_limit'], 1), e(lim['m2B_minus_limit'], 1)),
    'R1': '$%s$' % e(max(r['r2_rest'] for r in fin if r['r'] == 1e-2), 1),
    'R2': '$%s$' % e(max(r['r2_rest'] for r in fin if r['r'] == 1e-3), 1),
    'T1': ', '.join('$%.3f$' % v for v in C['large_r_tail'][0]['tail_R_to_10R']),
    'T2': ', '.join('$%s$' % (('%.3f' % v) if abs(v) > 1e-3 else e(v, 1)) for v in C['large_r_tail'][1]['tail_R_to_10R']),
    'DIFF': '; '.join(r'$q=%.1f$, $Y=%.1f$: $%.6f$ vs.\ $%.6f$' % (r['q'], r['Y'], r['numeric'], r['closed']) for r in C['difference_of_counterterms_2D']),
}
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'Part2b_UV_phase.tex'), 'w').write(s)
print('wrote Part2b_UV_phase.tex')
