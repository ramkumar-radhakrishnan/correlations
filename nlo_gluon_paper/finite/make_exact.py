"""Writes Finite_terms_exact.tex from exact_static.tex, checks.json, checks_exact.json (with simp_exact.json)."""
import json

C = json.load(open('checks.json'))
X = json.load(open('checks_exact.json'))
S = X['simp_exact']
s = open('exact_static.tex').read()


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'$%s\times10^{%d}$' % (m, int(p))


def worst(names):
    return max(max(max(r) for r in S[n]) for n in names)


v = C['verify_spread']
closed = ['G1-III', 'G1-IV', 'G2-I(1)', 'G2-I(2)', 'G3-I', 'G3-II']
plus = ['G1-V(1)', 'G3-III(1)', 'G3-III(2)', 'G3-V', 'G3-VI']
four = X['four_rho_symmetrized_max_sum_and_scale']
fa = four['all lines']
fv = four['upper limit vee-k+ (L3, L5-L10)']
n2 = [r for rows in S['N2_network_rel'].values() for r in rows]
sub = {
    'G1III': e(v['G1-III']), 'G1IV': e(v['G1-IV']), 'G1V1': e(v['G1-V(1)']), 'G2I1': e(v['G2-I(1)']),
    'G2I2': e(v['G2-I(2)']), 'G3I': e(v['G3-I']), 'G3II': e(v['G3-II']),
    'G3III1': e(max(v['G3-III(1)'], v['G3-V'], v['G3-VI'])), 'G3III2': e(v['G3-III(2)']),
    'CLOSEDMAX': e(max(v[k] for k in closed)), 'PLUSMAX': e(max(v[k] for k in plus)),
    'EXCL': e(worst(closed)),
    'EXINT': e(worst(['G1-V(1)', 'G3-III(1)', 'G3-III(2)', 'G3-V', 'G3-VI'])),
    'EXV2': e(worst(['G1-V(2a)', 'G1-V(2b fin)'])),
    'MASTERS': e(S['masters']), 'NOTEFORMS': e(S['note_forms']),
    'POLE': '$' + X['pole_delta_coefficient_units_as2pi_m2eps_Nc'].replace('log', r'\log') + '$',
    'B0': 'difference $%s$' % X['beta0_in_rows_II_plus_V'],
    'FIN': 'difference $%s$' % X['finite_minus_expected'],
    'FOURALL': r'max.\ $|\Sigma|$ %s (terms up to %.0f)' % (e(max(r[0] for r in fa)), max(r[1] for r in fa)),
    'FOURVK': r'max.\ $|\Sigma|$ = %.1f--%.1f (terms up to %.0f)' % (min(r[0] for r in fv), max(r[0] for r in fv),
                                                                    max(r[1] for r in fv)),
    'IVI': e(C['G1IV_plus_G3I_integrand_max_rel']), 'COMB': e(C['G2I1A_plus_G3II_combined_max_rel']),
    'N2EX': '%.2f--%.1f' % (min(n2), max(n2)),
}
for k, val in sub.items():
    s = s.replace('@@%s@@' % k, val)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open('Finite_terms_exact.tex', 'w').write(s)
print('wrote Finite_terms_exact.tex')
