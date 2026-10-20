"""Writes Part2_check.tex from part2check_static.tex and part2_user_check.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'part2_user_check.json')))
s = open(os.path.join(HERE, 'part2check_static.tex')).read()


def e(x, d=0):
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


mx = lambda key, lab: max(r[key] for r in C[lab])
mn = lambda key, lab: min(r[key] for r in C[lab])
labs = ['Lambda/k+=0.05', 'Lambda/k+=1e-3', 'Lambda/k+=1e-5']
names = [('rel_2a_user', 'your 2a (plus prescription)'),
         ('rel_2b_correct', r'2b finite, corrected, exact limits (Eq.~\eqref{eq:corr})'),
         ('rel_2b_user', 'your 2b finite, final form'),
         ('rel_2b_fix_colourB2', r'\quad with E1 fixed'),
         ('rel_2b_fix_colourB2_E4', r'\quad with E1, E2 fixed'),
         ('rel_2b_fix_all', r'\quad with E1, E2, E3 fixed')]
rows = [r'%s & %s\\' % (lab, ' & '.join('$%s$' % e(mx(k, l), 1) for l in labs)) for k, lab in names]
sub = {
    'A3': '$%s$' % e(mx('rel_2a_user', 'Lambda/k+=1e-3'), 1),
    'A5': '$%s$' % e(mx('rel_2a_user', 'Lambda/k+=1e-5'), 1),
    'U': r'$%.1f$--$%.1f$' % (min(mn('rel_2b_user', l) for l in labs), max(mx('rel_2b_user', l) for l in labs)),
    'F5': '$%s$' % e(mx('rel_2b_fix_all', 'Lambda/k+=1e-5'), 1),
    'C': '$%s$' % e(max(mx('rel_2b_correct', l) for l in labs), 0),
    'ROWS': '\n'.join(rows),
}
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
assert '@@' not in s, [l for l in s.splitlines() if '@@' in l]
open(os.path.join(HERE, 'Part2_check.tex'), 'w').write(s)
print('wrote Part2_check.tex')
