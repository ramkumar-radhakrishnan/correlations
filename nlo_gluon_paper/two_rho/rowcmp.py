"""Row-by-row comparison of the two-rho content: notes vs region-T reference (pointwise)."""
import pickle, re
import numpy as np
from collections import defaultdict
from core import Config, accumulate
import users as US
import reference as RF


def sub_of(tag):
    m = re.match(r'^(G[IV.]+[0-9]?|GI\.[IV]+[0-9]?|GII\.I[12]|GIII\.[IV]+|A[12]|B[12]|D)\.(\d)', tag)
    return int(m.group(2)) if m else None


def user_sets():
    G = US.genuine()
    r3 = US.R3()
    r4 = US.R4()
    S = defaultdict(list)
    for t in r3:
        key = t['tag'].split('.')[0] + '.' + t['tag'].split('.')[1]
        S[key].append(t)
    for t in r4:
        _, row, lab = t['tag'].split(':')
        line, k1, k2 = lab.rstrip('*').split('.')
        S['4:%s:%s' % (row, line) + ('.k2=%s' % k2 if row == 'G1' else '')].append(t)
    for k, v in G.items():
        S['G:' + k] = v
    S['D'] = US.D_terms()
    S['-RowIII2'] = US.neg(US.RowIII_part2())
    S['-RowIIaU'] = US.neg(US.RowII_alpha_U())
    return S


def subs(ts, ks):
    out = []
    for t in ts:
        n = int(t['tag'].split('.')[2][0])
        if n in ks:
            out.append(t)
    return out


if __name__ == '__main__':
    S = user_sets()
    R = RF.rows()
    print(sorted(S.keys()))
    pickle.dump(dict(S=dict(S)), open('user_sets.pkl', 'wb'))
