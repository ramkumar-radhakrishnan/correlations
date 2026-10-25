"""Which kernels occur: is the measured gluon (w, w') emitted by a valence charge in both amplitudes?
For every term: does its kernel contain a factor K(c-w).K(c'-w') with c, c' charges (not z)?"""
import json, os
import env  # noqa
import sets as S2
from parse_v4 import equations
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
CH = {'x', "x'"}


def lo_factor(t):
    for f in t['kern']:
        a, b, c, d = f[1:]
        pts = {(a, b), (b, a)}, {(c, d), (d, c)}
        if any(('w' in p and (set(p) - {'w'}) <= CH) for p in [(a, b), (c, d)]) and \
           any(("w'" in p and (set(p) - {"w'"}) <= CH) for p in [(a, b), (c, d)]):
            if {a, b} != {c, d}:
                return True
    return False


def soft_kern(t):
    """the factor that is not the LO one, as a readable string"""
    return ' '.join('K(%s-%s).K(%s-%s)' % f[1:] for f in t['kern'])


E = equations()
A = S2.A_sets()
out = {}
sets = {'Secs 1+2': [t for e in E if e['sec'].startswith(('Genuine', 'Symmetrizing')) for t in e['terms']],
        '4rho.tex lists (P1)': A['4rho'],
        'Sec 3 (3.3)': [e for e in E if e['sec'].startswith('JIMWLK') and e['kind'] == 'terms'][2]['terms'],
        'Sec 4 (4.2)': [e for e in E if e['sec'].startswith('BFKL') and e['kind'] == 'terms'][1]['terms']}
for k, ts in sets.items():
    n_lo = sum(lo_factor(t) for t in ts)
    no = Counter(soft_kern(t) for t in ts if not lo_factor(t))
    out[k] = dict(terms=len(ts), with_charge_charge_LO_factor=n_lo, without=dict(no))
    print(k, len(ts), 'with LO factor:', n_lo, 'without:', dict(no))
json.dump(out, open(os.path.join(HERE, 'structure.json'), 'w'), indent=1)
