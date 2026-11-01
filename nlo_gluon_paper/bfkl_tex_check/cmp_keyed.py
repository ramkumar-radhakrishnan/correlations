r"""Row by row, the 'written as' lists of the three-rho rows against the reference B(R) = B_old(R) + J(R)
(see cmp_written.py), compared key by key with ../two_rho/core.py: every kernel monomial at random points, the colour
numerically; terms whose colour does not depend on U(z) are integrated over z (dim. reg.), as in
../bfkl_final/two_rho_engine.py.  Also each row of Group I and Group II Row I alone and summed.  -> cmp_keyed.json"""
import json, os, sys
from collections import OrderedDict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cmp_written as W
TT = W.TT
import parse_v4 as P
sys.path.insert(0, os.path.join(HERE, '..', 'bfkl_final'))
from two_rho_engine import keyed, residual, neg


def J(EV, sub, subsub):
    return [t for e in EV if e['kind'] == 'terms' and e['sub'].startswith(sub) and e['subsub'].strip() == subsub
            for t in e['terms']]


if __name__ == '__main__':
    E = TT.equations()
    EV = [e for e in P.equations(W.EVOL) if e['sec'].startswith('Symmetrizing')]
    out = OrderedDict()
    groups = OrderedDict()
    for name, (user, eng, jrow) in W.MAP.items():
        ref = W.engine(eng) + J(EV, *jrow)
        u = W.pick(E, *user) if user else []
        groups[name] = (u, ref)
    groups['Group I (II-V) + Group II Row I'] = tuple(
        sum((groups[k][i] for k in ('G1 II', 'G1 III', 'G1 IV', 'G1 V(1)', 'G2 I first part')), []) for i in (0, 1))
    for name, (u, ref) in groups.items():
        res = []
        for seed in (1, 2):
            bad = residual(keyed({'user': u, 'ref': neg(ref)}, seed))
            res.append(len(bad))
            if seed == 1 and bad:
                worst = sorted(bad, key=lambda b: -b[1])[:3]
                out.setdefault(name + ' worst', [str(b[0]) + ' %.3e' % b[1] for b in worst])
        out[name] = res
        print(name, res, flush=True)
    json.dump(out, open(os.path.join(HERE, 'cmp_keyed.json'), 'w'), indent=1)
