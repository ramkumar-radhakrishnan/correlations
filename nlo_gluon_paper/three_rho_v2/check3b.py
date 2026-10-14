"""Three-rho cancellation: as sent, after the typo fixes, after F4, after F3 (cumulative)."""
import json
from collections import defaultdict
import env
from core import symmetrize
import acc3 as A3
import parse as P
import fixes as FX


def nonzero(entries, seed):
    cfg = A3.Config3(seed, True)
    acc = defaultdict(lambda: defaultdict(complex))
    for e in entries:
        A3.accumulate3(symmetrize(P.as_term(e), (3,)), cfg, acc, label=e['block'])
    bad = A3.report3(acc, show=False)
    return len(acc), len(bad), max([b[0] for b in bad] + [0.0])


states = []
E0 = P.parse()
T = FX.typo_fixed()
S4 = [e for e in T if (e['block'], e['n']) not in FX.F4_DROP]
S34 = FX.corrected()
S3 = [e for e in T if not (e['block'] == 'I.III' and e['cls'] in ('gamma', "gamma'"))] + FX.d_entries()
for name, ents in (('as sent', E0), ('typos T1-T3 fixed', T), ('typos + F4', S4), ('typos + F3', S3), ('typos + F3 + F4', S34)):
    r = [nonzero(ents, s) for s in (5, 17, 23)]
    states.append(dict(state=name, terms=len(ents), keys=r[0][0], bad=[x[1] for x in r], maxres=max(x[2] for x in r)))
    print('%-22s terms %3d  kernel keys %2d  non-cancelling %s  largest residual %.2e' % (name, len(ents), r[0][0], [x[1] for x in r], max(x[2] for x in r)))
json.dump(states, open('states3.json', 'w'), indent=1)
