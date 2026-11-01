"""Group I Rows II-V and Group II Row I (first part) of the 'written as' lists, summed, against the same reference."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cmp_written as W
TT = W.TT
import parse_v4 as P
E = TT.equations()
EV = [e for e in P.equations(W.EVOL) if e['sec'].startswith('Symmetrizing')]
def J(sub, subsub):
    return [t for e in EV if e['kind'] == 'terms' and e['sub'].startswith(sub) and e['subsub'].strip() == subsub for t in e['terms']]
uG1 = W.pick(E, 'Group I', 'Row II', '') + W.pick(E, 'Group I', 'Row III', '') + W.pick(E, 'Group I', 'Row IV', '')
rG1 = W.engine(['G1 II 2a', 'G1 III', 'G1 IV', 'G1 V(1)']) + J('Group I', 'Row II') + J('Group I', 'Row III') + J('Group I', 'Row IV') + J('Group I', 'Row V')
uG2 = W.pick(E, 'Group II', 'Row I', 'First')
rG2 = W.engine(['G2 I(1)']) + J('Group II', 'Row I')
out = {'Group I (II-V)': TT.both(uG1, rG1), 'Group II Row I first part': TT.both(uG2, rG2),
       'Group I + Group II Row I': TT.both(uG1 + uG2, rG1 + rG2)}
for k, v in out.items():
    print(k, v)
json.dump(out, open(os.path.join(HERE, 'cmp_written_groups.json'), 'w'), indent=1, default=float)
