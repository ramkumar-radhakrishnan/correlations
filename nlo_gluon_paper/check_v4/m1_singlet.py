"""Size of M1 in the singlet projection: Secs. 1+2 (as written) vs Sec. 3 collected form, U = U^T and general U."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'singlet')))
import env  # noqa
sys.path.append(HERE)
from core import Config  # noqa
from check_singlet import projected  # noqa
from parse_v4 import equations  # noqa
E = equations()
S12 = [t for e in E if e['sec'].startswith(('Genuine', 'Symmetrizing')) for t in e['terms']]
Jc = [e for e in E if e['sec'].startswith('JIMWLK') and e['kind'] == 'terms'][2]['terms']
out = {}
for sym in (True, False):
    r = []
    for s in (1, 2, 3):
        cfg = Config(s, sym)
        a, b = projected(S12, cfg), projected(Jc, cfg)
        r.append(abs(a - b) / abs(b))
    out['U=U^T' if sym else 'general U'] = r
print(out)
json.dump(out, open(os.path.join(HERE, 'm1_singlet.json'), 'w'), indent=1)
