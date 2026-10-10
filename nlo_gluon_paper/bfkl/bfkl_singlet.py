"""BFKL x LO (operator form, units N_P) in the colour-singlet projection vs Eq. (LOBFKL) of Appendix E."""
import env  # noqa
import json
from core import Config
import appF as AF
import refP
B = refP.bfkl()
out = []
for sym in (True, False):
    for seed in (1, 2, 3):
        cfg = Config(seed, sym)
        op = AF.singlet_value(B, cfg)
        eq = -AF.NC * 0.5 * (AF.lobfkl(cfg, 'p2', 'p1') + AF.lobfkl(cfg, 'p1', 'p2'))
        out.append(dict(sym=sym, seed=seed, op=op.real, eq=eq.real))
        print('symmetric' if sym else 'general  ', seed, ' operator form %.10f   Eq.(LOBFKL) %.10f   ratio %.12f' % (op.real, eq.real, op.real / eq.real))
json.dump(out, open('bfkl_singlet.json', 'w'), indent=1)
