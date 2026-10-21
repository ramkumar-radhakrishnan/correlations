"""Normalization of the targets against the formulas sent:
   sum of the log(vee/Lambda) two-rho terms (units N)  <->  JIMWLK of App. F (corrected), prefactor alpha_s^2 dY/(8pi^6(Nc^2-1));
   sum of the log(vee) two-rho terms (units N)          <->  Eq. (LOBFKL) of App. E.
Singlet projection <rho^a rho^b> = delta^{ab} mu/(Nc^2-1); one point (p1, p2, w, w', z), symmetrized in p1 <-> p2."""
import json
import env  # noqa
import sys, os
sys.path.insert(0, env.TWO)
from core import Config
import target as TG
import refP
from appF import singlet_value, appF, lobfkl

NC = 3
out = []
for sym in (True, False):
    for seed in (1, 2, 3):
        cfg = Config(seed, sym)
        J = singlet_value(TG.jimwlk(), cfg).real
        B = singlet_value(refP.bfkl(), cfg).real
        Fc = (0.5 * (appF(cfg, 'p2', 'p1') + appF(cfg, 'p1', 'p2'))).sum().real
        Fw = (0.5 * (appF(cfg, 'p2', 'p1', corrected=False) + appF(cfg, 'p1', 'p2', corrected=False))).sum().real
        E = -NC * (0.5 * (lobfkl(cfg, 'p2', 'p1') + lobfkl(cfg, 'p1', 'p2'))).real
        out.append(dict(sym=sym, seed=seed, JIMWLKxLO=J, AppF_corrected=Fc, AppF_as_sent=Fw, BFKLxLO=B, AppE=E))
        print('sym' if sym else 'gen', seed, 'J %.6f F(corr) %.6f F(sent) %.6f | B %.6f AppE %.6f' % (J, Fc, Fw, B, E))
json.dump(out, open('targets_check.json', 'w'), indent=1)
