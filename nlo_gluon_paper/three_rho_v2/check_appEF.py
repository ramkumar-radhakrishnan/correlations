"""Compare Appendices E and F of main.tex (as pasted) with JIMWLK x LO and BFKL x LO of Evolution_from_final.tex.
Colour-singlet projection <rho^a(p) rho^b(q)> = delta^{ab} mu(p,q)/(Nc^2-1); units P mu, P = alpha_s^2 dY/(8 pi^6 (Nc^2-1)).
Appendix F: d/dY = -H, dN = sum of dN^(1..5) (bracket); Eq. (LOBFKL): -Nc x bracket."""
import json
import sys
import env
from core import Config
import target as TG
import refP
sys.path.insert(0, env.TWO)
import appF as AF
import evo_check as EC

out = []
for sym in (True, False):
    for seed in (1, 2, 3, 4):
        cfg = Config(seed, sym)
        J_note = AF.singlet_value(EC.Jw, cfg).real          # Eq. (3.2) of the note, as written there
        B_note = AF.singlet_value(EC.Bw, cfg).real          # Eq. (6.2) of the note, as written there
        J_op = AF.singlet_value(TG.jimwlk(), cfg).real      # operator form, Eq. (3.1)
        B_op = AF.singlet_value(refP.bfkl(), cfg).real      # operator form, Eq. (6.1)
        Fw = 0.5 * (AF.appF(cfg, 'p2', 'p1', corrected=False) + AF.appF(cfg, 'p1', 'p2', corrected=False))
        Fc = 0.5 * (AF.appF(cfg, 'p2', 'p1') + AF.appF(cfg, 'p1', 'p2'))
        E = -3 * 0.5 * (AF.lobfkl(cfg, 'p2', 'p1') + AF.lobfkl(cfg, 'p1', 'p2')).real
        r = dict(sym=sym, seed=seed, J_note=J_note, J_op=J_op, F_written=Fw.sum().real, F_corrected=Fc.sum().real,
                 B_note=B_note, B_op=B_op, E=E,
                 F_pieces_written=[float(v.real) for v in Fw], F_pieces_corrected=[float(v.real) for v in Fc])
        out.append(r)
        print('%s seed %d | JIMWLK: note %.6f  App.F as pasted %.6f  App.F corrected %.6f | BFKL: note %.6f  App.E %.6f'
              % ('sym U ' if sym else 'gen U ', seed, J_note, r['F_written'], r['F_corrected'], B_note, E))
json.dump(out, open('appEF.json', 'w'), indent=1)
