# s13_V2_collinear.py -- Group I Row V second part (Sec. 6, eq:D-I.V-2o), operator-level colour: behaviour at the
# collinear point z -> x on circles |z-x| = lam: lam^2 * |angular mean| and lam^2 * max, with and without the 11/6 line.
# Run: python3 -I s13_V2_collinear.py -> s13_V2_collinear.log
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s10_uv_sec6 as S
field = S.lib.SmoothField(seed=7); out = []
for t in S.TERMS:
    if t.row == 'I.V' and '2nd part' in t.name:
        for seed in (3, 4):
            r = S.circle_mean(t, field, ['x', 'z'], nang=64, seed=seed)
            out.append('%-36s base %d: ' % (t.name, seed) + '  '.join('lam=%.0e: lam^2|mean|=%.3e lam^2 max=%.3e' % q for q in r))
            print(out[-1], flush=True)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 's13_V2_collinear.log'), 'w').write('\n'.join(out) + '\n')
