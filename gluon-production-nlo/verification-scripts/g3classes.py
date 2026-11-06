"""Per-class diagnosis of the 3rho cancellation: group collected terms by the
z-kernel class (alpha: K(x-z).K(y-z), alpha': K(x'-z).K(y-z), beta: K(w-z).K(y-z),
beta': K(w'-z).K(y-z)) and report sector-3 norms per class, with and without h.c."""
import numpy as np
from engine import Config, accumulate, total_norm
from g3cancel import get_block_terms

def zclass(term):
    for kf in term.kern:
        if kf[0] == 'K':
            pair = {kf[1], kf[2]}, {kf[3], kf[4]}
            for A, B in (pair, pair[::-1]):
                if 'z' in A and 'z' in B:
                    a = (A - {'z'}).pop(); b = (B - {'z'}).pop()
                    s = {a, b}
                    if s == {'x', 'y'}: return 'alpha'
                    if s == {'xp', 'y'}: return 'alphap'
                    if s == {'w', 'y'}: return 'beta'
                    if s == {'wp', 'y'}: return 'betap'
                    return f'other({a},{b})'
    return 'none'

if __name__ == "__main__":
    per = get_block_terms()
    allterms = [t for ts in per.values() for t in ts]
    classes = {}
    for t in allterms:
        classes.setdefault(zclass(t), []).append(t)
    print("class sizes:", {k: len(v) for k, v in classes.items()})
    cfg = Config(6, "su2", seed=11)
    for cname, ts in sorted(classes.items()):
        tA = accumulate(ts, cfg, sector_filter={3})
        tB = accumulate(ts + [t.hc() for t in ts], cfg, sector_filter={3})
        # reference: largest single block contribution in this class
        ref = max(total_norm(accumulate([t], cfg, sector_filter={3})).get(3, 0) for t in ts)
        print(f"  {cname}: ref={ref:.3e}  |sum|={total_norm(tA).get(3,0):.3e}  |sum+hc|={total_norm(tB).get(3,0):.3e}")
