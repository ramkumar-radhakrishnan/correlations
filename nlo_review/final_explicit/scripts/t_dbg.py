import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 's4_operator.py')).read()
exec(src[:src.index("worst = 0; worst6")])
U = randUs(S, rng)
for name, Cfn, typ, pts in STRUCT[-2:]:
    for sites in [(0,1,2),(0,0,1),(0,1,1),(0,1,0)]:
        R = {p: rho_site[{pts[i]: sites[i] for i in range(3)}.get(p, 0)] for p in PTS}
        simp, old = rowop(S, U, R, name)
        o2 = ordered(Cfn(S, U), typ, [rho_site[s] for s in sites])
        print(name, sites, np.abs(simp-o2).max(), np.abs(old-o2).max(), np.abs(simp-old).max())
