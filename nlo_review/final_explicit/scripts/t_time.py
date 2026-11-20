import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 's4_operator.py')).read()
src = src[:src.index("worst = 0; worst6")]
exec(src)
U = randUs(S, rng)
for name, Cfn, typ, pts in STRUCT[:3]:
    sites = (0, 1, 2); R = {p: rho_site[{pts[i]: sites[i] for i in range(3)}.get(p, 0)] for p in PTS}
    t = time.time(); lhs = rowop(S, U, R, name); t1 = time.time()
    T = Cfn(S, U); rp = [rho_site[s] for s in sites]; rhs = decomp(S, T, typ, rp, sites); t2 = time.time()
    print(name, '%.2f %.2f' % (t1-t, t2-t1), np.abs(lhs-rhs).max())
for name, Cfn, typ, pts in STRUCT:
    sites = (0, 0, 1); R = {p: rho_site[{pts[i]: sites[i] for i in range(3)}.get(p, 0)] for p in PTS}
    t = time.time(); lhs = rowop(S, U, R, name); lhs = lhs[0] if isinstance(lhs, tuple) else lhs; t1 = time.time()
    T = Cfn(S, U); rp = [rho_site[s] for s in sites]; rhs = decomp(S, T, typ, rp, sites); o2 = ordered(T, typ, rp); t2 = time.time()
    print(name, '%.2f %.2f' % (t1-t, t2-t1), np.abs(lhs-rhs).max(), np.abs(lhs-o2).max())
