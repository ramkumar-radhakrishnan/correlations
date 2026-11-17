# Extra (internal to Limiting_Cases, Group I Row III): do the parts N1+N2+N3a+N3b add up to the "Adding them together" formula?
import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ltx import *
from t1_simple import setsym
from t1_special import addd, evt, ev
for N, seed in ((2, 1), (2, 5), (3, 2), (2, 9)):
    env = make_env(N, seed); setsym(env)
    d = lines("L", 116, 116); parts = d.split("& =")
    a1 = parse(parts[1].rsplit(r"\notag", 1)[0], env); a2 = parse(parts[2], env)
    print(f"--- N={N} seed={seed}:  A0 (arctan form) = {a1.real:.12g}   A0 (arccos form) = {a2.real:.12g}")
    env["sym"]["A0"] = a1
    report(" L N1 first form (57-60, three logs) vs N1 after 'adding first and second' (68-70)", ev("L", 57, 60, env), ev("L", 68, 70, env))
    S = addd(ev("L", 68, 70, env), ev("L", 73, 78, env), ev("L", 81, 102, env), ev("L", 105, 113, env))
    T = ev("L", 120, 127, env)
    report(" L N1+N2+N3a+N3b (68-113) vs L 'Adding them together' (120-127)", S, T)
    # same, but leaving out the last term of each side (the U(w) term, kernel K(x'-w').K(y-w))
    t3b = rhs(lines("L", 105, 113)); m = r"& -\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\int_{x,y,z,x'}\log{\left[\frac{(z-w)^{2}}{(x-w)^{2}}\right]}\,\frac{(x'-w')\cdot(y-w)"
    assert t3b.count(m) == 1; t3b_head, t3b_last = t3b.split(m); t3b_last = m[1:] + t3b_last
    tS = rhs(lines("L", 120, 127)); m2 = r"&-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{abd}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{w,w'}e^{-ik\cdot(w'-w)}\int_{x,y,z,x'}\log\left(\frac{(z-w)^{2}}{(x-w)^{2}}\right)"
    assert tS.count(m2) == 1; tS_head, tS_last = tS.split(m2); tS_last = m2[1:] + tS_last
    S0 = addd(ev("L", 68, 70, env), ev("L", 73, 78, env), ev("L", 81, 102, env), evt(t3b_head + r"\Bigg\}", env))
    report(" ... all terms except the last one of each side", S0, evt(tS_head, env))
    pre = r"\frac{1}{(2\pi)^{3}}\int_{w,w'}e^{-ik(w'-w)}\Bigg\{"
    report(" ... last term only: L113 (colour [U(z)-U(y)]) vs L126-127 (colour [U(w)-U(y)])", evt(pre + t3b_last, env), evt(tS_last, env))
