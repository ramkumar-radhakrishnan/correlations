# Task 3: two-rho Group III Row III (2nd part): Group_III.tex 89-91 (G3) vs Other_terms 215-222 (O). Not a pointwise identity; partial checks.
import sys, re; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ltx import *
from t1_special import addd, scale, evt
def merge(d):
    out = {}
    for (r, tg), v in d.items(): out[(r, ())] = out.get((r, ()), 0) + v
    return out
tG = rhs(lines("G3", 89, 91)); mk = r"&  -\, \frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}\,{\color{blue}"
assert tG.count(mk) == 1; G1t, GLt = tG.split(mk); GLt = mk[1:] + GLt
logbr = r"{\color{blue}\left[\log{\left(\frac{k^{+}}{\Lambda}\right)}+\log{\left(\frac{\vee}{k^{+}}-1\right)}\right]}"; assert GLt.count(logbr) == 1
tO = rhs(lines("O", 215, 222))
parts = re.split(r"\\notag\s*\\\\\s*&\s*(?=[+-]\\,\s*\\frac\{1\}\{\(2\\pi\)\^\{3\}\})", tO); assert len(parts) == 4, len(parts)
O1, O2, O3, O4 = parts
a = O2.index(r"\notag \\ & -\frac{(x'-z)^{k'}}{(x'-z)^{2}}"); b = O2.index(r"\Bigg]\times", a); O2y = O2[:a] + O2[b:]      # keep only the K(x'-y') half
ologbr = r"\log{\left(1-\frac{k^{+}}{\vee}\right)}"; assert O4.count(ologbr) == 1
for N, seed in ((2, 41), (3, 42)):
    env = make_env(N, seed); env["vee"] = 6.3 * env["kp"]; env["Lam"] = 1e-9 * env["kp"]; print(f"--- N={N} seed={seed} p+/k+={env['pp']/env['kp']:.3f} vee/k+=6.3")
    report("(i) p+ integrand, plus brackets read as plain functions: G3:89-90[1st block]  vs  O215-216[1st block] + O216-218[2nd block, K(x'-y') half only]",
           merge(evt(G1t, env)), addd(merge(evt(O1, env)), merge(evt(O2y, env))))
    report("(ii) coefficient of the log bracket: G3:90-91 [log(k+/Lam)+log(vee/k+ -1)]  vs  O221-222 [log(1-k+/vee)]", evt(GLt.replace(logbr, "1"), env), evt(O4.replace(ologbr, "1"), env))
    e0 = dict(env); e0["pp"] = 1e-9 * env["kp"]
    T = {(k[0], ()): v * e0["pp"] for k, v in evt(G1t, e0).items() if k[1]}
    report("(iii) p+ * ([1/p+]_+ term of G3:89) at p+ -> 0   vs   coefficient of the log bracket in G3:90-91", T, evt(GLt.replace(logbr, "1"), env))
