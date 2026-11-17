# Task 3: DGLAP and running-coupling sections of Final_terms (F) vs Group_II.tex Row IV (G2 lines 217-266)
import sys, re; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ltx import *
from t1_special import addd, scale, evt
def prep(t):
    t = re.sub(r"\\mathcal\s*P_\{gg\}\(\\xi\)", r"\\Pgg ", t)
    return t.replace(r"g^{2}(\mu)", "g^{2}").replace(r"\alpha_{s}(\mu)", r"\frac{g^{2}}{4\pi}").rstrip().rstrip(",").rstrip(".")
for N, seed in ((2, 31), (3, 32)):
    env = make_env(N, seed); env["sym"].update(musq=1.7, xi=0.41, beta0=11 * N / 3, eps=1e300); print(f"--- N={N} seed={seed}, xi={env['sym']['xi']}")
    F_D = evt(prep(rhs(lines("F", 136, 137))), env)
    G_D = evt(prep(rhs(lines("G2", 251, 252)).replace(r"\Bigg\{", "")), env)          # first term in the braces of G2:251-253, 1/eps -> 0
    report("DGLAP coefficient of P_gg: F136-137 vs G2:251-252 (2/eps dropped; lower limit 0 vs k+/vee not tested)", F_D, G_D)
    report("DGLAP: F136-137 vs G2:222-223", F_D, evt(prep(rhs(lines("G2", 222, 223))), env))
    F_RC = evt(prep(rhs(lines("F", 144, 145))), env)
    report("RC: F144-145 vs G2:256-257 (2/eps dropped)", F_RC, evt(prep(rhs(lines("G2", 256, 257)).rsplit(r"\\[6pt]", 1)[0]), env))
    report("RC: F144-145 vs G2:227-228", F_RC, evt(prep(rhs(lines("G2", 227, 228))), env))
    # the (beta0/2) delta(1-xi) part of DGLAP: coefficient of P at xi=1 times beta0/2, compared with -RC
    e1 = dict(env); e1["sym"] = dict(env["sym"]); e1["sym"]["xi"] = 1.0
    dl = evt(prep(rhs(lines("F", 136, 137))), e1)
    dl = {(k[0], ()): v * env["sym"]["beta0"] / 2 for k, v in dl.items()}
    report("(beta0/2) delta(1-xi) part of DGLAP (F136-141) vs  -RC (F144-145)", dl, scale(F_RC, -1.0))
    # LO+RC formula F149-150:  minus RC must be the same expression with the square bracket -> 1
    t = prep(rhs(lines("F", 149, 150))); i = t.index(r"\left[1+"); j = t.index(r"\right]", i) + len(r"\right]")
    report("LO+RC (F149-150) - RC (F144-145)  vs  F149-150 with [1+...] -> 1", addd(evt(t, env), scale(F_RC, -1.0)), evt(t[:i] + t[j:], env))
