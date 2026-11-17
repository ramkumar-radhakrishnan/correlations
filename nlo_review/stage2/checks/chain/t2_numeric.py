# Task 2 (numerical confirmation): Other_terms (O) evaluated at vee = 1e300 (k+/vee -> 0, log(1-k+/vee) -> 0) vs Limiting_Cases (L)
import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ltx import *
from t2_diff import PAIRS
from t1_special import evt
def fixL(name, t):
    if name == "GII RowI 1a": t = t.replace(r"\Bigg]\Bigg\}", r"\Bigg]")
    return t
for N, seed in ((2, 3), (3, 4)):
    env = make_env(N, seed); env["sym"]["xi"] = 0.37; print(f"--- N={N} seed={seed}")
    dL = lines("L", 116, 116).split("& =")[1].rsplit(r"\notag", 1)[0]; dO = lines("O", 114, 114).split("& =")[1].rsplit(r"\notag", 1)[0]
    for name, a, b, oa, ob in PAIRS:
        if name in ("GI RowII N2b", "GI RowIII A0"): continue
        try:
            env["sym"]["A0"] = parse(dL, env); A = evt(fixL(name, rhs(lines("L", a, b))), env)
            env["sym"]["A0"] = parse(dO, env); B = evt(rhs(lines("O", oa, ob)), env)
            report(f"{name:16s} L{a}-{b} vs O{oa}-{ob}", A, B)
        except Exception as ex: print(name, "ERROR", type(ex).__name__, str(ex)[:200])
