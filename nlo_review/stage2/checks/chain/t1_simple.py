# Task 1: Final_terms (F) vs Limiting_Cases (L), rows where both sides can be parsed directly from the LaTeX.
import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ltx import *
def ev(f, a, b, env, tagmap=None, pre=None):
    t = rhs(lines(f, a, b))
    if pre: t = pre(t)
    return contract(parse(t, env), tagmap)
def setsym(env):
    d = lines("L", 131, 132)
    X, Z = d.split(r",\qquad")
    env["sym"]["Xperp"] = parse(rhs(X), env).v; env["sym"]["Zperp"] = parse(rhs(Z).rstrip().rstrip("."), env).v
def fix166(t):
    assert t.count(r"\Bigg]\Bigg\}") == 1
    return t.replace(r"\Bigg]\Bigg\}", r"\Bigg]")
Q2P = {"plus[q@0]": "plus[p@0]"}
ROWS = [
 ("3rho GI  RowIII (F172-178 vs L120-127 'Adding them together')", ("F", 172, 178), ("L", 120, 127), {}),
 ("3rho GI  RowIV  (F182 vs L137-138)", ("F", 182, 182), ("L", 137, 138), {}),
 ("3rho GI  RowV   (F186-190 vs L143; L uses q+ in 2nd term, q+ -> p+)", ("F", 186, 190), ("L", 143, 143), dict(tagmap=Q2P)),
 ("3rho GII RowI 1a (F196-197 vs L166)", ("F", 196, 197), ("L", 166, 166), dict(pre=fix166)),
 ("3rho GII RowI 1b (F200-201 vs L169)", ("F", 200, 201), ("L", 169, 169), {}),
 ("3rho GIII RowI  (F212-214 vs L207)", ("F", 212, 214), ("L", 207, 207), {}),
 ("3rho GIII RowII (F218-219 vs L211-212)", ("F", 218, 219), ("L", 211, 212), {}),
 ("3rho GIII RowIII 1st part (F223-226 vs L217-218)", ("F", 223, 226), ("L", 217, 218), {}),
 ("3rho GIII RowV  (F232-234 vs L233-234)", ("F", 232, 234), ("L", 233, 234), {}),
 ("3rho GIII RowVI (F238 vs L238)", ("F", 238, 238), ("L", 238, 238), {}),
 ("2rho GI  RowII N1 (F246-247 vs L41, vee->inf)", ("F", 246, 247), ("L", 41, 41), {}),
 ("2rho GI  RowV 2b (F267-276 vs L150-159)", ("F", 267, 276), ("L", 150, 159), {}),
]
if __name__ == "__main__":
    for N, seed in ((2, 1), (3, 2), (2, 7)):
        print(f"--- N={N} seed={seed}")
        for pp in (None, 0.05, 1.3):
            env = make_env(N, seed); setsym(env)
            if pp is not None: env["pp"] = env["qp"] = pp * env["kp"]
            if pp is not None: print(f"  p+/k+ = {pp}")
            for name, fa, la, kw in ROWS:
                if pp is not None and "plus" not in str(kw) and not any(s in name for s in ("RowV", "RowIII 1st", "RowVI", "2b")): continue
                if pp == 1.3 and "GI " in name: continue   # p+ in (0,k+) there
                try:
                    A = ev(*fa, env); B = ev(*la, env, **kw)
                    report(name, A, B)
                except Exception as ex:
                    print(name, "ERROR", type(ex).__name__, str(ex)[:300])
