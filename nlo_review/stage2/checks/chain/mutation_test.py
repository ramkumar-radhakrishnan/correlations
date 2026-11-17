# Sensitivity check of the comparison: random single-token mutations of the Final_terms side must be detected.
import sys, random; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ltx import *
from t1_simple import ROWS, setsym
random.seed(5)
SW = {"+": "-", "-": "+", "x": "z", "z": "w", "w": "x", "y": "x", "x'": "w'", "w'": "x'", "y'": "x'", "2": "3", "4": "8", "8": "4", "16": "8"}
def mutator(log):
    def m(t):
        t = list(t); d = 0; cand = []
        for j, tok in enumerate(t):
            if tok == r"\int": d = -1          # skip integral-sign decorations
            if d == -1:
                if tok == "}" and t[j + 1] not in ("^", "_"): d = 0
                continue
            if tok in SW: cand.append(j)
        j = random.choice(cand); log.append((j, t[j], " ".join(t[max(0, j - 6):j + 6]))); t[j] = SW[t[j]]; return t
    return m
env = make_env(2, 3); setsym(env)
tot = miss = 0
for name, fa, la, kw in ROWS:
    B = contract(parse(kw["pre"](rhs(lines(*la))) if "pre" in kw else rhs(lines(*la)), env), kw.get("tagmap"))
    for rep in range(40):
        log = []
        try: A = contract(parse(rhs(lines(*fa)), env, tokmut=mutator(log)))
        except Exception as ex: tot += 1; continue      # mutation produced an invalid formula: counted as detected
        md, ma, mb, _ = compare(A, B); tot += 1
        if np.isfinite(md) and md / max(ma, mb) < 1e-9: miss += 1; print("UNDETECTED:", name[:22], log[0])
print(f"mutations tried {tot}, undetected {miss}")
