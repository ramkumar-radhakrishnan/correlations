# character-level diff of LaTeX line ranges after stripping spacing (adapted from nlo_review/stage2/checks/tdiff.py)
import re, sys, difflib
UP = "/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
F = {"I": "16e62d79-Group_I.tex", "II": "24c35375-Group_II.tex", "III": "96b5b029-Group_III.tex", "O": "f10fa593-Other_terms.tex",
     "L": "56ca25e0-Limiting_Cases.tex", "F": "690ba62f-Final_terms.tex"}
def get(g, a, b):
    L = open(UP + F[g]).read().split("\n"); return "\n".join(L[a - 1:b])
def norm(s):
    s = re.sub(r"\\notag|\\nonumber", "", s); s = re.sub(r"\\\\(\[[^\]]*\])?", "", s); s = s.replace("&", "")
    s = s.replace(r"^{\prime}", "'").replace(r"^\prime", "'")
    s = re.sub(r"\\,|\;|\\!|\s+", "", s)
    s = s.replace(r"\left", "").replace(r"\right", "")
    for b in (r"\Bigg", r"\bigg", r"\Big", r"\big"): s = s.replace(b, "")
    s = re.sub(r"\{\\color\{(blue|red|green)\}", "{", s)
    return s
def show(g1, a, b, g2, oa, ob, ctx=40, maxlen=260):
    A = norm(get(g1, a, b)); B = norm(get(g2, oa, ob))
    sm = difflib.SequenceMatcher(None, A, B, autojunk=False)
    ops = [o for o in sm.get_opcodes() if o[0] != "equal"]
    print(f"##### {g1} {a}-{b} ({len(A)} ch) vs {g2} {oa}-{ob} ({len(B)} ch): {len(ops)} differences")
    for tag, i1, i2, j1, j2 in ops:
        print(f"  [{tag}] ...{A[max(0,i1-ctx):i1]}[[{A[i1:i2][:maxlen]}]] => [[{B[j1:j2][:maxlen]}]]{B[j2:j2+ctx]}...")
    return ops
if __name__ == "__main__":
    args = sys.argv[1:]
    for k in range(0, len(args), 6):
        g1, a, b, g2, oa, ob = args[k:k + 6]; show(g1, int(a), int(b), g2, int(oa), int(ob))
