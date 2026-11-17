# textual cross-check: does the text following each blue log(k+/Lambda) in Evolution_two_rho.tex / Evolution_three_rho.tex
# (bracketed entries) occur verbatim (after whitespace/markup normalisation) in the Group files?
import sys, re; sys.path.insert(0, ".")
from nc import UP, align_blocks
def norm(s):
    s = s.replace(r"^{\prime}", "'")
    s = "\n".join(l for l in s.split("\n") if not l.lstrip().startswith("%"))
    for x in (r"\notag", r"\nonumber", "\\\\", "&", r"\,", r"\;", r"\!", r"\left", r"\right", r"\times", r"\Bigg", r"\bigg", r"\Big", r"\big", " ", "\n", "\t"):
        s = s.replace(x, "")
    s = s.replace(r"\mathbf{x}", "x").replace(r"\mathbf{y}", "y")
    return s
G = {k: norm(open(UP + f).read()) for k, f in (("I", "16e62d79-Group_I.tex"), ("II", "24c35375-Group_II.tex"), ("III", "96b5b029-Group_III.tex"))}
BLUE = re.compile(r"\{\\color\{blue\}\\log\{?\(?\[?\\frac\{k\^\{\+\}\}\{\\Lambda\}\]?\)?\}?\}")
def chunks(fn, lo, hi):
    out = []
    for l, t in align_blocks(fn):
        if not (lo <= l <= hi): continue
        n = norm(t)
        parts = BLUE.split(n)
        for i, p in enumerate(parts[1:], 1):
            # text after the log up to the prefactor of the next blue term (cut at the last '\frac{1}{(2\pi)^{3}}' or 'ig^{4}' start)
            seg = p
            if i < len(parts) - 1:
                cut = max(seg.rfind(r"-\frac{1}{(2\pi)^{3}}"), seg.rfind(r"+\frac{1}{(2\pi)^{3}}"), seg.rfind(r"-\frac{ig^{4}"), seg.rfind(r"+\frac{ig^{4}"))
                if cut > 0: seg = seg[:cut]
            pre = parts[i - 1][-60:]
            out.append((l, i, pre, seg))
    return out
def report(fn, lo, hi, grp, label):
    for l, i, pre, seg in chunks(fn, lo, hi):
        seg0 = seg.rstrip("}")
        found = seg0 in G[grp]
        k = len(seg0)
        if not found:
            # longest prefix found
            a, b = 0, len(seg0)
            while a < b:
                m = (a + b + 1) // 2
                if seg0[:m] in G[grp]: a = m
                else: b = m - 1
            k = a
        print("%-22s line %3d term %d: %s  (%d/%d chars)%s" % (label, l, i, "FOUND verbatim in Group_" + grp if found else "NOT verbatim", k, len(seg0),
              "" if found else "\n      first differing text: ..." + seg0[max(0, k - 40):k] + " || " + seg0[k:k + 70]))
T2 = "e661e1ae-Evolution_two_rho.tex"; T3 = "7583871c-Evolution_three_rho.tex"
if __name__ == "__main__":
    sel = sys.argv[1] if len(sys.argv) > 1 else "2"
    if sel == "2":
        report(T2, 40, 47, "I", "2rho GI RII"); report(T2, 52, 66, "I", "2rho GI RV")
        report(T2, 67, 75, "II", "2rho GII RI"); report(T2, 80, 86, "II", "2rho GII RIV"); report(T2, 92, 105, "III", "2rho GIII RIII-2")
    else:
        report(T3, 39, 47, "I", "3rho GI RII"); report(T3, 48, 61, "I", "3rho GI RIII"); report(T3, 69, 77, "I", "3rho GI RIV"); report(T3, 82, 88, "I", "3rho GI RV")
        report(T3, 97, 103, "II", "3rho GII RI"); report(T3, 115, 122, "III", "3rho GIII RI"); report(T3, 127, 133, "III", "3rho GIII RII")
        report(T3, 138, 146, "III", "3rho GIII RIII-1"); report(T3, 155, 161, "III", "3rho GIII RV"); report(T3, 166, 172, "III", "3rho GIII RVI")
