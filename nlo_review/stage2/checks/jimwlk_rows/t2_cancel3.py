# (3) cancellation of the symmetrised (commuting) three-rho log terms
#   A: from the monomial lists ("From four rho's" + "Genuine terms" of all rows), grouped by tag and number of Wilson lines
#   B: the 12 "= 0" blocks of the section "Cancellation of all three rho's" as written (kernel and f^{abc} factored out)
import sys, numpy as np
sys.path.insert(0, ".")
from nc import *
from t1_sec2 import ROWS, B3, b3
PH = ("P0", "P1", "P2")
worlds = [World(s) for s in (21, 22)]
allm = []; perrow = {}
for name, r3, rj in ROWS:
    ms = []
    for l in b3(*r3): ms += list_monos_h(B3[l], True)
    perrow[name] = ms; allm += ms
print("A) total number of ordered three-rho monomials:", len(allm), " untagged:", sum(1 for m in allm if m.tag is None))
for W in worlds[:1]:
    print("   per row |sym coefficient| :", "  ".join("%s %.2f" % (n.strip(), np.abs(canon_sym(ms, W, PH)).max()) for n, ms in perrow.items()))
tags = ["alpha", "alpha'", "beta", "beta'"]
for tg in tags:
    for nU in (2, 3, 4):
        sel = [m for m in allm if m.tag == tg and len(m.Us) == nU]
        if not sel: print("   %-7s %dU: no terms" % (tg, nU)); continue
        v = [np.abs(canon_sym(sel, W, PH)).max() for W in worlds]
        sz = max(np.abs(canon_sym([m], worlds[0], PH)).max() for m in sel)
        print("   %-7s %dU: %3d monomials, max|sum of symmetric coefficient| = %.2e / %.2e   (single-term size %.2g)" % (tg, nU, len(sel), v[0], v[1], sz))
v = [np.abs(canon_sym(allm, W, PH)).max() for W in worlds]
print("   ALL: max|sum| = %.2e / %.2e" % tuple(v))

# B) the cancellation blocks as written
cb = [(l, t) for l, t in align_blocks("7583871c-Evolution_three_rho.tex") if l >= 254]
def split_terms(body):
    """terms are top-level {\\color{c} ... } groups"""
    out = []; i = 0
    while True:
        m = re.search(r"\{\\color\{\w+\}", body[i:])
        if not m: break
        st = i + m.start(); j = i + m.end(); depth = 1
        while depth:
            c = body[j]
            if c == '{' and body[j - 1] != '\\': depth += 1
            elif c == '}' and body[j - 1] != '\\': depth -= 1
            j += 1
        out.append(body[st:j]); i = j
    return out
print("B) cancellation blocks as written (f^{abc} and the header kernel restored):")
allc = []
for l, t in cb:
    a = t.index(r"\Bigg[") + 6; b = t.rindex(r"\Bigg]")
    head = t[:a]; body = t[a:b]
    hk = re.findall(r"\\mathcal\{K\}\(([^()]*)\)", head.split("e^{-ik(w'-w)}")[1])
    rest = body
    terms = split_terms(body)
    for x in terms: rest = rest.replace(x, "")
    rest = rest.replace(r"\notag", "").replace("\\\\", "").replace("&", "").replace(r"\,", "").strip()
    assert rest == "", ("unparsed text in block", l, rest)
    ms = []
    for x in terms:
        ms += parse_mono(x + r" f^{abc} " + r"\mathcal{K}(%s)\cdot\mathcal{K}(%s)\mathcal{K}(%s)\cdot\mathcal{K}(%s)" % tuple(hk))
    allc.append((l, hk, ms))
    v = [np.abs(canon_sym(ms, W, PH)).max() for W in worlds]
    sz = max(np.abs(canon_sym([m], worlds[0], PH)).max() for m in ms)
    print("   line %d  header kernel K(%s).K(%s) K(%s).K(%s): %2d terms (%d ordered monomials), nU=%s : max|sym sum| = %.2e / %.2e (term size %.2g)" % (
        l, *hk, len(terms), len(ms), sorted(set(len(m.Us) for m in ms)), v[0], v[1], sz))
import pickle
pickle.dump([(l, hk, len(ms)) for l, hk, ms in allc], open("cancel_blocks.pkl", "wb"))
