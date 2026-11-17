# three-rho rows: bracketed forms of Evolution_three_rho.tex (hand-copied operator structure, + c.c. where the row has it)
# versus the "Genuine terms" monomial lists, as ORDERED operator expressions (free algebra, modulo relabelling of x,y,x').
import sys, numpy as np, itertools
sys.path.insert(0, ".")
from nc import *
from t1_sec2 import ROWS, B3, b3
FB = r"[U^{ab'}(x')\rho^{b'}(x')-U^{ab'}(w')\rho^{b'}(x')]"
al = [("x'-w'", "x-w"), ("x-z", "y-z")]      # alpha kernel
be = [("x'-w'", "x-w"), ("w-z", "y-z")]      # beta kernel
# each entry: (coefficient relative to i g^4/(16 pi^5), kernel, operator string)
BR = {}
BR["GI  RII "] = (True, [
    (-0.5, be, FB + r"[f^{hbc}U^{ah}(w)\{\rho^{b}(x),\rho^{c}(y)\}]"),
    (+0.5, [("x'-w'", "y-w"), ("w-z", "x-z")], FB + r"[f^{hbc}U^{ah}(w)\{\rho^{b}(x),\rho^{c}(y)\}]"),
    # log term rewritten with  log((w-y)^2/(x-y)^2) = (1/pi) int_z [K(x-z).K(y-z) - K(w-z).K(y-z)]
    (-0.5, al, FB + r"[f^{abc}U^{bd}(x)U^{ce}(y)\{\rho^{d}(x),\rho^{e}(y)\}]"),
    (+0.5, be, FB + r"[f^{abc}U^{bd}(x)U^{ce}(y)\{\rho^{d}(x),\rho^{e}(y)\}]"),
    (+0.5, [("x'-w'", "y-w"), ("y-z", "x-z")], FB + r"[f^{abc}U^{bd}(x)U^{ce}(y)\{\rho^{d}(x),\rho^{e}(y)\}]"),
    (-0.5, [("x'-w'", "y-w"), ("w-z", "x-z")], FB + r"[f^{abc}U^{bd}(x)U^{ce}(y)\{\rho^{d}(x),\rho^{e}(y)\}]")])
BR["GI  RIII"] = (True, [
    (+1.0, al, r"f^{abd}" + FB + r"[U^{de}(x)U^{bc}(z)\rho^{e}(x)\rho^{c}(y)-U^{de}(x)U^{bc}(y)\rho^{e}(x)\rho^{c}(y)]"),
    (-2.0, [("x'-w'", "y-w"), ("x-z", "w-z")], r"f^{abd}" + FB + r"[U^{de}(x)U^{bc}(z)\rho^{e}(x)\rho^{c}(y)-U^{de}(x)U^{bc}(y)\rho^{e}(x)\rho^{c}(y)]")])
S4 = FB + r"[f^{cdc'}U^{ac}(w)U^{bd}(z)U^{be}(y)\rho^{e}(y)\rho^{c'}(x)-f^{abc}U^{be}(y)U^{cd}(x)\rho^{e}(y)\rho^{d}(x)]"
BR["GI  RIV "] = (True, [(-2.0, be, S4), (+1.0, al, S4)])
BR["GI  RV  "] = (True, [
    (-0.5, [("x'-w'", "y-w"), ("w-z", "y'-z")], FB + r"[f^{adc}U^{db}(z)U^{ce}(w)\{\rho^{e}(y),\rho^{b}(y')\}-f^{abc}U^{cd}(y)U^{be}(y')\{\rho^{d}(y),\rho^{e}(y')\}]"),
    (-0.5, [("x'-w'", "y'-w"), ("z-w", "y-z")], FB + r"[f^{adc}U^{db}(w)U^{ce}(z)\{\rho^{e}(y),\rho^{b}(y')\}-f^{abc}U^{cd}(y)U^{be}(y')\{\rho^{d}(y),\rho^{e}(y')\}]")])
G2a = r"[U^{ac'}(w')\{\rho^{c'}(x'),\rho^{b}(y)\}-U^{bc'}(z)U^{ad'}(x')U^{c'e'}(y)\{\rho^{d'}(x'),\rho^{e'}(y)\}][f^{cbd}U^{ac}(w)\rho^{d}(x)-f^{acd}U^{bc}(z)U^{de}(x)\rho^{e}(x)]"
G2b = r"[f^{c'bd'}U^{ac'}(w')\rho^{d'}(x)-f^{ac'd'}U^{bc'}(z)U^{d'e'}(x)\rho^{e'}(x)][U^{ac}(w)\{\rho^{c}(x'),\rho^{b}(y)\}-U^{bc}(z)U^{ad}(x')U^{ce}(y)\{\rho^{d}(x'),\rho^{e}(y)\}]"
BR["GII RI  "] = (False, [
    (-1.0, [("x-w", "x'-w'"), ("z-w", "y-z")], G2a), (-0.5, [("x-z", "y-z"), ("x-w", "x'-w'")], G2a),
    (+1.0, [("x-w'", "x'-w"), ("z-w'", "y-z")], G2b), (+0.5, [("x-z", "y-z"), ("x-w'", "x'-w")], G2b)])
SB = r"[f^{cbd}U^{ac}(w)\rho^{d}(x)-f^{acd}U^{bc}(z)U^{de}(x)\rho^{e}(x)]"
G31 = r"[U^{ab'}(x')U^{c'd'}(y)U^{c'b}(z)\rho^{b'}(x')\rho^{d'}(y)-U^{bd'}(z)U^{d'e'}(y)U^{ac'}(w')\rho^{c'}(x')\rho^{e'}(y)]" + SB
BR["GIII RI "] = (True, [(+2.0, be, G31), (-1.0, al, G31)])
G32 = r"[U^{c'd'}(y)U^{ab'}(x')U^{bc'}(z)\rho^{d'}(y)\rho^{b'}(x')-U^{ab'}(x')\rho^{b}(y)\rho^{b'}(x')]" + SB
BR["GIII RII"] = (True, [(+2.0, be, G32), (-1.0, al, G32)])
bp = [("x'-w'", "x-w"), ("w'-z", "y-z")]     # beta' kernel
FC = r"f^{c'd'a}[U^{e'c'}(w')U^{d'b}(z)\rho^{e'}(x')-U^{bd'}(z)U^{c'e'}(x')\rho^{e'}(x')]"
BR["GIII RIII"] = (True, [(-1.0, bp, FC + r"[U^{ac}(w)\{\rho^{c}(x),\rho^{b}(y)\}-U^{bc}(z)U^{ad}(x)U^{ce}(y)\{\rho^{d}(x),\rho^{e}(y)\}]")])
BR["GIII RV "] = (True, [(-2.0, bp, r"f^{c'ca}[U^{e'c'}(w')\rho^{e'}(x')-U^{c'e'}(x')\rho^{e'}(x')][U^{ce}(y)U^{ad}(x)\rho^{e}(y)\rho^{d}(x)-U^{ad}(w)U^{ce}(y)\rho^{e}(y)\rho^{d}(x)]")])
BR["GIII RVI"] = (True, [(-2.0, bp, FC + r"[U^{bc}(z)U^{ad}(x)U^{ce}(y)\rho^{d}(x)\rho^{e}(y)-U^{ad}(x)\rho^{d}(x)\rho^{b}(y)]")])

def build(row, flips=()):
    """Monos of the row with prefactor i (so that c.c. is -i * reversed): returns (monos, monos_cc);
       the caller combines  i*M + (-i)*Mcc ; to compare with lists (prefactor i) we use M - Mcc."""
    withcc, ents = BR[row]
    ms = []
    for c, ker, s in ents:
        ms += monos_from(s, coef=c, dots=ker)
    for m in ms:
        m.Us = [((j, i, p) if (i + j, p) in flips else (i, j, p)) for (i, j, p) in m.Us]
    if withcc:
        mc = cc(ms)
        for m in mc: m.coef = -m.coef
        return ms + mc
    return ms
def odiff(A, B, W, PH=("P0", "P1", "P2")):
    DA = canon_ordered(A, W, PH); DB = canon_ordered(B, W, PH)
    return max(np.abs(DA.get(k, 0) - DB.get(k, 0)).max() for k in set(DA) | set(DB)), max(np.abs(v).max() for v in DB.values())
def utokens(row):
    t = set()
    for m in build(row):
        for (i, j, p) in m.Us: t.add((i + j, p))
    return sorted(t)
if __name__ == "__main__":
    W = World(41); W2 = World(42)
    for name, r3, rj in ROWS[:-1]:
        lst = []
        for l in b3(*r3): lst += list_monos_h(B3[l], True)
        d, sz = odiff(build(name), lst, W)
        msg = "%-9s literal bracket form (+cc) vs list, ordered: max|diff| = %.2e (list size %.2g)" % (name, d, sz)
        if d > 1e-9:
            toks = utokens(name); best = None
            for r in range(1, 5):
                for comb in itertools.combinations(toks, r):
                    dd, _ = odiff(build(name, set(comb)), lst, W)
                    if dd < 1e-9: best = comb; break
                if best: break
            if best:
                d2, _ = odiff(build(name, set(best)), lst, W2)
                msg += "\n          -> agrees (%.1e) after transposing in the bracket form: %s" % (d2, ", ".join("U^{%s}(%s)" % t for t in best))
            else:
                msg += "\n          -> no agreement found by transposing up to 4 Wilson lines"
        print(msg)
