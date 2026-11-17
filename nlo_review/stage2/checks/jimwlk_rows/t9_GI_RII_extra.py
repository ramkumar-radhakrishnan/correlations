# Group_I.tex Row II, N_2a (l.66-68): common kernel  K^i(x'-w')[K^i(x-w) log((w-y)^2/(x-y)^2) - K^i(y-w) log((w-x)^2/(x-y)^2)]
# multiplies BOTH colour structures.  Evolution_three_rho.tex (l.44-46) keeps for the U(w) structure only the K(w-z).K(.-z) part.
# Difference ("extra") = +1/2 [K(x'-w').K(x-w) - K(x'-w').K(y-w)] K(x-z).K(y-z) FB f^{hbc} U^{ah}(w){rho^b(x),rho^c(y)}  (+cc), units i g^4/16pi^5
import sys; sys.path.insert(0,".")
from t4_brackets import *
from t1_sec2 import BJ, bj
PH=("P0","P1","P2")
BR["x"] = (True, [(+0.5, al, FB + r"[f^{hbc}U^{ah}(w)\{\rho^{b}(x),\rho^{c}(y)\}]"),
                  (-0.5, [("x'-w'", "y-w"), ("y-z", "x-z")], FB + r"[f^{hbc}U^{ah}(w)\{\rho^{b}(x),\rho^{c}(y)\}]")])
ex = build("x")
grp = (True, BR["GI  RII "][1] + BR["x"][1]); BR["g"] = grp
allm=[]
for name, r3, rj in ROWS:
    for l in b3(*r3): allm += list_monos_h(B3[l], True)
lst=[]
for l in b3(192,206): lst += list_monos_h(B3[l], True)
for seed in (71,72):
    W=World(seed)
    a2=[m for m in allm if m.tag=="alpha" and len(m.Us)==2]; ap2=[m for m in allm if m.tag=="alpha'" and len(m.Us)==2]
    print("seed %d: |sym 3-rho coeff of extra term| = %.3g ; all rows + extra: %.3g (all rows alone %.1e)" % (seed, np.abs(canon_sym(ex,W,PH)).max(), np.abs(canon_sym(allm+ex,W,PH)).max(), np.abs(canon_sym(allm,W,PH)).max()))
    D = 1j*two_from_three(ex, W); mj=[]
    for l in bj(289,313): mj += list_monos_h(BJ[l], False)
    A = two_list(mj, W)
    print("         two-rho contact part of the extra term: max %.3g (imag %.1e) ; J Row II list max %.3g" % (np.abs(D).max(), np.abs(D.imag).max(), np.abs(A).max()))
