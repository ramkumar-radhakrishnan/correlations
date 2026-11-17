# colour-rotation test: for U(all points) = V (constant) each row must vanish (it is a product of [U A - Abar] differences).
# For rows without four-rho terms the symmetric three-rho coefficient must vanish by itself at U=V.
import sys; sys.path.insert(0,".")
from t4_brackets import *
PH=("P0","P1","P2")
W=World(51); V=W.U["z"]
Wc=World(51); 
for k in Wc.U: Wc.U[k]=V
extra_GI_RII = (True, [  # the piece of the Group-I Row II N_2a formula that is absent in Evolution_three_rho.tex:
    (+0.5, al, FB + r"[f^{hbc}U^{ah}(w)\{\rho^{b}(x),\rho^{c}(y)\}]"),
    (-0.5, [("x'-w'", "y-w"), ("y-z", "x-z")], FB + r"[f^{hbc}U^{ah}(w)\{\rho^{b}(x),\rho^{c}(y)\}]")])
print("row: max|symmetric 3-rho coefficient| generic U  /  U=const")
for name, r3, rj in ROWS:
    lst=[]
    for l in b3(*r3): lst += list_monos_h(B3[l], True)
    print("  %-26s %.3g   %.2e" % (name, np.abs(canon_sym(lst,W,PH)).max(), np.abs(canon_sym(lst,Wc,PH)).max()))
BR["x"]=extra_GI_RII
ex=build("x")
lst=[]
for l in b3(192,206): lst += list_monos_h(B3[l], True)
print("GI RII list + extra Group-file piece: generic %.3g ; U=const %.2e ; extra piece alone generic %.3g" % (np.abs(canon_sym(lst+ex,W,PH)).max(), np.abs(canon_sym(lst+ex,Wc,PH)).max(), np.abs(canon_sym(ex,W,PH)).max()))
# rows whose U=const value is non-zero: sum of them?
grp=[]
for name, r3, rj in ROWS:
    for l in b3(*r3): grp += list_monos_h(B3[l], True)
print("all rows, U=const: %.2e" % np.abs(canon_sym(grp,Wc,PH)).max())
