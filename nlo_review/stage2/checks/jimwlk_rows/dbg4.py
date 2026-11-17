import sys; sys.path.insert(0,".")
from t3_4rho import *
W=World(5); PH=("P0","P1","P2","P3")
for nm,p in (("A",parts[0]),("B",parts[1]),("C",parts[2])):
    print(nm, [ (len(m.Us), m.coef, m.pts(), [u[2] for u in m.Us]) for m in p if len(m.Us)==2])
two=[m for m in M4 if len(m.Us)==2]
print(len(two), np.abs(canon_sym(two,W,PH)).max())
A2=[m for m in parts[0] if len(m.Us)==2]; B2=[m for m in parts[1] if len(m.Us)==2]; C2=[m for m in parts[2] if len(m.Us)==2]
for nm,s in (("A",A2),("Acc",cc(A2)),("B",B2),("C",C2),("Ccc",cc(C2)),("A+Acc",A2+cc(A2)),("C+Ccc",C2+cc(C2)),("B+C+Ccc",B2+C2+cc(C2)) ,("all",A2+cc(A2)+B2+C2+cc(C2))):
    print(nm, np.abs(canon_sym(s,W,PH)).max())
