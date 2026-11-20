# check the three reordering identities with random matrices (operator level, generic)
import numpy as np, itertools
rng=np.random.default_rng(1)
def R(): return rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
A,B,C=R(),R(),R()
cm=lambda X,Y:X@Y-Y@X
S2=lambda X,Y:(X@Y+Y@X)/2
def S3(X,Y,Z): return sum(np.linalg.multi_dot(p) for p in itertools.permutations([X,Y,Z]))/6
ord_=S3(A,B,C)+0.5*(S2(cm(A,B),C)+S2(cm(A,C),B)+S2(cm(B,C),A))+(cm(cm(A,B),C)-cm(cm(B,C),A))/6
print('ord', abs(A@B@C-ord_).max())
a23=2*S3(A,B,C)+S2(cm(A,B),C)+S2(cm(A,C),B)+(cm(cm(A,B),C)+cm(cm(A,C),B))/6
print('A23', abs(A@(B@C+C@B)-a23).max())
a12=2*S3(A,B,C)+S2(cm(A,C),B)+S2(cm(B,C),A)-(cm(cm(B,C),A)+cm(cm(A,C),B))/6
print('A12', abs((A@B+B@A)@C-a12).max())
