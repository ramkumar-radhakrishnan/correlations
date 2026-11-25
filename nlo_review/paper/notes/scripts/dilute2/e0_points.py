#!/usr/bin/env python3
"""dilute2 / e0: kinematic points for the F_P(q_a)+F_T(q_b) test.
Grid {q_a1,q_a2,q_a3} x {q_b1,q_b2}; each corner (q_a,q_b) has k = q_a+q_b and is mapped by a rotation and a
rescaling to the frame of the ref3_dilute scripts, k=(1,0), l=-q_b (so q=k+l=q_a), mu=|k|=1."""
import numpy as np
QA = {'a1': np.array([1.3, 0.7]), 'a2': np.array([0.5, 0.4]), 'a3': np.array([-0.4, 0.9])}
QB = {'b1': np.array([-0.3, -0.7]), 'b2': np.array([0.5, -0.4])}
def corner(a, b):
    qa, qb = QA[a], QB[b]; k = qa + qb; kn = np.linalg.norm(k)
    c, s = k / kn; R = np.array([[c, s], [-s, c]])       # rotates k to (kn,0)
    return dict(qa=qa, qb=qb, k=k, kn=kn, l=-(R @ qb) / kn, q=(R @ qa) / kn)
PTS = {a + b: corner(a, b) for a in QA for b in QB}
# additional points with |q_a| = |q_b| (k = (1,0), q_a = (1/2, t), q_b = (1/2, -t)); a2b2 is of this type (t = 0.4)
for _n, _t in (('s1', 0.8660254037844386), ('s2', 0.1339745962155614), ('s3', 1.8660254037844386)):
    PTS[_n] = dict(qa=np.array([0.5, _t]), qb=np.array([0.5, -_t]), k=np.array([1.0, 0.0]), kn=1.0, l=np.array([-0.5, _t]), q=np.array([0.5, _t]))
if __name__ == '__main__':
    for n, P in PTS.items():
        print(n, 'k=', P['k'], '|k|=%.6f' % P['kn'], 'l_scaled= %.12f %.12f' % tuple(P['l']), 'q_scaled=', P['q'])
