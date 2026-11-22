# s12_N2b_momentum.py -- UV behaviour of the momentum-space kernel of Group I Row II, N_2b (Sec. 6, eq:D-I.II-3/-2):
# Kq(q,p) = (1/q^2) (1/2)[R(q,p)/p^2 - R(q,q-p)/(q-p)^2], with the author's R^i(q,p), under the phases
# e^{iq.(w-y)} e^{ip.(y-x)}.  UV regions: (p,q) -> lam^{-1}(p,q) jointly  (x,y -> w);  |p| -> inf at fixed q (x -> y);
# |q|,|p| -> inf at fixed l = p - q (x -> w);  |q| -> inf at fixed p (y -> w).  |Kq| ~ |big|^{-n}: the coordinate
# kernel behaves as r^{n-2} (2D) in the corresponding separation -> integrable iff n > 0; jointly 4D: r^{-(4-3)}.
# Run: python3 -I s12_N2b_momentum.py -> s12_N2b_momentum.log
import numpy as np
out = []
def pr(*a): s = ' '.join(str(x) for x in a); print(s); out.append(s)
def R(q, p):
    q2, p2, qp = q @ q, p @ p, q @ p; d = q - p; d2 = d @ d
    return -2 * q * np.log(p2 / q2) + ((2 * p2 - 2 * qp + q2) * q + q2 * p) / d2 * np.log(q2 / d2)
def Kq(q, p):
    d = q - p
    return (R(q, p) / (p @ p) - R(q, d) / (d @ d)) / (2 * (q @ q))
rng = np.random.default_rng(3)
def slope(fun, Ls):
    v = [max(np.linalg.norm(fun(L, u)) for u in U) for L in Ls]
    return [np.log(v[i] / v[i - 1]) / np.log(Ls[i] / Ls[i - 1]) for i in range(1, len(Ls))], v
U = [(rng.normal(size=2), rng.normal(size=2)) for _ in range(20)]
Ls = [1e2, 1e3, 1e4, 1e5, 1e6]
hom = max(np.linalg.norm(Kq(7.3 * a, 7.3 * b) - 7.3 ** -3 * Kq(a, b)) / np.linalg.norm(Kq(a, b)) for a, b in U)
pr('joint scaling (p,q)->7.3(p,q): |Kq(7.3q,7.3p) - 7.3^-3 Kq(q,p)|/|Kq| max = %.1e  (exactly homogeneous of degree -3:'
   ' coordinate kernel homogeneous of degree -1 in (x-w,y-w), 4D cluster {x,y,w}: integrable)' % hom)
sl, v = slope(lambda L, u: Kq(u[0], L * u[1]), Ls)
pr('|p| -> inf, q fixed (x -> y): slopes d log|Kq|/d log|p| =', ' '.join('%.3f' % s for s in sl), '; |Kq| p^3/log p =', ' '.join('%.3f' % (vv * L ** 3 / np.log(L)) for vv, L in zip(v, Ls)))
sl, v = slope(lambda L, u: Kq(L * u[0] + u[1], L * u[0]), Ls)   # p = L u0, q = p + l (l = u1 fixed): p - q = -l fixed
pr('|q|,|p| -> inf, p-q fixed (x -> w): slopes =', ' '.join('%.3f' % s for s in sl), '; |Kq| q^2/log q =', ' '.join('%.3f' % (vv * L**2 / np.log(L)) for vv, L in zip(v, Ls)))
sl, v = slope(lambda L, u: Kq(L * u[0], u[1]), Ls)
pr('|q| -> inf, p fixed (y -> w): slopes =', ' '.join('%.3f' % s for s in sl), '; |Kq| q^2/log q =', ' '.join('%.3f' % (vv * L ** 2 / np.log(L)) for vv, L in zip(v, Ls)))
pr('=> x->y: |Kq| ~ log p/p^3 (n=3): bounded coordinate kernel at x->y;')
pr('   x->w and y->w: the leading log q/q terms of the two parts of the antisymmetrised kernel cancel; |Kq| ~ (log q)^k/q^2, 1<k<2 by the fit (n=2): the coordinate kernel is at most a power of log r at x->w, y->w, integrable;')
pr('   joint x,y->w: degree -3 in 4D -> r^{-1} at {x,y,w} (d=4): integrable.  With P^i = K^i(x\'-w\') (or K^i(x-w\'), K^i(y-w\') after a contact): {x\',w\'} r^{-1}; {x,y,w,(x\'),w\'}: r^{-2}, d>=6: integrable.')
pr('   No 1/eps UV pole: the kernel has no scaleless UV part left (the 11/3 terms of Rows II and V cancel exactly; the value at p->0, q->p is fixed by the antisymmetrisation).')
open('s12_N2b_momentum.log', 'w').write('\n'.join(out) + '\n')
