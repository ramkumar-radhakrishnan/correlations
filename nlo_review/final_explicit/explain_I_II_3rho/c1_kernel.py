# c1_kernel.py -- checks of the kernel of Group I Row II (three-rho part) for Explain_Row_I_II_3rho.tex
# Run: python3 -I c1_kernel.py  -> c1_kernel.log
# K1  xi integral of Psi^{6b} (source of the first logarithm of R^i), sympy
# K2  xi integral of the remainder of Psi^{4+5} (source of the second logarithm), mpmath quadrature
# K3  assembly of the Lambda-divergent pieces: Psi^{2a} + Psi^{3}  (k^+/Lambda terms cancel)
# K4  full bracket of [A^(3)]_{rho rho} = R_vee^i - log(k+/Lambda) * J^i ; R_vee -> R for vee -> infinity
# K5  prefactor: (LO conjugate) x (two-rho part of A^(3)) = - g^4/(64 pi^6 k^+)
# K6  R^i odd, homogeneous of degree 1; kernel homogeneous of degree -3
# K7  p -> q and p -> 0: each part ~ 2 q^i log(q^2/r^2)/r^2, antisymmetrised kernel ~ log(1/r)/r   (table as in s15)
# K8  large momenta (as in s12) and q -> 0
import sympy as sp, mpmath as mp, numpy as np
out = []
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)
mp.mp.dps = 30

# ---------------- K1
xi, k2, p2, L, V, kp = sp.symbols('xi k2 p2 Lambda vee kplus', positive=True)
F = sp.log(xi / (xi * k2 + p2)) / p2
assert sp.simplify(sp.diff(F, xi) - 1 / (xi * (xi * k2 + p2))) == 0
val = F.subs(xi, V / kp) - F.subs(xi, L / kp)
claim = (sp.log(V / L) + sp.log((L * k2 + kp * p2) / (V * k2 + kp * p2))) / p2
num = [abs(sp.N((val - claim).subs({k2: a, p2: b, L: c, V: d, kp: e}))) for a, b, c, d, e in
       [(1.3, 0.7, 0.01, 50., 2.), (0.2, 3.1, 1e-3, 1e3, 0.5)]]
pr('K1  int_{L/k+}^{V/k+} dxi/(xi (xi k^2 + p^2)) = [log(V/L) + log((L k^2 + k+ p^2)/(V k^2 + k+ p^2))]/p^2 : max dev %.1e' % max(num))

# ---------------- K2
def vec(a, b): return mp.matrix([a, b])
def dot(a, b): return a[0] * b[0] + a[1] * b[1]
def K2check(k, p, kpl, Lam, Vee):
    kk, pp, kdp = dot(k, k), dot(p, p), dot(k, p); d2 = kk - 2 * kdp + pp
    res = []
    for i in range(2):
        f = lambda x: -x / ((x * d2 + (1 - x) * kk) * (1 - x)) * ((2 / x * (pp - kdp) - pp + 2 * kdp) * k[i] + kk * p[i] / x)
        lhs = mp.quad(f, mp.linspace(kpl / Vee, kpl / (kpl + Lam), 12))
        C = (2 * pp - 2 * kdp + kk) * k[i] + kk * p[i]; B = pp * k[i] + kk * p[i]
        rhs = C / d2 * mp.log(1 + (Vee / kpl - 1) * kk / d2) - B / d2 * mp.log((Vee - kpl) / Lam) - k[i] * mp.log(Vee / kpl)
        res.append((lhs, rhs))
    return res
dev = 0
for (k, p, kpl, Lam, Vee) in [(vec(0.7, -0.4), vec(0.3, 0.9), mp.mpf(2), mp.mpf('1e-5'), mp.mpf(300)),
                             (vec(-1.1, 0.2), vec(0.5, 0.6), mp.mpf(1), mp.mpf('1e-6'), mp.mpf(4000))]:
    for lhs, rhs in K2check(k, p, kpl, Lam, Vee):
        dev = max(dev, abs(lhs - rhs) / abs(rhs))
pr('K2  -int_{k+/V}^{k+/(k+ + L)} dxi xi/((xi (k-p)^2 + (1-xi) k^2)(1-xi)) [...]  =  C^i/(k-p)^2 log[1+(V/k+ -1)k^2/(k-p)^2]'
   ' - B^i/(k-p)^2 log((V-k+)/L) - k^i log(V/k+) : max rel dev %.1e  (O(Lambda/k+) corrections; Lambda/k+ = 5e-6, 1e-6)' % dev)

# ---------------- K3 (symbolic, component form)
kx, ky, px, py, lam = sp.symbols('kx ky px py lam', real=True)   # lam = k+/Lambda - 1 - log(Lambda/k+)
lg = sp.Symbol('lg', real=True)                                  # lg = log(Lambda/k+)
K = sp.Matrix([kx, ky]); P = sp.Matrix([px, py]); kk = K.dot(K); pp = P.dot(P); kdp = K.dot(P); d2 = (K - P).dot(K - P)
# units: i g^3 f^{abc}{rho^b(-p),rho^c(-k+p)} / (32 pi^{9/2} p^2 k^2)
t2a = -2 * P * lam                                               # Eq. (grho3), after a <-> b
t3 = -(lg * kdp * K - (lg * kk + 2 * lam * d2) * P) / d2         # Eq. (tworho2)
target = lg * (kk * P - kdp * K) / d2
pr('K3  Psi^{2a} + Psi^{3} = log(Lambda/k+) [k^2 p^i - (k.p) k^i]/(k-p)^2 (k+/Lambda terms cancel):',
   sp.simplify(t2a + t3 - target) == sp.zeros(2, 1))

# ---------------- K4
def Rvec(q, p):
    q2, p2_, qp = q @ q, p @ p, q @ p; d = q - p; dd = d @ d
    return -2 * q * np.log(p2_ / q2) + ((2 * p2_ - 2 * qp + q2) * q + q2 * p) / dd * np.log(q2 / dd)
def bracket(k, p, kpl, Lam, Vee):
    kk, pp, kdp = k @ k, p @ p, k @ p; dd = kk - 2 * kdp + pp
    A = (kk * p - kdp * k) / dd; B = pp * k + kk * p; C = (2 * pp - 2 * kdp + kk) * k + kk * p
    return (np.log(Lam / kpl) * A - 2 * k * np.log(Vee / Lam) - 2 * k * np.log(kpl * pp / (Vee * kk + kpl * pp))
            + C / dd * np.log(1 + (Vee / kpl - 1) * kk / dd) - B / dd * np.log((Vee - kpl) / Lam) - k * np.log(Vee / kpl))
def Rvee(k, p, kpl, Vee):
    kk, pp, kdp = k @ k, p @ p, k @ p; dd = kk - 2 * kdp + pp
    B = pp * k + kk * p; C = (2 * pp - 2 * kdp + kk) * k + kk * p; e = kpl / Vee
    return -2 * k * np.log(pp / (kk + e * pp)) + C / dd * np.log(kk / dd + e * (1 - kk / dd)) - B / dd * np.log(1 - e)
def J(k, p):
    kk, pp, kdp = k @ k, p @ p, k @ p; dd = kk - 2 * kdp + pp
    return (kk * p - kdp * k) / dd + 2 * k + (pp * k + kk * p) / dd
rng = np.random.default_rng(11)
dev = 0; dev2 = []
for _ in range(200):
    k, p = rng.normal(size=2), rng.normal(size=2); kpl = np.exp(rng.normal()); Lam = kpl * 10 ** rng.uniform(-6, -1); Vee = kpl * 10 ** rng.uniform(1, 6)
    lhs = bracket(k, p, kpl, Lam, Vee); rhs = Rvee(k, p, kpl, Vee) - np.log(kpl / Lam) * J(k, p)
    dev = max(dev, np.linalg.norm(lhs - rhs) / np.linalg.norm(lhs))
pr('K4a bracket of [A^(3)]_{rho rho} = R_vee^i - log(k+/Lambda) [ (k^2 p^i - k.p k^i)/(k-p)^2 + 2 k^i + (p^2 k^i + k^2 p^i)/(k-p)^2 ]'
   ' (no log(vee/k+) left): max rel dev %.1e (200 random points)' % dev)
k, p = np.array([0.7, -0.4]), np.array([0.3, 0.9])
for Vee in [1e2, 1e4, 1e6, 1e8]:
    dev2.append(np.linalg.norm(Rvee(k, p, 1.0, Vee) - Rvec(k, p)) / np.linalg.norm(Rvec(k, p)))
pr('K4b |R_vee - R|/|R| at k+/vee = 1e-2, 1e-4, 1e-6, 1e-8:', ' '.join('%.1e' % d for d in dev2))

# ---------------- K5
g, kplus = sp.symbols('g kplus', positive=True)
LOc = sp.I * g / (sp.sqrt(2) * sp.pi * sp.sqrt(kplus))               # conjugate LO bracket: (i g/(sqrt2 pi sqrt k+)) K^i(x'-w') [U(x')-U(w')] rho(x')
A3 = sp.I * g ** 3 / (32 * sp.sqrt(2 * kplus) * sp.pi ** 5)          # two-rho part of A^(3) (mixed representation)
pr('K5  prefactor (LO^dagger) x ([A^(3)]_{rho rho}) =', sp.simplify(LOc * A3), ' (expected -g^4/(64 pi^6 k+))')

# ---------------- K6
def Kq(q, p):
    d = q - p
    return (Rvec(q, p) / (p @ p) - Rvec(q, d) / (d @ d)) / (2 * (q @ q))
U = [(rng.normal(size=2), rng.normal(size=2)) for _ in range(50)]
odd = max(np.linalg.norm(Rvec(-a, -b) + Rvec(a, b)) / np.linalg.norm(Rvec(a, b)) for a, b in U)
hom1 = max(np.linalg.norm(Rvec(3.7 * a, 3.7 * b) - 3.7 * Rvec(a, b)) / np.linalg.norm(Rvec(a, b)) for a, b in U)
hom3 = max(np.linalg.norm(Kq(3.7 * a, 3.7 * b) - 3.7 ** -3 * Kq(a, b)) / np.linalg.norm(Kq(a, b)) for a, b in U)
anti = max(np.linalg.norm(Kq(a, a - b) + Kq(a, b)) / np.linalg.norm(Kq(a, b)) for a, b in U)
pr('K6  R(-q,-p) = -R(q,p): %.1e ; R homogeneous of degree 1: %.1e ; kernel (1/q^2)(1/2)[R(q,p)/p^2 - R(q,q-p)/(q-p)^2]'
   ' homogeneous of degree -3: %.1e ; kernel odd under p -> q-p: %.1e' % (odd, hom1, hom3, anti))

# ---------------- K7 (mpmath, same point q as scripts/s15_N2b_p_to_q.py)
mp.mp.dps = 40
def Rm(q, p):
    q2 = q[0] ** 2 + q[1] ** 2; pq = p[0] ** 2 + p[1] ** 2; qp = q[0] * p[0] + q[1] * p[1]
    d = (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2; c = 2 * pq - 2 * qp + q2
    return [-2 * q[i] * mp.log(pq / q2) + (c * q[i] + q2 * p[i]) / d * mp.log(q2 / d) for i in range(2)]
n2 = lambda v: v[0] ** 2 + v[1] ** 2
q = [mp.mpf('0.7'), mp.mpf('-0.4')]; q2 = n2(q)
def terms(p):
    qmp = [q[0] - p[0], q[1] - p[1]]
    A = [x / n2(p) for x in Rm(q, p)]; B = [x / n2(qmp) for x in Rm(q, qmp)]
    return A, B, [(A[i] - B[i]) / 2 for i in range(2)]
for name, centre in [('p -> q (r = |p - q|)', q), ('p -> 0 (r = |p|)', [mp.mpf(0), mp.mpf(0)])]:
    pr('K7  ' + name + ':   r | max|A_x| r^2 | max|A_x| r^2/(2|q_x| log(q^2/r^2)) | max|B_x| r^2/(2|q_x| log(q^2/r^2)) | max|F_x| r | max|F_x| r/log(1/r)')
    for rs in ['1e-2', '1e-4', '1e-6', '1e-8']:
        r = mp.mpf(rs); mA = mB = mF = 0; N = 64
        for j in range(N):
            th = 2 * mp.pi * j / N; p = [centre[0] + r * mp.cos(th), centre[1] + r * mp.sin(th)]
            A, B, Fk = terms(p); mA = max(mA, abs(A[0]) * r ** 2); mB = max(mB, abs(B[0]) * r ** 2); mF = max(mF, abs(Fk[0]) * r)
        lead = 2 * abs(q[0]) * mp.log(q2 / r ** 2)
        pr('     %s  %s  %s  %s  %s  %s' % (rs, mp.nstr(mA, 6), mp.nstr(mA / lead, 6), mp.nstr(mB / lead, 6), mp.nstr(mF, 6), mp.nstr(mF / mp.log(1 / r), 6)))
pr('     (A = R(q,p)/p^2, B = R(q,q-p)/(q-p)^2, F = (A-B)/2;  q = (0.7,-0.4))')

# ---------------- K8
def slope(fun, Ls):
    v = [max(np.linalg.norm(fun(Lx, u)) for u in U[:20]) for Lx in Ls]
    return [np.log(v[i] / v[i - 1]) / np.log(Ls[i] / Ls[i - 1]) for i in range(1, len(Ls))], v
Ls = [1e2, 1e3, 1e4, 1e5, 1e6]
sl, v = slope(lambda Lx, u: Kq(u[0], Lx * u[1]), Ls)
pr('K8a |p| -> inf, q fixed (x -> y): slopes', ' '.join('%.3f' % s for s in sl), '; |K| p^3/log p =', ' '.join('%.2f' % (vv * Lx ** 3 / np.log(Lx)) for vv, Lx in zip(v, Ls)))
sl, v = slope(lambda Lx, u: Kq(Lx * u[0] + u[1], Lx * u[0]), Ls)
pr('K8b |q|,|p| -> inf, q-p fixed (x -> w): slopes', ' '.join('%.3f' % s for s in sl), '; |K| q^2/log q =', ' '.join('%.2f' % (vv * Lx ** 2 / np.log(Lx)) for vv, Lx in zip(v, Ls)))
sl, v = slope(lambda Lx, u: Kq(Lx * u[0], u[1]), Ls)
pr('K8c |q| -> inf, p fixed (y -> w): slopes', ' '.join('%.3f' % s for s in sl), '; |K| q^2/log q =', ' '.join('%.2f' % (vv * Lx ** 2 / np.log(Lx)) for vv, Lx in zip(v, Ls)))
es = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
sl, v = slope(lambda e, u: Kq(e * u[0], u[1]), es)
pr('K8d q -> 0, p fixed: slopes d log|K|/d log|q| =', ' '.join('%.3f' % s for s in sl), '; |K| q =', ' '.join('%.3f' % (vv * e) for vv, e in zip(v, es)))
sl, v = slope(lambda e, u: Rvec(e * u[0], u[1]) / (u[1] @ u[1]) / (e * e * (u[0] @ u[0])), es)
pr('K8e same for the un-antisymmetrised part R(q,p)/(p^2 q^2): slopes', ' '.join('%.3f' % s for s in sl))
open(__file__.replace('.py', '.log'), 'w').write('\n'.join(out) + '\n')
