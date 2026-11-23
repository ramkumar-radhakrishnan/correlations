# Check J: Section 4 of DIS_NLO.pdf (momentum space, elastic part) and the limits used in the paper.
#   (1) changes of variables (4.3)-(4.7) and the shifted arguments (4.14)                      [symbolic]
#   (2) f1 - f2 identity (4.10) and the energy-denominator identity (1.9)                      [symbolic]
#   (3) helicity sums (4.19), (4.21), (4.23) from the amplitudes (4.17), random complex T      [exact rationals]
#   (4) quasi-collinear factorisation of the after x after term -> massive splitting function  [symbolic]
#   (5) soft-gluon weights: one common LO factor in the three channels                         [symbolic]
import sympy as sp, itertools, random
I = sp.I; Rt = sp.Rational
ok = lambda b: "OK" if b else "FAIL"

# ---------------------------------------------------------------- (1) kinematics
z, xi = sp.symbols('zeta xi', positive=True); eta = 1 - z - xi
def v2(n): return sp.Matrix(sp.symbols(f'{n}_1 {n}_2', real=True))
c, R, Y, Rq, X, k1, k2, k3, l1, l2 = [v2(n) for n in ('c', 'R', 'Y', 'Rq', 'X', 'k1', 'k2', 'k3', 'l1', 'l2')]
# antiquark channel (4.3)
x_ = c + (1 - z)*R; y_ = c - z*R + xi/(1 - z)*Y; z_ = c - z*R - eta/(1 - z)*Y; yp = c - z*R
kt = (1 - z)*k1 - z*(k2 + k3); pt = (eta*k3 - xi*k2)/(1 - z)
chk = [sp.simplify(z*x_ + eta*y_ + xi*z_ - c), sp.simplify(y_ - z_ - Y), sp.simplify((eta*y_ + xi*z_)/(1 - z) - yp),
       sp.simplify(x_ - yp - R), sp.simplify(x_ - z_ - (R + eta/(1 - z)*Y)),
       sp.Matrix([sp.simplify(k1.dot(x_) + k2.dot(y_) + k3.dot(z_) - ((k1 + k2 + k3).dot(c) + kt.dot(R) - pt.dot(Y)))])]
print("(1a) antiquark variables (4.3), (4.6), (4.7):", ok(all(m == sp.zeros(*m.shape) for m in chk)))
# quark channel (4.4)
yq = c - (1 - eta)*Rq; xq = c + eta*Rq + xi/(1 - eta)*X; zq = c + eta*Rq - z/(1 - eta)*X; xp = c + eta*Rq
ktq = eta*(k1 + k3) - (1 - eta)*k2; ptq = (z*k3 - xi*k1)/(1 - eta)
chk = [sp.simplify(z*xq + eta*yq + xi*zq - c), sp.simplify(xq - zq - X), sp.simplify((z*xq + xi*zq)/(z + xi) - xp),
       sp.simplify(xp - yq - Rq), sp.simplify(zq - yq - (Rq - z/(1 - eta)*X)),
       sp.Matrix([sp.simplify(k1.dot(xq) + k2.dot(yq) + k3.dot(zq) - ((k1 + k2 + k3).dot(c) + ktq.dot(Rq) - ptq.dot(X)))])]
print("(1b) quark variables (4.4)-(4.7):           ", ok(all(m == sp.zeros(*m.shape) for m in chk)))
# shifted arguments (4.14) on the support of delta(k1+k2+k3)
sub = {k2[0]: -k1[0] - k3[0], k2[1]: -k1[1] - k3[1]}
kap = (kt - l1).subs(sub); rho = (pt + eta/(1 - z)*l1 - l2).subs(sub)
kapq = (ktq - l2).subs(sub); rhoq = (ptq + l1 - z/(1 - eta)*l2).subs(sub)
K2 = -k1 - k3
chk = [sp.simplify(kap - (k1 - l1)), sp.simplify(rho - ((k3 + l1 - l2) + xi/(1 - z)*(k1 - l1))),
       sp.simplify(kapq + (K2 + l2)), sp.simplify(rhoq - ((k3 + l1 - l2) + xi/(1 - eta)*(K2 + l2)))]
print("(1c) shifted momenta (4.14):                ", ok(all(m == sp.zeros(2, 1) for m in chk)))

# ---------------------------------------------------------------- (2) scalar identities
m, Q, kk, pp = sp.symbols('m Q k2 p2', positive=True)
eps2 = z*(1 - z)*Q**2 + m**2; om = xi*m/(1 - z); E2 = Q**2 + m**2*(1/z + 1/eta)
den3 = z*(1 - z)**2*(pp + om**2) + xi*eta*(kk + eps2)
f1 = (1 - z)/((kk + eps2)*den3); f2 = 1/(z*(1 - z)*(kk + eps2)*(pp + om**2))
print("(2a) f1 - f2 of (4.10):                     ", ok(sp.simplify(f1 - f2 + xi*eta/(z*(1 - z)*(pp + om**2)*den3)) == 0))
print("(2b) (1.9): zeta zb^2 om^2 + xi eta eps^2 = zeta zb xi eta E_f^2:",
      ok(sp.simplify(z*(1 - z)**2*om**2 + xi*eta*eps2 - z*(1 - z)*xi*eta*E2) == 0))

# ---------------------------------------------------------------- (3) helicity sums in momentum space
def e2(a, b): return [[0, 1], [-1, 0]][a][b]
def d(a, b): return 1 if a == b else 0
def phi(i, j, a, b, zq): return d(a, b)*((2*zq - 1)*d(i, j) + 2*I*e2(i, j)*a)
def tau(mm, n, a, b, x_, par): return d(a, b)*((2*par + x_)*d(mm, n) + 2*I*x_*e2(mm, n)*a)
def mu(p, a, b): return d(b, -a)*(2*a*d(p, 0) - I*d(p, 1))
H = (Rt(1, 2), Rt(-1, 2))
random.seed(11)
def rc(): return Rt(random.randint(-9, 9), random.randint(1, 5)) + I*Rt(random.randint(-9, 9), random.randint(1, 5))
def conj(a): return sp.conjugate(a)
worst = 0
for (zv, xv, mv) in [(Rt(1, 3), Rt(1, 5), Rt(2, 3)), (Rt(2, 7), Rt(3, 11), Rt(5, 4)), (Rt(3, 5), Rt(1, 10), Rt(0))]:
    ev = 1 - zv - xv; zb = 1 - zv; omv = xv*mv/zb; omq = xv*mv/(1 - ev)
    T = [[rc(), rc()], [rc(), rc()]]; Tk = [rc(), rc()]; Tp = [rc(), rc()]; T0 = rc(); TI = rc()
    U = [[rc(), rc()], [rc(), rc()]]; Uk = [rc(), rc()]; Up = [rc(), rc()]; U0 = rc(); UI = rc()
    def A(a1, a2, l, r):      # (4.17), antiquark emitter
        s = 0
        for lam in H:
            s += sum(phi(l, j, a1, lam, zv)*tau(r, n, lam, a2, xv, ev)*T[j][n] for j in (0, 1) for n in (0, 1))
            s += xv*omv*sum(phi(l, j, a1, lam, zv)*Tk[j] for j in (0, 1))*mu(r, lam, a2)
            s += mv*mu(l, a1, lam)*sum(tau(r, n, lam, a2, xv, ev)*Tp[n] for n in (0, 1))
            s += mv*xv*omv*mu(l, a1, lam)*mu(r, lam, a2)*T0
        return s - xv*ev*(d(l, r) + 2*I*e2(l, r)*a1)*d(a1, a2)*TI
    def At(a1, a2, l, r):     # (4.17), quark emitter
        s = 0
        for lam in H:
            s += sum(tau(n, r, a1, lam, xv, zv)*phi(l, j, lam, a2, zv + xv)*U[j][n] for j in (0, 1) for n in (0, 1))
            s += xv*omq*mu(r, a1, lam)*sum(phi(l, j, lam, a2, zv + xv)*Uk[j] for j in (0, 1))
            s += mv*sum(tau(n, r, a1, lam, xv, zv)*Up[n] for n in (0, 1))*mu(l, lam, a2)
            s += mv*xv*omq*mu(r, a1, lam)*mu(l, lam, a2)*U0
        return s - xv*zv*(d(l, r) - 2*I*e2(l, r)*a1)*d(a1, a2)*UI
    def hs(F, G): return sp.expand(sum(F(*h)*conj(G(*h)) for h in itertools.product(H, H, (0, 1), (0, 1))))
    n2 = lambda M_, N_: sum(M_[j][n]*conj(N_[j][n]) for j in (0, 1) for n in (0, 1))
    tr = lambda M_: M_[0][0] + M_[1][1]
    dot = lambda a, b: a[0]*conj(b[0]) + a[1]*conj(b[1])
    ee = lambda M_, N_: M_[0][0]*conj(N_[1][1]) + M_[1][1]*conj(N_[0][0]) - M_[0][1]*conj(N_[1][0]) - M_[1][0]*conj(N_[0][1])
    re = sp.re
    # (4.19)
    Hml = 8*(zv**2 + zb**2)*(ev**2 + zb**2)*n2(T, T) - 8*xv*(2*zv - 1)*(2*zb - xv)*ee(T, T) \
          - 16*zv*xv*ev**2*re(tr(T)*conj(TI)) + 8*xv**2*ev**2*TI*conj(TI)
    Hmv = 8*mv**2*(ev**2 + zb**2)*dot(Tp, Tp) + 8*xv**2*omv**2*(zv**2 + zb**2)*dot(Tk, Tk) + 8*mv**2*xv**2*omv**2*T0*conj(T0) \
          - 16*zv*ev*mv*xv*omv*re(tr(T)*conj(T0) - dot(Tk, Tp)) + 16*xv*ev*mv*xv*omv*re(T0*conj(TI))
    d1 = sp.simplify(hs(A, A) - Hml - Hmv)
    # (4.21)
    Hml = 8*(zv**2 + (1 - ev)**2)*(ev**2 + (1 - ev)**2)*n2(U, U) - 8*xv*(2*ev - 1)*(2*(1 - ev) - xv)*ee(U, U) \
          + 16*zv**2*xv*ev*re(tr(U)*conj(UI)) + 8*xv**2*zv**2*UI*conj(UI)
    Hmv = 8*mv**2*(zv**2 + (1 - ev)**2)*dot(Up, Up) + 8*xv**2*omq**2*(ev**2 + (1 - ev)**2)*dot(Uk, Uk) + 8*mv**2*xv**2*omq**2*U0*conj(U0) \
          + 16*zv*ev*mv*xv*omq*re(tr(U)*conj(U0) - dot(Uk, Up)) + 16*xv*zv*mv*xv*omq*re(U0*conj(UI))
    d2 = sp.simplify(hs(At, At) - Hml - Hmv)
    # (4.23)
    Hml = 8*(1 - xv - 2*zv*ev)*(xv*(1 - xv) + 2*zv*ev)*n2(U, T) - 8*xv*(1 - 2*zv - xv)**2*ee(U, T) \
          - 8*xv*ev*(zv + xv)**2*tr(U)*conj(TI) + 8*xv*zv*zb**2*UI*conj(tr(T))
    Hmv = 8*mv**2*(xv + 2*zv*ev)*dot(Up, Tp) - 8*(xv + 2*zv*ev)*xv**2*omq*omv*dot(Uk, Tk) \
          - 8*(zv + xv)**2*mv*xv*omv*tr(U)*conj(T0) + 8*zb**2*mv*xv*omq*U0*conj(tr(T)) \
          - 8*ev**2*mv*xv*omq*dot(Uk, Tp) + 8*zv**2*mv*xv*omv*dot(Up, Tk)
    d3 = sp.simplify(hs(At, A) - Hml - Hmv)
    print(f"(3) zeta={zv}, xi={xv}, m={mv}:  (4.19) {ok(d1 == 0)}   (4.21) {ok(d2 == 0)}   (4.23) {ok(d3 == 0)}")

# ---------------------------------------------------------------- (4) quasi-collinear limit (after x after, antiquark)
zg = sp.symbols('z_g', positive=True)           # gluon splitting fraction xi/(1-zeta)
k_, p_ = v2('k'), v2('p')
def Aaa(a1, a2, l, r, zv, xv):
    ev = 1 - zv - xv; omv = xv*m/(1 - zv)
    s = 0
    for lam in H:
        s += (sum(phi(l, j, a1, lam, zv)*k_[j] for j in (0, 1)) + m*mu(l, a1, lam)) * \
             (sum(tau(r, n, lam, a2, xv, ev)*p_[n] for n in (0, 1)) + xv*omv*mu(r, lam, a2))
    return s
zs = sp.symbols('zeta', positive=True); xs = zg*(1 - zs)
tot = sp.expand(sum(Aaa(*h, zs, xs)*sp.conjugate(Aaa(*h, zs, xs)) for h in itertools.product(H, H, (0, 1), (0, 1))))
k2s, p2s = k_.dot(k_), p_.dot(p_)
target = 8*((zs**2 + (1 - zs)**2)*k2s + m**2)*(1 - zs)**2*((1 + (1 - zg)**2)*p2s + zg**4*m**2)
print("(4a) sum|V_gam V_g|^2 = 8 [(zeta^2+zb^2) k^2 + m^2] zb^2 [(1+(1-z)^2) p^2 + z^4 m^2]:", ok(sp.simplify(tot - target) == 0))
lhs = ((1 + (1 - zg)**2)*p2s + zg**4*m**2)/(zg*(p2s + zg**2*m**2)**2)
rhs = ((1 + (1 - zg)**2)/zg - 2*zg*(1 - zg)*m**2/(p2s + zg**2*m**2))/(p2s + zg**2*m**2)
print("(4b) = massive splitting function [ (1+(1-z)^2)/z - 2 z (1-z) m^2/(p^2+z^2 m^2) ] / (p^2+z^2 m^2):", ok(sp.simplify(lhs - rhs) == 0))

# ---------------------------------------------------------------- (5) soft-gluon weights
zb = 1 - z
cq_b = sp.limit(8*(z**2 + zb**2)*(eta**2 + zb**2), xi, 0), sp.limit(8*(eta**2 + zb**2), xi, 0)
cq = sp.limit(8*(z**2 + (1 - eta)**2)*(eta**2 + (1 - eta)**2), xi, 0), sp.limit(8*(z**2 + (1 - eta)**2), xi, 0)
cx = sp.limit(8*(1 - xi - 2*z*eta)*(xi*(1 - xi) + 2*z*eta), xi, 0), sp.limit(8*(xi + 2*z*eta), xi, 0)
P = z**2 + zb**2
# explicit prefactors zeta^2, eta^2, zeta*eta and 1/(zeta zb)^2 from the two radial functions
w = [sp.simplify(pref*c_[0]/(z*zb)**2/(16*P)) for pref, c_ in ((z**2, cq_b), (zb**2, cq), (z*zb, cx))]
wm = [sp.simplify(pref*c_[1]/(z*zb)**2/16) for pref, c_ in ((z**2, cq_b), (zb**2, cq), (z*zb, cx))]
print("(5) soft limit, weights of [zeta^2+zb^2] eps^2 K1 K1 and of m^2 K0 K0 (antiquark, quark, interference):", w, wm,
      ok(w == [1, 1, 1] and wm == [1, 1, 1]))
