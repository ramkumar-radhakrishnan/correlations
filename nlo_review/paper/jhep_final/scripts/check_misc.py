#!/usr/bin/env python3
"""Further checks for Sections 2 and 3.
(1) charge algebra on a periodic lattice (su(2) spins): [rho^a(k),rho^b(p)] = i f^{abc} rho^c(k+p);
(2) commutator [N_2, A^(1)], Eq. (commN2A):  [sum_p w(p) rho^b(p) rho^b(-p), rho^a(-k)]
        = - i f^{abc} sum_p w(p) {rho^b(-p), rho^c(-k+p)}   for an even weight w;
(3) decomposition of the two-rho terms, Eqs. (psi1rr)+(psi2rr) -> (psi_anti)+(psi_comm);
(4) cancellation of the two-rho overlaps in the norm of the dressed one-gluon state, Eqs. (ov1)-(ov3);
(5) printed forms (psi1), (psi2) against the xi-forms (psi1_xi), (psi2_xi);
(6) matching factors: (2pi)^{3/2}/(2 pi^{3/2}) = sqrt2, N_2 from A^(1), three-rho part of A^(3).
Run:  python3 -I check_misc.py
"""
import numpy as np, itertools
rng = np.random.default_rng(3)
ok = True
# ---- lattice charges -------------------------------------------------------------------------------
Nsite = 3
sig = [np.array([[0, 1], [1, 0]], complex)/2, np.array([[0, -1j], [1j, 0]])/2, np.array([[1, 0], [0, -1]], complex)/2]
def site_op(a, x):
    ops = [np.eye(2, dtype=complex)]*Nsite; ops[x] = sig[a]
    out = ops[0]
    for o in ops[1:]: out = np.kron(out, o)
    return out
eps = np.zeros((3, 3, 3))
for a, b, c in itertools.permutations(range(3)): eps[a, b, c] = np.linalg.det(np.eye(3)[[a, b, c]])
def rho(a, n): return sum(np.exp(2j*np.pi*n*x/Nsite)*site_op(a, x) for x in range(Nsite))   # rho^a(k_n), f(k)=sum_x e^{ikx} f(x)
comm = lambda A, B: A@B-B@A
acomm = lambda A, B: A@B+B@A
e1 = max(np.abs(comm(rho(a, n), rho(b, m))-1j*sum(eps[a, b, c]*rho(c, n+m) for c in range(3))).max() for a in range(3) for b in range(3) for n in range(3) for m in range(3))
print("(1) charge algebra:", f"{e1:.1e}"); ok &= e1 < 1e-12
w = {0: 0.7, 1: 1.9, 2: 1.9}                        # even weight: w(p)=w(-p)  (modes 1 and 2 = -1)
e2 = 0
for a in range(3):
    for n in range(3):
        lhs = comm(sum(w[p]*rho(b, p)@rho(b, -p) for p in range(3) for b in range(3)), rho(a, -n))
        rhs = -1j*sum(eps[a, b, c]*w[p]*acomm(rho(b, -p), rho(c, -n+p)) for p in range(3) for b in range(3) for c in range(3))
        e2 = max(e2, np.abs(lhs-rhs).max())
print("(2) [N_2, A^(1)] identity:", f"{e2:.1e}", " (norm of the commutator:", f"{np.abs(lhs).max():.2f})"); ok &= e2 < 1e-12
# ---- (3),(4): two-rho terms with matrices for rho^a(k), rho^b(-p) ------------------------------------
a, b, nk, npp = 0, 1, 1, 2
Ra, Rb = rho(a, nk), rho(b, -npp)                   # rho^a(k), rho^b(-p)
kp, pp = 1.3, 0.8; k = rng.normal(size=2); p = rng.normal(size=2); k2, p2 = k@k, p@p
D = pp*k2-kp*p2
t1 = np.sqrt(kp/pp)*(Rb@Ra)/(k2*D)                  # (psi1rr) / (g^2 k^i p^j /4pi^3)
t2 = -np.sqrt(pp/kp)*(Ra@Rb)/(p2*D)                 # (psi2rr)
anti = -0.5/np.sqrt(kp*pp)/(k2*p2)*acomm(Ra, Rb)    # (psi_anti): -g^2/8pi^3 -> -1/2 in these units
com = 0.5/np.sqrt(kp*pp)/(k2*p2)*(kp*p2+pp*k2)/(kp*p2-pp*k2)*comm(Ra, Rb)   # (psi_comm): i f^{abc} rho^c(k-p) = [rho^a(k),rho^b(-p)]
e3 = np.abs(t1+t2-anti-com).max()/np.abs(t1).max()
print("(3) decomposition of the two-rho terms:", f"{e3:.1e}"); ok &= e3 < 1e-12
Rap = rho(b, -npp)                                   # rho^{a'}(-k') ; (ov1)+(ov2, exchange term)+(ov3) in units g^2 k k'/(4 pi^3 sqrt k^2 k'^2)
e4 = np.abs(Rap@Ra+Ra@Rap-acomm(Ra, Rap)).max()
print("(4) two-rho overlaps cancel:", f"{e4:.1e}"); ok &= e4 < 1e-12
# ---- (5) printed forms of the triple-vertex terms -----------------------------------------------------
d2 = np.eye(2); pi = np.pi
def psi1_printed(k, kp, p, pp):
    T = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            for l in range(2):
                br = (2*k[j]-p[j]+(pp-2*kp)/pp*p[j])*d2[i, l]+(2*p[i]-k[i]-(2*pp-kp)/kp*k[i])*d2[j, l]+((pp+kp)/(pp-kp)*(p-k)[l]-k[l]-p[l])*d2[i, j]
                T[i, j] += (p-k)[l]*br
    return 1j/(8*pi**3)*T/(np.sqrt(pp*kp)*(pp-kp)**2)/((k@k/kp-p@p/pp)*((k-p)@(k-p))/(kp-pp))
def psi2_printed(k, kp, p, pp):
    T = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            for l in range(2):
                br = (2*p[i]-k[i]+(kp-2*pp)/kp*k[i])*d2[j, l]+(2*k[j]-p[j]-(2*kp-pp)/pp*p[j])*d2[i, l]+((kp+pp)/(kp-pp)*(k-p)[l]-p[l]-k[l])*d2[i, j]
                T[i, j] += (k-p)[l]*br
    return 1j/(8*pi**3)*T/(np.sqrt(pp*kp)*(kp-pp)**2)/((k@k/kp-((k-p)@(k-p))/(kp-pp)-p@p/pp)*(k@k/kp-p@p/pp))
def bracket(a, b, xi): return np.outer(b, a)-np.outer(a, b)/xi-d2*(a@b)/(1-xi)     # a^m b^l [d_il d_jm - d_jl d_im/xi - d_ij d_lm/(1-xi)]
def psi1_xi(k, kp, p, pp):
    xi = pp/kp; return -1j/(4*pi**3)*np.sqrt(xi)/(1-xi)*bracket(p-k, p-xi*k, xi)/(((p-k)@(p-k))*(p@p-xi*k@k))/kp
def psi2_xi(k, kp, p, pp):
    xi = pp/kp; b_ = p-xi*k; return 1j/(4*pi**3)*xi**1.5/(1-xi)*bracket(k-p, b_, xi)/((p@p-xi*k@k)*(b_@b_))/kp
e5 = max(np.abs(psi1_printed(k, 1.0, p, 1.7)-psi1_xi(k, 1.0, p, 1.7)).max()/np.abs(psi1_xi(k, 1.0, p, 1.7)).max(),
         np.abs(psi2_printed(k, 1.0, p, 0.4)-psi2_xi(k, 1.0, p, 0.4)).max()/np.abs(psi2_xi(k, 1.0, p, 0.4)).max())
print("(5) Eqs.(psi1),(psi2) vs xi-forms:", f"{e5:.1e}"); ok &= e5 < 1e-12
# ---- (6) constants -----------------------------------------------------------------------------------
import sympy as sp
c = [sp.simplify((2*sp.pi)**sp.Rational(3, 2)/(2*sp.pi**sp.Rational(3, 2))-sp.sqrt(2)),          # A^(1) from Psi_LO
     sp.simplify(sp.Rational(1, 2)*2/(2*sp.pi)**3-1/(8*sp.pi**3)),                                 # N_2 = -1/2 int d3p/(2pi)^3 A^dag A, |A|^2 = 2 g^2 rho rho/(p^+ p^2)
     sp.simplify((2*sp.pi)**sp.Rational(3, 2)/(16*sp.pi**sp.Rational(9, 2))-sp.sqrt(2)/(8*sp.pi**3)),   # three-rho part of A^(3)
     sp.simplify((2*sp.pi)**sp.Rational(3, 2)/(32*sp.pi**sp.Rational(9, 2))-1/(8*sp.sqrt(2)*sp.pi**3)), # two-rho part of A^(3)
     sp.simplify((2*sp.pi)**sp.Rational(3, 2)/(32*sp.pi**sp.Rational(7, 2))-1/(8*sp.sqrt(2)*sp.pi**2))] # one-rho part of A^(3)
print("(6) matching constants (all must be 0):", c); ok &= all(x == 0 for x in c)
print("ALL CHECKS PASSED" if ok else "FAILED")
