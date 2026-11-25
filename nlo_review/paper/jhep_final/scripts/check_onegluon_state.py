#!/usr/bin/env python3
"""Checks for Section 3.2 (single-gluon incoming state): the components of the dressed one-gluon state
are rebuilt from the matrix elements (m1)-(m5) and the second-order formula (pt2), and compared with
the expressions printed in the paper (instantaneous terms, triple-vertex terms incl. overall sign,
xi-forms, coefficient C_1 and B_3, cancellation of the double pole at p^+ -> k^+).
Colour: the coefficient of  f^{abc} rho^c(k-p)  is compared (f^{bca}=f^{abc}, f^{acb}=-f^{abc}).
Run:  python3 -I check_onegluon_state.py
"""
import numpy as np
rng = np.random.default_rng(11)
g = 1.0; pi = np.pi
d2 = np.eye(2)
def E(v, vp): return v@v/(2*vp)

def m3(k, kp, p, pp, q, qp):
    """<g_l(q) g_j(p)|H_ggg|g_i(k)> without i g f^{abc} delta3; array [i,j,l]."""
    T = np.zeros((2, 2, 2))
    A = p - q + (qp-pp)/kp*k; B = k + q - (kp+qp)/pp*p; C = (kp+pp)/qp*q - p - k
    for i in range(2):
        for j in range(2):
            for l in range(2):
                T[i, j, l] = A[i]*d2[j, l] + B[j]*d2[i, l] + C[l]*d2[i, j]
    return T/(8*pi**1.5*np.sqrt(kp*pp*qp))
def emit(r, rp):   # <..g_m(r)|H_g|..> / rho^d(-r) ; vector index m
    return g*r/(4*pi**1.5*rp**1.5)

def psi1(k, kp, p, pp):
    """coefficient of f^{abc} rho^c(k-p) |g_j^b(p)>, p^+>k^+: H_g emits r=p-k, then H_ggg merges (k,r)->p. [i,j]"""
    r, rp = p-k, pp-kp
    # <g_j(p)|H_ggg|g_l(q=k) g_m(r)> = conj of <g_l(k) g_m(r)|H_ggg|g_j(p)> = -i f^{b d a}(...) ; m3 args: parent p, 'p'-slot=r (index m), 'q'-slot=k (index l=i)
    T = m3(p, pp, r, rp, k, kp)          # [j, m, l]
    amp = np.einsum('jmi,m->ij', T, emit(r, rp))
    # colour: (-i) f^{bda} -> (-i) f^{abc};  factor 1/2 * 2 (two equal terms of m2)
    return (-1j)*amp/((E(k, kp)-E(k, kp)-E(r, rp))*(E(k, kp)-E(p, pp)))
def psi2(k, kp, p, pp):
    """p^+<k^+: H_ggg splits k->(p, r=k-p), then H_g absorbs r. [i,j]"""
    r, rp = k-p, kp-pp
    T = m3(k, kp, r, rp, p, pp)          # [i, m, l=j]  colour i f^{a d b} = -i f^{abc}
    amp = np.einsum('imj,m->ij', T, emit(r, rp))   # absorption: g r^m rho^d(+r)/(...)  (hermitian conjugate of emission)
    return (-1j)*amp/((E(k, kp)-E(p, pp)-E(r, rp))*(E(k, kp)-E(p, pp)))
def inst(k, kp, p, pp):
    return 1j*g**2*(kp+pp)*d2/(2*(2*pi)**3*np.sqrt(kp*pp)*(kp-pp)**2)/(E(k, kp)-E(p, pp))

# ---- expressions printed in the paper -------------------------------------------------
def paper_inst(k, kp, p, pp):       # Eq. (psi_inst)
    return 1j*g**2/(2*pi)**3*np.sqrt(pp*kp)*(kp+pp)/(pp-kp)**2/(pp*k@k-kp*p@p)*d2
def bracket(k, p, xi):              # (p-k)^m (p-xi k)^l [d_il d_jm - d_jl d_im/xi - d_ij d_lm/(1-xi)]  -> [i,j]
    a, b = p-k, p-xi*k
    return np.outer(b, a) - np.outer(a, b)/xi - d2*(a@b)/(1-xi)
def paper_psi1_xi(k, kp, p, pp):    # Eq. (psi1_xi): includes dp^+ = k^+ dxi  -> divide by k^+
    xi = pp/kp
    return -1j*g**2/(4*pi**3)*np.sqrt(xi)/(1-xi)*bracket(k, p, xi)/((p-k)@(p-k)*(p@p-xi*k@k))/kp
def paper_psi2_xi(k, kp, p, pp):
    xi = pp/kp
    a, b = k-p, p-xi*k
    br = np.outer(b, a) - np.outer(a, b)/xi - d2*(a@b)/(1-xi)
    return 1j*g**2/(4*pi**3)*xi**1.5/(1-xi)*br/((p@p-xi*k@k)*(b@b))/kp
def paper_B3_third(k, kp, p, pp):   # third term of Eq. (B3mom), coefficient of f^{abc} rho^c(k-p)
    kmp = k-p; P = k*pp-kp*p
    br = np.zeros((2, 2), complex)
    for i in range(2):
        for j in range(2):
            for l in range(2):
                for m in range(2):
                    br[i, j] += kmp[l]*P[m]*(d2[j, m]*d2[i, l]/pp - d2[i, m]*d2[j, l]/kp - d2[l, m]*d2[i, j]/(pp-kp))
    pref = 2j*g**2*np.sqrt(kp*pp)/(kp-pp)**2/(k@k*pp-p@p*kp)
    if pp > kp: last = (pp-kp)/(kmp@kmp)
    else: last = 1/(k@k/kp-(kmp@kmp)/(kp-pp)-p@p/pp)
    return pref*br*last
def paper_B3_inst(k, kp, p, pp):
    return 1j*g**2*np.sqrt(pp*kp)*(kp+pp)/(pp-kp)**2/(pp*k@k-kp*p@p)*d2

ok = True
for trial in range(4):
    k, p = rng.normal(size=2), rng.normal(size=2); kp = 1.0+rng.random()
    for pp in (kp*(1.3+rng.random()), kp*(0.2+0.6*rng.random())):
        full = psi1(k, kp, p, pp) if pp > kp else psi2(k, kp, p, pp)
        xi_form = paper_psi1_xi(k, kp, p, pp) if pp > kp else paper_psi2_xi(k, kp, p, pp)
        e1 = np.abs(full-xi_form).max()/np.abs(full).max()
        e2 = np.abs((2*pi)**3*full-paper_B3_third(k, kp, p, pp)).max()/np.abs(full).max()/(2*pi)**3
        e3 = np.abs(inst(k, kp, p, pp)-paper_inst(k, kp, p, pp)).max()/np.abs(inst(k, kp, p, pp)).max()
        ok &= max(e1, e2, e3) < 1e-12
        print(f"p+/k+={pp/kp:5.2f}: triple-vertex term vs xi-form {e1:.1e}, vs B3 (third term) {e2:.1e}; instantaneous {e3:.1e}")
# double pole at p^+ -> k^+ : (inst + triple) * (p^+-k^+)^2 -> 0
k, p = rng.normal(size=2), rng.normal(size=2); kp = 1.0
print("cancellation of 1/(p+-k+)^2:  eps, |inst+triple|*eps^2, |inst|*eps^2")
for eps in (1e-2, 1e-3, 1e-4):
    for sgn in (+1, -1):
        pp = kp+sgn*eps
        tot = paper_B3_inst(k, kp, p, pp)+paper_B3_third(k, kp, p, pp)
        print(f"  p+-k+={sgn*eps:+.0e}: {np.abs(tot).max()*eps**2:.2e}   {np.abs(paper_B3_inst(k, kp, p, pp)).max()*eps**2:.2e}")

# ---- C_1 from <gg|H_ggg|g>/(2 (E_g(k)-E_gg)) * (2pi)^{3/2}, coefficient of i... f^{abc}  -------------
def C1_paper(k, kp, p, pp):
    P = p*kp-pp*k; T = np.zeros((2, 2, 2), complex)
    for i in range(2):
        for j in range(2):
            for l in range(2):
                for m in range(2):
                    T[i, j, l] += P[m]*(d2[i, m]*d2[j, l]/kp - d2[i, l]*d2[m, j]/pp - d2[l, m]*d2[i, j]/(kp-pp))
    return -1j*g*np.sqrt(kp*pp*(kp-pp))/(np.sqrt(2)*(P@P))*T
k, p = rng.normal(size=2), rng.normal(size=2); kp, pp = 1.7, 0.6
# the LCPT state is  int d3p d3q |g(q) g(p)> M/(2 dE): compare (2pi)^{3/2} M/(2 dE) with C_1
lhs = (2*pi)**1.5*1j*g*m3(k, kp, p, pp, k-p, kp-pp)/(2*(E(k, kp)-E(p, pp)-E(k-p, kp-pp)))
e = np.abs(lhs-C1_paper(k, kp, p, pp)).max()/np.abs(lhs).max(); ok &= e < 1e-12
print("C_1 from the matrix element vs Eq.(C1mom):", f"{e:.1e}")

# ---- a^dag a unitarity relation, triple-vertex part (fixes the sign of the p^+>k^+ term of B_3) -----
# for q^+ > p^+ :  B3(q->p)|_triple + [B3(p->q)|_triple]^dag + 2 A2^dag(q-p) C1(q,p) = 0   (coefficient of f^{cbe} rho^e(q-p), indices [k (in q), j (out p)])
q, p = rng.normal(size=2), rng.normal(size=2); qp, pp = 1.9, 0.7
r, rp = q-p, qp-pp
T1 = paper_B3_third(q, qp, p, pp)                 # [k, j]
T2 = paper_B3_third(p, pp, q, qp)                 # [j, k], colour f^{bce} rho^e(p-q)
P = p*qp-pp*q; AC = np.zeros((2, 2), complex)
for kk in range(2):
    for j in range(2):
        for l in range(2):
            for m in range(2):
                C1 = -1j*g*P[m]*np.sqrt(qp*pp*rp)/(np.sqrt(2)*(P@P))*(d2[kk, m]*d2[j, l]/qp-d2[kk, l]*d2[m, j]/pp-d2[l, m]*d2[kk, j]/rp)
                AC[kk, j] += 2*(-np.sqrt(2)*g*r[l]/(np.sqrt(rp)*(r@r)))*C1
res = T1-np.conj(T2).T+AC
e = np.abs(res).max()/np.abs(T1).max(); ok &= e < 1e-12
print("a^dag a unitarity relation, triple-vertex part: relative residual", f"{e:.1e}", " (with the opposite sign of the p^+>k^+ term:", f"{np.abs(T1+np.conj(T2).T+AC).max()/np.abs(T1).max():.2f})")
# ---- norm of |Psi_gg>: sum_{j,l} T T  and the xi integral -------------------------------------------
import sympy as sp
xi, dl = sp.symbols('xi delta', positive=True)
def Tt(i, j, l, m): return sp.KroneckerDelta(i, m)*sp.KroneckerDelta(j, l)-sp.KroneckerDelta(j, m)*sp.KroneckerDelta(i, l)/xi-sp.KroneckerDelta(l, m)*sp.KroneckerDelta(i, j)/(1-xi)
tot = sp.simplify(sum(Tt(i, j, l, m)**2 for i in range(2) for j in range(2) for l in range(2) for m in range(2)))
print("sum_{ijlm} T^2 - 4(1+1/xi^2+1/(1-xi)^2) =", sp.simplify(tot-4*(1+1/xi**2+1/(1-xi)**2)))
I = sp.integrate(xi*(1-xi)+(1-xi)/xi+xi/(1-xi), (xi, dl, 1-dl))
print("int_delta^{1-delta} [xi(1-xi)+(1-xi)/xi+xi/(1-xi)] = -2 log(delta) - 11/6 + O(delta):", sp.series(sp.simplify(I+2*sp.log(dl)+sp.Rational(11, 6)), dl, 0, 2))
print("ALL NUMERICAL CHECKS PASSED" if ok else "FAILED")
