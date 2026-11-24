# c3_expansion.py -- expansion of the kernel of Row I.II at p -> q and p -> 0 (Explain_Row_I_II_3rho.tex)
# Run: python3 -I c3_expansion.py -> c3_expansion.log
# E1  sympy series at p = q + r, r = eps*(cos phi, sin phi), L = log(q^2/r^2) kept as a symbol:
#     A = R(q,p)/p^2,  B = R(q,q-p)/(q-p)^2,  F = (A-B)/2
# E2  same at p -> 0 (r = p)
# E3  numerical check of the remainders (mpmath)
# E4  integral over the annulus eps < |p-q| < a: single term ~ log^2, antisymmetrised kernel converges
import sympy as sp, mpmath as mp
out = []
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)
eps, phi, L, Q = sp.symbols('epsilon phi L Q', positive=True)   # q = Q*(1,0) without loss of generality (rotation covariance)
q = sp.Matrix([Q, 0]); n = sp.Matrix([sp.cos(phi), sp.sin(phi)]); r = eps * n
def Rsym(qv, pv, logs):
    """R^i(q,p) with the two logarithms supplied (so that log(q^2/r^2) can be kept as the symbol L)."""
    q2 = qv.dot(qv); p2 = pv.dot(pv); qp = qv.dot(pv); d2 = (qv - pv).dot(qv - pv)
    return -2 * qv * logs[0] + ((2 * p2 - 2 * qp + q2) * qv + q2 * pv) / d2 * logs[1]
def ser(expr, order=1):
    return sp.simplify(sp.series(expr, eps, 0, order).removeO())
# ---- E1: p = q + r
p = q + r; p2 = p.dot(p)
A = Rsym(q, p, [sp.log(p2 / Q**2), L]) / p2                       # log(q^2/(q-p)^2) = L
s = -r                                                            # q - p = -r ;  q - s = p
B = Rsym(q, s, [-L, sp.log(Q**2 / p2)]) / (eps**2)                # log(s^2/q^2) = -L
A_ser = sp.Matrix([ser(A[i]) for i in range(2)]); B_ser = sp.Matrix([ser(B[i]) for i in range(2)])
qh = q / Q
A_claim = 2 * q * L / eps**2 + L * (r - 2 * qh.dot(r) * qh) / eps**2
B_claim = 2 * q * L / eps**2 - 2 * qh.dot(r) * qh / eps**2
F_claim = L * (r - 2 * qh.dot(r) * qh) / (2 * eps**2) + qh.dot(r) * qh / eps**2
def sing(M):   # keep eps^-2 and eps^-1
    return sp.Matrix([sp.simplify(sum(sp.expand(M[i]).coeff(eps, k) * eps**k for k in (-2, -1))) for i in range(2)])
pr('E1  p -> q.  terms of order eps^-2, eps^-1 (L = log(q^2/r^2), q = Q(1,0), r = eps(cos phi, sin phi)):')
pr('    A =', list(sing(A_ser))); pr('    B =', list(sing(B_ser)))
pr('    A - [2 q L/r^2 + L (r - 2 (qhat.r) qhat)/r^2] :', list(sp.simplify(sing(A_ser) - A_claim)))
pr('    B - [2 q L/r^2 - 2 (qhat.r) qhat/r^2]         :', list(sp.simplify(sing(B_ser) - B_claim)))
pr('    F - [L (r - 2 (qhat.r) qhat)/(2 r^2) + (qhat.r) qhat/r^2] :', list(sp.simplify(sing((A_ser - B_ser) / 2) - F_claim)))
F0 = sp.Matrix([sp.simplify(sp.expand((A_ser[i] - B_ser[i]) / 2).coeff(eps, 0)) for i in range(2)])
pr('    order eps^0 of F:', list(F0))
pr('    angular average of the eps^0 term of F:', [sp.simplify(sp.integrate(F0[i], (phi, 0, 2 * sp.pi)) / (2 * sp.pi)) for i in range(2)])
# ---- E2: p -> 0, r = p ; now L = log(q^2/p^2)
p = r; d2 = (q - p).dot(q - p)
A0 = Rsym(q, p, [-L, sp.log(Q**2 / d2)]) / eps**2
s = q - p
B0 = Rsym(q, s, [sp.log(d2 / Q**2), L]) / d2                      # q - s = p: log(q^2/p^2) = L
A0s = sp.Matrix([ser(A0[i]) for i in range(2)]); B0s = sp.Matrix([ser(B0[i]) for i in range(2)])
pr('E2  p -> 0 (r = p).')
pr('    A - [2 q L/r^2 + 2 (qhat.r) qhat/r^2]          :', list(sp.simplify(sing(A0s) - (2 * q * L / eps**2 + 2 * qh.dot(r) * qh / eps**2))))
pr('    B - [2 q L/r^2 - L (r - 2 (qhat.r) qhat)/r^2]  :', list(sp.simplify(sing(B0s) - (2 * q * L / eps**2 - L * (r - 2 * qh.dot(r) * qh) / eps**2))))
pr('    F - [L (r - 2 (qhat.r) qhat)/(2 r^2) + (qhat.r) qhat/r^2] :', list(sp.simplify(sing((A0s - B0s) / 2) - F_claim)))

# ---- E3, E4 numerics with a generic q
mp.mp.dps = 30
qv = [mp.mpf('0.7'), mp.mpf('-0.4')]; q2 = qv[0]**2 + qv[1]**2
def Rm(qq, pp):
    a2 = qq[0]**2 + qq[1]**2; b2 = pp[0]**2 + pp[1]**2; ab = qq[0] * pp[0] + qq[1] * pp[1]
    d = (qq[0] - pp[0])**2 + (qq[1] - pp[1])**2; c = 2 * b2 - 2 * ab + a2
    return [-2 * qq[i] * mp.log(b2 / a2) + (c * qq[i] + a2 * pp[i]) / d * mp.log(a2 / d) for i in range(2)]
def ABF(pp):
    s = [qv[0] - pp[0], qv[1] - pp[1]]; b2 = pp[0]**2 + pp[1]**2; s2 = s[0]**2 + s[1]**2
    A = [x / b2 for x in Rm(qv, pp)]; B = [x / s2 for x in Rm(qv, s)]
    return A, B, [(A[i] - B[i]) / 2 for i in range(2)]
def lead(rv):
    r2 = rv[0]**2 + rv[1]**2; Lr = mp.log(q2 / r2); qr = (qv[0] * rv[0] + qv[1] * rv[1]) / q2
    return [Lr * (rv[i] - 2 * qr * qv[i]) / (2 * r2) + qr * qv[i] / r2 for i in range(2)]
pr('E3  generic q = (0.7,-0.4): max over 64 angles of |F - F_lead| / log(q^2/r^2)  (F_lead = L (r - 2(qhat.r)qhat)/(2r^2) + (qhat.r)qhat/r^2)')
for name, c, sgn in [('p = q + r', qv, 1), ('p = r    ', [0, 0], 1)]:
    row = []
    for rs in ['1e-2', '1e-4', '1e-6', '1e-8']:
        rr = mp.mpf(rs); m = 0
        for j in range(64):
            th = 2 * mp.pi * j / 64; rv = [rr * mp.cos(th), rr * mp.sin(th)]
            F = ABF([c[0] + rv[0], c[1] + rv[1]])[2]; Fl = lead(rv)
            m = max(m, max(abs(F[i] - Fl[i]) for i in range(2)) / mp.log(q2 / rr**2))
        row.append(mp.nstr(m, 5))
    pr('    ' + name + ' : r = 1e-2, 1e-4, 1e-6, 1e-8 :', ' '.join(row), ' (bounded: the remainder is O(log(1/r)))')
a = mp.mpf('0.1')
def annulus(comp, which, e):
    f = lambda t, th: mp.exp(2 * t) * ABF([qv[0] + mp.exp(t) * mp.cos(th), qv[1] + mp.exp(t) * mp.sin(th)])[which][comp]
    return mp.quad(f, [mp.log(e), mp.log(a)], mp.linspace(0, 2 * mp.pi, 5))
pr('E4  int_{eps<|p-q|<0.1} d^2p (x component):  eps | single term A | pi q_x [log^2(q^2/eps^2) - log^2(q^2/a^2)] | antisymmetrised F')
for es in ['1e-2', '1e-4', '1e-6', '1e-8']:
    e = mp.mpf(es)
    pr('     %s  %s  %s  %s' % (es, mp.nstr(annulus(0, 0, e), 8), mp.nstr(mp.pi * qv[0] * (mp.log(q2 / e**2)**2 - mp.log(q2 / a**2)**2), 8), mp.nstr(annulus(0, 2, e), 10)))
open(__file__.replace('.py', '.log'), 'w').write('\n'.join(out) + '\n')
