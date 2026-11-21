# s8_plus.py -- checks of the plus-prescription forms (author's request of 2026-11-21).
# Every closed-form function of the previous version (E, J, Ehat, Psi, pi delta - i P 1/s) was replaced by
# its p+ integral with the author's plus prescriptions; here each replacement is checked as an identity
# (with the prefactors as printed), the p+ integrals being computed by quadrature directly from the
# definition [1/p+]_+ f = (f(p+) - f(0))/p+ (no closed forms on the "new" side).  Run: python3 -I s8_plus.py
import numpy as np
from scipy import integrate, special
pi = np.pi; gE = np.euler_gamma
Q = dict(limit=400, epsabs=1e-13, epsrel=1e-12)

# ---------------- p+ integrals from the definitions (quadrature) ----------------
def _osc(fun, a, b, om):  # int_a^b fun(p) e^{-i om p} dp, b may be inf (QAWO / QAWF)
    if om == 0:
        return integrate.quad(fun, a, b, **Q)[0]+0j
    if np.isinf(b):
        c = integrate.quad(lambda p: fun(p+a), 0, np.inf, weight='cos', wvar=abs(om), limlst=200)[0]
        s = integrate.quad(lambda p: fun(p+a), 0, np.inf, weight='sin', wvar=abs(om), limlst=200)[0]
        return np.exp(-1j*om*a)*(c-1j*np.sign(om)*s)
    c = integrate.quad(fun, a, b, weight='cos', wvar=om, limit=2000)[0]
    s = integrate.quad(fun, a, b, weight='sin', wvar=om, limit=2000)[0]
    return c-1j*s
def _small(fun, a, b):
    re = integrate.quad(lambda p: fun(p).real, a, b, **Q)[0]; im = integrate.quad(lambda p: fun(p).imag, a, b, **Q)[0]
    return re+1j*im
def plus0(w, s, kp, upper):
    """int_0^upper dp w(p) [e^{-i s p/kp} - 1]  (plus prescription subtracting at p=0; w ~ 1/p at p->0)"""
    om = s/kp
    a = min(upper, kp/max(1.0, abs(s)))   # small-p piece with at most one oscillation
    r = _small(lambda p: w(p)*(np.exp(-1j*om*p)-1.0) if p > 0 else -1j*om*(w(1e-300)*1e-300), 0, a)
    if a < upper:
        r += _osc(w, a, upper, om)-integrate.quad(w, a, upper, **Q)[0]
    return r
def IE(s, kp):     # int_0^kp dp [1/p]_+ e^{-isp/kp}
    return plus0(lambda p: 1.0/p, s, kp, kp)
def IEk(s, kp):    # (1/kp) int_0^kp dp [kp/p]_+ e^{-isp/kp}
    return plus0(lambda p: kp/p, s, kp, kp)/kp
def IEq(s, kp):    # int_0^kp dp [kp/(kp-p)]_+ e^{-isp/kp}  (subtraction at p = kp), direct quadrature
    om = s/kp
    f = lambda p: kp/(kp-p)*(np.exp(-1j*om*p)-np.exp(-1j*om*kp)) if p < kp else 1j*om*kp*np.exp(-1j*om*kp)
    n = int(abs(s)/(2*pi))+2; pts = np.linspace(0, kp, n+1)
    return sum(_small(f, pts[i], pts[i+1]) for i in range(n))
def IEshift(s, kp):  # int_0^kp dp [kp/p]_+ (e^{-is} e^{+isp/kp})   (Row V, third term; subtraction of the whole phase at p=0)
    om = s/kp
    f = lambda p: kp/p*(np.exp(-1j*s)*np.exp(1j*om*p)-np.exp(-1j*s)) if p > 0 else 1j*om*kp*np.exp(-1j*s)
    n = int(abs(s)/(2*pi))+2; pts = np.linspace(0, kp, n+1)
    return sum(_small(f, pts[i], pts[i+1]) for i in range(n))
def IJ(s, kp):     # int_0^kp dp e^{-isp/kp}
    return _osc(lambda p: 1.0, 0, kp, s/kp)
def IH(s, kp):     # int_0^inf dp e^{-isp/kp}/(p+kp)
    return _osc(lambda p: 1.0/(p+kp), 0, np.inf, s/kp)
def IP(s, kp):     # int_0^inf dp [kp/(p(p+kp))]_+ e^{-isp/kp}
    return plus0(lambda p: kp/(p*(p+kp)), s, kp, np.inf)
# ---------------- closed forms of the previous version ----------------
def E(s):
    si, ci = special.sici(abs(s)); return ci-gE-np.log(abs(s))-1j*np.sign(s)*si
def J(s): return (1-np.exp(-1j*s))/(1j*s)
def Eh(s): return np.exp(1j*s)*special.exp1(1j*s)
def Psi(s): return -gE-np.log(abs(s))-1j*pi/2*np.sign(s)-Eh(s)

out = []
def rep(name, dev): out.append((name, dev)); print('%-78s %.2e' % (name, dev))
S = [-23.0, -5.2, -1.7, -0.3, 0.3, 1.7, 5.2, 23.0]
KP = [0.37, 1.0, 2.9]
rel = lambda a, b: abs(a-b)/max(abs(b), 1e-300)
# (1) function level
rep('E(s) = int_0^k+ dp+ [1/p+]_+ e^{-isp+/k+}', max(rel(IE(s, k), E(s)) for s in S for k in KP))
rep('E(s) = (1/k+) int_0^k+ dp+ [k+/p+]_+ e^{-isp+/k+}', max(rel(IEk(s, k), E(s)) for s in S for k in KP))
rep('k+ e^{-is} E(-s) = int_0^k+ dp+ [k+/(k+-p+)]_+ e^{-isp+/k+}', max(rel(IEq(s, k), k*np.exp(-1j*s)*E(-s)) for s in S for k in KP))
rep('k+ e^{-is} E(-s) = int_0^k+ dp+ [k+/p+]_+ e^{-is} e^{+isp+/k+}', max(rel(IEshift(s, k), k*np.exp(-1j*s)*E(-s)) for s in S for k in KP))
rep('k+ J(s) = int_0^k+ dp+ e^{-isp+/k+}', max(rel(IJ(s, k), k*J(s)) for s in S for k in KP))
rep('Ehat(s) = int_0^inf dp+ e^{-isp+/k+}/(p++k+)', max(rel(IH(s, k), Eh(s)) for s in S for k in KP))
rep('Psi(s) = int_0^inf dp+ [k+/(p+(p++k+))]_+ e^{-isp+/k+}', max(rel(IP(s, k), Psi(s)) for s in S for k in KP))
# distribution: int ds g(s) int_0^inf dp+ e^{-isp+/k+} = k+ [pi g(0) - i P int g(s)/s ds], g(s) = exp(-(s-a)^2)
dmax = 0
for a in [-1.3, 0.0, 0.4, 2.2]:
    for k in KP:
        gh = lambda u: np.sqrt(pi)*np.exp(-1j*a*u-u*u/4)   # int ds g(s) e^{-isu}
        lhs = _small(lambda p: gh(p/k), 0, 60*k)
        pv = integrate.quad(lambda s: np.exp(-(s-a)**2), a-14, a+14, weight='cauchy', wvar=0.0)[0]
        rhs = k*(pi*np.exp(-a*a)-1j*pv)
        dmax = max(dmax, rel(lhs, rhs))
rep('int_0^inf dp+ e^{-isp+/k+} = k+[pi delta(s) - i P 1/s] (Gaussian test functions)', dmax)
# asymptotics, Eq. (asym) of the conventions
for SS in [150.0, 400.0, 2000.0]:
    da = 0
    for s in [-SS, SS]:
        k = 1.3
        da = max(da, abs(IE(s, k)-(-gE-np.log(abs(s))-1j*pi/2*np.sign(s)-np.exp(-1j*s)/(1j*s)))*s*s,
                 abs(IP(s, k)-(-gE-np.log(abs(s))-1j*pi/2*np.sign(s)-1/(1j*s)))*s*s, abs(IH(s, k)-1/(1j*s))*s*s)
    print('asymptotic forms Eq.(asym): max |remainder| * s^2 at |s| = %g: %.3f (bounded: remainder O(s^-2))' % (SS, da))
    out.append(('asym |rem| s^2 at |s|=%g (not a deviation)' % SS, 0.0))

# (2) prefactors as printed: old coefficient x old function  vs  new coefficient x (1/2pi) x p+ integral
rng = np.random.default_rng(7)
def chk(name, old, new):
    d = 0
    for s in S:
        for k in KP:
            d = max(d, rel(new(s, k), old(s, k)))
    rep(name, d)
P = lambda n, m: 1.0/(n*pi**m)
chk('I.III last term N3b (B-I3, pairs A, D Group I): ig^4f/(8pi^5)(1/k+)E -> ig^4f/(4pi^4)(1/k+)int dp+/2pi [1/p+]_+',
    lambda s, k: P(8, 5)/k*E(s), lambda s, k: P(4, 4)/k*IE(s, k)/(2*pi))
chk('mixed brackets (I.III in C4, D, E; Eq. C-pairform): E(s) -> int_0^k+ dp+ [1/p+]_+ e',
    lambda s, k: E(s), lambda s, k: IE(s, k))
A_, B_, C_ = rng.normal(size=3)+1j*rng.normal(size=3)
chk('I.V first part bracket (B-I5, D/E, pair B): ig^4/(32pi^5)(1/k+)[J A - E B - e^{-is}E(-s) C]',
    lambda s, k: P(32, 5)/k*(J(s)*A_-E(s)*B_-np.exp(-1j*s)*E(-s)*C_),
    lambda s, k: P(16, 4)/k**2*(IJ(s, k)*A_-IEk(s, k)*k*B_-IEshift(s, k)*C_)/(2*pi))
chk('pair B common form: ig^4/(16pi^5)(1/k+)E -> ig^4/(8pi^4)(1/k+)int dp+/2pi [1/p+]_+', lambda s, k: P(16, 5)/k*E(s), lambda s, k: P(8, 4)/k*IE(s, k)/(2*pi))
chk('pair C common form: g^4/(16pi^5)(1/k+)E -> g^4/(8pi^4)(1/k+)int dp+/2pi [1/p+]_+', lambda s, k: P(16, 5)/k*E(s), lambda s, k: P(8, 4)/k*IE(s, k)/(2*pi))
K1, K2 = rng.normal(size=2)
for nm, (n, m), (n2, m2) in [('III.III first: 16pi^5 -> 8pi^4', (16, 5), (8, 4)), ('III.V, III.VI: 8pi^5 -> 4pi^4', (8, 5), (4, 4))]:
    chk('Group III %s: (1/k+)[Ehat(K1-K2) - Psi K2] -> (1/k+)int_0^inf dp+/2pi[..]' % nm,
        lambda s, k: P(n, m)/k*(Eh(s)*(K1-K2)-Psi(s)*K2), lambda s, k: P(n2, m2)/k*(IH(s, k)*(K1-K2)-IP(s, k)*K2)/(2*pi))
d = 0
for a in [-1.3, 0.4, 2.2]:
    for k in KP:
        gh = lambda u: np.sqrt(pi)*np.exp(-1j*a*u-u*u/4)
        new = -P(8, 4)/k*_small(lambda p: gh(p/k), 0, 60*k)/k/(2*pi)
        pv = integrate.quad(lambda s: np.exp(-(s-a)**2), a-14, a+14, weight='cauchy', wvar=0.0)[0]
        old = -P(16, 5)/k*(pi*np.exp(-a*a)-1j*pv)
        d = max(d, rel(new, old))
rep('Group III instantaneous: -(1/16pi^5)(1/k+)[pi delta - iP/s] -> -(1/8pi^4)(1/k+)int dp+/2pi (1/k+) e (smeared)', d)
chk('III.III second part soft: -g^4/(8pi^5)(1/k+)Psi -> -g^4/(4pi^4)(1/k+)int dp+/2pi[..]_+', lambda s, k: -P(8, 5)/k*Psi(s), lambda s, k: -P(4, 4)/k*IP(s, k)/(2*pi))

# (3) kernels at random transverse configurations
def Kv(v): return v/np.dot(v, v)
def dot(a, b): return float(np.dot(a, b))
kvec = np.array([1.3, -0.4])
def IV2_parts(x, y, z, xp, wp):
    Kx = Kv(xp-wp); Kzx = Kv(z-x); D2 = dot(y-x, y-x)-dot(y-z, y-z); xz2 = dot(x-z, x-z)
    yx2, yz2 = dot(y-x, y-x), dot(y-z, y-z)
    B2 = -dot(y-z, Kx)/xz2+0.5*dot(y-z, Kx)*dot(Kzx, Kv(y-x))
    B1 = -dot(y-x, Kx)/xz2-0.5*dot(y-x, Kx)*dot(Kzx, Kv(y-z))
    def R(xi):
        return (xi*(1-xi)*D2/xz2*dot(Kx, Kzx)
                + (1-xi)*dot(Kx, Kzx)*(dot(y-z, Kv(x-z))+0.5*dot(y-z, Kv(y-x))-D2/(2*xz2))
                + xi*dot(Kx, Kzx)*(dot(y-x, Kv(x-z))-0.5*dot(y-x, Kv(y-z))-D2/(2*xz2))
                + yx2/yz2*B2+yz2/yx2*B1
                - dot(y-x, Kzx)*dot(Kv(x-z), Kx)+0.5*dot(y-x, Kzx)*dot(Kv(y-z), Kx)
                - dot(y-z, Kzx)*dot(Kv(x-z), Kx)-0.5*dot(Kv(y-x), Kx)*dot(y-z, Kzx))
    Dxi = lambda xi: xi*yx2+(1-xi)*yz2
    # author's bracket (Row V, N2a): poly part, plus-prescribed coefficients Z, V, extra bracket
    def Rauth(xi):
        return (xi*(1-xi)*D2/xz2*dot(Kx, Kzx)
                + (1-xi)*dot(Kx, Kzx)*(dot(y-z, Kv(x-z))+0.5*dot(y-z, Kv(y-x))-D2/(2*xz2))
                + xi*dot(Kx, Kzx)*(dot(y-x, Kv(x-z))-0.5*dot(y-x, Kv(y-z))-D2/(2*xz2))
                - dot(y-x, Kzx)*dot(Kv(x-z), Kx)+0.5*dot(y-x, Kzx)*dot(Kv(y-z), Kx)
                - dot(y-z, Kzx)*dot(Kv(x-z), Kx)-0.5*dot(Kv(y-x), Kx)*dot(y-z, Kzx)
                + (-dot(y-z, Kx)/xz2+0.5*dot(Kv(y-x), Kzx)*dot(y-z, Kx)-dot(y-x, Kx)/xz2-0.5*dot(y-x, Kx)*dot(Kv(y-z), Kzx)))
    Z = dot(y-z, Kx)/xz2-0.5*dot(Kv(y-x), Kzx)*dot(y-z, Kx)
    V = dot(y-x, Kx)/xz2+0.5*dot(y-x, Kx)*dot(Kv(y-z), Kzx)
    return R, Dxi, B1, B2, yx2, yz2, Rauth, Z, V
d_old = d_auth = 0
for trial in range(6):
    x, y, z, xp, wp = rng.uniform(-1, 1, size=(5, 2))
    s = dot(kvec, x-z); k = [0.6, 1.0, 2.3][trial % 3]
    R, Dxi, B1, B2, yx2, yz2, Rauth, Z, V = IV2_parts(x, y, z, xp, wp)
    reg = _small(lambda xi: np.exp(-1j*xi*s)*R(xi)/Dxi(xi), 0, 1)
    old = P(16, 5)/k*(reg-E(s)*B2/yz2-np.exp(-1j*s)*E(-s)*B1/yx2)
    regp = _small(lambda p: np.exp(-1j*s*p/k)*R(p/k)/Dxi(p/k), 0, k)
    new = P(8, 4)/k**2*(regp-k*IEk(s, k)*B2/yz2-IEq(s, k)*B1/yx2)/(2*pi)
    d_old = max(d_old, rel(new, old))
    # author's form: (1/k+) int dp+/2pi { f(p) Rauth + [k/p]_+ f Z + [k/(k-p)]_+ f V }, f = e/[p (y-x)^2 + (k-p)(y-z)^2]
    f = lambda p: np.exp(-1j*s*p/k)/(p*yx2+(k-p)*yz2)
    main = _small(lambda p: f(p)*Rauth(p/k), 0, k)
    pz = _small(lambda p: k/p*(f(p)-f(0.0))*Z if p > 0 else 0j, 0, k)
    pv_ = _small(lambda p: k/(k-p)*(f(p)-f(k))*V if p < k else 0j, 0, k)
    auth = P(8, 4)/k*(main+pz+pv_)/(2*pi)
    d_auth = max(d_auth, rel(new, auth))
rep('I.V second part brace, actual kernel (6 configs): new p+ form vs previous xi + E form', d_old)
rep('I.V second part brace, actual kernel: new p+ form vs author\'s form (plus on e/denominator)', d_auth)
# pair sums (A; B, C, D have the same p+ structure with other kappa, V) and the identity (C-pairform)
d_pair = d_id = 0
for trial in range(6):
    x, y, z, w, wp, xp = rng.uniform(-1, 1, size=(6, 2)); k = [0.6, 1.0, 2.3][trial % 3]
    Vf = dot(Kv(xp-wp), Kv(y-w))
    kap = lambda u: dot(Kv(x-z), Kv(u-z))
    s1, s3 = dot(kvec, z-w), dot(kvec, z-wp)
    old = -P(8, 5)/k*Vf*((kap(w)-kap(wp))*E(s1)+kap(wp)*(E(s1)-Psi(s3)))
    I1, I3 = IE(s1, k), IP(s3, k)
    new = -P(4, 4)/k*Vf*((kap(w)-kap(wp))*I1/(2*pi)+kap(wp)*(I1-I3)/(2*pi))
    d_pair = max(d_pair, rel(new, old))
    d_id = max(d_id, abs((kap(w)*I1-kap(wp)*I3)-((kap(w)-kap(wp))*I1+kap(wp)*(I1-I3))))
rep('pair sum (Eq. C-A form), actual kernels: new vs previous', d_pair)
rep('identity Eq. (C-pairform) with the p+ integrals (absolute)', d_id)
# KZ_LE + KZ_E = ell (C4 split)
d = 0
for trial in range(4):
    x, z, w = rng.uniform(-1, 1, size=(3, 2)); k = 1.1
    ell = -0.5*np.log(dot(z-w, z-w)/dot(x-w, x-w)); I1 = IE(dot(kvec, z-w), k)
    d = max(d, abs((ell-I1)+I1-ell), rel(ell-I1, ell-E(dot(kvec, z-w))))
rep('C4 split: (ell - int[1/p+]_+ e) + int[1/p+]_+ e = ell; blue part vs previous ell - E', d)
# (4) tail of pair A with the p+ integrals (colour stripped, units of the prefactor): z^2 x members and z^2 x sum
print('pair A tail, z = R(cos0.7, sin0.7), x,y,x\',w,w\' fixed:')
x, y, xp, w, wp = np.array([0.3, -0.2]), np.array([-0.5, 0.4]), np.array([0.1, 0.6]), np.array([0.2, 0.1]), np.array([-0.3, -0.4])
k = 1.0; tails = []
for Rr in [1e2, 1e3, 1e4, 1e5]:
    z = Rr*np.array([np.cos(0.7), np.sin(0.7)])
    Vf = dot(Kv(xp-wp), Kv(y-w)); kap = lambda u: dot(Kv(x-z), Kv(u-z))
    m1 = Vf*kap(w)*IE(dot(kvec, z-w), k); m3 = -Vf*kap(wp)*IP(dot(kvec, z-wp), k)
    tails.append((Rr, Rr**2*m1, Rr**2*m3, Rr**2*(m1+m3)))
    print('  |z|=%.0e  z^2*GroupI = %s  z^2*GroupIII = %s  |z^2*sum| = %.3e' % (Rr, np.round(Rr**2*m1, 3), np.round(Rr**2*m3, 3), abs(Rr**2*(m1+m3))))
ratio = tuple(abs(tails[i+1][3])/abs(tails[i][3]) for i in range(3))
print('  ratios |z^2 sum|(10R)/|z^2 sum|(R): %.3f, %.3f, %.3f (1/|z| fall-off: 0.1)' % ratio)
open('s8_plus.log', 'w').write('\n'.join('%-78s %.2e' % t for t in out)+'\n'+'pair A tail: '+'; '.join('|z|=%.0e: z^2 I=%.3f%+.3fi, z^2 III=%.3f%+.3fi, |z^2 sum|=%.2e' % (t[0], t[1].real, t[1].imag, t[2].real, t[2].imag, abs(t[3])) for t in tails)+'\n')
print('max deviation', max(v for _, v in out))
