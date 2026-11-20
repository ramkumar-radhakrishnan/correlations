# Kernel checks for Final_terms_explicit.tex (paper phase convention e^{-ik(w-w')}).
# (1) Group I Row III: author's combined block == paper M^{ij} + l K(y,w,z) P.K(x-w)
# (2) last term of N_3b: author's l-form vs verified -E form (not equal; tails)
# (3) Group I Row V second part: author's plus prescriptions == F_reg - E B2 - e^{-is}E(-s) B1 (quadrature)
# (4) Group III Row III second part: integrand identity author == regular + soft (Psi) integrand
# (5) Group II Row I second part: author's N2a..N2e == paper K2a..K2f
# (6) recoil functions: closed forms vs quadrature; vee-cut convergence of Ehat
# (7) Group II Row IV: p+ -> xi substitution
import numpy as np
from scipy import integrate, special
rng = np.random.default_rng(7)
dot = lambda a, b: a @ b
sq = lambda a: a @ a
Kv = lambda v: v/sq(v)
Kab = lambda a, b, z: dot(Kv(a-z), Kv(b-z))
gE = np.euler_gamma
out = []
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

# ---------- recoil functions ----------
def Ecal(s):  # int_0^1 dt/t (e^{-ist}-1)
    si, ci = special.sici(abs(s))
    return ci - gE - np.log(abs(s)) - 1j*np.sign(s)*si
def Ehat(s):  # int_0^inf dt e^{-ist}/(1+t) = e^{is} E1(is)
    return np.exp(1j*s)*special.exp1(1j*s)
def Psi(s):
    return -gE - np.log(abs(s)) - 0.5j*np.pi*np.sign(s) - Ehat(s)
def Jf(s): return (1-np.exp(-1j*s))/(1j*s)
def quadc(f, a, b, **kw):
    r = integrate.quad(lambda t: f(t).real, a, b, limit=2000, **kw)[0]
    i = integrate.quad(lambda t: f(t).imag, a, b, limit=2000, **kw)[0]
    return r+1j*i

# (6) closed forms
worst = 0
for s in [0.37, -1.3, 4.2, -9.7]:
    e1 = quadc(lambda t: (np.exp(-1j*s*t)-1)/t, 0, 1)
    worst = max(worst, abs(e1-Ecal(s)))
    # Psi: split [0,1] + [1,inf) with oscillatory weights
    p1 = quadc(lambda t: (np.exp(-1j*s*t)-1)/(t*(1+t)), 0, 1)
    # int_1^inf (e^{-ist}-1)/(t(1+t)) = int_1^inf cos(st)/(t(1+t)) - i sgn int sin(|s|t)/(t(1+t)) - log2
    c = integrate.quad(lambda t: 1/(t*(1+t)), 1, np.inf, weight='cos', wvar=abs(s))[0]
    sn = integrate.quad(lambda t: 1/(t*(1+t)), 1, np.inf, weight='sin', wvar=abs(s))[0]
    p2 = c - 1j*np.sign(s)*sn - np.log(2)
    worst = max(worst, abs(p1+p2-Psi(s)))
    # Ehat
    c = integrate.quad(lambda t: 1/(1+t), 0, np.inf, weight='cos', wvar=abs(s))[0]
    sn = integrate.quad(lambda t: 1/(1+t), 0, np.inf, weight='sin', wvar=abs(s))[0]
    worst = max(worst, abs(c-1j*np.sign(s)*sn-Ehat(s)))
pr('(6) recoil functions closed form vs quadrature, max |diff| =', '%.1e' % worst)
# vee-cut: int_0^T dt e^{-ist}/(1+t) -> Ehat(s)
s = 1.7
for T in [1e1, 1e2, 1e3, 1e4]:
    c = integrate.quad(lambda t: np.cos(s*t)/(1+t), 0, T, limit=20000)[0]
    sn = integrate.quad(lambda t: np.sin(s*t)/(1+t), 0, T, limit=20000)[0]
    pr('(6) vee-cut T=%g: |int_0^T e^{-ist}/(1+t) - Ehat| = %.2e' % (T, abs(c-1j*sn-Ehat(s))))

# ---------- (1) Group I Row III ----------
def M_paper(x, w, z):
    a = x-w; b = x-z; c = z-w; A = sq(a); B = sq(b); C = sq(c)
    A0 = np.arccos(dot(a, c)/np.sqrt(A*C))/np.sqrt(A*B-dot(a, b)**2)
    d = np.eye(2); O = np.outer
    M = np.log(B/A)*(d/(2*C)-O(a, b)/(A*B)+1/(C*(A-B))*(O(a, c)*B/A+O(c, b)*A/B-O(c, a)-O(b, c)))
    M += A0/C*(O(c, a)+O(b, c)-dot(a, b)*(O(c, b)/B+O(a, c)/A))
    M += -1/(2*C)*np.log(C/A)*O(c, b)/B
    return M
def ell(z, w, x): return -0.5*np.log(sq(z-w)/sq(x-w))
def rowIII_paper(x, y, z, w, xp, wp):
    P = Kv(xp-wp)
    return P@M_paper(x, w, z)@Kv(y-z) + ell(z, w, x)*Kab(y, w, z)*dot(P, Kv(x-w))
def rowIII_author(x, y, z, w, xp, wp):
    P = Kv(xp-wp); a = x-w; b = x-z; c = z-w; A = sq(a); B = sq(b); C = sq(c); KY = Kv(y-z)
    t = (dot(P, KY)/(2*C)+dot(P, Kv(z-w))*dot(KY, Kv(x-z))-dot(P, Kv(x-w))*dot(KY, Kv(z-w))
         - dot(P, Kv(x-w))*dot(KY, Kv(x-z)))*np.log(B/A)
    ac = np.arccos((A-dot(a, b))/np.sqrt(A*C))
    t += ac/(np.sqrt(A*B-dot(a, b)**2)*C)*(dot(P, c)*dot(KY, a)-dot(a, b)/B*dot(P, c)*dot(KY, b)
                                          + dot(P, b)*dot(KY, c)-dot(a, b)/A*dot(P, a)*dot(KY, c))
    t += 1/(2*C)*np.log(C/A)*(dot(P, Kv(x-w))*dot(KY, c)-dot(P, c)*dot(KY, Kv(x-z)))
    return t
rel = 0
for _ in range(200):
    pts = [rng.normal(size=2)*2 for _ in range(6)]
    p, q = rowIII_paper(*pts), rowIII_author(*pts)
    rel = max(rel, abs(p-q)/abs(p))
pr('(1) Row I.III author combined block vs paper M^ij + l K P.K: max rel diff = %.1e (200 random points)' % rel)

# ---------- (2) last term of N_3b ----------
kvec = np.array([1.0, 0.0])
x, y, w, wp, xp = [rng.normal(size=2) for _ in range(5)]
n = np.array([np.cos(0.7), np.sin(0.7)])
pr('(2) last term of N_3b, z = R n, units 2iC_f, colour C_A stripped:')
for R in [1e2, 1e3, 1e4]:
    z = R*n
    P = Kv(xp-wp)
    TA = ell(z, w, x)*Kab(x, w, z)*dot(P, Kv(y-w))                     # author
    TV = -Ecal(kvec@(z-w))*Kab(x, w, z)*dot(P, Kv(y-w))                # verified (Group I member of pair A)
    T3 = +Psi(kvec@(z-wp))*Kab(x, wp, z)*dot(P, Kv(y-w))               # Group III Row V member (after y'->w', colour id.)
    pr('    |z|=%g: z^2*author=%.4f%+.4fi  z^2*verified=%.4f%+.4fi  z^2*(verified+III.V)=%.2e  z^2*(author+III.V)=%.3f%+.3fi'
       % (R, (R*R*TA).real, (R*R*TA).imag, (R*R*TV).real, (R*R*TV).imag, abs(R*R*(TV+T3)), (R*R*(TA+T3)).real, (R*R*(TA+T3)).imag))
d = []
for _ in range(5):
    x, y, z, w, xp, wp = [rng.normal(size=2) for _ in range(6)]
    P = Kv(xp-wp)
    d.append(abs(ell(z, w, x)+Ecal(kvec@(z-w))))
pr('(2) pointwise |l + E| at 5 random points (zero iff equal):', ' '.join('%.3f' % v for v in d))

# ---------- (3) Group I Row V second part ----------
def rowV2_author(s, x, y, z, xp, wp):
    P = Kv(xp-wp); KZ = Kv(z-x); DY = sq(y-x)-sq(y-z)
    D = lambda t: t*sq(y-x)+(1-t)*sq(y-z)
    T1 = dot(y-z, Kv(x-z))+0.5*dot(y-z, Kv(y-x))-DY/(2*sq(x-z))
    T2 = dot(y-x, Kv(x-z))-0.5*dot(y-x, Kv(y-z))-DY/(2*sq(x-z))
    reg = lambda t: (t*(1-t)*DY/sq(z-x)*dot(P, KZ)+(1-t)*dot(P, KZ)*T1+t*dot(P, KZ)*T2
                     - dot(y-x, KZ)*dot(Kv(x-z), P)+0.5*dot(y-x, KZ)*dot(Kv(y-z), P)
                     - dot(y-z, KZ)*dot(Kv(x-z), P)-0.5*dot(Kv(y-x), P)*dot(y-z, KZ)
                     + (-dot(y-z, P)/sq(x-z)+0.5*dot(Kv(y-x), KZ)*dot(y-z, P)-dot(y-x, P)/sq(x-z)
                        - 0.5*dot(y-x, P)*dot(Kv(y-z), KZ)))
    c1 = dot(y-z, P)/sq(x-z)-0.5*dot(Kv(y-x), KZ)*dot(y-z, P)       # coefficient of [1/t]_+
    c2 = dot(y-x, P)/sq(x-z)+0.5*dot(y-x, P)*dot(Kv(y-z), KZ)       # coefficient of [1/(1-t)]_+
    h = lambda t: np.exp(-1j*s*t)/D(t)
    I = quadc(lambda t: h(t)*reg(t), 0, 1)
    I += c1*quadc(lambda t: (h(t)-h(0))/t, 0, 1)
    I += c2*quadc(lambda t: (h(t)-h(1))/(1-t), 0, 1)
    return I
def rowV2_paper(s, x, y, z, xp, wp):
    P = Kv(xp-wp); KZ = Kv(z-x); O = np.outer; d = np.eye(2)
    YX = sq(y-x); YZ = sq(y-z); XZ = sq(x-z); DY = YX-YZ
    B1 = O(y-x, x-z)/XZ-O(y-x, y-z)/(2*YZ)
    B2 = O(y-z, x-z)/XZ+O(y-z, y-x)/(2*YX)          # B2^{im} = (y-z)^i(x-z)^m/XZ + (y-x)^m(y-z)^i/(2YX)
    T1 = dot(y-z, x-z)/XZ+dot(y-x, y-z)/(2*YX)-DY/(2*XZ)
    T2 = dot(y-x, x-z)/XZ-dot(y-x, y-z)/(2*YZ)-DY/(2*XZ)
    Q = O(x-z, y-x)/XZ-O(y-z, y-x)/(2*YZ)+O(x-z, y-z)/XZ+O(y-x, y-z)/(2*YX)   # Q^{im}, i first index
    N = lambda xi: d*xi*(1-xi)*DY/XZ+d*((1-xi)*T1+xi*T2)+YX/YZ*B2+YZ/YX*B1-Q
    F = quadc(lambda xi: np.exp(-1j*xi*s)*(P@N(xi)@KZ)/(xi*YX+(1-xi)*YZ), 0, 1)
    F += -Ecal(s)*(P@B2@KZ)/YZ-np.exp(-1j*s)*Ecal(-s)*(P@B1@KZ)/YX
    return F
rel = 0
for _ in range(20):
    x, y, z, xp, wp = [rng.normal(size=2) for _ in range(5)]
    s = kvec@(x-z)*1.3
    a, b = rowV2_author(s, x, y, z, xp, wp), rowV2_paper(s, x, y, z, xp, wp)
    rel = max(rel, abs(a-b)/abs(b))
pr('(3) Row I.V 2nd part: author plus prescriptions vs F_reg - E B2/(y-z)^2 - e^{-is}E(-s)B1/(y-x)^2: max rel diff = %.1e (20 points)' % rel)
pr('    11/6: int_0^1 dt [t(1-t)-2] =', integrate.quad(lambda t: t*(1-t)-2, 0, 1)[0], '(= -11/6 = %.6f)' % (-11/6))

# ---------- (4) Group III Row III second part, integrand identity ----------
def III3_author_integrand(t, s, x, xp, yp, z, w):
    a = x-w; b = x-z; c = z-w; A = sq(a); B = sq(b); C = sq(c)
    KK = dot(Kv(yp-z), Kv(xp-yp))
    I_ = dot(a, Kv(xp-yp))*dot(Kv(z-w), Kv(yp-z))+0.5*dot(Kv(x-z), Kv(yp-z))*dot(a, Kv(xp-yp))
    V = (dot(b, Kv(yp-z))*dot(Kv(z-w), Kv(xp-yp))-0.5*dot(b, Kv(yp-z))*dot(Kv(x-w), Kv(xp-yp))
         + dot(a, Kv(yp-z))*dot(Kv(z-w), Kv(xp-yp))+0.5*dot(Kv(x-z), Kv(xp-yp))*dot(a, Kv(yp-z)))
    h = lambda tt: np.exp(-1j*s*tt)/(tt*B+A)
    br = KK*(t*(B-A)/(C*(1+t)**2)+t/(1+t)*(dot(b, Kv(z-w))-0.5*dot(b, Kv(x-w))-(B-A)/(2*C))
             + 1/(1+t)*(dot(a, Kv(z-w))+0.5*dot(a, Kv(x-z))-(B-A)/(2*C)))-V
    # -([1/(t(1+t))]_+ + 1/(1+t)) I acting on h
    return h(t)*br-I_*((h(t)-h(0))/(t*(1+t))+h(t)/(1+t))
def III3_paper_integrand(t, s, x, xp, yp, z, w):
    a = x-w; b = x-z; c = z-w; A = sq(a); B = sq(b); C = sq(c)
    KK = dot(Kv(yp-z), Kv(xp-yp))
    I_ = dot(a, Kv(xp-yp))*c@Kv(yp-z)/C+dot(b, Kv(yp-z))*dot(a, Kv(xp-yp))/(2*B)
    V = (dot(b, Kv(yp-z))*dot(c, Kv(xp-yp))/C-dot(b, Kv(yp-z))*dot(a, Kv(xp-yp))/(2*A)
         + dot(a, Kv(yp-z))*dot(c, Kv(xp-yp))/C+dot(b, Kv(xp-yp))*dot(a, Kv(yp-z))/(2*B))
    reg = np.exp(-1j*s*t)/(t*B+A)*(KK*(t*(B-A)/(C*(1+t)**2)+t/(1+t)*(dot(b, c)/C-dot(a, b)/(2*A)-(B-A)/(2*C))
                                       + 1/(1+t)*(dot(a, c)/C+dot(a, b)/(2*B)-(B-A)/(2*C)))+(B-A)/((1+t)*A)*I_-V)
    soft = -I_/A*(np.exp(-1j*s*t)-1)/(t*(1+t))   # integrates to -Psi(s) I/A
    return reg+soft
rel = 0
for _ in range(200):
    x, xp, yp, z, w = [rng.normal(size=2) for _ in range(5)]
    t = rng.uniform(0.01, 30); s = rng.normal()*2
    a, b = III3_author_integrand(t, s, x, xp, yp, z, w), III3_paper_integrand(t, s, x, xp, yp, z, w)
    rel = max(rel, abs(a-b)/abs(b))
pr('(4) Group III Row III 2nd part: author integrand vs regular + soft integrand: max rel diff = %.1e (200 points)' % rel)

# ---------- (5) Group II Row I second part ----------
def II1_author(x, xp, z, w, wp):
    a = x-w; ap = xp-wp; b = x-z; bp = xp-z; c = z-w; cp = z-wp
    A, Ap, B, Bp, C, Cp = map(sq, (a, ap, b, bp, c, cp))
    D = Ap*B-A*Bp
    br8 = (dot(a, Kv(z-wp))*dot(bp, Kv(z-w))-0.5*dot(a, Kv(xp-wp))*dot(Kv(z-w), bp)
           + 0.5*dot(Kv(x-z), bp)*dot(a, Kv(z-wp))-0.25*dot(bp, Kv(x-z))*dot(Kv(xp-wp), a)
           + dot(ap, Kv(z-w))*dot(b, Kv(z-wp))+0.5*dot(ap, Kv(z-w))*dot(b, Kv(xp-z))
           - 0.5*dot(b, Kv(z-wp))*dot(ap, Kv(x-w))-0.25*dot(ap, Kv(x-w))*dot(b, Kv(xp-z)))
    N2a = br8*np.log(B/Bp)/D
    N2b = -br8*np.log(A/Ap)/D
    N2c = np.log(B/A)/D*(-(B/A)*(dot(ap, a)*dot(c, cp)/(C*Cp)+dot(ap, a)*dot(bp, c)/(2*C*Bp)+dot(a, ap)*dot(b, cp)/(2*B*Cp)
                                 + dot(b, bp)*dot(ap, a)/(4*B*Bp))
                         - (A/B)*(dot(bp, b)*dot(c, cp)/(C*Cp)-dot(bp, b)*dot(ap, c)/(2*C*Ap)-dot(a, cp)*dot(b, bp)/(2*A*Cp)
                                  + dot(b, bp)*dot(ap, a)/(4*A*Ap))
                         + 1/(B-A)*(B*(dot(a, c)*(Bp-Ap)/(2*C*Cp)+dot(b, a)*(Bp-Ap)/(4*B*Cp)+dot(ap, cp)*(B-A)/(2*C*Cp)
                                       + dot(bp, ap)*(B-A)/(4*C*Bp))
                                    - A*(dot(b, c)*(Bp-Ap)/(2*C*Cp)-dot(b, a)*(Bp-Ap)/(4*A*Cp)+dot(bp, cp)*(B-A)/(2*C*Cp)
                                         - dot(bp, ap)*(B-A)/(4*C*Ap)))
                         - 0.5*B*A*(Bp-Ap)/(C*Cp*(B-A)))
    N2d = np.log(Bp/Ap)/D*((Bp/Ap)*(dot(ap, a)*dot(Kv(z-w), Kv(z-wp))+0.5*dot(ap, a)*dot(Kv(xp-z), Kv(z-w))
                                    + 0.5*dot(a, ap)*dot(Kv(x-z), Kv(z-wp))+0.25*dot(Kv(x-z), Kv(xp-z))*dot(ap, a))
                           + (Ap/Bp)*(dot(bp, b)*dot(Kv(z-w), Kv(z-wp))-0.5*dot(bp, b)*dot(Kv(xp-wp), Kv(z-w))
                                      - 0.5*dot(Kv(x-w), Kv(z-wp))*dot(b, bp)+0.25*dot(b, bp)*dot(Kv(xp-wp), Kv(x-w)))
                           + 1/(Ap-Bp)*(Bp*(dot(a, Kv(z-w))*(Bp-Ap)/(2*Cp)+dot(a, Kv(x-z))*(Bp-Ap)/(4*Cp)
                                            + dot(ap, Kv(z-wp))*(B-A)/(2*C)+dot(ap, Kv(xp-z))*(B-A)/(4*C))
                                        - Ap*(dot(b, Kv(z-w))*(Bp-Ap)/(2*Cp)-dot(b, Kv(x-w))*(Bp-Ap)/(4*Cp)
                                              + dot(bp, Kv(z-wp))*(B-A)/(2*C)-dot(bp, Kv(xp-wp))*(B-A)/(4*C)))
                           + 0.5*Ap*Bp*(B-A)/(C*Cp*(Bp-Ap)))
    N2e = -0.5/(C*Cp)
    return N2a+N2b+N2c+N2d+N2e, (N2a, N2b, N2c, N2d)
def II1_paper(x, xp, z, w, wp):
    a = x-w; ap = xp-wp; b = x-z; bp = xp-z; c = z-w; cp = z-wp
    A, Ap, B, Bp, C, Cp = map(sq, (a, ap, b, bp, c, cp))
    D = Ap*B-A*Bp; l = np.log(B/A); lp = np.log(Bp/Ap)
    K2a = (l-lp)/D*(dot(a, cp)*dot(bp, c)/(C*Cp)-dot(a, ap)*dot(c, bp)/(2*C*Ap)+dot(b, bp)*dot(a, cp)/(2*B*Cp)
                    - dot(b, bp)*dot(a, ap)/(4*B*Ap)+dot(ap, c)*dot(b, cp)/(C*Cp)+dot(ap, c)*dot(b, bp)/(2*C*Bp)
                    - dot(b, cp)*dot(a, ap)/(2*A*Cp)-dot(a, ap)*dot(b, bp)/(4*A*Bp))
    K2b = 1/D*(Bp/Ap*lp-B/A*l)*(dot(a, ap)*dot(c, cp)/(C*Cp)+dot(a, ap)*dot(bp, c)/(2*C*Bp)+dot(a, ap)*dot(b, cp)/(2*B*Cp)
                                + dot(b, bp)*dot(a, ap)/(4*B*Bp))
    K2c = 1/D*(Ap/Bp*lp-A/B*l)*(dot(b, bp)*dot(c, cp)/(C*Cp)-dot(b, bp)*dot(ap, c)/(2*C*Ap)-dot(a, cp)*dot(b, bp)/(2*A*Cp)
                                + dot(b, bp)*dot(a, ap)/(4*A*Ap))
    K2d = 1/D*(B*l/(B-A)+Bp*lp/(Ap-Bp))*(dot(a, c)*(Bp-Ap)/(2*C*Cp)+dot(a, b)*(Bp-Ap)/(4*B*Cp)+dot(ap, cp)*(B-A)/(2*C*Cp)
                                         + dot(ap, bp)*(B-A)/(4*C*Bp))
    K2e = -1/D*(A*l/(B-A)+Ap*lp/(Ap-Bp))*(dot(b, c)*(Bp-Ap)/(2*C*Cp)-dot(a, b)*(Bp-Ap)/(4*A*Cp)+dot(bp, cp)*(B-A)/(2*C*Cp)
                                           - dot(ap, bp)*(B-A)/(4*C*Ap))
    K2f = 1/(2*C*Cp)*(-A*B*(Bp-Ap)/((B-A)*D)*l-1+Ap*Bp*(B-A)/((Bp-Ap)*D)*lp)
    return K2a+K2b+K2c+K2d+K2e+K2f
rel = 0
for _ in range(200):
    x, xp, z, w, wp = [rng.normal(size=2) for _ in range(5)]
    a, b = II1_author(x, xp, z, w, wp)[0], II1_paper(x, xp, z, w, wp)
    rel = max(rel, abs(a-b)/abs(b))
pr('(5) Group II Row I 2nd part: author N2a+..+N2e vs paper K2a+..+K2f: max rel diff = %.1e (200 points)' % rel)
# pole on D=0: approach the surface
x, xp, z, w = [rng.normal(size=2) for _ in range(4)]
# choose wp such that D = A'B - AB' = 0 + eps: solve along a ray for |xp-wp|^2
bp = xp-z; A = sq(x-w); B = sq(x-z); Bp = sq(bp)
Ap0 = A*Bp/B
u = rng.normal(size=2); u /= np.linalg.norm(u)
for eps in [1e-3, 1e-5]:
    vals = []
    for sgn in [1, -1]:
        wp = xp-u*np.sqrt(Ap0*(1+sgn*eps))
        tot, parts = II1_author(x, xp, z, w, wp)
        vals.append((tot, parts[0]+parts[1], parts[2]+parts[3], parts[0]))
    pr('(5) near D=0 (rel. distance %g): total %.6f / %.6f ; N2a alone %.3e / %.3e ; N2a+N2b %.4f / %.4f'
       % (eps, vals[0][0], vals[1][0], vals[0][3], vals[1][3], vals[0][1], vals[1][1]))

# ---------- (7) Group II Row IV substitution ----------
g = lambda u: np.exp(-0.3*u)*np.cos(u)
for T in [50.0]:
    lhs = integrate.quad(lambda p: p/(p+1)**2*g(1+p), 0, T, limit=500)[0]
    rhs = integrate.quad(lambda xi: (1-xi)/xi*g(1/xi), 1/(1+T), 1, limit=500)[0]
    pr('(7) Row II.IV: int_0^T dp+ p+/(p++k+)^2 g(1+p+/k+) - int_{1/(1+T)}^1 dxi (1-xi)/xi g(1/xi) = %.1e' % (lhs-rhs))
open('s2_kernels.log', 'w').write('\n'.join(out)+'\n')
