"""Row 8: exhaustive UV sweep. Every collapsing subset of the integration points,
on (a) the transverse integrand, (b) the z-integrated form, (c) the p+ integrated form."""
import numpy as np, itertools, mpmath as mp
mp.mp.dps = 25

Kp = 1.0; Pv = 30.0; Lam = 1e-8; kt = np.array([0.7,-0.4]); Rcut = 1e4

def K(a,b):
    d = a-b; return d/(d@d)

# ---------------------------------------------------------------- (a) transverse integrand
def F_transverse(pt, s=0.37):
    x,xp,y,yp,z = pt
    m,mp_,kk,kp_ = K(y,z),K(yp,z),K(x,y),K(xp,yp)
    P = s; S = Kp+P; dperp = 2.0
    t1 = dperp*P/S**2 * (kk@m)*(kp_@mp_)
    t2 = 2/Kp * ((m@kp_)*(mp_@kk) - (mp_@kp_)*(kk@m))
    t3 = (P/Kp**2 + 1/P) * (m@mp_)*(kk@kp_)
    return np.exp(-1j*(kt@(yp-y))*(1+s))*(t1+t2+t3)

# ---------------------------------------------------------------- (b)/(c) after z, after p+
def I1(k):
    a = abs(k)*Pv
    Cin = float(mp.euler + mp.log(a) - mp.ci(a)) if a > 1e-12 else a*a/4
    Si  = float(mp.si(a)) if a > 1e-12 else a
    return np.log(Pv/Lam) - Cin - 1j*np.sign(k)*Si
def I2(k):
    if abs(k)*Pv < 1e-8: return Pv**2/2 + 0j
    return complex((1 - mp.e**(-1j*k*Pv)*(1+1j*k*Pv))/(-k**2))
def I3(k):
    if abs(k) < 1e-10:
        return complex(mp.quad(lambda t: t/(Kp+t)**2, [0, Pv]))
    V = Kp + Pv
    return complex(mp.e**(1j*k*Kp)*((1+1j*k*Kp)*(mp.ei(-1j*k*V)-mp.ei(-1j*k*Kp))
                                    + Kp*(mp.e**(-1j*k*V)/V - mp.e**(-1j*k*Kp)/Kp)))

def F_zdone(pt, s=0.37):
    """z already integrated (cutoff R), p+ NOT yet."""
    x,xp,y,yp,_ = pt
    r = yp-y; r2 = r@r
    X = K(x,y)@K(xp,yp)
    R_ = (r@K(x,y))*(r@K(xp,yp))/r2
    L = np.log(Rcut**2/r2) + 1.0
    P = s
    val = L*X/P + P/Kp**2*L*X + P/(Kp+P)**2*((L+1)*X - 2*R_)
    return np.exp(-1j*(kt@r))*np.exp(-1j*(kt@r)*P/Kp)*val

def F_full(pt):
    """z and p+ both integrated."""
    x,xp,y,yp,_ = pt
    r = yp-y; r2 = r@r; kap = (kt@r)/Kp
    X = K(x,y)@K(xp,yp)
    R_ = (r@K(x,y))*(r@K(xp,yp))/r2
    L = np.log(Rcut**2/r2) + 1.0
    return np.exp(-1j*(kt@r))*(L*X*(I1(kap)+I2(kap)/Kp**2) + ((L+1)*X - 2*R_)*I3(kap))

# ---------------------------------------------------------------- the sweep
names = ['x',"x'",'y',"y'",'z']
base  = [np.array([0.31,-0.77]), np.array([-0.52,0.19]), np.array([1.13,0.42]),
         np.array([-0.31,0.88]), np.array([0.05,-0.23])]
u = lambda t: np.array([np.cos(t), np.sin(t)])
OFF = [0.0, 2.1, 4.0, 1.05, 5.2]          # distinct angular offsets, no accidental coincidence

def sweep(F, idxs, rhos=(1e-2,1e-4,1e-6), nth=32):
    """collapse the points in idxs onto base[idxs[0]] at scale rho."""
    n = len(idxs); out = []
    for rho in rhos:
        acc = 0.0
        for j,th in enumerate(np.linspace(0, 2*np.pi, nth, endpoint=False)):
            pt = [b.copy() for b in base]
            c  = base[idxs[0]]
            for a,i in enumerate(idxs):
                pt[i] = c + rho*u(th + OFF[a])
            try: acc += abs(F(pt))
            except Exception: return None
        out.append(acc/nth * rho**(2*(n-1)))
    return out

def verdict(v):
    if v is None or not np.isfinite(v).all(): return "n/a"
    a,b,c = v
    if a == 0: return "zero"
    r1, r2 = b/a, c/b
    if r1 > 3 or r2 > 3:    return "*** POWER DIVERGENT ***"
    if 0.3 < r1 < 3 and 0.3 < r2 < 3: return "*** LOG DIVERGENT ***"
    return "convergent"

for label, F, pts in [("(a) TRANSVERSE INTEGRAND  (all 5 points)", F_transverse, [0,1,2,3,4]),
                      ("(b) AFTER the z integral  (x,x',y,y')",    F_zdone,      [0,1,2,3]),
                      ("(c) AFTER z AND p+        (x,x',y,y')",    F_full,       [0,1,2,3])]:
    print("="*78); print(label); print("="*78)
    print("%-22s %11s %11s %11s   %s" % ("collapsing subset","rho=1e-2","1e-4","1e-6","verdict"))
    flag = False
    for n in range(2, len(pts)+1):
        for sub in itertools.combinations(pts, n):
            v = sweep(F, list(sub))
            vd = verdict(v)
            if vd.startswith("***"): flag = True
            s = "{"+",".join(names[i] for i in sub)+"}"
            if v is None: print("%-22s %s" % (s, "n/a")); continue
            print("%-22s %11.3e %11.3e %11.3e   %s" % (s, *v, vd))
    print("  -->", "UV DIVERGENCE FOUND" if flag else "NO UV DIVERGENCE IN ANY REGION")
    print()
