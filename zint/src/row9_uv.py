"""Row 9 (|B_2|^2): exhaustive transverse UV scan of N1 (3 rho) and N2 (2 rho)."""
import numpy as np, itertools
P, K = 0.37, 1.0; S = P+K
kt = np.array([0.7,-0.4])

def Wji(x, z, w):
    """the B_2 bracket as a matrix W[j,i] (S-multiplied form)."""
    xz, xw, zw = x-z, x-w, z-w
    s = lambda v: v@v
    t1 = np.eye(2)*(s(xz)-s(xw))/(2*s(zw))
    t2 = (S/K)*( np.outer(xz, zw)/s(zw) - np.outer(xz, xw)/(2*s(xw)) )
    t3 = (S/P)*( np.outer(zw, xw)/s(zw) + np.outer(xz, xw)/(2*s(xz)) )
    return t1 + t2 + t3
def Dden(x, z, w):
    return P*((x-z)@(x-z)) + K*((x-w)@(x-w))
def Kv(a, b):
    d = a-b; return d/(d@d)

def N1_int(pt):
    """pt = (x, y, z, w, xp, wp)"""
    x,y,z,w,xp,wp = pt
    W = Wji(x,z,w)/S
    val = np.einsum('ji,j,i->', W, Kv(y,z), Kv(xp,wp))/Dden(x,z,w)
    return np.exp(-1j*(kt@(wp-w)))*val
def N2_int(pt):
    """pt = (x, xp, z, w, wp)"""
    x,xp,z,w,wp = pt
    val = np.einsum('ji,ji->', Wji(x,z,w), Wji(xp,z,wp))
    return np.exp(-1j*(kt@(wp-w)))*val*(P*K/S**2)/(Dden(x,z,w)*Dden(xp,z,wp))

u  = lambda t: np.array([np.cos(t), np.sin(t)])
OFF = [0.0, 2.1, 4.0, 1.05, 5.2, 3.3]

def sweep(F, base, names, rhos=(1e-2,1e-4,1e-6), nth=28):
    n_pts = len(base)
    print("%-26s %11s %11s %11s   %s" % ("collapsing subset","rho=1e-2","1e-4","1e-6","verdict"))
    found = []
    for n in range(2, n_pts+1):
        for sub in itertools.combinations(range(n_pts), n):
            vals = []
            for rho in rhos:
                acc = 0.0
                for th in np.linspace(0, 2*np.pi, nth, endpoint=False):
                    pt = [b.copy() for b in base]; c = base[sub[0]]
                    for a,i in enumerate(sub): pt[i] = c + rho*u(th+OFF[a])
                    try: acc += abs(F(pt))
                    except Exception: acc = np.nan
                vals.append(acc/nth * rho**(2*(n-1)))
            if not np.isfinite(vals).all(): continue
            r1, r2 = vals[1]/vals[0], vals[2]/vals[1]
            if r1 > 3 or r2 > 3:      v = "*** POWER DIVERGENT ***"
            elif r1 > 0.3 and r2 > 0.3: v = "*** LOG DIVERGENT ***"
            else:                      v = "convergent"
            if v.startswith("*"): found.append(sub)
            lab = "{"+",".join(names[i] for i in sub)+"}"
            if v.startswith("*") or n <= 3:
                print("%-26s %11.3e %11.3e %11.3e   %s" % (lab, *vals, v))
    return found

print("="*80); print("N2  (two rho)   points: x, x', z, w, w'"); print("="*80)
b2 = [np.array([0.31,-0.77]), np.array([-0.52,0.19]), np.array([0.05,-0.23]),
      np.array([1.13,0.42]),  np.array([-0.31,0.88])]
f2 = sweep(N2_int, b2, ['x',"x'",'z','w',"w'"])
print("  ==>", ("UV FOUND: "+", ".join("{"+",".join(['x',"x'",'z','w',"w'"][i] for i in s)+"}"
        for s in f2)) if f2 else "no UV")

print(); print("="*80); print("N1  (three rho)  points: x, y, z, w, x', w'"); print("="*80)
b1 = [np.array([0.31,-0.77]), np.array([0.62,0.41]), np.array([0.05,-0.23]),
      np.array([1.13,0.42]),  np.array([-0.52,0.19]), np.array([-0.31,0.88])]
f1 = sweep(N1_int, b1, ['x','y','z','w',"x'","w'"])
print("  ==>", ("UV FOUND: "+", ".join("{"+",".join(['x','y','z','w',"x'","w'"][i] for i in s)+"}"
        for s in f1)) if f1 else "no UV")
