"""Row 11: transverse scan of the five terms, and the p+ pole structure."""
import numpy as np, itertools
K = 1.0
def build(P):
    S = P+K
    def Bt(i,j,kk,m):
        d = lambda a,b: 1.0 if a==b else 0.0
        return d(i,kk)*d(j,m) - (K/P)*d(j,kk)*d(i,m) - (K/(K-P))*d(i,j)*d(kk,m)
    B = np.array([[[[Bt(i,j,kk,m) for m in range(2)] for kk in range(2)]
                   for j in range(2)] for i in range(2)])
    def terms(x,y,z,xp,w,wp):
        s = lambda v: v@v
        Kw  = (xp-wp)/s(xp-wp); Kyz = (y-z)/s(y-z)
        A   = s(x-z); Bq = s(x-w)
        Dm  = P*A - K*Bq                      # p+(x-z)^2 - k+(x-w)^2
        W   = K*Bq - P*A                      # = -Dm
        V   = K*(x-w) - P*(x-z)
        zw  = z-w
        t_I   = (Kw@Kyz) * S/(K-P)**2 / Dm
        t_III = -(Kw@(x-w)/Bq)*((x-z)@Kyz/A) * (K*Bq + P*A)/Dm /(P*K)
        core  = np.einsum('ijkm,i,j,k,m->', B, Kw, Kyz, zw, V/W)
        t_IIa = core/(K*(K-P))/s(zw)
        t_IIb = -P/(K-P)*np.einsum('ijkm,i,j,k,m->', B, Kw, Kyz, zw, V/W)/s(V)
        return np.array([t_I, t_III, t_IIa, t_IIb])
    return terms

names=['x','y','z',"x'",'w',"w'"]
base=[np.array([0.31,-0.77]),np.array([0.62,0.41]),np.array([0.05,-0.23]),
      np.array([-0.52,0.19]),np.array([1.13,0.42]),np.array([-0.31,0.88])]
u=lambda t:np.array([np.cos(t),np.sin(t)]); OFF=[0,2.1,4.0,1.05,5.2,3.3]
F = build(0.37)
lbl=["I (terms 1,2)","III (term 3)","IIa (term 4)","IIb (term 5)"]
print("="*80); print("TRANSVERSE COLLAPSE SCAN  (integrand x rho^{2(n-1)}); * = divergent"); print("="*80)
print("%-18s %s" % ("subset", "  ".join("%-16s"%l for l in lbl)))
bad=[]
for n in range(2,7):
    for sub in itertools.combinations(range(6),n):
        vals=[]
        for rho in (1e-2,1e-4,1e-6):
            acc=np.zeros(4)
            for th in np.linspace(0,2*np.pi,24,endpoint=False):
                pt=[b.copy() for b in base]; c=base[sub[0]]
                for a,i in enumerate(sub): pt[i]=c+rho*u(th+OFF[a])
                acc+=np.abs(F(*pt))
            vals.append(acc/24*rho**(2*(n-1)))
        vals=np.array(vals)
        r1=vals[1]/np.maximum(vals[0],1e-300); r2=vals[2]/np.maximum(vals[1],1e-300)
        flag=(r1>0.3)&(r2>0.3)
        if flag.any():
            bad.append((sub,flag))
            print("%-18s %s   <-- %s" % ("{"+",".join(names[i] for i in sub)+"}",
                  "  ".join("%-16.4e"%v for v in vals[2]),
                  ", ".join(lbl[i] for i in range(4) if flag[i])))
print("  divergent subsets:", "none" if not bad else len(bad))

print(); print("="*80); print("p+ SINGULARITIES  (the real story in this row)"); print("="*80)
x,y,z,xp,w,wp = base
A = (x-z)@(x-z); Bq=(x-w)@(x-w)
print("   (a) p+ = k+   :  double pole 1/(k+-p+)^2 in terms 1,2,4,5.")
print("       residues at p+ -> k+, coefficient of 1/(k+-p+)^2, in units of")
print("          [(x'-w').(y-z)/((x'-w')^2(y-z)^2)] / [(x-w)^2-(x-z)^2] :")
for eps in (1e-3,1e-4,1e-5):
    for sgn,tag in ((+1,"p+ > k+"),(-1,"p+ < k+")):
        Fe = build(K+sgn*eps); t = Fe(x,y,z,xp,w,wp)
        pref = np.array([1/(8*np.pi**4), -1/(8*np.pi**4), 1/(4*np.pi**4), -1/(4*np.pi**4)])
        norm = eps**2/(((xp-wp)/((xp-wp)@(xp-wp)))@((y-z)/((y-z)@(y-z))))*((Bq-A))
        r = t*pref*norm*(8*np.pi**4)
        if eps==1e-5:
            print("       eps=%.0e %s : I %8.4f   IIa %8.4f   IIb %8.4f   (III %8.4f)"
                  % (eps,tag,r[0],r[2],r[3],r[1]))
print()
print("   (b) W = 0  i.e.  p+ (x-z)^2 = k+ (x-w)^2  ->  p+_* = k+ (x-w)^2/(x-z)^2 = %.6f k+"
      % (Bq/A))
print("       this lies INSIDE [Lambda, V] for generic transverse points, in EVERY term.")
print("       It is a simple pole ON the integration contour.")
print("   (c) p+ = 0 : term 3 has 1/(p+ k+); term 4's bracket has -k+/p+ .")
print("   (d) V = 0  i.e. k+(x-w) = p+(x-z) as VECTORS, at x_* = (k+ w - p+ z)/(k+-p+):")
print("       term 5 has 1/V^2 there, but its numerator V^m vanishes linearly,")
print("       so the integrand is ~1/rho against rho drho -- integrable.")
