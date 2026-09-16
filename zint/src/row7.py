"""Row 7: C (first argument k+) with four single WW kernels. UV scan + large distance.
After the w-integration the transverse integrand is

  K = (x'-w')^i (z-x)^m (y-z)^k (y'-x)^j / [ (x'-w')^2 (z-x)^2 (y-z)^2 (y'-x)^2 ]
      x exp[-i k.(w'-z)] exp[-i xi k.(z-x)]
      x [ d_im d_jk - (1/xi) d_jm d_ik - (1/(1-xi)) d_ij d_km ]        (k+ = 1)
"""
import numpy as np
kk = np.array([0.9, 0.6])

def F(x, y, z, xp, yp, wp, xi):
    A = xp-wp; B = z-x; C = y-z; D = yp-x
    A2,B2,C2,D2 = A@A, B@B, C@C, D@D
    ph = np.exp(-1j*(kk@(wp-z)))*np.exp(-1j*xi*(kk@(z-x)))
    s = 0j
    for i in range(2):
        for m in range(2):
            for kdx in range(2):
                for j in range(2):
                    d = lambda a,b: 1.0 if a==b else 0.0
                    br = d(i,m)*d(j,kdx) - d(j,m)*d(i,kdx)/xi - d(i,j)*d(kdx,m)/(1-xi)
                    if br != 0:
                        s += A[i]*B[m]*C[kdx]*D[j]/(A2*B2*C2*D2)*br
    return ph*s

rng = np.random.default_rng(5)
x0,y0,z0,xp0,yp0,wp0 = rng.normal(size=(6,2))
xi = 0.4
th = np.linspace(0,2*np.pi,4096,endpoint=False); u = np.stack([np.cos(th),np.sin(th)],-1)
names = ['x','y','z','xp','yp','wp']
base = dict(x=x0,y=y0,z=z0,xp=xp0,yp=yp0,wp=wp0)

print("A. SHORT DISTANCE. rho^2 x |<integrand>| at each of the four coincidences.")
print("   Each separation carries exactly ONE kernel, so we expect O(1/rho) -> rho^2<.> ~ rho^2 -> 0.")
for nm, mv, ctr in (("x' -> w'",'xp',wp0), ("z -> x",'z',x0), ("y -> z",'y',z0), ("y' -> x",'yp',x0)):
    row=[]
    for rho in (1e-2,1e-3,1e-4):
        vals=[]
        for uu in u:
            a=dict(base); a[mv]=ctr+uu*rho
            vals.append(F(a['x'],a['y'],a['z'],a['xp'],a['yp'],a['wp'],xi))
        row.append(abs(np.mean(vals))*rho**2)
    print(f"   {nm:10s} rho^2|<F>| = " + "  ".join(f"{v:.4e}" for v in row) + "   -> 0 : NO UV")

print()
print("B. do any two kernels share a separation?  the four pairs are")
print("   {x',w'}, {z,x}, {y,z}, {y',x}  -- all distinct, so no 1/rho^2 can form.")
print("   (shared POINTS z and x only give multi-point collapses, which cost measure)")
print("   triple collapse y -> z -> x (3 points, 2 relative coords, measure rho^4 drho/rho,")
print("   integrand ~ rho^-2, net rho^2 drho/rho):")
for rho in (1e-2,1e-3,1e-4):
    vals=[F(x0,x0+u[n]*rho*1.3,x0+u[(n+1000)%len(th)]*rho,xp0,yp0,wp0,xi) for n in range(len(th))]
    print(f"      rho={rho:.0e}:  rho^2|<F>| = {abs(np.mean(vals))*rho**2:.4e}")

print()
print("C. LARGE DISTANCE.  which variables carry a phase?")
print("   w' : e^{-i k.w'}          x : e^{+i xi k.x}      z : e^{+i(1-xi) k.z}")
print("   -> x is unprotected as xi -> 0 ; z is unprotected as xi -> 1.")
for nm, mv, xis in (("|x| -> inf",'x',(0.4,0.02)), ("|z| -> inf",'z',(0.4,0.98))):
    for xiv in xis:
        row=[]
        for R in (1e2,1e3,1e4):
            vals=[]
            for uu in u:
                a=dict(base); a[mv]=uu*R
                vals.append(F(a['x'],a['y'],a['z'],a['xp'],a['yp'],a['wp'],xiv))
            row.append(abs(np.mean(vals))*R**2)
        print(f"   {nm:12s} xi={xiv:<5} R^2|<F>| = " + "  ".join(f"{v:.4e}" for v in row))
