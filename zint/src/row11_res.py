"""Row 11: the p+ = k+ double-pole residues, and the large-distance transverse end."""
import numpy as np
K=1.0
def terms(P, x,y,z,xp,w,wp):
    """the four structures WITH their numerical prefactors (g^4 f/pi^4 stripped)."""
    S=P+K; s=lambda v: v@v
    d=lambda a,b: 1.0 if a==b else 0.0
    B=np.array([[[[d(i,kk)*d(j,m)-(K/P)*d(j,kk)*d(i,m)-(K/(K-P))*d(i,j)*d(kk,m)
                   for m in range(2)] for kk in range(2)] for j in range(2)] for i in range(2)])
    Kw=(xp-wp)/s(xp-wp); Kyz=(y-z)/s(y-z)
    A=s(x-z); Bq=s(x-w); Dm=P*A-K*Bq; W=-Dm; V=K*(x-w)-P*(x-z); zw=z-w
    core=np.einsum('ijkm,i,j,k,m->',B,Kw,Kyz,zw,V/W)
    return np.array([
        ( 1/(8*np.pi**4))*(Kw@Kyz)*S/(K-P)**2/Dm,                                   # I
        (-1/(8*np.pi**4))*(Kw@(x-w)/Bq)*((x-z)@Kyz/A)*(K*Bq+P*A)/Dm/(P*K),          # III
        ( 1/(4*np.pi**4))*core/(K*(K-P))/s(zw),                                     # IIa
        (-1/(4*np.pi**4))*(P/(K-P))*np.einsum('ijkm,i,j,k,m->',B,Kw,Kyz,zw,V/W)/s(V)# IIb
    ])
x=np.array([0.31,-0.77]); y=np.array([0.62,0.41]); z=np.array([0.05,-0.23])
xp=np.array([-0.52,0.19]); w=np.array([1.13,0.42]); wp=np.array([-0.31,0.88])
print("="*78); print("DOUBLE POLE AT p+ = k+ :  lim (k+-p+)^2 x term"); print("="*78)
print("  %-8s %13s %13s %13s %13s" % ("eps","I","III","IIa","IIb"))
for eps in (1e-3,1e-4,1e-5,1e-6):
    for sgn,tag in ((+1,">"),(-1,"<")):
        r = eps**2*terms(K+sgn*eps,x,y,z,xp,w,wp)
        print("  %-8s %13.8f %13.8f %13.8f %13.8f" % (("%+.0e"%(sgn*eps)),*r))
Kw=(xp-wp)/((xp-wp)@(xp-wp)); Kyz=(y-z)/((y-z)@(y-z))
A=(x-z)@(x-z); Bq=(x-w)@(x-w)
pred = (Kw@Kyz)/(4*np.pi**4*(A-Bq))
print()
print("  predicted residue for I :  K/(4 pi^4 (A-B)) = %13.8f   with A=(x-z)^2=%.4f, B=(x-w)^2=%.4f"
      % (pred,A,Bq))
r=1e-6**2*terms(K+1e-6,x,y,z,xp,w,wp)
print("  sum of the four residues  = %13.8f" % r.sum())
print("  I + IIa                   = %13.8f      (both present for p+ < k+ in your ranges)" % (r[0]+r[2]))
print("  I + IIb                   = %13.8f      (both present for p+ > k+ in your ranges)" % (r[0]+r[3]))
print()
print("  => IIa and IIb have EQUAL AND OPPOSITE residues; whichever one accompanies I")
print("     decides whether the double pole cancels.  With your ranges it cancels ABOVE")
print("     k+ (terms 1 and 5) and does NOT cancel BELOW (terms 2 and 4).")

print(); print("="*78); print("LARGE-DISTANCE TRANSVERSE END  (|A(R)| R^2 flat = log divergent)"); print("="*78)
base=[x,y,z,xp,w,wp]; names=['x','y','z',"x'",'w',"w'"]
print("  %-5s %s" % ("var", "  ".join("%-14s"%l for l in ["I","III","IIa","IIb"])))
for i,nm in enumerate(names):
    out=[]
    for R in (1e2,1e4,1e6):
        acc=np.zeros(4)
        for th in np.linspace(0,2*np.pi,512,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*np.array([np.cos(th),np.sin(th)])
            acc+=terms(0.37,*pt)
        out.append(np.abs(acc/512)*R**2)
    out=np.array(out)
    verd=["log-div" if (out[2,j]>0.2*out[0,j]) else "conv" for j in range(4)]
    print("  %-5s %s   %s" % (nm, "  ".join("%-14.4e"%v for v in out[2]), ", ".join(verd)))
