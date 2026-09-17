"""Row 11: large-distance transverse end of the BARE integrand (phase set to 1)."""
import numpy as np
K=1.0; P=0.37
d=lambda a,b:1.0 if a==b else 0.0
B4=np.array([[[[d(i,kk)*d(j,m)-(K/P)*d(j,kk)*d(i,m)-(K/(K-P))*d(i,j)*d(kk,m)
    for m in range(2)] for kk in range(2)] for j in range(2)] for i in range(2)])
def terms(x,y,z,xp,w,wp):
    s=lambda v:v@v
    Kw=(xp-wp)/s(xp-wp); Kyz=(y-z)/s(y-z)
    A=s(x-z); Bq=s(x-w); Dm=P*A-K*Bq; W=-Dm; V=K*(x-w)-P*(x-z); zw=z-w
    core=np.einsum('ijkm,i,j,k,m->',B4,Kw,Kyz,zw,V/W)
    return np.array([ (Kw@Kyz)*(P+K)/(K-P)**2/Dm,
        -(Kw@(x-w)/Bq)*((x-z)@Kyz/A)*(K*Bq+P*A)/Dm/(P*K),
        core/(K*(K-P))/s(zw),
        -(P/(K-P))*core/s(V) ])
base=[np.array([0.31,-0.77]),np.array([0.62,0.41]),np.array([0.05,-0.23]),
      np.array([-0.52,0.19]),np.array([1.13,0.42]),np.array([-0.31,0.88])]
names=['x','y','z',"x'",'w',"w'"]; lbl=["I","III","IIa","IIb"]
scale=np.abs(terms(*base))
print("="*78); print("LARGE-DISTANCE END of the BARE integrand (phase set to 1)"); print("="*78)
print("  angular average x R^2;  flat = log divergent.  noise floor ~ 1e-8 x scale")
print("  %-4s %-6s %13s %13s %13s   %s"%("var","term","R=1e2","R=1e4","R=1e6","verdict"))
for i,nm in enumerate(names):
    rows=[]
    for R in (1e2,1e4,1e6):
        acc=np.zeros(4)
        for th in np.linspace(0,2*np.pi,4096,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*np.array([np.cos(th),np.sin(th)])
            acc+=terms(*pt)
        rows.append(acc/4096*R**2)
    rows=np.array(rows)
    for j in range(4):
        v=abs(rows[2,j])
        flat = v>1e-6*scale[j] and v>0.2*abs(rows[0,j])
        print("  %-4s %-6s %13.4e %13.4e %13.4e   %s"
              %(nm if j==0 else "", lbl[j], rows[0,j],rows[1,j],rows[2,j],
                "*** LOG DIVERGENT ***" if flat else "convergent"))
