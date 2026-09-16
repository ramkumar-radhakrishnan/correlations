"""Row 8: the large-distance end, done with enough angular resolution to be meaningful."""
import numpy as np
from scipy.linalg import expm
rng = np.random.default_rng(5)
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
Tg=lam/2; A_=rng.normal(size=(8,2)); B_=rng.normal(size=8)
def Uadj(p):
    al=np.exp(-0.25*(p@p))*(A_@p+B_); H=sum(al[a]*Tg[a] for a in range(8)); V=expm(1j*H)
    return np.array([[(2*np.trace(Tg[a]@V@Tg[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
kt=np.array([0.7,-0.4]); Kp=1.0; s=0.37; S=Kp+s
def Kv(p,q): d=p-q; return d/(d@d)
def F(pt):
    x,xp,y,yp,z=pt
    m,mp_,kk,kp_=Kv(y,z),Kv(yp,z),Kv(x,y),Kv(xp,yp)
    t1=2*s/S**2*(kk@m)*(kp_@mp_)
    t2=2/Kp*((m@kp_)*(mp_@kk)-(mp_@kp_)*(kk@m))
    t3=(s/Kp**2+1/s)*(m@mp_)*(kk@kp_)
    col=np.einsum('Ec,ec->Ee',Uadj(yp)-Uadj(xp),Uadj(y)-Uadj(x))
    return np.exp(-1j*(kt@(yp-y))*S/Kp)*(t1+t2+t3)*col[0,0]
names=['x',"x'",'y',"y'",'z']
base=[np.array([0.31,-0.77]),np.array([-0.52,0.19]),np.array([1.13,0.42]),
      np.array([-0.31,0.88]),np.array([0.05,-0.23])]

print("="*78)
print("LARGE-DISTANCE END.  angular average A(R), resolution scaled to the oscillation.")
print("The radial integral converges iff |A(R)| R^2 falls faster than R^0.")
print("="*78)
print("%-5s %9s %13s %13s %13s %13s   %s" % ("var","R:","1e2","1e3","1e4","1e5","fitted power of |A|R^2"))
for i,nm in enumerate(names):
    vals=[]
    for R in (1e2,1e3,1e4,1e5):
        # resolve the phase: need >> k*R/(2 pi) periods; z carries no phase so few points suffice
        osc = 0 if i==4 else abs(kt)@np.ones(2)*R*S/Kp
        nth = int(min(max(4096, 40*osc), 4_000_000))
        th=np.linspace(0,2*np.pi,nth,endpoint=False)
        pts=R*np.stack([np.cos(th),np.sin(th)],-1)
        acc=0j
        for c in range(0,nth,20000):
            blk=pts[c:c+20000]
            for p in blk:
                pt=[b.copy() for b in base]; pt[i]=p
                acc+=F(pt)
        vals.append(abs(acc/nth)*R**2)
    p=np.polyfit(np.log([1e2,1e3,1e4,1e5]),np.log(vals),1)[0]
    print("%-5s %9s %13.5e %13.5e %13.5e %13.5e   %+.2f  %s"
          % (nm,"",*vals,p,"CONVERGENT" if p<-0.2 else "*** LOG DIVERGENT ***"))
