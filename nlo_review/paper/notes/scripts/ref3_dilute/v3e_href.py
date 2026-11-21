import sys; sys.path.insert(0,'/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import numpy as np, v3c_fast as f
from scipy.integrate import quad
def href(p,k,l,Lam,a=0.5):
    """reference: PV by symmetric pairing around s0, 1/|s-1| subtracted around s=1, all with adaptive quad."""
    k2=k@k; s0=(p@p)/k2
    F=lambda s: f.Fvec(np.array([s]),p,k,l)[0]
    e=1e-8; c=0.5*(F(1+e)+F(1-e))*e
    Ainf=F(1e10)*1e10
    pts=[Lam,1-a,1-Lam,1+Lam,1+a]
    S=max(20.0,8*s0)
    def integ(g,x0,x1,extra=()):
        brk=sorted(set([x0,x1]+[x for x in extra if x0<x<x1]))
        return sum(quad(g,u,v,limit=500,epsabs=1e-13,epsrel=1e-11)[0] for u,v in zip(brk[:-1],brk[1:]))
    s1=np.sqrt(p@p/k2); sst=(p@k)/k2
    extra=[s1,sst,s1/2,2*s1,s0/2,2*s0]
    tot=0.0
    # neighbourhood of s0 (PV) unless near 1 or Lam
    d=min(0.25*abs(s0-1),0.5*(s0-Lam) if s0>Lam else 0,0.1*s0)
    regs=[(Lam,1-a),(1+a,S)]
    def pvsum(x0,x1):
        if d>0 and x0<s0-d and s0+d<x1:
            v=integ(lambda s:F(s),x0,s0-d,extra)+integ(lambda s:F(s),s0+d,x1,extra)
            v+=quad(lambda t:F(s0+t)+F(s0-t),0,d,limit=500,epsabs=1e-13,epsrel=1e-11)[0]
            return v
        return integ(lambda s:F(s),x0,x1,extra)
    tot+=pvsum(Lam,1-a)
    fm=lambda s:F(s)-c/abs(s-1)
    # near 1, s0 possibly inside
    def pvsum_m(x0,x1):
        if d>0 and x0<s0-d and s0+d<x1:
            v=integ(fm,x0,s0-d,extra)+integ(fm,s0+d,x1,extra)
            v+=quad(lambda t:fm(s0+t)+fm(s0-t),0,d,limit=500,epsabs=1e-13,epsrel=1e-11)[0]
            return v
        return integ(fm,x0,x1,extra)
    tot+=pvsum_m(1-a,1-1e-13)+pvsum_m(1+1e-13,1+a)+2*c*np.log(a/Lam)
    tot+=pvsum(1+a,S)
    tot+=quad(lambda s:F(s)-Ainf/s,S,np.inf,limit=500)[0]-Ainf*np.log(S)
    return tot
if __name__=='__main__':
    pass
    k=np.array([1.,0.]); l=np.array([.3,.7])
    for p,Lam in ((np.array([0.59149491,0.58847801]),1e-5),(np.array([-1.15011636,-0.17971109]),3e-6),(np.array([-0.73354855,-0.27733323]),3e-6),(np.array([1.00001034e+00,-4.61054906e-06]),1e-5),(np.array([4.31539186e-06,-2.23393508e-08]),3e-6)):
        print(p, 'ref',href(p,k,l,Lam),'fast',f.h(p,k,l,Lam)[0],'quad',pr.hfull(p,Lam)[0],flush=True)
