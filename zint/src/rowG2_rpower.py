"""Settle the power of (y'-y)^2 in the second term of G^{mm'} by direct 2d quadrature."""
import numpy as np
def G_cut(r,R,n=3000,m=1024):
    """int_{|u|<R} d^2u  u^m (u+r)^m' / [u^2 (u+r)^2]"""
    lo,hi=np.log(1e-8*np.linalg.norm(r)),np.log(R)
    xg,wg=np.polynomial.legendre.leggauss(n)
    lr=0.5*(hi-lo)*xg+0.5*(hi+lo); wr=0.5*(hi-lo)*wg
    rho=np.exp(lr); ph=(np.arange(m)+0.5)*(2*np.pi/m)
    ux=np.outer(rho,np.cos(ph)); uy=np.outer(rho,np.sin(ph))
    vx=ux+r[0]; vy=uy+r[1]
    den=(ux**2+uy**2)*(vx**2+vy**2)
    w=(wr*rho*rho)[:,None]*(2*np.pi/m); U=[ux,uy]; V=[vx,vy]
    return np.array([[float((U[a]*V[b]/den*w).sum()) for b in range(2)] for a in range(2)])
print("cutoff master:  G^{mm'} = (pi/2)(log(R^2/r^2)+1) delta^{mm'} - pi r^m r^m' / r^POWER")
print("  (at |r| = 1 the two candidates coincide -- use |r| != 1 to tell them apart)")
print("  %-18s %14s %14s"%("r","maxdiff, 1/r^2","maxdiff, 1/r^4"))
dd=np.eye(2); R=400.0
for rv in ([1.0,0.0],[0.6,-0.8],[2.0,1.5],[0.4,0.2],[3.0,-1.0]):
    r=np.array(rv); r2=r@r; num=G_cut(r,R)
    base=np.pi/2*(np.log(R**2/r2)+1)*dd
    d2=np.abs(num-(base-np.pi*np.outer(r,r)/r2)).max()
    d4=np.abs(num-(base-np.pi*np.outer(r,r)/r2**2)).max()
    print("  %-18s %14.2e %14.2e"%(rv,d2,d4))
print()
print("  => the second term is  - pi r^m r^m' / r^2  =  - pi rhat^m rhat^m' ,")
print("     one power of (y'-y)^2 in the denominator, not two.")
print("     (dimensions agree: G is dimensionless at eps=0, and so is rhat rhat.)")
