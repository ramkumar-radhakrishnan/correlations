"""Row 25 (vectorised): the large-|z| log, the phase cutoff, and p+ -> 0."""
import numpy as np
s=lambda v:v@v; KT=np.array([0.6,-0.35]); KP=1.0; kabs=np.sqrt(KT@KT)
x =np.array([0.41,-0.63]); xp=np.array([-0.58,0.22]); yp=np.array([0.09,0.51]); w=np.array([-0.25,0.86])
X =(xp-yp)/s(xp-yp)
def val(Z,P,K=KP,dperp=2.0,phase=True,only=None):
    """Z : (n,2) array of z points.  Returns the integrand at each."""
    zx,zy=Z[:,0],Z[:,1]
    xz=np.stack([x[0]-zx,x[1]-zy],1); xw=(x-w)[None,:]+0*xz
    zw=np.stack([zx-w[0],zy-w[1]],1); yz=np.stack([yp[0]-zx,yp[1]-zy],1)
    A=(xz**2).sum(1); B=s(x-w); C=(zw**2).sum(1); Dv=P*A+K*B
    N=yz/((yz**2).sum(1))[:,None]
    S=P+K
    NX=(N*X).sum(1)
    def out(a,b): return (a*N).sum(1)*(b*X).sum(1)     # contracts R=outer(a,b) with N^m X^k'
    g={}
    g['d_perp']  = NX*dperp*P/(2*C*S**2)*(A-B)
    g['S/k+']    = NX*(P/(K*S))*(((xz*zw).sum(1))/C-((xz*xw).sum(1))/(2*B)-(A-B)/(2*C))
    g['S/p+']    = NX*(1/S)*(((xw*zw).sum(1))/C+((xz*xw).sum(1))/(2*A)-(A-B)/(2*C))
    g['-1/p+']   = -(1/P)*( out(zw,xw)/C + out(xz,xw)/(2*A) )
    g['-p+/k+^2']= -(P/K**2)*( out(zw,xz)/C - out(xw,xz)/(2*B) )
    g['-1/k+']   = -(1/K)*( out(xz,zw)/C - out(xz,xw)/(2*B) + out(xw,zw)/C + out(xw,xz)/(2*A) )
    tot=g[only] if only else sum(g.values())
    tot=tot/Dv
    if phase: tot=tot*np.exp(-1j*(P/K)*(yz@KT))
    return tot
