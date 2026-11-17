# Net finite remainder of the pair  (Group I Row III last block) + (Group III Row V soft term)^dagger
# Claim:  int_w e^{ik.w} I2[F](w) = int_z e^{ik.z} F(z) [ 2 pi lam(x-z)^2 - D(k,x-z) ],
#         lam(b) = ln(2 e^{-gamma}/(|k||b|)),   D(k,b) = int d^2v e^{-ik.v} ln((b+v)^2/b^2)/v^2
# so the ln^2 of Row V cancels against the Row III block and only D is left.
import mpmath as mp, numpy as np
from scipy import integrate, special
g=float(mp.euler)
# (a) the radial identity used:  2 pi int_delta^inf dv/v J0(kv) 2 ln(v/b) + 2 pi ln^2(b/delta) = 2 pi lam^2
k,b=1.3,0.7
lam=np.log(2*np.exp(-g)/(k*b))
for delta in (1e-3,1e-5):
    f=lambda t: special.j0(k*np.exp(t))*2*(t-np.log(b))
    val=0; edges=np.concatenate([np.linspace(np.log(delta),np.log(1.0/k),40),np.log(np.linspace(1.0/k,4000/k,6000))[1:]])
    for lo,hi in zip(edges[:-1],edges[1:]): val+=integrate.quad(f,lo,hi,limit=200)[0]
    print(f"delta={delta:g}:  lhs = {2*np.pi*val+2*np.pi*np.log(b/delta)**2:.5f}   2 pi lam^2 = {2*np.pi*lam**2:.5f}")
# (b) D(k,b) numerically (polar coordinates around v=0, Gaussian damping removed by extrapolation)
def D(kv,bv,eta):
    kv=np.array(kv); bv=np.array(bv); b2=bv@bv
    def inner(r,part):
        def gph(ph):
            v=np.array([r*np.cos(ph),r*np.sin(ph)]); w=bv+v
            val=np.exp(-1j*kv@v)*np.log((w@w)/b2)/(r*r)*np.exp(-eta*r*r)
            return val.real if part==0 else val.imag
        pts=[np.arctan2(-bv[1],-bv[0])%(2*np.pi)]
        return integrate.quad(gph,0,2*np.pi,limit=400,points=pts)[0]*r
    bn=np.sqrt(b2); out=0
    for part in (0,1):
        s=0
        for lo,hi in ((0,0.5*bn),(0.5*bn,bn),(bn,2*bn),(2*bn,8),(8,25),(25,60)):
            s+=integrate.quad(inner,lo,hi,args=(part,),limit=400)[0]
        out+=s*(1 if part==0 else 1j)
    return out
print("   |k||b|   angle     D (eta->0 extrapolated)          2 pi lam^2")
for kk,bb,th in ((0.2,0.5,0.6),(1.0,0.7,0.6),(1.0,0.7,2.0),(3.0,1.0,0.6)):
    kv=(kk,0.0); bv=(bb*np.cos(th),bb*np.sin(th))
    d1=D(kv,bv,0.004); d2=D(kv,bv,0.002); d=2*d2-d1
    lam=np.log(2*np.exp(-g)/(kk*bb))
    print(f"   {kk*bb:5.2f}   {th:4.1f}    {d.real:+9.4f} {d.imag:+9.4f} i        {2*np.pi*lam**2:9.4f}")
