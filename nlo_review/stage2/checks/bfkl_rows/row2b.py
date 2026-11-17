import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from tworho import *
np.seterr(all="ignore")
B=blocks("9ea88190-BFKL.tex")
# equivalent bracket after using  int_z [K(x,y,z)-K(w,y,z)] = pi log[(w-y)^2/(x-y)^2]  on the log term of lines 109-110
eq=r"""-\frac{1}{k^{+}}\log{\left(\frac{\vee}{k^{+}}\right)}\frac{1}{(2\pi)^{3}}\,\frac{i\,g^{4}\,f^{hbc}}{32\pi^{5}}\,\frac{(x'-w')^{i}}{(x'-w')^{2}}\,\Bigg[\frac{(x-w)^{i}}{(x-w)^{2}} - \frac{(y-w)^{i}}{(y-w)^{2}}\Bigg]\frac{(x-z)\cdot (y-z)}{(x-z)^{2}\,(y-z)^{2}}
\times\Big[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\,\rho^{b'}(x')\Big]\times\, \Big[U^{ah}(w)\,\{\rho^{b}(x),\rho^{c}(y)\}\Big]"""
if __name__=="__main__":
  for seed in (1,2):
    env=setup(3,seed); rng=np.random.default_rng(100+seed); u=rng.normal(size=2); v=rng.normal(size=2)
    L=Finv(env,[rhs(B[113])],u,v)
    for cc in (True,False):
        a=Finv(env,[eq],u,v,cc=cc); a2=Finv(env,[eq],u,v,cc=cc,eps=1e-6)
        print(seed,"equivalent z-form, cc",cc,"|2rho|=%.4f |list|=%.4f |2rho-list|=%.3e  eps-dep %.1e"%(np.max(np.abs(a)),np.max(np.abs(L)),np.max(np.abs(a-L)),np.max(np.abs(a-a2))))
# numerical check of the integral identity
from scipy import integrate
rng=np.random.default_rng(5); x,y,w=rng.normal(size=(3,2))
def K(a,b,z): return ((a-z)@(b-z))/(((a-z)@(a-z))*((b-z)@(b-z)))
def integrand(r,th):
    z=y+r*np.array([np.cos(th),np.sin(th)]); return r*(K(x,y,z)-K(w,y,z))
angs=sorted([np.arctan2(*(x-y)[::-1])%(2*np.pi),np.arctan2(*(w-y)[::-1])%(2*np.pi)])
rs=sorted([np.linalg.norm(x-y),np.linalg.norm(w-y)])
tot=0
tb=[0]+angs+[2*np.pi]; rb=[0]+rs+[50.]
for i in range(3):
    for j in range(3):
        val,err=integrate.dblquad(integrand,tb[i],tb[i+1],rb[j],rb[j+1],epsabs=1e-7,epsrel=1e-7); tot+=val
# tail r>50: integrand ~ r*( (x-w).zhat /r^2 )/r -> angular average 0 at O(1/r^2); remaining O(1/r^3)*r -> negligible ~1e-3
print("int_z [K(x,y,z)-K(w,y,z)] = %.5f   pi*log((w-y)^2/(x-y)^2) = %.5f"%(tot,np.pi*np.log(((w-y)@(w-y))/((x-y)@(x-y)))))
