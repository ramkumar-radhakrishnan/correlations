# Check M: closed forms of the helicity sums h_AB(a;b) = sum_{l1 l2 l r} Phi_A(a) Phi_B^*(b), A,B in {1,2},
# for the antiquark, quark and interference channels, against brute-force sums of the momentum-space vertex functions
#   Phi_1 = [(1-zeta) N/(k^2+eps^2) - xi eta I]/D,  Phi_2 = N/[zeta(1-zeta)(k^2+eps^2)(p^2+om^2)],  D = zeta(1-zeta)^2(p^2+om^2)+xi eta(k^2+eps^2)
#   (quark: zeta<->eta, N->Nt, I->It, eps->epst, om->omt)
import numpy as np, itertools
rng=np.random.default_rng(5)
H=(0.5,-0.5); E2=np.array([[0,1],[-1,0]])
def d(a,b): return 1.0 if a==b else 0.0
def varphi(i,j,l1,l2,zq): return d(l1,l2)*((2*zq-1)*d(i,j)+2j*E2[i,j]*l1)
def tau(mm,n,l,l2,xi,pd): return d(l,l2)*((2*pd+xi)*d(mm,n)+2j*xi*E2[mm,n]*l)
def mu(p,a,b): return d(b,-a)*(2*a*d(p,0)-1j*d(p,1))
def run(z,x,Q,m):
    e=1-z-x; zb=1-z; om=x*m/zb; omt=x*m/(1-e); eps2=z*zb*Q*Q+m*m; epst2=e*(1-e)*Q*Q+m*m
    hl=list(itertools.product(H,H,(0,1),(0,1)))
    def N(k,p):  return np.array([sum((sum(varphi(l,j,l1,lam,z)*k[j] for j in (0,1))+m*mu(l,l1,lam))*(sum(tau(r,n,lam,l2,x,e)*p[n] for n in (0,1))+x*om*mu(r,lam,l2)) for lam in H) for l1,l2,l,r in hl])
    def Nt(K,P): return np.array([sum((sum(tau(n,r,l1,lam,x,z)*P[n] for n in (0,1))+x*omt*mu(r,l1,lam))*(sum(varphi(l,j,lam,l2,z+x)*K[j] for j in (0,1))+m*mu(l,lam,l2)) for lam in H) for l1,l2,l,r in hl])
    Iq =np.array([(d(l,r)+2j*E2[l,r]*l1)*d(l1,l2) for l1,l2,l,r in hl])
    Iqt=np.array([(d(l,r)-2j*E2[l,r]*l1)*d(l1,l2) for l1,l2,l,r in hl])
    D =lambda k,p: z*zb**2*(p@p+om**2)+x*e*(k@k+eps2)
    Dt=lambda K,P: e*(1-e)**2*(P@P+omt**2)+x*z*(K@K+epst2)
    Ph1=lambda k,p: (zb*N(k,p)/(k@k+eps2)-x*e*Iq)/D(k,p)
    Ph2=lambda k,p: N(k,p)/(z*zb*(k@k+eps2)*(p@p+om**2))
    Pt1=lambda K,P: ((1-e)*Nt(K,P)/(K@K+epst2)-x*z*Iqt)/Dt(K,P)
    Pt2=lambda K,P: Nt(K,P)/(e*(1-e)*(K@K+epst2)*(P@P+omt**2))
    # closed-form building blocks (Check L)
    dot=lambda u,v: u@v; crs=lambda u,v: u[0]*v[1]-u[1]*v[0]
    def Kq(ka,pa,kb,pb):  return (8*(z*z+zb*zb)*(e*e+zb*zb)*dot(ka,kb)*dot(pa,pb)-8*x*(2*z-1)*(2*zb-x)*crs(ka,kb)*crs(pa,pb)
                                  +8*m*m*(e*e+zb*zb)*dot(pa,pb)+8*x*x*om*om*(z*z+zb*zb)*dot(ka,kb)+8*m*m*x*x*om*om
                                  -8*z*e*m*x*om*dot(ka-kb,pa-pb))
    Jq=lambda k,p: 8*z*e*dot(k,p)-8*m*x*om
    def Kt(Ka,Pa,Kb,Pb):  return (8*(z*z+(1-e)**2)*(e*e+(1-e)**2)*dot(Ka,Kb)*dot(Pa,Pb)-8*x*(2*e-1)*(2*(1-e)-x)*crs(Ka,Kb)*crs(Pa,Pb)
                                  +8*m*m*(z*z+(1-e)**2)*dot(Pa,Pb)+8*x*x*omt*omt*(e*e+(1-e)**2)*dot(Ka,Kb)+8*m*m*x*x*omt*omt
                                  +8*z*e*m*x*omt*dot(Ka-Kb,Pa-Pb))
    Jt=lambda K,P: -8*z*e*dot(K,P)-8*m*x*omt
    def Kx(Ka,Pa,kb,pb):  return (8*(1-x-2*z*e)*(x*(1-x)+2*z*e)*dot(Ka,kb)*dot(Pa,pb)-8*x*(1-2*z-x)**2*crs(Ka,kb)*crs(Pa,pb)
                                  +8*m*m*(x+2*z*e)*dot(Pa,pb)-8*(x+2*z*e)*x*x*omt*om*dot(Ka,kb)
                                  -8*(z+x)**2*m*x*om*dot(Ka,Pa)+8*zb**2*m*x*omt*dot(kb,pb)
                                  -8*e*e*m*x*omt*dot(Ka,pb)+8*z*z*m*x*om*dot(Pa,kb))
    Jx =lambda K,P: 8*(z+x)**2*dot(K,P)
    Jxb=lambda k,p: -8*zb**2*dot(k,p)
    A=lambda k: k@k+eps2; B=lambda p: p@p+om**2; At=lambda K: K@K+epst2; Bt=lambda P: P@P+omt**2
    cf={}
    cf['qbar h11']=lambda a,b: (zb**2*Kq(*a,*b)/(A(a[0])*A(b[0]))-zb*x*e*Jq(*a)/A(a[0])-zb*x*e*Jq(*b)/A(b[0])+8*x*x*e*e)/(D(*a)*D(*b))
    cf['qbar h12']=lambda a,b: (zb*Kq(*a,*b)/A(a[0])-x*e*Jq(*b))/(z*zb*D(*a)*A(b[0])*B(b[1]))
    cf['qbar h21']=lambda a,b: (zb*Kq(*a,*b)/A(b[0])-x*e*Jq(*a))/(z*zb*D(*b)*A(a[0])*B(a[1]))
    cf['qbar h22']=lambda a,b: Kq(*a,*b)/(z*z*zb*zb*A(a[0])*B(a[1])*A(b[0])*B(b[1]))
    cf['q h11']=lambda a,b: ((1-e)**2*Kt(*a,*b)/(At(a[0])*At(b[0]))-(1-e)*x*z*Jt(*a)/At(a[0])-(1-e)*x*z*Jt(*b)/At(b[0])+8*x*x*z*z)/(Dt(*a)*Dt(*b))
    cf['q h12']=lambda a,b: ((1-e)*Kt(*a,*b)/At(a[0])-x*z*Jt(*b))/(e*(1-e)*Dt(*a)*At(b[0])*Bt(b[1]))
    cf['q h21']=lambda a,b: ((1-e)*Kt(*a,*b)/At(b[0])-x*z*Jt(*a))/(e*(1-e)*Dt(*b)*At(a[0])*Bt(a[1]))
    cf['q h22']=lambda a,b: Kt(*a,*b)/(e*e*(1-e)**2*At(a[0])*Bt(a[1])*At(b[0])*Bt(b[1]))
    cf['x h11']=lambda a,b: ((1-e)*zb*Kx(*a,*b)/(At(a[0])*A(b[0]))-(1-e)*x*e*Jx(*a)/At(a[0])-zb*x*z*Jxb(*b)/A(b[0]))/(Dt(*a)*D(*b))
    cf['x h12']=lambda a,b: ((1-e)*Kx(*a,*b)/At(a[0])-x*z*Jxb(*b))/(z*zb*Dt(*a)*A(b[0])*B(b[1]))
    cf['x h21']=lambda a,b: (zb*Kx(*a,*b)/A(b[0])-x*e*Jx(*a))/(e*(1-e)*D(*b)*At(a[0])*Bt(a[1]))
    cf['x h22']=lambda a,b: Kx(*a,*b)/(e*(1-e)*z*zb*At(a[0])*Bt(a[1])*A(b[0])*B(b[1]))
    bf={'qbar h11':(Ph1,Ph1),'qbar h12':(Ph1,Ph2),'qbar h21':(Ph2,Ph1),'qbar h22':(Ph2,Ph2),
        'q h11':(Pt1,Pt1),'q h12':(Pt1,Pt2),'q h21':(Pt2,Pt1),'q h22':(Pt2,Pt2),
        'x h11':(Pt1,Ph1),'x h12':(Pt1,Ph2),'x h21':(Pt2,Ph1),'x h22':(Pt2,Ph2)}
    worst={}
    for t in range(4):
        a=(rng.normal(size=2),rng.normal(size=2)); b=(rng.normal(size=2),rng.normal(size=2))
        for k,(F,G) in bf.items():
            brute=np.sum(F(*a)*np.conj(G(*b))); clos=cf[k](a,b)
            worst[k]=max(worst.get(k,0),abs(brute-clos)/abs(brute))
    return worst
for pars in [(0.3,0.2,1.1,0.8),(0.45,0.15,0.7,1.6),(0.2,0.5,1.4,0.0)]:
    w=run(*pars); print("zeta,xi,Q,m =",pars,"  max rel. diff over 4 random (a,b):", {k:f"{v:.0e}" for k,v in w.items()})
