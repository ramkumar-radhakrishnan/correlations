# Check L: helicity/polarisation sums of the momentum-space vertex numerators, for two DIFFERENT momentum arguments
#   N^{lr}_{l1 l2}(k,p)  = sum_lam [phi^{lj}_{l1 lam}(zeta) k^j + m mu^l_{l1 lam}] [tau^{rn}_{lam l2}(xi,eta) p^n + xi om mu^r_{lam l2}]   (antiquark emission)
#   Nt^{lr}_{l1 l2}(K,P) = sum_lam [tau^{nr}_{l1 lam}(xi,zeta) P^n + xi omt mu^r_{l1 lam}] [phi^{lj}_{lam l2}(zeta+xi) K^j + m mu^l_{lam l2}] (quark emission)
#   I^{lr} = (d^{lr} + 2i eps^{lr} l1) d_{l1 l2},   It^{lr} = (d^{lr} - 2i eps^{lr} l1) d_{l1 l2}
import sympy as sp, itertools
I=sp.I; Rat=sp.Rational
z,x,m=sp.symbols('zeta xi m',positive=True); e=1-z-x; zb=1-z
om=x*m/zb; omt=x*m/(1-e)
def eps2(a,b): return [[0,1],[-1,0]][a][b]
def d(a,b): return 1 if a==b else 0
def varphi(i,j,l1,l2,zq): return d(l1,l2)*((2*zq-1)*d(i,j) + 2*I*eps2(i,j)*l1)
def tau(mm,n,l,l2,xi_,par_d): return d(l,l2)*((2*par_d+xi_)*d(mm,n) + 2*I*xi_*eps2(mm,n)*l)
def mu(p,a,b): return d(b,-a)*(2*a*d(p,0) - I*d(p,1))
H=(Rat(1,2),Rat(-1,2))
def V(n): return sp.symbols(f'{n}1 {n}2',real=True)
ka,pa,kb,pb=V('ka'),V('pa'),V('kb'),V('pb')
def N(l1,l2,l,r,k,p):
    return sum((sum(varphi(l,j,l1,lam,z)*k[j] for j in (0,1))+m*mu(l,l1,lam))*(sum(tau(r,n,lam,l2,x,e)*p[n] for n in (0,1))+x*om*mu(r,lam,l2)) for lam in H)
def Nt(l1,l2,l,r,K,P):
    return sum((sum(tau(n,r,l1,lam,x,z)*P[n] for n in (0,1))+x*omt*mu(r,l1,lam))*(sum(varphi(l,j,lam,l2,z+x)*K[j] for j in (0,1))+m*mu(l,lam,l2)) for lam in H)
def Iq(l1,l2,l,r): return (d(l,r)+2*I*eps2(l,r)*l1)*d(l1,l2)
def Iqt(l1,l2,l,r): return (d(l,r)-2*I*eps2(l,r)*l1)*d(l1,l2)
def hs(A,B): return sp.expand(sum(A(*h)*sp.conjugate(B(*h)) for h in itertools.product(H,H,(0,1),(0,1))))
dot=lambda u,v: u[0]*v[0]+u[1]*v[1]; crs=lambda u,v: u[0]*v[1]-u[1]*v[0]
def decompose(poly):
    cands={'(ka.kb)(pa.pb)':dot(ka,kb)*dot(pa,pb),'(ka x kb)(pa x pb)':crs(ka,kb)*crs(pa,pb),'(ka.pb)(pa.kb)':dot(ka,pb)*dot(pa,kb),
           'i(ka x kb)(pa.pb)':I*crs(ka,kb)*dot(pa,pb),'i(ka.kb)(pa x pb)':I*dot(ka,kb)*crs(pa,pb),'i(ka x pb)(pa.kb)':I*crs(ka,pb)*dot(pa,kb)}
    vs={'ka':ka,'pa':pa,'kb':kb,'pb':pb}
    for (n1,v1),(n2,v2) in itertools.combinations(vs.items(),2):
        cands[f'({n1}.{n2})']=dot(v1,v2); cands[f'i({n1} x {n2})']=I*crs(v1,v2)
    for n1,v1 in vs.items(): cands[f'{n1}^2']=dot(v1,v1)
    cands['1']=sp.Integer(1)
    cs=sp.symbols(f'c0:{len(cands)}')
    ans=sum(c_*v_ for c_,v_ in zip(cs,cands.values()))
    comps=[*ka,*pa,*kb,*pb]
    sol=sp.solve(sp.Poly(sp.expand(poly-ans),*comps).coeffs(),cs,dict=True)[0]
    return {k_:sp.factor(sp.simplify(sol.get(c_,0))) for k_,c_ in zip(cands,cs) if sp.simplify(sol.get(c_,0))!=0}
def show(title,poly):
    print(f"== {title}")
    for k_,v_ in decompose(poly).items(): print(f"     {k_:22s}: {v_}")
show("K_qbar(a;b) = sum N(ka,pa) N*(kb,pb)", hs(lambda *h: N(*h,ka,pa), lambda *h: N(*h,kb,pb)))
show("J_qbar(a)   = sum N(ka,pa) I*",          hs(lambda *h: N(*h,ka,pa), Iq))
show("sum I I*",                               hs(Iq,Iq))
show("K_q(a;b)    = sum Nt(ka,pa) Nt*(kb,pb)", hs(lambda *h: Nt(*h,ka,pa), lambda *h: Nt(*h,kb,pb)))
show("J_q(a)      = sum Nt(ka,pa) It*",         hs(lambda *h: Nt(*h,ka,pa), Iqt))
show("sum It It*",                             hs(Iqt,Iqt))
show("K_x(a;b)    = sum Nt(ka,pa) N*(kb,pb)  [a: quark-emission momenta, b: antiquark-emission momenta]", hs(lambda *h: Nt(*h,ka,pa), lambda *h: N(*h,kb,pb)))
show("J_x(a)      = sum Nt(ka,pa) I*",          hs(lambda *h: Nt(*h,ka,pa), Iq))
show("Jbar_x(b)   = sum It N*(kb,pb)",          hs(Iqt, lambda *h: N(*h,kb,pb)))
show("sum It I*",                              hs(Iqt,Iq))
