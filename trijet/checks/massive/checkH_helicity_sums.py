# Check H: helicity/polarisation sums  sum_{lam1 lam2 l r} M(x) M*(xbar)  for the massive amplitudes, fully symbolic in zeta, xi.
#   antiquark channel (before SW):
#     M^{lr}_{l1 l2} = sum_lam [ phi^{lj}_{l1 lam}(z) tau^{rn}_{lam l2}(xi,eta) R^j Y^n G11 + i xi om phi^{lj} mu^r R^j G10
#                               - i m mu^l tau^{rn} Y^n G01 + m xi om mu^l mu^r G00 ] - (xi eta/zb)(d^{lr}+2i eps^{lr} l1) d_{l1 l2} GE
#   quark channel (before SW):
#     Mq^{lr}_{l1 l2} = sum_lam [ tau^{nr}_{l1 lam}(xi,z) phi^{lj}_{lam l2}(z+xi) Rt^j X^n Gq11 + i xi omq mu^r_{l1 lam} phi^{lj} Rt^j Gq10
#                               - i m tau^{nr} mu^l_{lam l2} X^n Gq01 + m xi omq mu^r mu^l Gq00 ] - (xi z/(z+xi))(d^{lr}-2i eps^{lr} l1) d_{l1 l2} GE
#   results are decomposed on 2D rotation-invariant tensors; index convention of the note: i<->R, m<->Y, j<->Rbar, n<->Ybar.
import sympy as sp, itertools
I=sp.I; Rat=sp.Rational
z,x,m=sp.symbols('zeta xi m',positive=True); e=1-z-x; zb=1-z
om=x*m/zb; omq=x*m/(1-e)
def eps2(a,b): return [[0,1],[-1,0]][a][b]
def d(a,b): return 1 if a==b else 0
def varphi(i,j,l1,l2,zq): return d(l1,l2)*((2*zq-1)*d(i,j) + 2*I*eps2(i,j)*l1)
def tau(mm,n,l,l2,xi_,par_d): return d(l,l2)*((2*par_d+xi_)*d(mm,n) + 2*I*xi_*eps2(mm,n)*l)
def mu(p,a,b): return d(b,-a)*(2*a*d(p,0) - I*d(p,1))
H=(Rat(1,2),Rat(-1,2))
def vec(n): return sp.symbols(f'{n}1 {n}2',real=True)
Rv,Yv,Rb,Yb,Rt,Xv,Rtb,Xb = [vec(n) for n in ('R','Y','Rb','Yb','Rt','X','Rtb','Xb')]
Gs=sp.symbols('G11 G10 G01 G00 GE',real=True); Gbs=sp.symbols('Gb11 Gb10 Gb01 Gb00 GbE',real=True)
Gq=sp.symbols('Gq11 Gq10 Gq01 Gq00 GqE',real=True); Gqb=sp.symbols('Gqb11 Gqb10 Gqb01 Gqb00 GqbE',real=True)
def Mqb(l1,l2,l,r,Rv,Yv,G,inst=True):
    s=0
    for lam in H:
        s+= sum(varphi(l,j,l1,lam,z)*tau(r,n,lam,l2,x,e)*Rv[j]*Yv[n] for j in (0,1) for n in (0,1))*G[0]
        s+= I*x*om*sum(varphi(l,j,l1,lam,z)*Rv[j] for j in (0,1))*mu(r,lam,l2)*G[1]
        s+= -I*m*mu(l,l1,lam)*sum(tau(r,n,lam,l2,x,e)*Yv[n] for n in (0,1))*G[2]
        s+= m*x*om*mu(l,l1,lam)*mu(r,lam,l2)*G[3]
    if inst: s+= -(x*e/zb)*(d(l,r)+2*I*eps2(l,r)*l1)*d(l1,l2)*G[4]
    return s
def Mq(l1,l2,l,r,Rv,Xv,G,inst=True):
    s=0
    for lam in H:
        s+= sum(tau(n,r,l1,lam,x,z)*varphi(l,j,lam,l2,z+x)*Rv[j]*Xv[n] for j in (0,1) for n in (0,1))*G[0]
        s+= I*x*omq*mu(r,l1,lam)*sum(varphi(l,j,lam,l2,z+x)*Rv[j] for j in (0,1))*G[1]
        s+= -I*m*sum(tau(n,r,l1,lam,x,z)*Xv[n] for n in (0,1))*mu(l,lam,l2)*G[2]
        s+= m*x*omq*mu(r,l1,lam)*mu(l,lam,l2)*G[3]
    if inst: s+= -(x*z/(z+x))*(d(l,r)-2*I*eps2(l,r)*l1)*d(l1,l2)*G[4]
    return s
def hsum(A,B):
    tot=0
    for l1,l2,l,r in itertools.product(H,H,(0,1),(0,1)):
        tot+= sp.expand(A(l1,l2,l,r)*sp.conjugate(B(l1,l2,l,r)))
    return sp.expand(tot)
# ---------- tensor decomposition
def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def crs(a,b): return a[0]*b[1]-a[1]*b[0]
def basis(vs):
    n=len(vs)
    if n==4:
        A,B,C,D=vs   # indices i,m,j,n
        return {'d^{ij}d^{mn}':dot(A,C)*dot(B,D),'e^{ij}e^{mn}':crs(A,C)*crs(B,D),'d^{in}d^{jm}':dot(A,D)*dot(B,C),
                'i e^{ij}d^{mn}':I*crs(A,C)*dot(B,D),'i d^{ij}e^{mn}':I*dot(A,C)*crs(B,D),'i e^{in}d^{jm}':I*crs(A,D)*dot(B,C)}
    if n==2:
        A,B=vs; return {'d':dot(A,B),'i e':I*crs(A,B)}
    return {'1':sp.Integer(1)}
def decompose(poly,vs):
    b=basis(vs); cs=sp.symbols(f'c0:{len(b)}')
    ansatz=sum(c*v for c,v in zip(cs,b.values()))
    comps=[c for v in vs for c in v]
    eqs=sp.Poly(sp.expand(poly-ansatz),*comps).coeffs() if comps else [sp.expand(poly-ansatz)]
    sol=sp.solve(eqs,cs,dict=True)
    assert sol, "decomposition failed"
    sol=sol[0]; return {k:sp.factor(sp.simplify(sol.get(c,c))) for k,c in zip(b.keys(),cs) if sp.simplify(sol.get(c,c))!=0}
def report(title,S,Gl,Gr,vecmap):
    print(f"\n==== {title} ====")
    out={}
    for gi in Gl:
        for gj in Gr:
            c=sp.expand(S.coeff(gi).coeff(gj)) if True else 0
            if c==0: continue
            # which vectors appear
            vs=[v for v in vecmap[(str(gi)[-2:] if not str(gi).endswith('E') else 'E', str(gj)[-2:] if not str(gj).endswith('E') else 'E')]]
            dec=decompose(c,vs)
            print(f"  [{gi} x {gj}]  vectors {[v[0].name[:-1] for v in vs]}:")
            for k,v in dec.items(): print(f"        {k:16s}: {v}")
            out[f"{gi}*{gj}"]={k:str(v) for k,v in dec.items()}
    return out
def vmap(A,B,Ab,Bb):
    left={'11':[A,B],'10':[A],'01':[B],'00':[],'E':[]}
    right={'11':[Ab,Bb],'10':[Ab],'01':[Bb],'00':[],'E':[]}
    return {(a,b):left[a]+right[b] for a in left for b in right}
res={}
S=hsum(lambda *h: Mqb(*h,Rv,Yv,Gs), lambda *h: Mqb(*h,Rb,Yb,Gbs))
res['direct_qbar']=report("DIRECT, antiquark emission: sum M(x) M*(xbar)",S,Gs,Gbs,vmap(Rv,Yv,Rb,Yb))
S=hsum(lambda *h: Mq(*h,Rt,Xv,Gq), lambda *h: Mq(*h,Rtb,Xb,Gqb))
res['direct_q']=report("DIRECT, quark emission: sum Mq(x) Mq*(xbar)",S,Gq,Gqb,vmap(Rt,Xv,Rtb,Xb))
S=hsum(lambda *h: Mq(*h,Rt,Xv,Gq), lambda *h: Mqb(*h,Rb,Yb,Gbs))
res['interference']=report("INTERFERENCE: sum Mq(x) M*(xbar)",S,Gq,Gbs,vmap(Rt,Xv,Rb,Yb))

