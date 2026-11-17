# Independent check of the JIMWLK expression (singlet and explicit forms of the message)
# against  -H_JIMWLK acting on  O = Tr{(U^dag(x)-U^dag(w))(U(y)-U(w'))},   H from Eq.(HLO) of Calculation_main:
#   H = (alpha_s/2pi^2) int_{t,u,v} K(u,v,t) [U^dag(u)-U^dag(t)]_{ce} J^c(u) [U(v)-U(t)]_{ed} J^d(v),
#   [J^a(z),U(x)] = delta(x-z) U(x) T^a ,  [J^a(z),U^dag(x)] = -delta(x-z) T^a U^dag(x).
import numpy as np, itertools, sys
from scipy.linalg import expm
def run(N,seed):
    rng=np.random.default_rng(seed)
    # generators of SU(N) fundamental
    gens=[]
    for i in range(N):
        for j in range(i+1,N):
            m=np.zeros((N,N),complex); m[i,j]=m[j,i]=0.5; gens.append(m)
            m=np.zeros((N,N),complex); m[i,j]=-0.5j; m[j,i]=0.5j; gens.append(m)
    for k in range(1,N):
        m=np.zeros((N,N),complex)
        for i in range(k): m[i,i]=1
        m[k,k]=-k; gens.append(m/np.sqrt(2*k*(k+1)))
    t=np.array(gens); n=len(t)
    f=np.zeros((n,n,n))
    for a in range(n):
        for b in range(n):
            c=t[a]@t[b]-t[b]@t[a]
            for d in range(n): f[a,b,d]=np.real(-2j*np.trace(c@t[d]))
    T=-1j*f   # (T^e)_{ab} = -i f^{eab}
    def adj():
        H=sum(rng.normal()*t[a] for a in range(n)); V=expm(1j*H)
        return np.real(np.array([[2*np.trace(t[a]@V@t[b]@V.conj().T) for b in range(n)] for a in range(n)]))
    pts=['x','y','w','wp','t']
    U={p:adj() for p in pts}
    Ud={p:U[p].T for p in pts}
    tr=np.trace
    # O = sum_{p in {x,w}, q in {y,wp}} s_p s_q Tr(Ud[p] U[q])
    terms=[('x','y',1),('x','wp',-1),('w','y',-1),('w','wp',1)]
    def D1(d,v):       # J^d(v) O  (number)
        s=0
        for p,q,sg in terms:
            if q==v: s+=sg*tr(Ud[p]@U[q]@T[d])
            if p==v: s+=-sg*tr(T[d]@Ud[p]@U[q])
        return s
    def D2(c,u,d,v):   # J^c(u) J^d(v) O
        s=0
        for p,q,sg in terms:
            if q==v:
                # Tr(Ud[p] U[q] T^d)
                if q==u: s+=sg*tr(Ud[p]@U[q]@T[c]@T[d])
                if p==u: s+=-sg*tr(T[c]@Ud[p]@U[q]@T[d])
            if p==v:
                # -Tr(T^d Ud[p] U[q])
                if p==u: s+=sg*tr(T[d]@T[c]@Ud[p]@U[q])
                if q==u: s+=-sg*tr(T[d]@Ud[p]@U[q]@T[c])
        return s
    act=['x','y','w','wp']
    HO={}   # coefficient of K(u,v,t) in  sum_{u,v} K(u,v,t) [...]
    for u in act:
        for v in act:
            A=Ud[u]-Ud['t']; B=U[v]-U['t']
            val=0
            for c in range(n):
                for d in range(n):
                    M=(A@B)[c,d]
                    if abs(M)>0: val+=M*D2(c,u,d,v)
            if u==v:
                # J^c(u) acting on U(v) inside [U(v)-U(t)]:  [U^dag(u)-U^dag(t)]_{ce} [U(u) T^c]_{ed} J^d(u) O
                for d in range(n):
                    co=sum((A@U[u]@T[c])[c,d] for c in range(n))
                    val+=co*D1(d,u)
            key=tuple(sorted((u,v)))
            HO[key]=HO.get(key,0)+val
    # ---------- user's singlet form (z -> t) ----------
    z='t'
    def TT(A,B): return sum(tr(T[e]@A@T[e]@B) for e in range(n))
    def T1(A,B): return sum(tr(T[e]@A)*tr(T[e]@B) for e in range(n))
    P=lambda a,b: Ud[a]@U[b]
    S={}
    S[('w','wp')]=4*N*tr(P('w','wp'))-2*TT(P('w','wp'),P('wp',z))-2*TT(P('w','wp'),P(z,'w'))
    S[('wp','x')]=-4*N*tr(P('x','wp'))+2*TT(P('wp','x'),P('x',z))+2*TT(P('wp','x'),P(z,'wp'))
    S[('x','y')]=4*N*tr(P('x','y'))-2*TT(P('x','y'),P('y',z))-2*TT(P('x','y'),P(z,'x'))
    S[('w','y')]=-4*N*tr(P('w','y'))+2*TT(P('w','y'),P('y',z))+2*TT(P('w','y'),P(z,'w'))
    S[('wp','wp')]=(-2*N*tr(P('w','wp'))+2*N*tr(P('x','wp'))-TT(P('wp','x'),P('wp',z))-TT(P('wp','x'),P(z,'wp'))
                    +TT(P('w','wp'),P('wp',z))+TT(P('w','wp'),P(z,'wp'))+T1(P('wp','x'),P('wp',z))+T1(P('w','wp'),P('wp',z)))
    S[('w','w')]=(-2*N*tr(P('w','wp'))+2*N*tr(P('w','y'))+TT(P('w','wp'),P('w',z))+TT(P('w','wp'),P(z,'w'))
                  -TT(P('w','y'),P('w',z))-TT(P('w','y'),P(z,'w'))-T1(P('w','wp'),P('w',z))+T1(P('w','y'),P('w',z)))
    S[('x','x')]=(2*N*tr(P('x','wp'))-2*N*tr(P('x','y'))-TT(P('wp','x'),P('x',z))-TT(P('wp','x'),P(z,'x'))
                  +TT(P('x','y'),P('x',z))+TT(P('x','y'),P(z,'x'))-T1(P('wp','x'),P('x',z))-T1(P('x','y'),P('x',z)))
    S[('y','y')]=(2*N*tr(P('w','y'))-2*N*tr(P('x','y'))-TT(P('w','y'),P('y',z))-TT(P('w','y'),P(z,'y'))
                  +TT(P('x','y'),P('y',z))+TT(P('x','y'),P(z,'y'))-T1(P('w','y'),P('y',z))+T1(P('x','y'),P('y',z)))
    # ---------- user's explicit form: x' -> y, rho^b(x')rho^c(x) -> delta^{bc} ----------
    E=np.einsum; I=np.eye(n); Uz=U[z]; Ux=U['x']; Uy=U['y']; Uw=U['w']; Uwp=U['wp']
    X={}
    # index conventions: 'lc' = rho^l(x') rho^c(x) ; 'cl' = rho^c(x') rho^l(x)
    X[('x','x')]=2*N*tr(Uwp.T@Ux)-2*N*tr(Uy.T@Ux)+2*E('abc,deh,da,el,hb,lc->',f,f,Uz,Uwp,Ux,I)-2*E('abc,deh,da,el,hb,lc->',f,f,Uz,Uy,Ux,I)
    X[('y','y')]=2*N*tr(Uy.T@Uw)-2*N*tr(Uy.T@Ux)+2*E('abc,deh,da,el,hb,cl->',f,f,Uz,Uw,Uy,I)-2*E('abc,deh,da,el,hb,cl->',f,f,Uz,Ux,Uy,I)
    X[('x','y')]=(2*N*tr(Uy.T@Ux)+2*E('abc,ade,hb,hd,ce->',f,f,Uy,Ux,I)-2*E('abc,deh,da,eb,hl,lc->',f,f,Uz,Ux,Uy,I)
                  -2*E('abc,deh,da,eb,hl,cl->',f,f,Uz,Uy,Ux,I))
    X[('w','w')]=-2*N*tr(Uwp.T@Uw)+2*N*tr(Uy.T@Uw)-2*E('abc,deh,da,el,hb,lc->',f,f,Uz,Uwp,Uw,I)+2*E('abc,deh,da,el,hb,lc->',f,f,Uz,Uy,Uw,I)
    X[('wp','wp')]=-2*N*tr(Uwp.T@Uw)+2*N*tr(Uwp.T@Ux)-2*E('abc,deh,da,el,hb,cl->',f,f,Uz,Uw,Uwp,I)+2*E('abc,deh,da,el,hb,cl->',f,f,Uz,Ux,Uwp,I)
    X[('w','wp')]=(2*N*tr(Uwp.T@Uw)+2*E('abc,ade,hb,hd,ce->',f,f,Uwp,Uw,I)-2*E('abc,deh,da,eb,hl,lc->',f,f,Uz,Uw,Uwp,I)
                   -2*E('abc,deh,da,eb,hl,cl->',f,f,Uz,Uwp,Uw,I))
    X[('wp','x')]=(-2*N*tr(Uwp.T@Ux)-2*E('abc,ade,hb,hd,ce->',f,f,Uwp,Ux,I)+2*E('abc,deh,da,eb,hl,lc->',f,f,Uz,Ux,Uwp,I)
                   +2*E('abc,deh,da,eb,hl,cl->',f,f,Uz,Uwp,Ux,I))
    X[('w','y')]=(-2*N*tr(Uy.T@Uw)-2*E('abc,ade,hb,hd,ce->',f,f,Uy,Uw,I)+2*E('abc,deh,da,eb,hl,lc->',f,f,Uz,Uw,Uy,I)
                  +2*E('abc,deh,da,eb,hl,cl->',f,f,Uz,Uy,Uw,I))
    print(f"SU({N}) seed {seed}:  kernel   -H.O (independent)     singlet form       explicit form")
    ok=True
    for key in sorted(HO):
        a=-HO[key]; s=S.get(key); e=X.get(key)
        fl=lambda v: "   (absent)   " if v is None else f"{np.real(v):+14.8f}"
        good=(s is not None and abs(a-s)<1e-8 and abs(a-e)<1e-8) or (s is None and abs(a)<1e-8)
        ok&=good
        print(f"   K{key!s:14s} {np.real(a):+14.8f} (im {abs(np.imag(a)):.0e})  {fl(s)}   {fl(e)}   {'OK' if good else 'MISMATCH'}")
    return ok
r=[run(3,1),run(3,7),run(2,3)]
print("all kernels agree:",all(r))
