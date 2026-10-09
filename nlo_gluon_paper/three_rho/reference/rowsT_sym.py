"""Region T rows, symbolic.  Names: measured gluon at w (amp) / wb (conj), its source x / xb;
soft gluon at z; soft emitter y (amp) / yb (conj) in real rows; u, v = soft-loop emitters in virtual rows."""
from sym import *
I=1j
def V(p,q): return (p,q)
# ---- two-gluon amplitude pieces at W, measured source X, soft emitter Y. each term: (coef,Us,fs,rs,mvec,svec)
def amp2(W,X,Y):
    mv=(X,W)
    P={}
    P['BB']  =[(1,[('a','c',W)],[],[('h',Y),('c',X)],mv,(Y,'z'))]
    P['sfkB']=[(-I,[('a','c',W)],[('h','c','d')],[('d',X)],mv,(W,'z'))]
    P['AA']  =[(-1,[('m','h','z'),('m','n',Y),('a','c',X)],[],[('n',Y),('c',X)],mv,(Y,'z'))]
    P['sfkA']=[(I,[('m','h','z'),('c','d',X)],[('m','a','c')],[('d',X)],mv,(W,'z'))]
    P['X1']  =[(1,[('m','h','z'),('m','n',Y),('a','c',X)],[],[('n',Y),('c',X)],mv,(Y,'z')),
               (-1,[('m','h','z'),('m','n',Y),('a','c',W)],[],[('n',Y),('c',X)],mv,(Y,'z'))]
    P['X2']  =[(1,[('a','c',X),('m','h','z'),('m','n',Y)],[],[('c',X),('n',Y)],mv,(Y,'z')),
               (-1,[('a','c',X)],[],[('c',X),('h',Y)],mv,(Y,'z'))]
    P['C']   =[(I,[('m','h','z'),('c','d',W)],[('m','a','c')],[('d',X)],mv,(W,'z')),
               (-I,[('m','h','z'),('c','d',X)],[('m','a','c')],[('d',X)],mv,(W,'z'))]
    P['B']=P['BB']+P['sfkB']+P['AA']+P['sfkA']
    return P
def realrow(L1,L2,tag):
    """sum_{a,h} conj(L1 at wb)^dag * L2 at w"""
    out=[]
    for k1,(c1,U1,f1,r1,m1,s1) in enumerate(L1):
        for k2,(c2,U2,f2,r2,m2,s2) in enumerate(L2):
            t1=fresh(Term(c1,U1,f1,r1,()),'_c%d'%k1); t2=fresh(Term(c2,U2,f2,r2,()),'_%d'%k2)
            t=mul(conj(t1),t2,tag); t.kern=((m1,m2),(s1,s2)); out.append(t)
    return out
A=amp2('w','x','y'); B=amp2('wb','xb','yb')
def real_rows():
    R={}
    R['G2-I']=realrow(B['B'],A['B'],'G2-I'); R['G2-II']=realrow(B['X1'],A['X1'],'G2-II')
    R['G2-III']=realrow(B['X2'],A['X2'],'G2-III'); R['G2-IV']=realrow(B['C'],A['C'],'G2-IV')
    for nm,(p,q) in {'G3-I':('X1','B'),'G3-II':('X2','B'),'G3-III':('C','B'),'G3-IV':('X2','X1'),'G3-V':('C','X1'),'G3-VI':('C','X2')}.items():
        R[nm]=realrow(B[p],A[q],nm)+realrow(B[q],A[p],nm)
    return R
# ---- one-gluon (virtual) pieces at W with measured source X; soft loop emitters u,v at vertex z
NC=3
def amp1(W,X,u,v):
    K=lambda p,q: ((p,'z'),(q,'z'))
    mv=(X,W); P={}
    P['I_Nafter']=[(0.5,[('a','c',X)],[],[('c',X),('h',u),('h',v)],mv,K(u,v))]
    P['I_Nbarbefore']=[(-0.5,[('m','n',u),('m','e',v),('a','c',W)],[],[('n',u),('e',v),('c',X)],mv,K(u,v))]
    P['II_A3']=[(-0.5,[('a','c',W)],[],[('h',u),('h',v),('c',X)],mv,K(u,v)),
                (I,[('a','c',W)],[('h','c','d')],[('h',u),('d',X)],mv,K(u,W)),
                (-0.5*NC,[('a','c',W)],[],[('c',X)],mv,K(W,W))]
    P['II_A3bar']=[(0.5,[('a','c',X),('m','n',u),('m','e',v)],[],[('c',X),('n',u),('e',v)],mv,K(u,v))]
    P['III']=[(-1,[('a','c',X),('m','h','z'),('m','n',u)],[],[('c',X),('n',u),('h',v)],mv,K(u,v))]
    P['IV']=[(1,[('m','h','z'),('m','n',u),('a','c',W)],[],[('n',u),('h',v),('c',X)],mv,K(u,v)),
             (-I,[('m','h','z'),('m','n',u),('a','c',W)],[('h','c','d')],[('n',u),('d',X)],mv,K(u,W))]
    # V: K(W,v) (U^T(z)U(W))^{hn} U^{ac}(W) (T^n)_{cd} rho^h(v) rho^d(X) ;  K(W,W) U^{ac}(W) U^{mh}(z)U^{mn}(W) (T^n T^h)_{cd} rho^d(X)
    P['V']=[(-I,[('m','h','z'),('m','n',W),('a','c',W)],[('n','c','d')],[('h',v),('d',X)],mv,K(W,v)),
            (-1,[('a','c',W),('m','h','z'),('m','n',W)],[('n','c','e'),('h','e','d')],[('d',X)],mv,K(W,W))]
    # Ia: -K(u,W) U^{mg}(u) rho^g(u) U^{ac}(W) U^{mn}(W) (T^n)_{cd} rho^d(X) ; -1/2 Nc K(W,W) U^{ac}(W) rho^c(X)
    P['Ia']=[(I,[('m','g',u),('a','c',W),('m','n',W)],[('n','c','d')],[('g',u),('d',X)],mv,K(u,W)),
             (-0.5*NC,[('a','c',W)],[],[('c',X)],mv,K(W,W))]
    return P
def LOpiece(W,X): return [(1,[('a','b',W)],[],[('b',X)],(X,W),()),(-1,[('a','b',X)],[],[('b',X)],(X,W),())]
def virtrow(L,P,tag):
    out=[]
    for k1,(c1,U1,f1,r1,m1,s1) in enumerate(L):
        for k2,(c2,U2,f2,r2,m2,s2) in enumerate(P):
            t1=fresh(Term(c1,U1,f1,r1,()),'_c%d'%k1); t2=fresh(Term(c2,U2,f2,r2,()),'_%d'%k2)
            t=mul(conj(t1),t2,tag); t.kern=((m1,m2),)+tuple(s1)+tuple(s2) if False else ((m1,m2),)+((s1+s2),) if (s1+s2) else ((m1,m2),)
            out.append(t)
    return out
def virt_rows():
    PA=amp1('w','x','u','v'); PB=amp1('wb','xb','u','v'); LA=LOpiece('w','x'); LB=LOpiece('wb','xb')
    R={}
    for n in PA:
        # LO(conj)^dag piece(amp)  +  piece(conj)^dag LO(amp)
        R['G1-'+n]=virtrow(LB,PA[n],'G1-'+n)+virtrow(PB[n],LA,'G1-'+n)
    return R
def fix_kern(t):
    # kernel stored as ((m1,m2), (s1,s2)) where s=(K pair) ; normalize to tuple of dot-pairs
    out=[]
    for kp in t.kern:
        if len(kp)==2 and isinstance(kp[0],tuple) and isinstance(kp[0][0],str): out.append(kp)
        else:
            for x in kp: out.append(x)
    t.kern=tuple(out); return t
if __name__=='__main__':
    import rowsT
    pos,D,U,Xs,reps=rowsT.setup()
    c=rowsT.Ctx(pos,U,Xs,reps); N=Num(c)
    Rn=rowsT.rows(c)
    Rs={**real_rows(),**virt_rows()}
    for k,L in Rs.items():
        op=sum(N.op(t) for t in L)
        print('%-16s numeric %.4g  symbolic %.4g  diff %.2e  (%d terms)'%(k,np.abs(Rn[k]).max(),np.abs(op).max(),np.abs(op-Rn[k]).max(),len(L)))
