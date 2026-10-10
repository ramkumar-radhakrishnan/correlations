"""Region P (p harder than k): rows, symbolic.  Measured k at w/wb with valence source x/xb;
hard gluon p at z with source y (amp) / yb (conj); p-loop emitters u,v in virtual rows."""
from sym import *
from rowsT_sym import realrow, virtrow, LOpiece
I=1j
def amp2P(W,X,Y):
    mv=(X,W); kv=('z',W); pv=(Y,'z'); P={}
    P['BB']  =[(1,[('a','c',W)],[],[('c',X),('h',Y)],mv,pv)]
    P['kfpB']=[(-I,[('a','c',W)],[('c','h','e')],[('e',Y)],kv,pv)]
    P['AA']  =[(-1,[('m','h','z'),('a','c',X),('m','n',Y)],[],[('c',X),('n',Y)],mv,pv)]
    P['kfpA']=[(I,[('m','h','z'),('e','g',Y)],[('a','m','e')],[('g',Y)],kv,pv)]
    P['X1']  =[(1,[('m','h','z'),('m','n',Y),('a','c',X)],[],[('n',Y),('c',X)],mv,pv),
               (-1,[('m','h','z'),('m','n',Y),('a','c',W)],[],[('n',Y),('c',X)],mv,pv)]
    P['X2']  =[(1,[('a','c',X),('m','h','z'),('m','n',Y)],[],[('c',X),('n',Y)],mv,pv),
               (-1,[('a','c',X)],[],[('c',X),('h',Y)],mv,pv)]
    P['C']   =[(-I,[('m','h','z'),('n','g',Y)],[('a','m','n')],[('g',Y)],kv,pv),
               (I,[('m','h','z'),('n','g','z')],[('a','m','n')],[('g',Y)],kv,pv)]
    P['B']=P['BB']+P['kfpB']+P['AA']+P['kfpA']
    return P
def real_rows():
    A=amp2P('w','x','y'); B=amp2P('wb','xb','yb'); R={}
    R['G2-I']=realrow(B['B'],A['B'],'G2-I'); R['G2-II']=realrow(B['X1'],A['X1'],'G2-II')
    R['G2-III']=realrow(B['X2'],A['X2'],'G2-III'); R['G2-IV']=realrow(B['C'],A['C'],'G2-IV')
    for nm,(p,q) in {'G3-I':('X1','B'),'G3-II':('X2','B'),'G3-III':('C','B'),'G3-IV':('X2','X1'),'G3-V':('C','X1'),'G3-VI':('C','X2')}.items():
        R[nm]=realrow(B[p],A[q],nm)+realrow(B[q],A[p],nm)
    return R
def amp1P(W,X,u,v):
    K=lambda p,q: ((p,'z'),(q,'z'))
    mv=(X,W); kv=('z',W); P={}
    P['I_Nafter']=[(0.5,[('a','c',X)],[],[('c',X),('h',u),('h',v)],mv,K(u,v))]
    P['I_Nbarbefore']=[(-0.5,[('m','n',u),('m','e',v),('a','c',W)],[],[('n',u),('e',v),('c',X)],mv,K(u,v))]
    P['II_A3']=[(-0.5,[('a','c',W)],[],[('c',X),('h',u),('h',v)],mv,K(u,v))]
    P['II_A3bar']=[(0.5,[('m','n',u),('m','e',v),('a','c',X)],[],[('n',u),('e',v),('c',X)],mv,K(u,v))]
    P['III']=[(-1,[('m','h','z'),('m','n',u),('a','c',X)],[],[('n',u),('c',X),('h',v)],mv,K(u,v)),
              (I,[('m','hp','z'),('m','g',u),('a','n','z')],[('n','hp','h')],[('g',u),('h',v)],kv,K(u,v))]
    P['IV']=[(1,[('m','h','z'),('m','n',u),('a','c',W)],[],[('n',u),('c',X),('h',v)],mv,K(u,v)),
             (-I,[('m','h','z'),('m','n',u),('a','c',W)],[('c','h','e')],[('n',u),('e',v)],kv,K(u,v))]
    return P
def virt_rows():
    PA=amp1P('w','x','u','v'); PB=amp1P('wb','xb','u','v'); LA=LOpiece('w','x'); LB=LOpiece('wb','xb'); R={}
    for n in PA: R['G1-'+n]=virtrow(LB,PA[n],'G1-'+n)+virtrow(PB[n],LA,'G1-'+n)
    return R
if __name__=='__main__':
    import rowsT, pickle
    pos,D,U,Xs,reps=rowsT.setup(reps=('q','qb'))
    c=rowsT.Ctx(pos,U,Xs,reps); N=Num(c)
    r=pickle.load(open('semiP_2.pkl','rb'))
    LO=sum(N.op(t) for t in virtrow(LOpiece('wb','xb'),[(1,[('a','b','w')],[],[('b','x')],('x','w'),()),(-1,[('a','b','x')],[],[('b','x')],('x','w'),())],'LO')[:0]) if False else None
    V=virt_rows(); Rr=real_rows()
    opV=sum(N.op(t) for L in V.values() for t in L); opR=sum(N.op(t) for L in Rr.values() for t in L)
    print('virtual: model %.4g rows %.4g diff %.2e'%(np.abs(r['virtual']).max(),np.abs(opV).max(),np.abs(opV-r['virtual']).max()))
    print('real   : model %.4g rows %.4g diff %.2e'%(np.abs(r['real']).max(),np.abs(opR).max(),np.abs(opR-r['real']).max()))
