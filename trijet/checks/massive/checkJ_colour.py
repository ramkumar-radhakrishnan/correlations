# Check J: large-Nc reduction of every Wilson-line structure W used in the cross sections.
# Colour matrices of the amplitudes (adjoint index a):
#   B(x,y,z)^a  = U^{ab}(z) V(x) t^b V^dag(y) - t^a        (gluon emitted before the shock wave, either emitter)
#   Aqb(x,y')^a = V(x) V^dag(y') t^a - t^a                 (antiquark emits after the shock wave)
#   Aq(x',y)^a  = t^a V(x') V^dag(y) - t^a                 (quark emits after the shock wave)
#   W(A;B) = 1/(C_F N) sum_a Tr[A^a B^a dagger]
# U^{ab}(z) = 2 Tr[t^a V(z) t^b V^dag(z)].  S(a,b)=Tr[V_a V_b^dag]/N,  Q(a,b,c,d)=Tr[V_a V_b^dag V_c V_d^dag]/N.
# For each W we print |W_exact - W_largeNc| for growing N; it must fall like 1/N^2.
import numpy as np, scipy.linalg as sl
rng=np.random.default_rng(7)
def gens(N):
    T=[]
    for i in range(N):
        for j in range(i+1,N):
            M=np.zeros((N,N),complex); M[i,j]=M[j,i]=0.5; T.append(M)
            M=np.zeros((N,N),complex); M[i,j]=-0.5j; M[j,i]=0.5j; T.append(M)
    for k in range(1,N):
        M=np.zeros((N,N),complex)
        for i in range(k): M[i,i]=1
        M[k,k]=-k; M*=1/np.sqrt(2*k*(k+1)); T.append(M)
    return np.array(T)
def randV(N,eps=1.2):
    A=(rng.normal(size=(N,N))+1j*rng.normal(size=(N,N)))/np.sqrt(2*N); H=(A+A.conj().T)/2
    V=sl.expm(1j*eps*H); return V/np.linalg.det(V)**(1/N)
def run(N):
    T=gens(N); CF=(N*N-1)/(2*N); I=np.eye(N)
    pts=['x','y','z','yp','xb','yb','zb','ybp','xp','xbp']
    V={p:randV(N) for p in pts}; Vd={p:V[p].conj().T for p in pts}
    n=len(T); Tf=T.reshape(n,N*N)
    def U(p):
        M=np.matmul(np.matmul(V[p],T),Vd[p])            # M_b = V t^b V^dag
        return 2*(Tf@M.transpose(0,2,1).reshape(n,N*N).T).real   # U^{ab} = 2 Tr[t^a M_b]
    S=lambda a,b: np.trace(V[a]@Vd[b])/N
    Q=lambda a,b,c,d: np.trace(V[a]@Vd[b]@V[c]@Vd[d])/N
    def B(x,y,z):
        X=np.matmul(np.matmul(V[x],T),Vd[y]); return np.tensordot(U(z),X,axes=(1,0))-T
    def Aqb(x,yp): return np.matmul(V[x]@Vd[yp],T)-T
    def Aq(xp,y): return np.matmul(T,V[xp]@Vd[y])-T
    W=lambda A,Bm: np.einsum('aij,aij->',A,Bm.conj())/(CF*N)
    out={}
    # antiquark direct
    out['qbar W(x,y,z;xb,yb,zb)']   =(W(B('x','y','z'),B('xb','yb','zb')), Q('x','z','zb','xb')*Q('z','y','yb','zb')-S('x','z')*S('z','y')-S('zb','xb')*S('yb','zb')+1)
    out["qbar W(x,y',y';xb,yb,zb)"] =(W(Aqb('x','yp'),B('xb','yb','zb')), Q('x','yp','zb','xb')*S('yb','zb')-S('x','yp')-S('zb','xb')*S('yb','zb')+1)
    out["qbar W(x,y,z;xb,yb',yb')"] =(W(B('x','y','z'),Aqb('xb','ybp')), Q('x','z','ybp','xb')*S('z','y')-S('x','z')*S('z','y')-S('ybp','xb')+1)
    out["qbar W(x,y',y';xb,yb',yb')"]=(W(Aqb('x','yp'),Aqb('xb','ybp')), Q('x','yp','ybp','xb')-S('x','yp')-S('ybp','xb')+1)
    # quark direct: the note writes W(y,x,z;...) = complex conjugate of the physical structure -> compare conj
    out['q W(y,x,z;yb,xb,zb)* ']    =(W(B('x','y','z'),B('xb','yb','zb')), np.conj(Q('y','z','zb','yb')*Q('z','x','xb','zb')-S('y','z')*S('z','x')-S('zb','yb')*S('xb','zb')+1))
    out["q W(y,x',x';yb,xb,zb)* "]  =(W(Aq('xp','y'),B('xb','yb','zb')), np.conj(Q('y','xp','zb','yb')*S('xb','zb')-S('y','xp')-S('zb','yb')*S('xb','zb')+1))
    out["q W(y,x,z;yb,xb',xb')* "]  =(W(B('x','y','z'),Aq('xbp','yb')), np.conj(Q('y','z','xbp','yb')*S('z','x')-S('y','z')*S('z','x')-S('xbp','yb')+1))
    out["q W(y,x',x';yb,xb',xb')* "]=(W(Aq('xp','y'),Aq('xbp','yb')), np.conj(Q('y','xp','xbp','yb')-S('y','xp')-S('xbp','yb')+1))
    # interference (quark emission in the amplitude, antiquark emission in the conjugate)
    out['int W(x,y,z;xb,yb,zb)']    =(W(B('x','y','z'),B('xb','yb','zb')), Q('x','z','zb','xb')*Q('z','y','yb','zb')-S('x','z')*S('z','y')-S('zb','xb')*S('yb','zb')+1)
    out["int W(x',y,x';xb,yb,zb)"]  =(W(Aq('xp','y'),B('xb','yb','zb')), S('zb','xb')*Q('xp','y','yb','zb')-S('xp','y')-S('zb','xb')*S('yb','zb')+1)
    out["int W(x,y,z;xb,yb',yb')"]  =(W(B('x','y','z'),Aqb('xb','ybp')), Q('x','z','ybp','xb')*S('z','y')-S('x','z')*S('z','y')-S('ybp','xb')+1)
    out["int W(x',y,x';xb,yb',yb')"]=(W(Aq('xp','y'),Aqb('xb','ybp')), S('ybp','xb')*S('xp','y')-S('xp','y')-S('ybp','xb')+1)
    return out
Ns=(4,8,16,32); nsamp=6
err={N:{} for N in Ns}
for N in Ns:
    for s_ in range(nsamp):
        for k,(ex,ln) in run(N).items(): err[N].setdefault(k,[]).append(abs(ex-ln))
print("mean |W_exact - W_largeNc| over",nsamp,"configurations; last column N^2 x error at N=32 (stays O(1) => 1/N^2 corrections)")
for k in err[Ns[0]]:
    e=[np.mean(err[N][k]) for N in Ns]
    print(f"{k:32s}", "  ".join(f"N={N}: {v:.1e}" for N,v in zip(Ns,e)), f"   N^2*err(32) = {Ns[-1]**2*e[-1]:.2f}")
