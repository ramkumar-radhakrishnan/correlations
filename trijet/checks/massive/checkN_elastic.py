# Check N: Fourier transform of the elastic (no-F) terms.
#  (a) the elastic part of all 12 Gaussian Wilson-line structures factorises: [amplitude dipoles - 1][conjugate dipoles - 1]
#  (b) the elastic amplitude written with the integrals T^{jn}, T_k, T_p, T_0, T_I, and the helicity sums
#      H = H^(0) + H^(m) (antiquark), Ht (quark), Hx (interference), against brute-force sums with a discrete dipole measure
import sympy as sp, numpy as np, itertools
# ---------------- (a)
S=lambda a,b: sp.Symbol(f'S_{a}{b}')
def Qg(a,b,c,d): return S(a,b)*S(c,d)          # elastic (F-free) part of the Gaussian quadrupole
W={
 'qbar W(x,y,z;..)':   (Qg('x','z','zb','xb')*Qg('z','y','yb','zb')-S('x','z')*S('z','y')-S('zb','xb')*S('yb','zb')+1, (S('x','z')*S('z','y')-1)*(S('zb','xb')*S('yb','zb')-1)),
 "qbar W(x,y',y';..)": (Qg('x','yp','zb','xb')*S('yb','zb')-S('x','yp')-S('zb','xb')*S('yb','zb')+1, (S('x','yp')-1)*(S('zb','xb')*S('yb','zb')-1)),
 "qbar W(..;xb,yb',yb')":(Qg('x','z','ybp','xb')*S('z','y')-S('x','z')*S('z','y')-S('ybp','xb')+1, (S('x','z')*S('z','y')-1)*(S('ybp','xb')-1)),
 "qbar W(x,y',y';xb,yb',yb')":(Qg('x','yp','ybp','xb')-S('x','yp')-S('ybp','xb')+1, (S('x','yp')-1)*(S('ybp','xb')-1)),
 'q W(y,x,z;..)':      (Qg('y','z','zb','yb')*Qg('z','x','xb','zb')-S('y','z')*S('z','x')-S('zb','yb')*S('xb','zb')+1, (S('y','z')*S('z','x')-1)*(S('zb','yb')*S('xb','zb')-1)),
 "q W(y,x',x';..)":    (Qg('y','xp','zb','yb')*S('xb','zb')-S('y','xp')-S('zb','yb')*S('xb','zb')+1, (S('y','xp')-1)*(S('zb','yb')*S('xb','zb')-1)),
 "q W(..;yb,xb',xb')": (Qg('y','z','xbp','yb')*S('z','x')-S('y','z')*S('z','x')-S('xbp','yb')+1, (S('y','z')*S('z','x')-1)*(S('xbp','yb')-1)),
 "q W(y,x',x';yb,xb',xb')":(Qg('y','xp','xbp','yb')-S('y','xp')-S('xbp','yb')+1, (S('y','xp')-1)*(S('xbp','yb')-1)),
 'int W(x,y,z;..)':    (Qg('x','z','zb','xb')*Qg('z','y','yb','zb')-S('x','z')*S('z','y')-S('zb','xb')*S('yb','zb')+1, (S('x','z')*S('z','y')-1)*(S('zb','xb')*S('yb','zb')-1)),
 "int W(x',y,x';..)":  (S('zb','xb')*Qg('xp','y','yb','zb')-S('xp','y')-S('zb','xb')*S('yb','zb')+1, (S('xp','y')-1)*(S('zb','xb')*S('yb','zb')-1)),
 "int W(..;xb,yb',yb')":(Qg('x','z','ybp','xb')*S('z','y')-S('x','z')*S('z','y')-S('ybp','xb')+1, (S('x','z')*S('z','y')-1)*(S('ybp','xb')-1)),
 "int W(x',y,x';xb,yb',yb')":(S('ybp','xb')*S('xp','y')-S('xp','y')-S('ybp','xb')+1, (S('xp','y')-1)*(S('ybp','xb')-1)),
}
print("(a) elastic part factorises:", all(sp.expand(a-b)==0 for a,b in W.values()))
# ---------------- (b)
rng=np.random.default_rng(11); H=(0.5,-0.5); E2=np.array([[0,1],[-1,0]])
def d(a,b): return 1.0 if a==b else 0.0
def varphi(i,j,l1,l2,zq): return d(l1,l2)*((2*zq-1)*d(i,j)+2j*E2[i,j]*l1)
def tau(mm,n,l,l2,xi,pd): return d(l,l2)*((2*pd+xi)*d(mm,n)+2j*xi*E2[mm,n]*l)
def mu(p,a,b): return d(b,-a)*(2*a*d(p,0)-1j*d(p,1))
hl=list(itertools.product(H,H,(0,1),(0,1)))
def run(z,x,Q,m):
    e=1-z-x; zb=1-z; om=x*m/zb; omt=x*m/(1-e); eps2=z*zb*Q*Q+m*m; epst2=e*(1-e)*Q*Q+m*m
    # momentum-space vertex functions (antiquark: Phi1, Phi2; quark: Pt1, Pt2), as in Sec. 4
    def Nq(k,p): return np.array([sum((sum(varphi(l,j,l1,lam,z)*k[j] for j in (0,1))+m*mu(l,l1,lam))*(sum(tau(r,n,lam,l2,x,e)*p[n] for n in (0,1))+x*om*mu(r,lam,l2)) for lam in H) for l1,l2,l,r in hl])
    def Nt(K,P): return np.array([sum((sum(tau(n,r,l1,lam,x,z)*P[n] for n in (0,1))+x*omt*mu(r,l1,lam))*(sum(varphi(l,j,lam,l2,z+x)*K[j] for j in (0,1))+m*mu(l,lam,l2)) for lam in H) for l1,l2,l,r in hl])
    I =np.array([(d(l,r)+2j*E2[l,r]*l1)*d(l1,l2) for l1,l2,l,r in hl]); It=np.array([(d(l,r)-2j*E2[l,r]*l1)*d(l1,l2) for l1,l2,l,r in hl])
    D =lambda k,p: z*zb**2*(p@p+om**2)+x*e*(k@k+eps2);  Dt=lambda K,P: e*(1-e)**2*(P@P+omt**2)+x*z*(K@K+epst2)
    f1=lambda k,p: zb/((k@k+eps2)*D(k,p)); f2=lambda k,p: 1/(z*zb*(k@k+eps2)*(p@p+om**2)); g=lambda k,p: 1/D(k,p)
    ft1=lambda K,P: (1-e)/((K@K+epst2)*Dt(K,P)); ft2=lambda K,P: 1/(e*(1-e)*(K@K+epst2)*(P@P+omt**2)); gt=lambda K,P: 1/Dt(K,P)
    Phi1=lambda k,p: Nq(k,p)*f1(k,p)-x*e*I*g(k,p); Phi2=lambda k,p: Nq(k,p)*f2(k,p)
    Pt1=lambda K,P: Nt(K,P)*ft1(K,P)-x*z*It*gt(K,P); Pt2=lambda K,P: Nt(K,P)*ft2(K,P)
    # final momenta with k1+k2+k3=0; discrete dipole measure: D(l) -> weights w_i at points l_i with sum w = 1
    k1,k3=rng.normal(size=2),rng.normal(size=2); k2=-k1-k3
    kt=k1; pt=k3+x/zb*k1; Kq=-k2; Pq=k3+x/(1-e)*k2
    L=[rng.normal(size=2) for _ in range(3)]; w=rng.random(3); w/=w.sum()
    # elastic amplitudes, brute force
    Ab=sum(w[a]*w[b]*(Phi1(kt-L[a],pt+e/zb*L[a]-L[b])-Phi2(kt-L[a],pt)) for a in range(3) for b in range(3))-(Phi1(kt,pt)-Phi2(kt,pt))
    At=sum(w[a]*w[b]*(Pt1(Kq-L[b],Pq+L[a]-z/(1-e)*L[b])-Pt2(Kq-L[b],Pq)) for a in range(3) for b in range(3))-(Pt1(Kq,Pq)-Pt2(Kq,Pq))
    # the T integrals
    def Ts(kA,pA,kB,pB,F1,F2,G,kb,pb):
        T=np.zeros((2,2)); Tk=np.zeros(2); Tp=np.zeros(2); T0=0.; TI=0.
        for a in range(3):
            for b in range(3):
                ww=w[a]*w[b]; ka,pa=kA(a,b),pA(a,b); kB_,pB_=kB(a,b),pB(a,b)
                T =T +ww*(np.outer(ka,pa)*F1(ka,pa)-np.outer(kB_,pB_)*F2(kB_,pB_))
                Tk=Tk+ww*(ka*F1(ka,pa)-kB_*F2(kB_,pB_)); Tp=Tp+ww*(pa*F1(ka,pa)-pB_*F2(kB_,pB_))
                T0=T0+ww*(F1(ka,pa)-F2(kB_,pB_));        TI=TI+ww*G(ka,pa)
        d0=F1(kb,pb)-F2(kb,pb)
        return T-np.outer(kb,pb)*d0, Tk-kb*d0, Tp-pb*d0, T0-d0, TI-G(kb,pb)
    T,Tk,Tp,T0,TI=Ts(lambda a,b: kt-L[a], lambda a,b: pt+e/zb*L[a]-L[b], lambda a,b: kt-L[a], lambda a,b: pt, f1,f2,g,kt,pt)
    U,Uk,Up,U0,UI=Ts(lambda a,b: Kq-L[b], lambda a,b: Pq+L[a]-z/(1-e)*L[b], lambda a,b: Kq-L[b], lambda a,b: Pq, ft1,ft2,gt,Kq,Pq)
    ee=lambda A,B: A[0,0]*np.conj(B[1,1])+A[1,1]*np.conj(B[0,0])-A[0,1]*np.conj(B[1,0])-A[1,0]*np.conj(B[0,1])   # eps^{ij}eps^{mn}A^{im}B^{jn*}
    tt=lambda A,B: np.sum(A*np.conj(B)); tr=lambda A: A[0,0]+A[1,1]; vv=lambda a,b: a@np.conj(b)
    # antiquark: H = H0 + Hm
    H0=(8*(z*z+zb*zb)*(e*e+zb*zb)*tt(T,T)-8*x*(2*z-1)*(2*zb-x)*ee(T,T)-16*z*x*e*e*np.real(tr(T)*np.conj(TI))+8*x*x*e*e*abs(TI)**2)
    Hm=(8*m*m*(e*e+zb*zb)*vv(Tp,Tp)+8*x*x*om*om*(z*z+zb*zb)*vv(Tk,Tk)+8*m*m*x*x*om*om*abs(T0)**2
        -16*z*e*m*x*om*np.real(tr(T)*np.conj(T0)-vv(Tk,Tp))+16*x*e*m*x*om*np.real(T0*np.conj(TI)))
    # quark
    Ht0=(8*(z*z+(1-e)**2)*(e*e+(1-e)**2)*tt(U,U)-8*x*(2*e-1)*(2*(1-e)-x)*ee(U,U)+16*z*z*x*e*np.real(tr(U)*np.conj(UI))+8*x*x*z*z*abs(UI)**2)
    Htm=(8*m*m*(z*z+(1-e)**2)*vv(Up,Up)+8*x*x*omt*omt*(e*e+(1-e)**2)*vv(Uk,Uk)+8*m*m*x*x*omt*omt*abs(U0)**2
         +16*z*e*m*x*omt*np.real(tr(U)*np.conj(U0)-vv(Uk,Up))+16*x*z*m*x*omt*np.real(U0*np.conj(UI)))
    # interference (quark amplitude x antiquark conjugate)
    Hx0=(8*(1-x-2*z*e)*(x*(1-x)+2*z*e)*tt(U,T)-8*x*(1-2*z-x)**2*ee(U,T)-8*x*e*(z+x)**2*tr(U)*np.conj(TI)+8*x*z*zb**2*UI*np.conj(tr(T)))
    Hxm=(8*m*m*(x+2*z*e)*vv(Up,Tp)-8*(x+2*z*e)*x*x*omt*om*vv(Uk,Tk)-8*(z+x)**2*m*x*om*tr(U)*np.conj(T0)
         +8*zb**2*m*x*omt*U0*np.conj(tr(T))-8*e*e*m*x*omt*vv(Uk,Tp)+8*z*z*m*x*om*vv(Up,Tk))
    out={'antiquark |A|^2':(np.sum(abs(Ab)**2),H0+Hm),'quark |A|^2':(np.sum(abs(At)**2),Ht0+Htm),'interference A_q A_qbar*':(np.sum(At*np.conj(Ab)),Hx0+Hxm)}
    return {k:abs(a-b)/abs(a) for k,(a,b) in out.items()}, (abs(Hm/H0),abs(Htm/Ht0))
for pars in [(0.3,0.2,1.1,0.8),(0.45,0.15,0.7,1.6),(0.2,0.5,1.4,0.0)]:
    for t in range(2):
        r,ratio=run(*pars)
        print("(b) zeta,xi,Q,m =",pars," rel. diff brute force vs formula:",{k:f"{v:.0e}" for k,v in r.items()},
              f"  |H^(m)/H^(0)| = {ratio[0]:.2f}, {ratio[1]:.2f}")
