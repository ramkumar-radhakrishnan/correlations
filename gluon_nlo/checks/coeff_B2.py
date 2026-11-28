"""Check the momentum-space B2 (f^{abc} part) against (2pi)^3 x [symmetrised psi_ggrho amplitude]/K+.
Gluon a: (ka+, ka) with polarisation i; gluon b: (kb+, kb) with polarisation j.
Common factor i g^2 f^{abc} rho^c(-ka-kb) stripped."""
import numpy as np
rng=np.random.default_rng(3)
d=np.eye(2)
def ours(kap,ka,kbp,kb):
    Kp=kap+kbp; xi=kap/Kp; K=ka+kb; P=ka-xi*K
    # psi_ggrho bracket, gluon a=(xi K+, xi K+P) index i, gluon b=((1-xi)K+,(1-xi)K-P) index j
    br=np.zeros((2,2))
    A=xi*K+P; B=(1-xi)*K-P
    for i in range(2):
        for j in range(2):
            br[i,j]=(d[i,j]/2*((1-2*xi)*K@K-2*K@P)+K[j]*P[i]/xi+K[i]*P[j]/(1-xi)
                     +K@K*A[i]*B[j]*(xi*B@B-(1-xi)*A@A)/(2*xi*(1-xi)*(A@A)*(B@B)))
    return -np.sqrt(xi*(1-xi))/(Kp*(xi*(1-xi)*K@K+P@P)*(K@K))*br   # 8pi^3 * (-1/8pi^3)
def user(kp,k,pp,p):  # B^{ba}_{2ji}(k+,k;p+,p): (a,i)<->k, (b,j)<->p ; returns [i,j]
    out=np.zeros((2,2)); s=p+k
    for i in range(2):
        for j in range(2):
            T1=p[j]*k[i]*np.sqrt(pp*kp)/(2*(pp*k@k+kp*p@p))*(1/(kp*p@p)-1/(k@k*pp))
            T2=0
            for l in range(2):
                for m in range(2):
                    T2+=-(s[l]*(p[m]*kp-pp*k[m]))/((s@s)*np.sqrt(pp*kp)*(pp+kp)*(p@p/pp+k@k/kp))*(
                        d[l,m]*d[j,i]/(pp+kp)-d[i,m]*d[j,l]/kp-d[j,m]*d[i,l]/pp)
            T3=-(pp-kp)*d[j,i]/(2*np.sqrt(pp*kp)*(pp+kp)**2*(k@k/kp+p@p/pp))
            out[i,j]=T1+T2+T3
    return out
for t in range(3):
    kap,kbp=rng.uniform(0.2,2,2); ka,kb=rng.normal(size=(2,2))
    o=ours(kap,ka,kbp,kb); osw=ours(kbp,kb,kap,ka)
    sym=0.5*(o-osw.T)   # exchange: f^{bac}=-f^{abc}
    u=user(kap,ka,kbp,kb)
    print("raw  ratio user/ours:",np.round(u/o,4).ravel())
    print("sym  ratio user/sym :",np.round(u/sym,4).ravel())
