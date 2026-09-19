"""Row 18: check the fully expanded B_2-squared bracket, group by group."""
import numpy as np, itertools
from scipy.linalg import expm
s=lambda v:v@v; d=lambda a,b:a@b
def Wji(x,z,w,P,K):
    S=P+K; A,B,C=s(x-z),s(x-w),s(z-w)
    return (np.eye(2)*(A-B)/(2*C)
            +(S/K)*(np.outer(x-z,z-w)/C-np.outer(x-z,x-w)/(2*B))
            +(S/P)*(np.outer(z-w,x-w)/C+np.outer(x-z,x-w)/(2*A)))
rng=np.random.default_rng(31)
print("="*78); print("1.  THE TRANSVERSE BRACKET, GROUP BY GROUP"); print("="*78)
for trial in range(3):
    x,xp,z,w,wp=[rng.normal(size=2)*1.4 for _ in range(5)]
    P,K=(0.37,1.0) if trial==0 else (rng.uniform(.2,3),rng.uniform(.3,2))
    S=P+K; dperp=2.0
    A,B,C=s(x-z),s(x-w),s(z-w); Ap,Bp,Cp=s(xp-z),s(xp-wp),s(z-wp)
    g1=dperp*(A-B)*(Ap-Bp)/(4*C*Cp)
    g2=(S/K)*( d(x-z,z-w)*(Ap-Bp)/(2*C*Cp) - d(x-z,x-w)*(Ap-Bp)/(4*B*Cp)
              +d(xp-z,z-wp)*(A-B)/(2*C*Cp) - d(xp-z,xp-wp)*(A-B)/(4*C*Bp) )
    g3=(S/P)*( d(x-w,z-w)*(Ap-Bp)/(2*C*Cp) + d(x-z,x-w)*(Ap-Bp)/(4*A*Cp)
              +d(xp-wp,z-wp)*(A-B)/(2*C*Cp) + d(xp-z,xp-wp)*(A-B)/(4*C*Ap) )
    g4=(S/K)**2*( d(xp-z,x-z)*d(z-w,z-wp)/(C*Cp) - d(xp-z,x-z)*d(xp-wp,z-w)/(2*C*Bp)
                 -d(x-w,z-wp)*d(x-z,xp-z)/(2*B*Cp) + d(x-z,xp-z)*d(xp-wp,x-w)/(4*B*Bp) )
    g5=(S/P)**2*( d(xp-wp,x-w)*d(z-w,z-wp)/(C*Cp) + d(xp-wp,x-w)*d(xp-z,z-w)/(2*C*Ap)
                 +d(x-w,xp-wp)*d(x-z,z-wp)/(2*A*Cp) + d(x-z,xp-z)*d(xp-wp,x-w)/(4*A*Ap) )
    g6=(S*S/(P*K))*( d(x-w,z-wp)*d(xp-z,z-w)/(C*Cp) - d(x-w,xp-wp)*d(z-w,xp-z)/(2*C*Bp)
                    +d(x-z,xp-z)*d(x-w,z-wp)/(2*A*Cp) - d(xp-z,x-z)*d(xp-wp,x-w)/(4*A*Bp)
                    +d(xp-wp,z-w)*d(x-z,z-wp)/(C*Cp) + d(xp-wp,z-w)*d(x-z,xp-z)/(2*C*Ap)
                    -d(x-z,z-wp)*d(xp-wp,x-w)/(2*B*Cp) - d(xp-wp,x-w)*d(xp-z,x-z)/(4*B*Ap) )
    tot=g1+g2+g3+g4+g5+g6
    ref=np.einsum('ji,ji->',Wji(x,z,w,P,K),Wji(xp,z,wp,P,K))
    # group-by-group against the structure decomposition
    def taus(X,Z,W):
        AA,BB,CC=s(X-Z),s(X-W),s(Z-W)
        return (np.eye(2)*(AA-BB)/(2*CC), np.outer(X-Z,(Z-W)/CC-(X-W)/(2*BB)),
                np.outer((Z-W)/CC+(X-Z)/(2*AA),X-W))
    t1,t2,t3=taus(x,z,w); u1,u2,u3=taus(xp,z,wp); con=lambda a,b:np.einsum('ji,ji->',a,b)
    exact=[con(t1,u1), (S/K)*(con(t2,u1)+con(t1,u2)), (S/P)*(con(t3,u1)+con(t1,u3)),
           (S/K)**2*con(t2,u2), (S/P)**2*con(t3,u3), (S*S/(P*K))*(con(t2,u3)+con(t3,u2))]
    print("\n   trial %d   p+=%.4f k+=%.4f"%(trial+1,P,K))
    for i,(nm,gv,ev) in enumerate(zip(["d_perp  (T1T1')","S/k+    (T1T2'+T2T1')","S/p+    (T1T3'+T3T1')",
            "S^2/k+^2 (T2T2')","S^2/p+^2 (T3T3')","S^2/p+k+ (T2T3'+T3T2')"],[g1,g2,g3,g4,g5,g6],exact)):
        print("     %-24s yours %16.9f   exact %16.9f   diff %.2e"%(nm,gv,ev,abs(gv-ev)))
    print("     %-24s yours %16.9f   exact %16.9f   diff %.2e"%("TOTAL",tot,ref,abs(tot-ref)))

print(); print("="*78); print("2.  THE COLOUR BRACKET, AND WHICH N_c TERM IS RIGHT"); print("="*78)
Nc=3
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2; f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3): f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
def Uadj(al):
    V=expm(1j*sum(al[a]*T[a] for a in range(8)))
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
Ux,Uxp,Uz=Uadj(rng.normal(size=8)),Uadj(rng.normal(size=8)),Uadj(rng.normal(size=8))
# AS YOU WROTE IT: bra carries U^{e'd'}(x')  (rho index FIRST)
C22_asis=np.einsum('axy,azw,bx,fy,bz,we->fe',f,f,Uz,Uxp,Uz,Ux)
print("   as written  (bra U^{e'd'}(x')):  C22 - N_c U^{e'd}(x')U^{de}(x)      = %.2e   <-- your term 4"
      %np.abs(C22_asis-Nc*(Uxp@Ux)).max())
# CONJUGATE ORDERING: bra carries U^{d'e'}(x')  (f index FIRST, matching the ket)
C22_conj=np.einsum('axy,azw,bx,yf,bz,we->fe',f,f,Uz,Uxp,Uz,Ux)
print("   conjugate   (bra U^{d'e'}(x')):  C22 - N_c [U^dag(x')U(x)]^{e'e}     = %.2e"
      %np.abs(C22_conj-Nc*(Uxp.T@Ux)).max())
print("   || N_c U(x')U(x) - N_c U^dag(x')U(x) ||                              = %.4f"
      %np.abs(Nc*(Uxp@Ux)-Nc*(Uxp.T@Ux)).max())
print()
print("   -> your N_c term is the EXACT consequence of the bra ordering you wrote.")
print("      But that ordering is the transpose of the conjugate of the ket, and with")
print("      the ket's ordering the collapse gives the adjoint DIPOLE  U^dag(x')U(x).")
