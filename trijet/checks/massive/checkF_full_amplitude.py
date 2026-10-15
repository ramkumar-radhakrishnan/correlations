# Check F: complete LO amplitude gamma*_T -> q qbar g (before the shock wave, regular + instantaneous, BOTH emitters),
# computed from Feynman numerators with Lepage-Brodsky spinors, versus the two-component massive expressions.
#   physical:  A = ubar(k1)[ eps_g (K+m) eps_gam /(K^2-m^2) + eps_gam (-P+m) eps_g /(P^2-m^2) ] v(k2) / (q^- - sum k_i^-)
#              with K = q-k2, P = q-k1 off shell (minus components fixed by the photon vertex)
#   two-component (user's conventions + mass terms fitted in checkE):
#     T_qbar = - z zb   sum_lam V_gam V_g / [(kt^2+eps^2) den3]  + z xi eta (d^{lr} + 2i eps^{lr} lam1) d_{l1 l2} / den3
#     T_q    = + e (1-e) sum_lam W_g W_gam / [(K^2+eps_q^2) den3q] - z xi eta (d^{lr} - 2i eps^{lr} lam1) d_{l1 l2} / den3q
#   Requirement: A = C * (T_qbar + T_q) with ONE constant C for all helicities, polarisations, momenta, masses,
#   and separately A_qbar = C T_qbar, A_q = C T_q.
import sympy as sp, itertools, random
exec(open('lb.py').read())
def slash(ap,am,ax,ay): return R(1,2)*(gp*am+gm*ap) - g[0]*ax - g[1]*ay
def eps2(a,b): return [[0,1],[-1,0]][a][b]
def d(a,b): return 1 if a==b else 0
def varphi(i,j,l1,l2,zq): return d(l1,l2)*((2*zq-1)*d(i,j) + 2*I*eps2(i,j)*l1)
def tau(mm,n,l,l2,xi_,par_d): return d(l,l2)*((2*par_d+xi_)*d(mm,n) + 2*I*xi_*eps2(mm,n)*l)
def mu(p,a,b): return d(b,-a)*(2*a*d(p,0) - I*d(p,1))
H=(R(1,2),R(-1,2))
random.seed(3)
def rnd(): return R(random.randint(-9,9),random.randint(1,7))
for (z,x,qp,Q2,m) in [(R(1,3),R(1,5),R(7,4),R(3,2),R(2,3)), (R(2,7),R(3,11),R(5,3),R(1,4),R(5,4)),
                      (R(1,2),R(1,3),R(1,1),R(2,1),R(0)), (R(3,5),R(1,10),R(9,4),R(1,1),R(1,7))]:
    zb=1-z; e=1-z-x; consts=set()
    for trial in range(2):
        k1=[rnd(),rnd()]; k3=[rnd(),rnd()]; k2=[-k1[0]-k3[0],-k1[1]-k3[1]]
        mom={'k1':(z*qp,*k1),'k2':(e*qp,*k2),'k3':(x*qp,*k3)}
        minus=lambda p_,mass: (p_[1]**2+p_[2]**2+mass**2)/p_[0]
        qm=-Q2/qp
        E3 = qm - minus(mom['k1'],m) - minus(mom['k2'],m) - minus(mom['k3'],0)
        P=(zb*qp, qm-minus(mom['k1'],m), -k1[0], -k1[1]); Kv=((1-e)*qp, qm-minus(mom['k2'],m), -k2[0], -k2[1])
        PP=P[0]*P[1]-P[2]**2-P[3]**2-m**2; KK=Kv[0]*Kv[1]-Kv[2]**2-Kv[3]**2-m**2
        # two-component variables
        kt=k1; pt=[k3[0]+x/zb*k1[0], k3[1]+x/zb*k1[1]]
        Kt=[-k2[0],-k2[1]]; pq=[k3[0]+x/(1-e)*k2[0], k3[1]+x/(1-e)*k2[1]]
        epsq2=z*zb*Q2+m**2; epsq2q=e*(1-e)*Q2+m**2; om=x*m/zb; omq=x*m/(1-e)
        den3 = x*e*(kt[0]**2+kt[1]**2+epsq2) + z*zb**2*(pt[0]**2+pt[1]**2+om**2)
        den3q= x*z*(Kt[0]**2+Kt[1]**2+epsq2q) + e*(1-e)**2*(pq[0]**2+pq[1]**2+omq**2)
        for L1,L2,l,r in itertools.product((1,-1),(1,-1),(0,1),(0,1)):
            l1,l2=R(L1,2),R(-L2,2)
            ub=bar(u(*mom['k1'],m,L1)); vv=v(*mom['k2'],m,L2)
            eg=slash_eps(r,k3[0],k3[1],x*qp); ea=slash_eps(l,0,0,qp)
            Aq =(ub*eg*(slash(*Kv)+m*sp.eye(4))*ea*vv)[0]/KK/E3
            Aqb=(ub*ea*(-slash(*P)+m*sp.eye(4))*eg*vv)[0]/PP/E3
            Tqb = -z*zb*sum((sum(varphi(l,j,l1,lam,z)*kt[j] for j in (0,1))+m*mu(l,l1,lam))
                            *(sum(tau(r,n,lam,l2,x,e)*pt[n] for n in (0,1))+x*om*mu(r,lam,l2)) for lam in H) \
                  /((kt[0]**2+kt[1]**2+epsq2)*den3) \
                  + z*x*e*(d(l,r)+2*I*eps2(l,r)*l1)*d(l1,l2)/den3
            Tq  = +e*(1-e)*sum((sum(tau(n,r,l1,lam,x,z)*pq[n] for n in (0,1))+x*omq*mu(r,l1,lam))
                            *(sum(varphi(l,j,lam,l2,z+x)*Kt[j] for j in (0,1))+m*mu(l,lam,l2)) for lam in H) \
                  /((Kt[0]**2+Kt[1]**2+epsq2q)*den3q) \
                  - z*x*e*(d(l,r)-2*I*eps2(l,r)*l1)*d(l1,l2)/den3q
            for A_,T_ in ((Aqb,Tqb),(Aq,Tq),(Aqb+Aq,Tqb+Tq)):
                A_=sp.nsimplify(sp.simplify(A_)); T_=sp.simplify(T_)
                if T_==0:
                    consts.add('OK0' if A_==0 else 'MISMATCH(T=0)')
                else:
                    consts.add(sp.simplify(A_/T_/qp))
    consts.discard('OK0')
    print(f"zeta={z}, xi={x}, q+={qp}, Q2={Q2}, m={m}:  A/(q^+ T) over 2 momenta x 16 hel/pol x (qbar, q, sum) =", consts,
          " ; divided by sqrt(zeta*eta):", {sp.nsimplify(c/sp.sqrt(z*e)) for c in consts if c!='MISMATCH(T=0)'})
