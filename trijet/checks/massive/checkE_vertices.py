# Check E: massive two-component vertices in the user's conventions, fitted to Lepage-Brodsky spinors.
#  antiquark channel:  N_q̄(m) = sum_{lam'} [ubar(k1) eps_l v(P,lam')][vbar(P,lam') eps_r v(k2)]
#                         = c1 * sum_lam  V_gam^l_{lam1 lam}(k~) V_g^r_{lam lam2}(p~)
#       V_gam^l = varphi^{lj}(zeta) k~^j + m Gam^l ,   V_g^r = tau^{rn}(xi,eta) p~^n + m Gp^r
#  quark channel:      N_q(m) = sum_{lam'} [ubar(k1) eps_r u(K,lam')][ubar(K,lam') eps_l v(k2)]
#                         = c1q * sum_lam  W_g^r_{lam1 lam}(p~q) W_gam^l_{lam lam2}(K)
#       W_g^r = tau^{nr}(xi,zeta) p~q^n + m Hp^r ,  W_gam^l = varphi^{lj}(zeta+xi) K^j + m H^l
import sympy as sp, itertools
exec(open('lb.py').read())
z,x = R(1,3), R(1,5); zb=1-z; eta=1-z-x; qp=R(7,4)
kx,ky,px,py,m = sp.symbols('kx ky px py m', real=True)
def eps2(a,b): return [[0,1],[-1,0]][a][b]
def d(a,b): return 1 if a==b else 0
def varphi(i,j,l1,l2,zq): return d(l1,l2)*((2*zq-1)*d(i,j) + 2*I*eps2(i,j)*l1)
def tau(mm,n,l,l2,xi_,par_d): return d(l,l2)*((2*par_d+xi_)*d(mm,n) + 2*I*xi_*eps2(mm,n)*l)
H=(R(1,2),R(-1,2))
# --- antiquark channel kinematics (q_perp=0): k1=k~, P=-k~, k3 = p~ + (xi/zb) P, k2 = P-k3
k1=[kx,ky]; P=[-kx,-ky]; k3=[px+x/zb*P[0], py+x/zb*P[1]]; k2=[P[0]-k3[0],P[1]-k3[1]]
Nqb={}
for L1,L2,l,r in itertools.product((1,-1),(1,-1),(0,1),(0,1)):
    Nqb[(L1,L2,l,r)] = sp.expand(sum(
        (bar(u(z*qp,*k1,m,L1))*slash_eps(l,0,0,qp)*v(zb*qp,*P,m,Lp))[0] *
        (bar(v(zb*qp,*P,m,Lp))*slash_eps(r,*k3,x*qp)*v(eta*qp,*k2,m,L2))[0] for Lp in (1,-1)))
# --- quark channel kinematics: K=-k2 (parent quark, fraction 1-eta), k2 = antiquark;  k3 = pq + (xi/(1-eta)) K, k1 = K-k3
Kx,Ky=kx,ky   # reuse symbols: (kx,ky)=K, (px,py)=p~q
K=[Kx,Ky]; k2q=[-Kx,-Ky]; k3q=[px+x/(1-eta)*K[0], py+x/(1-eta)*K[1]]; k1q=[K[0]-k3q[0],K[1]-k3q[1]]
Nq={}
for L1,L2,l,r in itertools.product((1,-1),(1,-1),(0,1),(0,1)):
    Nq[(L1,L2,l,r)] = sp.expand(sum(
        (bar(u(z*qp,*k1q,m,L1))*slash_eps(r,*k3q,x*qp)*u((1-eta)*qp,*K,m,Lp))[0] *
        (bar(u((1-eta)*qp,*K,m,Lp))*slash_eps(l,0,0,qp)*v(eta*qp,*k2q,m,L2))[0] for Lp in (1,-1)))
# unknown mass vertices  (index: [pol][lamA][lamB])
def unk(name): return {(p,a,b): sp.Symbol(f'{name}_{p}_{"p" if a>0 else "m"}{"p" if b>0 else "m"}') for p in (0,1) for a in H for b in H}
Gam,Gp,Hm,Hp = unk('Gam'),unk('Gp'),unk('H'),unk('Hp')
c1=15*sp.sqrt(35)/7
def model_qb(l1,l2,l,r):
    return c1*sum((sum(varphi(l,j,l1,lam,z)*[kx,ky][j] for j in (0,1)) + m*Gam[(l,l1,lam)]) *
                  (sum(tau(r,n,lam,l2,x,eta)*[px,py][n] for n in (0,1)) + m*Gp[(r,lam,l2)]) for lam in H)
c1q=sp.Symbol('c1q')
def model_q(l1,l2,l,r):
    return c1q*sum((sum(tau(n,r,l1,lam,x,z)*[px,py][n] for n in (0,1)) + m*Hp[(r,l1,lam)]) *
                   (sum(varphi(l,j,lam,l2,z+x)*[kx,ky][j] for j in (0,1)) + m*Hm[(l,lam,l2)]) for lam in H)
def solve(N, model, extra=()):
    eqs=[]
    for (L1,L2,l,r),val in N.items():
        l1,l2 = R(L1,2), R(-L2,2)
        diff = sp.expand(val - model(l1,l2,l,r))
        eqs += sp.Poly(diff, kx,ky,px,py,m).coeffs()
    return eqs
# antiquark: linear (order m) equations first, then check order m^2
eqs=solve(Nqb, model_qb)
lin=[e for e in eqs if sp.Poly(e,*Gam.values(),*Gp.values()).total_degree()<=1]
sol=sp.solve(lin, list(Gam.values())+list(Gp.values()), dict=True)
print("antiquark-channel solutions:", len(sol))
s=sol[0]; free=[v_ for v_ in list(Gam.values())+list(Gp.values()) if v_ not in s]
print("  free params:", free)
rest=[sp.simplify(e.subs(s)) for e in eqs]; rest=[e for e in rest if e!=0]
print("  remaining (m^2) equations:", rest)
for k_,v_ in s.items():
    if v_!=0: print("   ",k_,"=",sp.nsimplify(sp.simplify(v_/m if False else v_)))
# quark channel: massless part first -> c1q
eqsq=solve(Nq, model_q)
m0=[sp.simplify(e) for e in eqsq if not any(e.has(v_) for v_ in list(Hm.values())+list(Hp.values()))]
print("\nquark channel massless eqs ->", sp.solve([e for e in m0 if e!=0], c1q))
eqsq=[sp.expand(e.subs(c1q,c1)) for e in eqsq]
unkq=list(Hm.values())+list(Hp.values())
linq=[e for e in eqsq if e!=0 and sp.Poly(e,*unkq).total_degree()<=1]
solq=sp.solve(linq, unkq, dict=True)
print("quark-channel solutions:", len(solq)); sq=solq[0]
print("  free params:", [v_ for v_ in unkq if v_ not in sq])
print("  remaining (m^2) equations:", [e for e in (sp.simplify(e.subs(sq)) for e in eqsq) if e!=0])
for k_,v_ in sq.items():
    if v_!=0: print("   ",k_,"=",v_)
# ---- compact forms:  mu^p_{ab} = delta_{b,-a} (2a delta^{p1} - i delta^{p2})
def mu(p,a,b): return d(b,-a)*(2*a*d(p,0) - I*d(p,1))
ok=True
for p in (0,1):
  for a in H:
    for b in H:
      ok &= sp.simplify(s[Gam[(p,a,b)]] - mu(p,a,b))==0
      ok &= sp.simplify(s[Gp[(p,a,b)]]  - (x**2/zb)*mu(p,a,b))==0
      ok &= sp.simplify(sq[Hm[(p,a,b)]] - mu(p,a,b))==0
      ok &= sp.simplify(sq[Hp[(p,a,b)]] - (x**2/(1-eta))*mu(p,a,b))==0
print("\ncompact forms: Gam = mu, Gp = (xi^2/(1-zeta)) mu, H = mu, Hp = (xi^2/(zeta+xi)) mu,"
      "  mu^p_{ab} = delta_{b,-a}(2a delta^{p1} - i delta^{p2}) :", ok)
