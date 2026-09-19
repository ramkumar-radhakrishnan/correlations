"""Row 17 = N_2, the (1rho)x(1rho) piece of |B_2|^2.  Prefactor and colour."""
import sympy as sp, numpy as np, itertools
from scipy.linalg import expm
g,pi=sp.symbols('g pi',positive=True); P,K=sp.symbols('pplus kplus',positive=True); S=P+K; I=sp.I
print("="*78); print("1.  PREFACTOR"); print("="*78)
B1=I*g**2*sp.sqrt(P*K)/(4*pi**2*S)
over=sp.Integer(4)/(2*pi)**3
got=sp.simplify(over*B1*sp.conjugate(B1).subs({sp.conjugate(g):g,sp.conjugate(pi):pi,
      sp.conjugate(P):P,sp.conjugate(K):K}))
want=1/(2*pi)**3*g**4/(4*pi**4)*P*K/S**2
print("   (4/(2pi)^3) B^(1) B^(1)dag =",sp.simplify(got))
print("   quoted (1/(2pi)^3)(g^4/4pi^4) p+k+/S^2  ->  ratio",sp.simplify(got/want))
print("   the vertex's 1/S is kept INSIDE the quoted brackets (they carry S/k+ and S/p+),")
print("   and the p+k+/S^2 out front -- consistent.   CORRECT.")

print(); print("="*78); print("2.  COLOUR"); print("="*78)
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2; f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3):
    f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
print("   B^{bc}_2 carries f^{cba} rho^a(x) (free b,c).  Relabel a -> d:")
print("     ket term 1 : U^{ac}(w) f^{cbd} rho^d(x)              -> as you wrote   OK")
print("     ket term 2 : -U^{bc}(z) f^{acd} U^{de}(x) rho^e(x)   -> as you wrote   OK")
print("   different structure constants f^{cbd} vs f^{acd}, so f is NOT an overall")
print("   factor, and you correctly left it inside.   CORRECT.")
print()
print("   BUT the bra's second term reads  f^{ac'd'} U^{bc'}(z) U^{e'd'}(x') rho^{e'}(x')")
print("   while the ket's reads             f^{acd}  U^{bc}(z)  U^{de}(x)  rho^{e}(x).")
print("   U^{e'd'}(x') vs U^{de}(x): the Wilson line indices are TRANSPOSED between them.")
rng=np.random.default_rng(21)
def U(al):
    H=sum(al[a]*T[a] for a in range(8)); V=expm(1j*H)
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
Ux=U(rng.normal(size=8))
print("   ||U(x) - U(x)^T|| = %.4f  -> they are different matrices; make the ordering uniform."
      %np.abs(Ux-Ux.T).max())

print(); print("="*78); print("3.  MEASURE"); print("="*78)
print("   you write int_{z,w,w'} but x and x' appear throughout (in the kernels and in")
print("   rho^{d}(x), rho^{d'}(x')).  It should be  int_{x,x',z,w,w'} .")
