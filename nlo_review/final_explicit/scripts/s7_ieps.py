# i-epsilon (out-state) prescription of Group I Row III (ref3_prescriptions.tex):
# (1) P^i B^{ij} K^j(y-z) in the cos/sin form == the bracket multiplying log((x-z)^2/(x-w)^2) in the row (author's bracket)
# (2) z->w limit of the Delta N kernel: -> P.K(y-w)/(2(z-w)^2) (angular mean)
# (3) the log(k+/Lambda) (Coulomb) term: T(w,w')^dagger = -T(w',w) as operators (SU(3) charges on 3 sites)
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fe_colour import SUN
rng = np.random.default_rng(11)
sq = lambda v: v@v; Kv = lambda v: v/sq(v); cr = lambda a, b: a[0]*b[1]-a[1]*b[0]
out = []
def pr(s): print(s); out.append(s)
w1 = 0
for _ in range(500):
    x, y, z, w, xp, wp = [rng.normal(size=2)*1.5 for _ in range(6)]
    P = Kv(xp-wp); KY = Kv(y-z); a = x-w; b = x-z; A = sq(a); B = sq(b); C = sq(z-w)
    c2 = ((a@b)**2-cr(a, b)**2)/(A*B); s2 = 2*(a@b)*cr(a, b)/(A*B)
    compact = (c2*(P@KY)+s2*cr(P, KY))/(2*C)
    author = P@KY/(2*C)+(P@Kv(z-w))*(KY@Kv(x-z))-(P@Kv(x-w))*(KY@Kv(z-w))-(P@Kv(x-w))*(KY@Kv(x-z))
    w1 = max(w1, abs(compact-author)/abs(author))
pr('(1) P^i B^{ij} K^j(y-z), cos/sin form vs the bracket multiplying the log in Group I Row III: max rel dev %.1e (500 points)' % w1)
x, y, w, xp, wp = [rng.normal(size=2) for _ in range(5)]
P = Kv(xp-wp)
for r in [1e-2, 1e-4]:
    vals = []
    for th in np.linspace(0, 2*np.pi, 400, endpoint=False):
        z = w+r*np.array([np.cos(th), np.sin(th)]); a = x-w; b = x-z
        c2 = ((a@b)**2-cr(a, b)**2)/(sq(a)*sq(b)); s2 = 2*(a@b)*cr(a, b)/(sq(a)*sq(b))
        vals.append((c2*(P@Kv(y-z))+s2*cr(P, Kv(y-z)))/(2*r*r))
    pr('(2) |z-w|=%g: angular mean of 2(z-w)^2 P B K(y-z) = %.8f, P.K(y-w) = %.8f' % (r, 2*r*r*np.mean(vals), P@Kv(y-w)))
# (3) operator identity for the Coulomb term
S = SUN(3); n = S.n; I3 = np.eye(3)
rho = []
for s in range(3):
    m = []
    for a_ in range(n):
        f = [I3, I3, I3]; f[s] = S.t[a_]; m.append(np.kron(np.kron(f[0], f[1]), f[2]))
    rho.append(np.array(m))
pos = [rng.normal(size=2) for _ in range(3)]; Us = [S.adj(S.randV(rng)) for _ in range(3)]
W = {'w': rng.normal(size=2), 'wp': rng.normal(size=2)}; UW = {'w': S.adj(S.randV(rng)), 'wp': S.adj(S.randV(rng))}
def T(wa, wb):   # sum_{x',x,y} K(x'-w_b).K(y-w_a) f^{abd}[U(x')-U(w_b)]^{ab'}rho^{b'}(x') U^{de}(x)rho^e(x) [U(w_a)-U(y)]^{bc} rho^c(y)
    tot = 0
    Q = sum(np.einsum('de,eij->dij', Us[s], rho[s]) for s in range(3))
    for sp in range(3):
        Lop = np.einsum('ab,bij->aij', Us[sp]-UW[wb], rho[sp])
        for sy in range(3):
            Rop = np.einsum('bc,cij->bij', UW[wa]-Us[sy], rho[sy])
            kk = Kv(pos[sp]-W[wb])@Kv(pos[sy]-W[wa])
            tot = tot+kk*np.einsum('abd,aij,djk,bkl->il', S.f, Lop, Q, Rop, optimize=True)
    return tot
T1 = T('w', 'wp'); T2 = T('wp', 'w')
pr('(3) Coulomb term: |T(w,w\')^dagger + T(w\',w)| / |T| = %.1e  (|T| = %.2f); i.e. the log(k+/Lambda) term cancels against its complex conjugate' % (np.abs(T1.conj().T+T2).max()/np.abs(T1).max(), np.abs(T1).max()))
open('s7_ieps.log', 'w').write('\n'.join(out)+'\n')
