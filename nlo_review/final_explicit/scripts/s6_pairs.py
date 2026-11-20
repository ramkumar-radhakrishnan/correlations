# Pairs B and C: the two Group I Row V recoil terms, written in the row variables, renamed to the common
# variables (soft z, hard w), equal the common form kernel x colour used in Sec. C (pointwise, random points,
# random SU(3) Wilson lines).  Also: Group III members after y'->w' (pairs A,B,C,D) kernel forms.
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fe_colour import SUN
from scipy import special
dot = lambda a, b: a @ b
sq = lambda a: a @ a
Kv = lambda v: v/sq(v)
Kab = lambda a, b, z: dot(Kv(a-z), Kv(b-z))
def Ecal(s):
    si, ci = special.sici(abs(s))
    return ci - np.euler_gamma - np.log(abs(s)) - 1j*np.sign(s)*si
E = lambda *a: np.einsum(*a, optimize=True)
S = SUN(3); rng = np.random.default_rng(5); f = S.f; N = 3
k = np.array([0.8, -0.5]); out = []
def C_I5(U, x, z, y, yp):   # tensor (x', y, y') of the Row V first-part colour, U at named points
    L = U['xp']-U['wp']
    return E('ga,gde,dc,eb->abc', L, f, U[x], U[z])-E('ga,gde,eb,dc->abc', L, f, U[y], U[yp])
def C_B(U):  # tensor (x', y, x)
    L = U['xp']-U['wp']
    return E('ga,gdc,db,ce->aeb', L, f, U['z'], U['w'])-E('ga,gbc,cd,be->ade', L, f, U['y'], U['x'])
sym = lambda T: 0.5*(T+T.transpose(0, 2, 1))
wB = wC = 0
for _ in range(20):
    P = {p: rng.normal(size=2)*1.5 for p in ['x', 'y', 'z', 'w', 'xp', 'wp']}
    U = {p: S.adj(S.randV(rng)) for p in P}
    KP = Kv(P['xp']-P['wp'])
    X, Y, Z, W = P['x'], P['y'], P['z'], P['w']
    # pair B. Term1 (row: soft x, hard z, charges y, y'): x->Z, z->W, y'->X ; colour U(x)->U(z), U(z)->U(w)
    s1 = k@(Z-W)
    K1 = Ecal(s1)*dot(KP, Kv(Y-W))*dot(Kv(W-Z), Kv(X-Z))
    T1 = C_I5(U, 'z', 'w', 'y', 'x')                      # charges (x', y, X)
    # Term2 (row: hard x, soft z): x->W, z->Z, y->X, y'->Y ; e^{-is}E(-s): E(k(z-x)) -> E(k(Z-W))
    K2 = Ecal(k@(Z-W))*dot(KP, Kv(Y-W))*dot(Kv(Z-W), Kv(X-Z))
    T2 = C_I5(U, 'w', 'z', 'x', 'y').transpose(0, 2, 1)   # charges (x', X, Y) -> reorder to (x', Y, X)
    lhs = -0.5j*(K1*sym(T1)+K2*sym(T2))
    rhs = -1j*Ecal(k@(Z-W))*Kab(X, W, Z)*dot(KP, Kv(Y-W))*sym(C_B(U))
    wB = max(wB, np.abs(lhs-rhs).max()/np.abs(rhs).max())
    # pair C. colour of row: [U(x')-U(w')]^{ab'}[-f^{adc}f^{fbe}U^{db}(x)U^{ce}(z)+N U^{af}(y)]
    L = U['xp']-U['wp']
    col = lambda a, b, c: E('aB,af->Bf', L, N*U[c])-E('aB,adc,fbe,db,ce->Bf', L, f, f, U[a], U[b])
    B2 = lambda x, y, z: np.outer(y-z, x-z)/sq(x-z)+np.outer(y-z, y-x)/(2*sq(y-x))
    B1 = lambda x, y, z: np.outer(y-x, x-z)/sq(x-z)-np.outer(y-x, y-z)/(2*sq(y-z))
    # Term1 (soft x, hard z, charge y): x->Z, z->W, y->X ; -E(s) P K(z-x) B2/(y-z)^2
    t1 = -Ecal(k@(Z-W))*(KP@B2(Z, X, W)@Kv(W-Z))/sq(X-W)*col('z', 'w', 'x')
    # Term2 (hard x, soft z): x->W, y->X ; -e^{-is}E(-s) P K(z-x) B1/(y-x)^2
    t2 = -Ecal(k@(Z-W))*(KP@B1(W, X, Z)@Kv(Z-W))/sq(X-W)*col('w', 'z', 'x')
    CC = E('aB,ah->Bh', L, N*U['x'])-E('aB,adc,hbe,db,ce->Bh', L, f, f, U['z'], U['w'])
    rhsC = -Ecal(k@(Z-W))*dot(KP, Kv(X-W))*(Kab(X, W, Z)-2/sq(Z-W))*CC
    wC = max(wC, np.abs(t1+t2-rhsC).max()/np.abs(rhsC).max())
s = 'pair B: Group I Row V recoil terms renamed to common variables vs common form: max rel dev %.1e; pair C: %.1e (20 random configurations, SU(3))' % (wB, wC)
print(s); open('s6_pairs.log', 'w').write(s+'\n')
