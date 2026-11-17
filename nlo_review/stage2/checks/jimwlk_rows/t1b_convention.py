# sensitivity of the row-by-row comparison to the reordering convention:
#  (S) full symmetrisation [used by the author]  vs  (O) reorder to the fixed order rho(x') rho(y) rho(x) (full commutator per inversion)
import sys; sys.path.insert(0,".")
from t1_sec2 import *
REF = {"x'": 0, "y": 1, "x": 2}
def two_from_three_ref(monos, world, phys=("P0","P1")):
    M = np.zeros((world.n, world.n), complex)
    for m in monos:
        pts = m.pts()
        for i, j in ((0,1),(0,2),(1,2)):
            if REF[pts[i]] < REF[pts[j]]: continue      # already in reference order
            k = 3 - i - j
            for pm, pn in (phys, phys[::-1]):
                a = {pts[i]: pm, pts[j]: pm, pts[k]: pn}
                Tm = np.moveaxis(ev(m, world, a), (i, j, k), (0, 1, 2))
                X = 1j * np.einsum("abc,abk->ck", world.f, Tm)
                M += X if pm == phys[0] else X.T
    return M
W = World(11); totS = 0; totO = 0; totJ = 0
for name, r3, rj in ROWS:
    m3 = []
    for l in b3(*r3): m3 += list_monos_h(B3[l], True)
    S = 1j * two_from_three(m3, W); O = 1j * two_from_three_ref(m3, W)
    totS = totS + S; totO = totO + O
    if rj:
        mj = []
        for l in bj(*rj): mj += list_monos_h(BJ[l], False)
        A = two_list(mj, W); totJ = totJ + A
        print("%-10s |list - S| = %.1e   |list - O| = %.2g" % (name, np.abs(A - S).max(), np.abs(A - O).max()))
print("sum over all rows incl. N_a,N_b:  |S - O| = %.2e  (|S| = %.3g)" % (np.abs(totS - totO).max(), np.abs(totS).max()))
