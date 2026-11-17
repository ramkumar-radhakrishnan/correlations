# (2) two-rho contact terms from full symmetrisation of the three-rho monomial lists of Evolution_three_rho.tex
#     ("Genuine terms", c.c. already included there) versus the section-2 lists of Evolution_JIMWLK__terms.tex
import sys, numpy as np
sys.path.insert(0, ".")
from nc import *
T3 = "7583871c-Evolution_three_rho.tex"; J = "9aa89920-Evolution_JIMWLK__terms.tex"
B3 = dict(align_blocks(T3)); BJ = dict(align_blocks(J))
def b3(lo, hi): return [l for l in sorted(B3) if lo <= l <= hi]
def bj(lo, hi): return [l for l in sorted(BJ) if lo <= l <= hi]
# row -> (line range of the three-rho list in T3 , line range of section-2 list in J)
ROWS = [("GI  RII ", (192, 206), (289, 313)), ("GI  RIII", (207, 211), (314, 364)), ("GI  RIV ", (212, 219), (365, 416)),
        ("GI  RV  ", (220, 228), (417, 435)), ("GII RI  ", (229, 232), (436, 470)), ("GIII RI ", (233, 236), (478, 528)),
        ("GIII RII", (237, 240), (529, 579)), ("GIII RIII", (241, 244), (580, 598)), ("GIII RV ", (245, 248), (601, 627)),
        ("GIII RVI", (249, 252), (628, 654)), ("from 4rho (3-rho N_a,N_b)", (178, 190), None)]
if __name__ == "__main__":
    worlds = [World(s) for s in (11, 12)]
    for name, r3, rj in ROWS:
        m3 = []
        for l in b3(*r3): m3 += list_monos_h(B3[l], True)
        res = []
        for W in worlds:
            D = 1j * two_from_three(m3, W)        # overall i of the three-rho prefactor
            if rj:
                mj = []
                for l in bj(*rj): mj += list_monos_h(BJ[l], False)
                A = two_list(mj, W)
                res.append((np.abs(D.imag).max(), np.abs(D).max(), np.abs(A).max(), np.abs(D - A).max()))
            else:
                res.append((np.abs(D.imag).max(), np.abs(D).max(), 0, 0))
        print("%-10s 3rho blocks %s (%d ordered monomials)  J blocks %s : |derived|=%.3g |listed|=%.3g  max|diff|=%.2e / %.2e" % (
            name, b3(*r3), len(m3), bj(*rj) if rj else None, res[0][1], res[0][2], res[0][3], res[1][3]))
