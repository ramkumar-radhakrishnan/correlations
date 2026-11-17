import sys; sys.path.insert(0,".")
from t3_4rho import *
from orient import fix_orientation
allflips = {}
for nm, p in zip("ABC", parts):
    for m in p:
        for fl in fix_orientation(m):
            allflips.setdefault((nm,) + fl, 0); allflips[(nm,) + fl] += 1
print("U's transposed relative to the literal text of eq.(4rho) (part, indices as typed, point): count of expanded monomials")
for k, v in sorted(allflips.items()): print("   ", k, v)
M4 = parts[0] + cc(parts[0]) + parts[1] + parts[2] + cc(parts[2])
PH4 = ("P0", "P1", "P2", "P3")
for seed in (31, 32):
    W = World(seed)
    print("seed %d: symmetric four-rho coefficient after orientation fix: max|.| = %.2e" % (seed, np.abs(canon_sym(M4, W, PH4)).max()))
    m3 = list_monos_h(B3[179], True) + list_monos_h(B3[185], True)
    s3_auth = 1j * canon_sym(m3, W, ("P0", "P1", "P2")); s3_der = three_from_four_sym(M4, W, ("P0", "P1", "P2"))
    print("    three-rho symmetric coefficient: derived %.3g  listed(N_a+N_b) %.3g   max|diff| = %.2e" % (np.abs(s3_der).max(), np.abs(s3_auth).max(), np.abs(s3_der - s3_auth).max()))
    full = two_from_four(M4, W)
    A1 = two_list(list_monos_h(BJ[230], False), W); A2 = two_list(list_monos_h(BJ[252], False), W)
    t3 = 1j * two_from_three(m3, W)
    print("    two-rho: |full|=%.3g |N_1|=%.3g |N_2|=%.3g ; full-(N_1+N_2): %.2e ; full-N_2: %.2e ; full-N_1: %.2e ; sym(Na+Nb)-N_1: %.2e" % (
        np.abs(full).max(), np.abs(A1).max(), np.abs(A2).max(), np.abs(full - A1 - A2).max(), np.abs(full - A2).max(), np.abs(full - A1).max(), np.abs(t3 - A1).max()))
