# (2) "From 4 rho's": full symmetrisation of eq.(4rho) of 4rho.tex versus the lists
import sys, numpy as np, itertools
sys.path.insert(0, ".")
from nc import *
F4 = "2d955658-4rho.tex"; J = "9aa89920-Evolution_JIMWLK__terms.tex"; T3 = "7583871c-Evolution_three_rho.tex"
B4 = dict(align_blocks(F4)); BJ = dict(align_blocks(J)); B3 = dict(align_blocks(T3))
eq = B4[64].replace(r"^{\prime}", "'")
op = [m.start() for m in re.finditer(r"\\Bigg\\\{", eq)]; cl = [m.start() for m in re.finditer(r"\\Bigg\\\}", eq)]
assert len(op) == 4 and len(cl) == 4
P = [eq[op[0] + 7:cl[0]], eq[op[1] + 7:cl[1]], eq[op[3] + 7:cl[2]]]
assert eq[cl[0]:cl[0] + 20].replace(" ", "").startswith(r"\Bigg\}+\rm{c.c}") and eq[cl[2]:cl[2] + 20].replace(" ", "").startswith(r"\Bigg\}+\rm{c.c}")
def conv(t):
    t = t.replace(r"\bigg[", "[").replace(r"\bigg]", "]").replace(r"\bigg(", "<(").replace(r"\bigg)", ")>")
    return t.replace(r"\times", " ").replace(r"\notag", " ").replace("\\\\", " ").replace("&", " ")
KER = [("x'-w'", "u-w"), ("x-z", "y-z")]
parts = [monos_from(conv(t), dots=KER) for t in P]
M4 = parts[0] + cc(parts[0]) + parts[1] + parts[2] + cc(parts[2])
print("eq.(4rho): ordered four-rho monomials  A:%d (+cc)  B:%d  C:%d (+cc)" % tuple(len(p) for p in parts))

def canon_sym4(monos, W, phys):
    return canon_sym(monos, W, phys)
def three_from_four_sym(monos, W, phys):
    """symmetric three-rho coefficient (axes ordered as phys) of the contact terms 1/2 sum_{i<j} S([Xi,Xj] Xk Xl); no overall factor"""
    tot = np.zeros((W.n,) * 3, complex)
    for m in monos:
        pts = m.pts()
        for i, j in itertools.combinations(range(4), 2):
            k, l = [q for q in range(4) if q not in (i, j)]
            for perm in itertools.permutations(phys):
                a = {pts[i]: perm[0], pts[j]: perm[0], pts[k]: perm[1], pts[l]: perm[2]}
                T = np.moveaxis(ev(m, W, a), (i, j, k, l), (0, 1, 2, 3))
                X = 0.5 * 1j * np.einsum("abc,abkl->ckl", W.f, T)       # axes at (perm0, perm1, perm2)
                tot += np.transpose(X, [perm.index(p) for p in phys])
    return tot
if __name__ == "__main__":
    for seed in (31, 32):
        W = World(seed)
        s4 = canon_sym4(M4, W, ("P0", "P1", "P2", "P3"))
        sz = max(np.abs(canon_sym4([m], W, ("P0", "P1", "P2", "P3"))).max() for m in M4[:40])
        print("seed %d: symmetric four-rho coefficient of eq.(4rho): max|.| = %.2e (single term %.2g)" % (seed, np.abs(s4).max(), sz))
        byU = {}
        for m in M4: byU.setdefault(len(m.Us), []).append(m)
        print("    by number of Wilson lines:", {k: "%.1e" % np.abs(canon_sym4(v, W, ("P0", "P1", "P2", "P3"))).max() for k, v in byU.items()})
        # three-rho
        m3 = list_monos_h(B3[179], True) + list_monos_h(B3[185], True)   # N_a, N_b of Evolution_three_rho (prefactor i)
        s3_auth = 1j * canon_sym(m3, W, ("P0", "P1", "P2"))
        s3_der = three_from_four_sym(M4, W, ("P0", "P1", "P2"))
        print("    three-rho (symmetric coefficient): derived %.3g  listed(N_a+N_b) %.3g   max|diff| = %.2e" % (np.abs(s3_der).max(), np.abs(s3_auth).max(), np.abs(s3_der - s3_auth).max()))
        # two-rho
        full = two_from_four(M4, W)
        mj1 = list_monos_h(BJ[230], False); mj2 = list_monos_h(BJ[252], False)
        A1 = two_list(mj1, W); A2 = two_list(mj2, W)
        t3 = 1j * two_from_three(m3, W)
        print("    two-rho: |full sym of eq.(4rho)| = %.3g (imag %.1e) ; |J N_1| = %.3g ; |J N_2| = %.3g ; |sym of N_a+N_b| = %.3g" % (
            np.abs(full).max(), np.abs(full.imag).max(), np.abs(A1).max(), np.abs(A2).max(), np.abs(t3).max()))
        for nm, X in (("full - (N_1+N_2)", full - A1 - A2), ("full - N_1", full - A1), ("full - N_2", full - A2),
                      ("sym(N_a+N_b) - N_1", t3 - A1), ("sym(N_a+N_b) - N_2", t3 - A2), ("full - sym(N_a+N_b) - N_1", full - t3 - A1),
                      ("full - sym(N_a+N_b) - N_2", full - t3 - A2)):
            print("       %-28s max|.| = %.2e" % (nm, np.abs(X).max()))
