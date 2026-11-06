"""JIMWLK: (a) flat adjoint N1+N2 (lines 657-707) vs colour-singlet trace form
(738-762, adjoint traces);  (b) singlet vs user's pasted JIMWLK dN^(1..5);
(c) BFKL singlet vs user's pasted LO-BFKL."""
import itertools, math, cmath
import numpy as np
from texparse import parse_block
from blocks import block_body, FJ, FB, PI
from singlet_bfkl import Cfg, Kdot, Ksc, eval_adjoint

def Tad(cfg):
    """adjoint generators (T^e)_{ab} = -i f^{eab}"""
    return [-1j * cfg.f[e] for e in range(cfg.d)]

def jimwlk_trace(cfg):
    N = cfg.N
    TA = Tad(cfg)
    tot = 0
    vs = ["x", "y", "z", "w", "wp"]
    for perm in itertools.permutations(range(cfg.nsites), 5):
        siteof = dict(zip(vs, perm))
        P = {v: cfg.pos[s] for v, s in siteof.items()}
        UA = {v: cfg.Uadj[siteof[v]] for v in vs}
        ph = cmath.exp(-1j * float(np.dot(cfg.kvec, P["wp"] - P["w"])))
        mu = cfg.mu2[siteof["x"], siteof["y"]]
        x, y, z, w, wp = (P[v] for v in vs)
        Ud = lambda v: UA[v].T
        Uf = lambda v: UA[v]
        tr = np.trace
        def T2(m1, m2):
            return sum(np.trace(TA[e] @ m1 @ TA[e] @ m2) for e in range(cfg.d))
        def TT(m1, m2):
            return sum(np.trace(TA[e] @ m1) * np.trace(TA[e] @ m2) for e in range(cfg.d))
        pre = Kdot(x - w, y - wp)
        val = 0
        val += pre * Ksc(w, wp, z) * (4 * N * tr(Ud("w") @ Uf("wp"))
               - 2 * T2(Ud("w") @ Uf("wp"), Ud("wp") @ Uf("z"))
               - 2 * T2(Ud("w") @ Uf("wp"), Ud("z") @ Uf("w")))
        val += pre * Ksc(x, wp, z) * (-4 * N * tr(Ud("x") @ Uf("wp"))
               + 2 * T2(Ud("wp") @ Uf("x"), Ud("x") @ Uf("z"))
               + 2 * T2(Ud("wp") @ Uf("x"), Ud("z") @ Uf("wp")))
        val += pre * Ksc(x, y, z) * (4 * N * tr(Ud("x") @ Uf("y"))
               - 2 * T2(Ud("x") @ Uf("y"), Ud("y") @ Uf("z"))
               - 2 * T2(Ud("x") @ Uf("y"), Ud("z") @ Uf("x")))
        val += pre * Ksc(y, w, z) * (-4 * N * tr(Ud("w") @ Uf("y"))
               + 2 * T2(Ud("w") @ Uf("y"), Ud("y") @ Uf("z"))
               + 2 * T2(Ud("w") @ Uf("y"), Ud("z") @ Uf("w")))
        val += pre * Ksc(wp, wp, z) * (-2 * N * tr(Ud("w") @ Uf("wp")) + 2 * N * tr(Ud("x") @ Uf("wp"))
               - T2(Ud("wp") @ Uf("x"), Ud("wp") @ Uf("z")) - T2(Ud("wp") @ Uf("x"), Ud("z") @ Uf("wp"))
               + T2(Ud("w") @ Uf("wp"), Ud("wp") @ Uf("z")) + T2(Ud("w") @ Uf("wp"), Ud("z") @ Uf("wp"))
               + TT(Ud("wp") @ Uf("x"), Ud("wp") @ Uf("z")) + TT(Ud("w") @ Uf("wp"), Ud("wp") @ Uf("z")))
        val += pre * Ksc(w, w, z) * (-2 * N * tr(Ud("w") @ Uf("wp")) + 2 * N * tr(Ud("w") @ Uf("y"))
               + T2(Ud("w") @ Uf("wp"), Ud("w") @ Uf("z")) + T2(Ud("w") @ Uf("wp"), Ud("z") @ Uf("w"))
               - T2(Ud("w") @ Uf("y"), Ud("w") @ Uf("z")) - T2(Ud("w") @ Uf("y"), Ud("z") @ Uf("w"))
               - TT(Ud("w") @ Uf("wp"), Ud("w") @ Uf("z")) + TT(Ud("w") @ Uf("y"), Ud("w") @ Uf("z")))
        val += pre * Ksc(x, x, z) * (2 * N * tr(Ud("x") @ Uf("wp")) - 2 * N * tr(Ud("x") @ Uf("y"))
               - T2(Ud("wp") @ Uf("x"), Ud("x") @ Uf("z")) - T2(Ud("wp") @ Uf("x"), Ud("z") @ Uf("x"))
               + T2(Ud("x") @ Uf("y"), Ud("x") @ Uf("z")) + T2(Ud("x") @ Uf("y"), Ud("z") @ Uf("x"))
               - TT(Ud("wp") @ Uf("x"), Ud("x") @ Uf("z")) - TT(Ud("x") @ Uf("y"), Ud("x") @ Uf("z")))
        val += pre * Ksc(y, y, z) * (2 * N * tr(Ud("w") @ Uf("y")) - 2 * N * tr(Ud("x") @ Uf("y"))
               - T2(Ud("w") @ Uf("y"), Ud("y") @ Uf("z")) - T2(Ud("w") @ Uf("y"), Ud("z") @ Uf("y"))
               + T2(Ud("x") @ Uf("y"), Ud("y") @ Uf("z")) + T2(Ud("x") @ Uf("y"), Ud("z") @ Uf("y"))
               - TT(Ud("w") @ Uf("y"), Ud("y") @ Uf("z")) + TT(Ud("x") @ Uf("y"), Ud("y") @ Uf("z")))
        tot += ph * mu / (N * N - 1) * val
    return tot / (16 * PI**5)

def jimwlk_target(cfg, s5=(1, -1, -1, 1)):
    """user's pasted dN^(1..5); t==z; overall C set to 1 (shape comparison);
    s5: sign pattern of the four dN^(5) terms (paste ambiguous)."""
    N = cfg.N
    TA = Tad(cfg)
    tot = 0
    vs = ["x", "y", "z", "w", "wp"]
    for perm in itertools.permutations(range(cfg.nsites), 5):
        siteof = dict(zip(vs, perm))
        P = {v: cfg.pos[s] for v, s in siteof.items()}
        UA = {v: cfg.Uadj[siteof[v]] for v in vs}
        ph = cmath.exp(-1j * float(np.dot(cfg.kvec, P["wp"] - P["w"])))
        mu = cfg.mu2[siteof["x"], siteof["y"]]
        x, y, z, w, wp = (P[v] for v in vs)
        Ud = lambda v: UA[v].T
        Uf = lambda v: UA[v]
        def T2(m1, m2):
            return sum(np.trace(TA[e] @ m1 @ TA[e] @ m2) for e in range(cfg.d))
        def TT(m1, m2):
            return sum(np.trace(TA[e] @ m1) * np.trace(TA[e] @ m2) for e in range(cfg.d))
        pre = Kdot(x - w, y - wp)
        D1 = Ud("x") - Ud("z"); D2 = Uf("y") - Uf("z")
        v1 = (Ksc(y, x, z) * T2(Ud("x") @ Uf("y"), D1 @ D2)
              - Ksc(wp, x, z) * T2(Ud("x") @ Uf("wp"), D1 @ D2)
              - Ksc(y, w, z) * T2(Ud("w") @ Uf("y"), D1 @ D2)
              + Ksc(wp, w, z) * T2(Ud("w") @ Uf("wp"), D1 @ D2))
        E1 = Ud("y") - Ud("z"); E2 = Uf("x") - Uf("z")
        v2 = (Ksc(x, y, z) * T2(Ud("x") @ Uf("y"), E1 @ E2)
              - Ksc(x, wp, z) * T2(Ud("x") @ Uf("wp"), E1 @ E2)
              - Ksc(w, y, z) * T2(Ud("w") @ Uf("y"), E1 @ E2)
              + Ksc(w, wp, z) * T2(Ud("w") @ Uf("wp"), E1 @ E2))
        v3 = -(Ksc(y, y, z) * T2((Ud("x") - Ud("w")) @ Uf("y"), D1 @ D2)
               - Ksc(wp, wp, z) * T2((Ud("x") - Ud("w")) @ Uf("wp"), D1 @ D2))
        v4 = -(Ksc(x, x, z) * T2(Ud("x") @ (Uf("y") - Uf("wp")), E1 @ E2)
               - Ksc(w, w, z) * T2(Ud("w") @ (Uf("y") - Uf("wp")), E1 @ E2))
        v5 = -(s5[0] * Ksc(y, y, z) * TT((Ud("y") - Ud("z")) @ Uf("y"), (Ud("x") - Ud("w")) @ Uf("y"))
               + s5[1] * Ksc(wp, wp, z) * TT((Ud("wp") - Ud("z")) @ Uf("wp"), (Ud("x") - Ud("w")) @ Uf("wp"))
               + s5[2] * Ksc(x, x, z) * TT((Ud("x") - Ud("z")) @ Uf("x"), Ud("x") @ (Uf("y") - Uf("wp")))
               + s5[3] * Ksc(w, w, z) * TT((Ud("w") - Ud("z")) @ Uf("w"), Ud("w") @ (Uf("y") - Uf("wp"))))
        tot += ph * mu / (N * N - 1) * pre * (v1 + v2 + v3 + v4 + v5)
    return tot

def bfkl_target(cfg):
    N = cfg.N
    tot = 0
    vs = ["x", "y", "z", "w", "wp"]
    for perm in itertools.permutations(range(cfg.nsites), 5):
        siteof = dict(zip(vs, perm))
        P = {v: cfg.pos[s] for v, s in siteof.items()}
        UA = {v: cfg.Uadj[siteof[v]] for v in vs}
        ph = cmath.exp(-1j * float(np.dot(cfg.kvec, P["wp"] - P["w"])))
        mu = cfg.mu2[siteof["x"], siteof["y"]]
        x, y, z, w, wp = (P[v] for v in vs)
        Ud = lambda v: UA[v].T
        Uf = lambda v: UA[v]
        tr = np.trace
        KD = 2 * Ksc(x, y, z) - Ksc(x, x, z) - Ksc(y, y, z)
        v = (Kdot(x - w, y - wp) * KD * tr((Ud("x") - Ud("w")) @ (Uf("y") - Uf("wp")))
             - 2 * Kdot(z - w, z - wp) * Ksc(x, y, z) * tr((Ud("z") - Ud("w")) @ (Uf("z") - Uf("wp")))
             + Kdot(z - w, y - wp) * (2 * Ksc(x, y, z) - Ksc(x, x, z)) * tr((Ud("z") - Ud("w")) @ (Uf("y") - Uf("wp")))
             + Kdot(x - w, z - wp) * (2 * Ksc(x, y, z) - Ksc(y, y, z)) * tr((Ud("x") - Ud("w")) @ (Uf("z") - Uf("wp"))))
        tot += ph * mu / (N * N - 1) * N * v
    return -tot   # overall minus of the target formula

if __name__ == "__main__":
    # (a) JIMWLK flat N1+N2 vs trace form
    b1 = block_body(FJ, 656, 686)
    b2 = block_body(FJ, 687, 707)
    terms = parse_block(b1, extra_coef=1 / (16 * PI**5), tag="J1") + \
            parse_block(b2, extra_coef=1 / (16 * PI**5), tag="J2")
    print("JIMWLK flat terms:", len(terms))
    for N, seed in ((2, 21), (3, 22)):
        cfg = Cfg(N, 5, seed)
        A = eval_adjoint(terms, cfg)
        S = jimwlk_trace(cfg)
        print(f"(a) N={N}: adjoint={A:.6e}  trace={S:.6e}  reldiff={abs(A - S)/max(abs(A), abs(S)):.3e}")
    # (b) singlet vs pasted target: ratio should be coordinate- and N-independent
    import itertools as it
    for s5 in [(1,-1,-1,1), (1,1,1,1), (1,-1,1,-1), (-1,1,1,-1)]:
        rats = []
        for N, seed in ((2, 21), (2, 31), (3, 22)):
            cfg = Cfg(N, 5, seed)
            S = jimwlk_trace(cfg)
            T = jimwlk_target(cfg, s5)
            rats.append((S / T).real if abs(T) > 0 else float('nan'))
        print(f"(b) s5={s5}: ratios J_user/J_target = {[f'{r:.6f}' for r in rats]}")
    # (c) BFKL singlet vs pasted target
    from singlet_bfkl import eval_adjoint as ea
    bodyB = block_body(FB, 570, 605, marker="e^{-ik\\cdot(w'-w)}")
    termsB = parse_block(bodyB, extra_coef=1 / (16 * PI**5), tag="Bflat")
    rats = []
    for N, seed in ((2, 21), (2, 31), (3, 22)):
        cfg = Cfg(N, 5, seed)
        A = ea(termsB, cfg)
        T = bfkl_target(cfg)
        rats.append((A / T).real if abs(T) > 0 else float('nan'))
    print(f"(c) ratios BFKL_user/BFKL_target = {[f'{r:.6f}' for r in rats]}")
