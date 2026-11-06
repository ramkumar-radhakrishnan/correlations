"""Decisive adjudication: apply the JIMWLK Hamiltonian to the LO cross-section and
compare with the user's JIMWLK final (flat N1+N2, already parsed & verified == trace form).

N_LO (unit prefactor) = sum_{s1,P1; s2,P2} s1 s2 K(x-w).K(x'-w')
      Sum_c U^{c e'}(P1) U^{c e}(P2) rho^{e'}(x') rho^{e}(x),
      P1 in {w' (+), x' (-)},  P2 in {w (+), x (-)}.

H = lam * Int_z sum_{p,q in U-positions} K(p,q,z) [ J_L J_L + J_R J_R - 2 J_L^a U^{ab}(z) J_R^b ]
with J_L^a U(P) = T^a_adj U(P), J_R^a U(P) = U(P) T^a_adj, (T^a)_{bc} = -i f^{abc}.

We generate H.N_LO as engine Terms (f-vertex insertions) and compare the 2rho sector
tensor against the parsed JIMWLK final using the engine, fitting one overall constant.
Convention variants scanned: J_L<->J_R swap, sign of the cross term.
"""
import math
import numpy as np
from engine import Term, Config, accumulate
from texparse import parse_block
from blocks import block_body, FJ, PI

K = lambda a, b, c, d: ('K', a, b, c, d)
I2 = lambda a, b: ('inv2', a, b)
U = lambda i, j, p: ('U', i, j, p)
F = lambda i, j, k: ('f', i, j, k)

def hn_lo_terms(swapLR=False, crosssign=-2.0):
    """Return list of engine Terms for H.N_LO (lam=1)."""
    out = []
    for s1, P1 in ((+1, "wp"), (-1, "xp")):
        for s2, P2 in ((+1, "w"), (-1, "x")):
            s = s1 * s2
            base_kern = [K("x", "w", "xp", "wp")]
            rhos = [("ep", "xp"), ("e", "x")]
            # the two U factors: U^{c e'}(P1), U^{c e}(P2) (first index contracted)
            # ---------- diagonal p=q=P (both factors) ----------
            for P, other, oidx in ((P1, P2, "e"), (P2, P1, "ep")):
                # J_L J_L -> T^aT^a U = Nc U ; J_R J_R -> U T^aT^a = Nc U : total 2Nc
                t = Term(s * 2.0, base_kern + [I2(P, "z")],
                         [U("c", "ep", P1), U("c", "e", P2)], rhos, tag=f"diagNc_{P}")
                t.Nc_pow = 1
                out.append(t)
                # cross on same factor: crosssign * U^{ab}(z) (T^a U(P) T^b)
                # (T^a U T^b)^{c m} = - f^{a c d} f^{b g m} U^{d g}(P)
                if P == P1:
                    cf = [F("a", "c", "d"), F("b", "g", "ep"), U("d", "g", P1), U("c", "e", P2), U("a", "b", "z")]
                else:
                    cf = [F("a", "c", "d"), F("b", "g", "e"), U("d", "g", P2), U("c", "ep", P1), U("a", "b", "z")]
                out.append(Term(s * crosssign * (-1.0), base_kern + [I2(P, "z")], cf, rhos,
                                tag=f"diagX_{P}"))
    # ---------- off-diagonal pairs (P1, P2): factor 2 from (p,q)+(q,p) in JJ terms ----------
    for s1, P1 in ((+1, "wp"), (-1, "xp")):
        for s2, P2 in ((+1, "w"), (-1, "x")):
            s = s1 * s2
            base_kern = [K("x", "w", "xp", "wp"), K(P1, "z", P2, "z")]
            rhos = [("ep", "xp"), ("e", "x")]
            # LL: sum_c (T^aU1)^{c e'} (T^aU2)^{c e} = -Nc [U1^T U2]^{e'e}
            t = Term(s * 2.0 * (-1.0), base_kern,
                     [U("c", "ep", P1), U("c", "e", P2)], rhos, tag="offLL")
            t.Nc_pow = 1
            out.append(t)
            # RR: sum_c (U1T^a)^{ce'} (U2T^a)^{ce} = [T^{aT} U1^T U2 T^a] = -[T^a M T^a]
            #  [T^a M T^a]^{e'e} = -f^{a e' d} f^{a g e} M^{d g}
            #  M^{dg} = sum_c U^{cd}(P1) U^{cg}(P2)
            out.append(Term(s * 2.0 * (+1.0), base_kern,
                            [F("a", "ep", "d"), F("a", "g", "e"), U("c", "d", P1), U("c", "g", P2)],
                            rhos, tag="offRR"))
            # cross: crosssign * U^{ab}(z) [ J_L^a(P1) J_R^b(P2) + J_L^a(P2) J_R^b(P1) ]
            # J_L^a(P1)J_R^b(P2): sum_c (T^aU1)^{ce'} (U2T^b)^{ce} = [U1^T (-T^a) U2 T^b]^{e'e}
            #   = -[T^a... ] in components: -(T^a)_{?}... compute:
            #   [U1^T T^a U2 T^b]^{e' e} = U^{c e'}(P1) (T^a)_{c d} U^{d g}(P2) (T^b)_{g e}
            #   = U^{ce'}(P1) (-i f^{acd}) U^{dg}(P2) (-i f^{bge}) = - f^{acd} f^{bge} U^{ce'}U^{dg}
            #   so J_L^a(P1)J_R^b(P2) = -1 * (that) = + f^{acd} f^{bge} U^{ce'}(P1) U^{dg}(P2)
            if not swapLR:
                out.append(Term(s * crosssign * (+1.0), base_kern,
                                [U("a", "b", "z"), F("a", "c", "d"), F("b", "g", "e"),
                                 U("c", "ep", P1), U("d", "g", P2)], rhos, tag="offX1"))
                # J_L^a(P2) J_R^b(P1): sum_c (U1 T^b)^{c e'} (T^a U2)^{c e}
                #   = [T^{bT} U1^T T^a U2]^{e' e} = -[T^b U1^T T^a U2]
                #   [T^b U1^T T^a U2]^{e'e} = (T^b)_{e'd} U^{cd}(P1)... careful:
                #   (T^b)_{e' m} (U1^T)^{m c'} (T^a)_{c' d} U^{d e}(P2)
                #   = (-i f^{b e' m}) U^{m ... }  -> implement with explicit f's:
                out.append(Term(s * crosssign * (+1.0), base_kern,
                                [U("a", "b", "z"), F("b", "ep2", "m"), F("a", "c2", "d2"),
                                 U("c2", "m", P1), U("d2", "e", P2), ('delta', 'ep2', 'ep')],
                                rhos, tag="offX2"))
            else:
                # swapped convention: U^{ab}(z) J_R^a J_L^b
                out.append(Term(s * crosssign * (+1.0), base_kern,
                                [U("b", "a", "z"), F("a", "c", "d"), F("b", "g", "e"),
                                 U("c", "ep", P1), U("d", "g", P2)], rhos, tag="offX1s"))
                out.append(Term(s * crosssign * (+1.0), base_kern,
                                [U("b", "a", "z"), F("b", "ep2", "m"), F("a", "c2", "d2"),
                                 U("c2", "m", P1), U("d2", "e", P2), ('delta', 'ep2', 'ep')],
                                rhos, tag="offX2s"))
    return out

def fit_ratio(tA, tB):
    """fit tA ~ lam * tB over the sector-2 arrays; return lam, relative residual"""
    A = tA[2].ravel()
    B = tB[2].ravel()
    lam = np.vdot(B, A) / np.vdot(B, B)
    res = np.linalg.norm(A - lam * B) / max(np.linalg.norm(A), 1e-30)
    return lam, res

if __name__ == "__main__":
    b1 = block_body(FJ, 656, 686)
    b2 = block_body(FJ, 687, 707)
    user = parse_block(b1, extra_coef=1.0, tag="J1") + parse_block(b2, extra_coef=1.0, tag="J2")
    print("user JIMWLK final terms:", len(user))
    for group in ("su2", "su3"):
        cfg = Config(5, group, seed=41)
        tU = accumulate(user, cfg, sector_filter={2})
        for swap in (False, True):
            mine = hn_lo_terms(swapLR=swap)
            tM = accumulate(mine, cfg, sector_filter={2})
            lam, res = fit_ratio(tU, tM)
            print(f"{group} swapLR={swap}: lam = {lam:.6f}  rel-residual = {res:.3e}")
