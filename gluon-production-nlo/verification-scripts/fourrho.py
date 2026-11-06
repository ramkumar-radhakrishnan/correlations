"""4-rho L_J sources, transcribed from 4rho_corrected.tex eq. (lines 65-68).
Kernel (common): K(x'-w').K(u-w) K(x-z).K(y-z); prefactor 1/(16 pi^5) per group;
blue log = L_J + log(vee/k+): L_J coefficient = 1.
Groups A and C carry '+c.c.'; group B does not.
"""
import math
from engine import Term

PI = math.pi
K = lambda a, b, c, d: ('K', a, b, c, d)
U = lambda i, j, p: ('U', i, j, p)
F = lambda i, j, k: ('f', i, j, k)
KERN = [K("xp", "wp", "u", "w"), K("x", "z", "y", "z")]

def dedup(color, rhos):
    """rename repeated rho indices, adding deltas; also indices repeated between
    two rhos only (U-indices repeat legitimately via contraction)."""
    seen = set()
    for c in color:
        pass
    out_r, extra = [], []
    used = {}
    for n, (ix, p) in enumerate(rhos):
        if ix in used:
            nix = ix + f"_{n}"
            extra.append(('delta', ix, nix))
            out_r.append((nix, p))
        else:
            used[ix] = True
            out_r.append((ix, p))
    return color + extra, out_r

def prod(coef, br1_list, br2_list, tag):
    """br*_list: list of (sign, colorfactors, rholist); product keeps br1 rhos left."""
    out = []
    for s1, c1, r1 in br1_list:
        for s2, c2, r2 in br2_list:
            color, rhos = dedup(c1 + c2, r1 + r2)
            out.append(Term(coef * s1 * s2, KERN, color, rhos, tag=tag))
    return out

def anti(i1, p1, i2, p2):
    """{rho^i1(p1), rho^i2(p2)} -> two orderings"""
    return [[(i1, p1), (i2, p2)], [(i2, p2), (i1, p1)]]

def terms_4rho():
    T = []
    c = 1 / (16 * PI**5)
    # Br1 (common to group A): [U^{ab'}(x') rho^{b'}(x') - rho^{b'}(x') U^{ab'}(w')]
    BR1 = [(+1, [U("a", "bp", "xp")], [("bp", "xp")]),
           (-1, [U("a", "bp", "wp")], [("bp", "xp")])]
    # ---- A1
    A1b = [(+1, [U("a", "b", "u")], [("b", "u"), ("c", "x"), ("c", "y")]),
           (-1, [U("a", "b", "w"), U("c", "d", "x"), U("c", "e", "y")], [("d", "x"), ("e", "y"), ("b", "u")])]
    T += prod(-c, BR1, A1b, "4r_A1")
    # ---- A2
    A2b = [(+1, [U("a", "b", "w")], [("c", "x"), ("c", "y"), ("b", "u")]),
           (-1, [U("c", "d", "x"), U("c", "e", "y"), U("a", "b", "u")], [("d", "x"), ("e", "y"), ("b", "u")])]
    T += prod(+c, BR1, A2b, "4r_A2")
    # ---- A3
    A3b = []
    for rr in anti("d", "x", "c", "u"):
        A3b.append((+1, [U("a", "c", "w"), U("b", "d", "z"), U("b", "e", "y")], [("e", "y")] + rr))
    A3b.append((-1, [U("b", "c", "y"), U("b", "d", "x"), U("a", "e", "u")], [("c", "y"), ("d", "x"), ("e", "u")]))
    A3b.append((-1, [U("b", "c", "y"), U("a", "e", "u"), U("b", "d", "x")], [("c", "y"), ("e", "u"), ("d", "x")]))
    T += prod(-c, BR1, A3b, "4r_A3")
    # ---- A4
    A4b = []
    for pe in ("y", "z"):
        sgn = +1 if pe == "y" else -1
        A4b.append((sgn, [U("b", "c", "x"), U("a", "d", "u"), U("b", "e", pe)], [("c", "x"), ("d", "u"), ("e", "y")]))
        A4b.append((sgn, [U("a", "d", "u"), U("b", "c", "x"), U("b", "e", pe)], [("d", "u"), ("c", "x"), ("e", "y")]))
    T += prod(-c, BR1, A4b, "4r_A4")

    # ---- B1 (coefficient 1/2)
    B1a = []
    for rr in anti("b", "x", "cp", "xp"):
        B1a.append((+1, [U("cp", "a", "wp")], rr))
    B1a.append((-1, [U("cp", "dp", "x"), U("a", "d", "xp"), U("b", "cp", "z")], [("dp", "x"), ("d", "xp")]))
    B1a.append((-1, [U("a", "d", "xp"), U("cp", "dp", "x"), U("b", "cp", "z")], [("d", "xp"), ("dp", "x")]))
    B1b = []
    for rr in anti("c", "u", "b", "y"):
        B1b.append((+1, [U("a", "c", "w")], rr))
    B1b.append((-1, [U("b", "c", "z"), U("c", "e", "y"), U("a", "bp", "u")], [("e", "y"), ("bp", "u")]))
    B1b.append((-1, [U("b", "c", "z"), U("a", "bp", "u"), U("c", "e", "y")], [("bp", "u"), ("e", "y")]))
    T += prod(+c * 0.5, B1a, B1b, "4r_B1")
    # ---- B2 (coefficient 2)
    B2a = [(+1, [U("a", "bp", "xp"), U("cp", "ap", "x"), U("cp", "b", "z")], [("bp", "xp"), ("ap", "x")]),
           (-1, [U("b", "dp", "z"), U("dp", "ap", "x"), U("cp", "a", "wp")], [("cp", "xp"), ("ap", "x")])]
    B2b = [(+1, [U("b", "c", "z"), U("c", "d", "y"), U("a", "e", "u")], [("d", "y"), ("e", "u")]),
           (-1, [U("a", "c", "w"), U("d", "e", "y"), U("d", "b", "z")], [("e", "y"), ("c", "u")])]
    T += prod(+c * 2, B2a, B2b, "4r_B2")
    # ---- B3 (coefficient 2)
    B3a = [(+1, [U("cp", "dp", "x"), U("a", "ep", "xp"), U("b", "cp", "z")], [("dp", "x"), ("ep", "xp")]),
           (-1, [U("a", "bp", "xp")], [("b", "x"), ("bp", "xp")])]
    B3b = [(+1, [U("b", "c", "z"), U("a", "d", "u"), U("c", "e", "y")], [("d", "u"), ("e", "y")]),
           (-1, [U("a", "c", "u")], [("c", "u"), ("b", "y")])]
    T += prod(+c * 2, B3a, B3b, "4r_B3")

    # ---- C1
    C1a = [(+1, [U("a", "bp", "xp"), U("cp", "ap", "x"), U("cp", "b", "z")], [("bp", "xp"), ("ap", "x")]),
           (-1, [U("dp", "b", "z"), U("dp", "ap", "x"), U("a", "cp", "wp")], [("cp", "xp"), ("ap", "x")])]
    C1b = []
    for rr in anti("c", "u", "b", "y"):
        C1b.append((+1, [U("a", "c", "w")], rr))
    C1b.append((-1, [U("c", "b", "z"), U("c", "e", "y"), U("a", "ep", "u")], [("e", "y"), ("ep", "u")]))
    C1b.append((-1, [U("c", "b", "z"), U("a", "ep", "u"), U("c", "e", "y")], [("ep", "u"), ("e", "y")]))
    T += prod(+c, C1a, C1b, "4r_C1")
    # ---- C2
    C2a = [(+1, [U("cp", "dp", "x"), U("a", "ep", "xp"), U("cp", "b", "z")], [("dp", "x"), ("ep", "xp")]),
           (-1, [U("a", "bp", "xp")], [("b", "x"), ("bp", "xp")])]
    C2b = []
    for rr in anti("c", "u", "b", "y"):
        C2b.append((+1, [U("a", "c", "w")], rr))
    C2b.append((-1, [U("c", "b", "z"), U("c", "e", "y"), U("a", "ap", "u")], [("e", "y"), ("ap", "u")]))
    C2b.append((-1, [U("c", "b", "z"), U("a", "ap", "u"), U("c", "e", "y")], [("ap", "u"), ("e", "y")]))
    T += prod(+c, C2a, C2b, "4r_C2")
    # ---- C3 (coefficient 2)
    C3a = C2a
    C3b = [(+1, [U("c", "b", "z"), U("c", "d", "y"), U("a", "e", "u")], [("d", "y"), ("e", "u")]),
           (-1, [U("a", "c", "w"), U("d", "e", "y"), U("d", "b", "z")], [("e", "y"), ("c", "u")])]
    T += prod(+c * 2, C3a, C3b, "4r_C3")

    # hc flags: A and C groups carry +c.c., B does not
    flagged = []
    for t in T:
        hc = t.tag.startswith("4r_A") or t.tag.startswith("4r_C")
        flagged.append((t, hc))
    return flagged
