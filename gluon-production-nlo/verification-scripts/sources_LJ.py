"""Independent source set: every log(k+/Lambda)-proportional term, transcribed from
the (previously verified) Group_I/II/III_corrected.tex results and 4rho_corrected.tex.

Conventions:
 - common factor (1/(2pi)^3) g^4 (1/k+) log(k+/Lambda) stripped; pi powers numeric.
 - positions renamed so the phase is always e^{-ik(wp-w)} ("standard frame").
 - [2/eps + 2gamma_E + log(r^2 mu^2/4)] = -(1/pi) * Sum_z K(a-z).K(b-z)  (dim-reg identity),
   applied with a fresh internal variable; log(r1^2/r2^2) likewise as a difference.
 - anticommutators expanded as AB + BA (exact).
 - U factors commute with everything; only rho-rho ordering matters.
 - 'hc_group' tags: Group I, Group III and 4rho groups 1,3 come 'with c.c.';
   Group II and 4rho group 2 are complete as written.
"""
import math
from engine import Term

PI = math.pi
K = lambda a, b, c, d: ('K', a, b, c, d)
I2 = lambda a, b: ('inv2', a, b)
U = lambda i, j, p: ('U', i, j, p)
F = lambda i, j, k: ('f', i, j, k)

SRC = []   # list of (term, needs_hc: bool)

def add(coef, kern, color, rhos, hc, tag, phase=("wp", "w")):
    SRC.append((Term(coef, kern, color, rhos, phase=phase, tag=tag), hc))

# ============================================================ GROUP I (with c.c.)
# ---- Row II, N1  [2rho]  (G-I eq (52): L_J part of the UV x dipole term)
# + (Nc/16 pi^4) [Ebar + l(x-w)] KK x [U(x')-U(w')]rho(x') [U(w)-U(x)]rho(x)
# -> -(Nc/16 pi^5) Sum_z K(x-z).K(w-z) K(x'-w').K(x-w) x brackets
for s1, p1 in ((+1, "xp"), (-1, "wp")):
    for s2, p2 in ((+1, "w"), (-1, "x")):
        add(-(1 / (16 * PI**5)) * s1 * s2,
            [K("xp", "wp", "x", "w"), K("x", "z", "w", "z")],
            [U("a", "bp", p1), U("a", "b", p2)],
            [("bp", "xp"), ("b", "x")], True, "GI_R2_N1", )
SRC[-4][0].Nc_pow = 1; SRC[-3][0].Nc_pow = 1; SRC[-2][0].Nc_pow = 1; SRC[-1][0].Nc_pow = 1

# ---- Row II, N2 2a  [3rho]  (G-I eq (64) 2a with FT logs -> z-kernel form)
# 2a = +(i/32 pi^4) L_J Int (x'-w')^i [ K_i(x-w) log((w-y)^2/(x-y)^2)
#                                      - K_i(y-w) log((w-x)^2/(x-y)^2) ] x Br1 x Br2
# log((w-y)^2/(x-y)^2) -> -(1/pi) Sum_z [K(w-z).K(y-z) - K(x-z).K(y-z)]
# log((w-x)^2/(x-y)^2) -> -(1/pi) Sum_z [K(w-z).K(x-z) - K(y-z).K(x-z)]
# Br1 = [U^{ab'}(x') - U^{ab'}(w')] rho^{b'}(x')
# Br2 = [f^{hbc} U^{ah}(w) {rho^b(x),rho^c(y)} - f^{abc} U^{bd}(x) U^{ce}(y) {rho^d(x),rho^e(y)}]
_kpieces = [
    (+1, [K("xp", "wp", "x", "w"), K("w", "z", "y", "z")]),
    (-1, [K("xp", "wp", "x", "w"), K("x", "z", "y", "z")]),
    (-1, [K("xp", "wp", "y", "w"), K("w", "z", "x", "z")]),
    (+1, [K("xp", "wp", "y", "w"), K("y", "z", "x", "z")]),
]
for sk, kern in _kpieces:
    for s1, p1 in ((+1, "xp"), (-1, "wp")):
        for s2, colr, r1, r2 in (
            (+1, [F("h", "b", "c"), U("a", "h", "w")], ("b", "x"), ("c", "y")),
            (-1, [F("a", "b", "c"), U("b", "d", "x"), U("c", "e", "y")], ("d", "x"), ("e", "y")),
        ):
            for rA, rB in ((r1, r2), (r2, r1)):   # anticommutator
                add(-(1j / (32 * PI**5)) * sk * s1 * s2,
                    kern, [U("a", "bp", p1)] + colr,
                    [("bp", "xp"), rA, rB], True, "GI_R2_2a")

# ---- Row III  [3rho]  (verified L_J part: +i/16pi^5 K2  - i/8pi^5 R2)
for s1, p1 in ((+1, "xp"), (-1, "wp")):
    for s2, pz in ((+1, "z"), (-1, "y")):
        base_color = [F("a", "b", "d"), U("a", "bp", p1), U("d", "e", "x"), U("b", "c", pz)]
        rhos = [("bp", "xp"), ("e", "x"), ("c", "y")]
        add(+(1j / (16 * PI**5)) * s1 * s2,
            [K("xp", "wp", "x", "w"), K("x", "z", "y", "z")],
            base_color, rhos, True, "GI_R3_K2")
        add(-(1j / (8 * PI**5)) * s1 * s2,
            [K("xp", "wp", "y", "z"), I2("z", "w")],
            base_color, rhos, True, "GI_R3_R2")

# ---- Row IV  [3rho]  (G-I eq (136-138) blue term; K12 kernel)
# +(i/8pi^5) L_J K12 x [U(x')-U(w')]rho(x')
#   x [f^{cdc'} U^{ac}(w) U^{bd}(z) U^{be}(y) rho^e(y) rho^{c'}(x)
#      - f^{abc} U^{be}(y) U^{cd}(x) rho^e(y) rho^d(x)]
_k12 = [(-1, [K("xp", "wp", "x", "w"), K("w", "z", "y", "z")]),
        (+0.5, [K("xp", "wp", "x", "w"), K("x", "z", "y", "z")])]
for sk, kern in _k12:
    for s1, p1 in ((+1, "xp"), (-1, "wp")):
        add(+(1j / (8 * PI**5)) * sk * s1,
            kern,
            [U("a", "bp", p1), F("c", "d", "cp"), U("a", "c", "w"), U("b", "d", "z"), U("b", "e", "y")],
            [("bp", "xp"), ("e", "y"), ("cp", "x")], True, "GI_R4a")
        add(-(1j / (8 * PI**5)) * sk * s1,
            kern,
            [U("a", "bp", p1), F("a", "b", "c"), U("b", "e", "y"), U("c", "d", "x")],
            [("bp", "xp"), ("e", "y"), ("d", "x")], True, "GI_R4b")

# ---- Row V part 1  [3rho]  (G-I eq (154)/(238) extraction terms)
# ext1 (phase (w',z->w)):  y'->z rename:
#   -(i/32pi^5) K(x'-w').K(y-w) K(w-x).K(z-x) x Br1 x
#     [f^{adc} U^{db}(x) U^{ce}(w) {rho^e(y),rho^b(z)} - f^{abc} U^{cd}(y) U^{be}(z) {rho^d(y),rho^e(z)}]
for s1, p1 in ((+1, "xp"), (-1, "wp")):
    for s2, colr, r1, r2 in (
        (+1, [F("a", "d", "c"), U("d", "b", "x"), U("c", "e", "w")], ("e", "y"), ("b", "z")),
        (-1, [F("a", "b", "c"), U("c", "d", "y"), U("b", "e", "z")], ("d", "y"), ("e", "z")),
    ):
        for rA, rB in ((r1, r2), (r2, r1)):
            add(-(1j / (32 * PI**5)) * s1 * s2,
                [K("xp", "wp", "y", "w"), K("w", "x", "z", "x")],
                [U("a", "bp", p1)] + colr,
                [("bp", "xp"), rA, rB], True, "GI_R5p1_e1")
# ext2 (phase (w',x->w)): renames x->w, y'->x:
#   kernel K(x'-w').K(x-w) K(z-w).K(y-z); colors with U(w),U(z),rho(y),rho(x)
for s1, p1 in ((+1, "xp"), (-1, "wp")):
    for s2, colr, r1, r2 in (
        (+1, [F("a", "d", "c"), U("d", "b", "w"), U("c", "e", "z")], ("e", "y"), ("b", "x")),
        (-1, [F("a", "b", "c"), U("c", "d", "y"), U("b", "e", "x")], ("d", "y"), ("e", "x")),
    ):
        for rA, rB in ((r1, r2), (r2, r1)):
            add(-(1j / (32 * PI**5)) * s1 * s2,
                [K("xp", "wp", "x", "w"), K("z", "w", "y", "z")],
                [U("a", "bp", p1)] + colr,
                [("bp", "xp"), rA, rB], True, "GI_R5p1_e2")

# ---- Row V part 2a  [2rho]  (G-I eq (182)/(242) extraction terms; real coefficient)
# u1 (phase (w', z->w); rename z->w):
#  -(1/16pi^5) (x'-w')^i (w-x)^m /((x'-w')^2 (w-x)^2 (y-w)^2)
#     [ (y-w)^i (x-w)^m/(x-w)^2 + (y-x)^m (y-w)^i /(2 (y-x)^2) ]
#  x [U(x')-U(w')] [ -f^{adc} f^{fbe} U^{db}(x) U^{ce}(w) + Nc U^{af}(x) ] rho^{b'}(x') rho^f(y)
# decompose kernel: piece1 = K(x'-w').K(y-w) * K(w-x).K(x-w)?? careful:
#  (x'-w')^i (y-w)^i /((x'-w')^2 (y-w)^2) * (w-x)^m (x-w)^m/((w-x)^2 (x-w)^2)
#  = K(x'-w').K(y-w) * ( -1/(x-w)^2 )     [since (w-x).(x-w) = -(x-w)^2]
#  piece2 = K(x'-w').K(y-w) * (w-x).(y-x) / ((w-x)^2 (y-x)^2) * (1/2)
#         = (1/2) K(x'-w').K(y-w) K(w-x).K(y-x)
_u1_k = [(-1, [K("xp", "wp", "y", "w"), I2("x", "w")]),
         (+0.5, [K("xp", "wp", "y", "w"), K("w", "x", "y", "x")])]
# u2 (phase (w', x->w); rename x->w, then kernel vars):
#  -(1/16pi^5) (x'-w')^i (z-w)^m /((x'-w')^2 (z-w)^2 (y-w)^2)
#     [ (y-w)^i (w-z)^m/(w-z)^2 - (y-w)^i (y-z)^m/(2 (y-z)^2) ]
#  piece1 = K(x'-w').K(y-w) * (z-w).(w-z)/((z-w)^2(w-z)^2) = -K(x'-w').K(y-w)/(z-w)^2
#  piece2 = -(1/2) K(x'-w').K(y-w) K(z-w).K(y-z)
_u2_k = [(-1, [K("xp", "wp", "y", "w"), I2("z", "w")]),
         (-0.5, [K("xp", "wp", "y", "w"), K("z", "w", "y", "z")])]
# colors: u1 after rename z->w:  [-ff U^{db}(x)U^{ce}(w) + Nc U^{af}(x)]
#         u2 after rename x->w:  [-ff U^{db}(w)U^{ce}(z) + Nc U^{af}(w)]
for kpieces, cpos, tag in ((_u1_k, ("x", "w", "x"), "GI_R5p2a_u1"),
                           (_u2_k, ("w", "z", "w"), "GI_R5p2a_u2")):
    pU1, pU2, pNc = cpos
    for sk, kern in kpieces:
        for s1, p1 in ((+1, "xp"), (-1, "wp")):
            add(-(1 / (16 * PI**5)) * sk * s1 * (-1),
                kern,
                [U("a", "bp", p1), F("a", "d", "c"), F("f", "b", "e"), U("d", "b", pU1), U("c", "e", pU2)],
                [("bp", "xp"), ("f", "y")], True, tag + "_ff")
            t = Term(-(1 / (16 * PI**5)) * sk * s1,
                     kern, [U("a", "bp", p1), U("a", "f", pNc)],
                     [("bp", "xp"), ("f", "y")], tag=tag + "_Nc")
            t.Nc_pow = 1
            SRC.append((t, True))

# ---- Row V 2bUV  [2rho]
# L_J part: +(Nc/8 pi^4)[Ebar+l(y-z)] e^{-ik(w'-z)} K(x'-w').K(y-z) [U(x')-U(w')][U(y)-U(z)] rr
# rename z->w (phase), y->x, dummy z for the Z-identity:
# -> -(Nc/8 pi^5) Sum_z K(x-z).K(w-z) K(x'-w').K(x-w) [U(x')-U(w')][U(x)-U(w)] rho(x')rho(x)
for s1, p1 in ((+1, "xp"), (-1, "wp")):
    for s2, p2 in ((+1, "x"), (-1, "w")):
        t = Term(-(1 / (8 * PI**5)) * s1 * s2,
                 [K("x", "z", "w", "z"), K("xp", "wp", "x", "w")],
                 [U("a", "bp", p1), U("a", "f", p2)],
                 [("bp", "xp"), ("f", "x")], tag="GI_R5_2bUV")
        t.Nc_pow = 1
        SRC.append((t, True))

# ---- Row V 2b finite  [2rho]  (G-I eq (233) L_J terms v1,v2,v3)
# v1 (phase (w',x)): rename x->w, then y stays, z stays, and old-w?? no w present.
#   kernel: [ (y-w).(x'-w') (w-z).(z-w) / ((w-z)^2 (x'-w')^2 (z-w)^2 (y-w)^2)
#             - (y-w).(x'-w') (y-z).(z-w) / (2 (y-z)^2 (x'-w')^2 (z-w)^2 (y-w)^2) ]
#   = -K(x'-w').K(y-w) /(z-w)^2  - (1/2) K(x'-w').K(y-w) K(y-z).K(z-w)
#   wait: first: (w-z).(z-w) = -(z-w)^2: num/(den) = -(y-w).(x'-w')/((x'-w')^2(y-w)^2 (w-z)^2)
#   second: -(1/2)(y-w).(x'-w')(y-z).(z-w)/((y-z)^2(z-w)^2 (x'-w')^2(y-w)^2)
#   colors: Nc[U(x')-U(w')][U^{af}(y)-U^{af}(w)] rho(x') rho(y)   (U(x)->U(w) after rename)
_v1_k = [(-1, [K("xp", "wp", "y", "w"), I2("z", "w")]),
         (-0.5, [K("xp", "wp", "y", "w"), K("y", "z", "z", "w")])]
for sk, kern in _v1_k:
    for s1, p1 in ((+1, "xp"), (-1, "wp")):
        for s2, p2 in ((+1, "y"), (-1, "w")):
            t = Term(-(1 / (16 * PI**5)) * sk * s1 * s2, kern,
                     [U("a", "bp", p1), U("a", "f", p2)],
                     [("bp", "xp"), ("f", "y")], tag="GI_R5_2bf_v1")
            t.Nc_pow = 1
            SRC.append((t, True))
# v2 (phase (w',z)): rename z->w, old w absent; x,y stay; old z->w:
#   kernel: [ (y-w).(x'-w') (x-w).(w-x)/((x-w)^2 (x'-w')^2 (w-x)^2 (y-w)^2)
#             + (y-x).(w-x)(y-w).(x'-w')/(2 (y-x)^2 (x'-w')^2 (w-x)^2 (y-w)^2) ]
#   = -K(x'-w').K(y-w)/(x-w)^2 + (1/2) K(x'-w').K(y-w) K(y-x).K(w-x)
#   colors: Nc[U(x')-U(w')][U(y)-U(x)] rho(x') rho(y)
_v2_k = [(-1, [K("xp", "wp", "y", "w"), I2("x", "w")]),
         (+0.5, [K("xp", "wp", "y", "w"), K("y", "x", "w", "x")])]
for sk, kern in _v2_k:
    for s1, p1 in ((+1, "xp"), (-1, "wp")):
        for s2, p2 in ((+1, "y"), (-1, "x")):
            t = Term(-(1 / (16 * PI**5)) * sk * s1 * s2, kern,
                     [U("a", "bp", p1), U("a", "f", p2)],
                     [("bp", "xp"), ("f", "y")], tag="GI_R5_2bf_v2")
            t.Nc_pow = 1
            SRC.append((t, True))
# v3 (phase (w',z)): rename z->w: coefficient -(2/16pi^5) = -(1/8pi^5)
#   kernel K(x'-w').K(y-w) * (y-x).(y-w)/((y-x)^2 (x-w)^2)   [non-K form -> 'fun']
def _v3fun(P):
    import numpy as np
    yx = P["y"] - P["x"]; yw = P["y"] - P["w"]; xw = P["x"] - P["w"]
    return float((yx @ yw) / ((yx @ yx) * (xw @ xw)))
for s1, p1 in ((+1, "xp"), (-1, "wp")):
    for s2, p2 in ((+1, "y"), (-1, "w")):
        t = Term(-(1 / (8 * PI**5)) * s1 * s2,
                 [K("xp", "wp", "y", "w"), ('fun', _v3fun, ["x", "y", "w"])],
                 [U("a", "bp", p1), U("a", "f", p2)],
                 [("bp", "xp"), ("f", "y")], tag="GI_R5_2bf_v3")
        t.Nc_pow = 1
        SRC.append((t, True))

# ============================================================ GROUP II (complete as written)
# ---- Row I part 1  [3rho]  (G-II eq (52) blue terms; K12 kernel both halves)
# half A: -(i/16pi^5) K12(x,w,z | x',w' spectators) x
#   [U^{ac'}(w'){rho^{c'}(x'),rho^b(y)} - U^{c'b}(z) U^{ad'}(x') U^{c'e'}(y) {rho^{d'}(x'),rho^{e'}(y)}]
#   x [f^{cbd} U^{ac}(w) rho^d(x) - f^{acd} U^{cb}(z) U^{de}(x) rho^e(x)]
for sk, kern in _k12:
    for s1, colr1, r1a, r1b in (
        (+1, [U("a", "cp", "wp")], ("cp", "xp"), ("b2", "y")),
        (-1, [U("cp", "b2", "z"), U("a", "dp", "xp"), U("cp", "ep", "y")], ("dp", "xp"), ("ep", "y")),
    ):
        for s2, colr2, r2 in (
            (+1, [F("c", "b2", "d"), U("a", "c", "w")], ("d", "x")),
            (-1, [F("a", "c", "d"), U("c", "b2", "z"), U("d", "e", "x")], ("e", "x")),
        ):
            for rA, rB in ((r1a, r1b), (r1b, r1a)):
                add(-(1j / (16 * PI**5)) * sk * s1 * s2,
                    kern, colr1 + colr2,
                    [rA, rB, r2], False, "GII_R1p1_A")
# half B: +(i/16pi^5) K12 with w<->w' swapped kernel (positions (x-w')... and (x'-w)):
#   kernel: (x'-w)^i (y-z)^j/((x'-w)^2(y-z)^2) [ (x-w')^i(z-w')^j/((z-w')^2(x-w')^2)
#            + (x-z)^j (x-w')^i/(2(x-z)^2 (x-w')^2) ]
#  = -K(x'-w).K(x-w') K(w'-z).K(y-z) + (1/2) K(x'-w).K(x-w') K(x-z).K(y-z)
_k12B = [(-1, [K("xp", "w", "x", "wp"), K("wp", "z", "y", "z")]),
         (+0.5, [K("xp", "w", "x", "wp"), K("x", "z", "y", "z")])]
for sk, kern in _k12B:
    for s1, colr1, r1 in (
        (+1, [F("cp", "b2", "dp"), U("a", "cp", "wp")], ("dp", "x")),
        (-1, [F("a", "cp", "dp"), U("cp", "b2", "z"), U("dp", "ep", "x")], ("ep", "x")),
    ):
        for s2, colr2, r2a, r2b in (
            (+1, [U("a", "c", "w")], ("c", "xp"), ("b2", "y")),
            (-1, [U("c", "b2", "z"), U("a", "d", "xp"), U("c", "e", "y")], ("d", "xp"), ("e", "y")),
        ):
            for rA, rB in ((r2a, r2b), (r2b, r2a)):
                add(+(1j / (16 * PI**5)) * sk * s1 * s2,
                    kern, colr1 + colr2,
                    [r1, rA, rB], False, "GII_R1p1_B")

# ---- Row I part 2, 2b1  [2rho]  (G-II eq (2b) first piece: L_J x Bil_b x bracket)
import numpy as np
def _2b1fun(P):
    xz = P["x"] - P["z"]; xw = P["x"] - P["w"]; zw = P["z"] - P["w"]; zwp = P["z"] - P["wp"]
    xpz = P["xp"] - P["z"]; xpwp = P["xp"] - P["wp"]
    A = xz @ xz; B = xw @ xw; C = zw @ zw; Cp = zwp @ zwp
    Ap = xpz @ xpz; Bp = xpwp @ xpwp
    Bil = ((xpwp @ xw) * (zw @ zwp) / (C * Cp)
           + (xpwp @ xw) * (xpz @ zw) / (2 * C * Ap)
           + (xw @ xpwp) * (xz @ zwp) / (2 * A * Cp)
           + (xz @ xpz) * (xpwp @ xw) / (4 * A * Ap))
    den = Bp * A * B - B**2 * Ap
    return float(Bil / den * (A - B * Ap / Bp))
_2b1_colors = [
    (+1, [F("cp", "b2", "dp"), F("c", "b2", "d"), U("a", "cp", "wp"), U("a", "c", "w")],
     [("dp", "xp"), ("d", "x")], 0),
    (-1, [F("cp", "b2", "dp"), F("a", "c", "d"), U("a", "cp", "wp"), U("c", "b2", "z"), U("d", "e", "x")],
     [("dp", "xp"), ("e", "x")], 0),
    (-1, [F("a", "cp", "dp"), F("c", "b2", "d"), U("a", "c", "w"), U("cp", "b2", "z"), U("dp", "ep", "xp")],
     [("ep", "xp"), ("d", "x")], 0),
    (+1, [U("d", "ep", "xp"), U("d", "e", "x")], [("ep", "xp"), ("e", "x")], 1),
]
for s, colr, rhos, ncp in _2b1_colors:
    t = Term(+(1 / (8 * PI**5)) * s, [('fun', _2b1fun, ["x", "xp", "z", "w", "wp"])],
             colr, rhos, tag="GII_R1p2_2b1")
    t.Nc_pow = ncp
    SRC.append((t, False))

# ---- Row IV  [2rho]
# +(Nc/8 pi^5) Sum_z K(w-z).K(w'-z) K(x-w).K(x'-w') [U(w')-U(x')][U(w)-U(x)] rho(x')rho(x)
for s1, p1 in ((+1, "wp"), (-1, "xp")):
    for s2, p2 in ((+1, "w"), (-1, "x")):
        t = Term(+(1 / (8 * PI**5)) * s1 * s2,
                 [K("w", "z", "wp", "z"), K("x", "w", "xp", "wp")],
                 [U("c", "ep", p1), U("c", "e", p2)],
                 [("ep", "xp"), ("e", "x")], tag="GII_R4")
        t.Nc_pow = 1
        SRC.append((t, False))

# ============================================================ GROUP III (with c.c.)
# ---- Row I  [3rho]:  -(i/8pi^5) K12 x
#  [U^{ab'}(x') U^{c'd'}(y) U^{c'b}(z) rho^{b'}(x') rho^{d'}(y)
#   - U^{d'b}(z) U^{d'e'}(y) U^{ac'}(w') rho^{c'}(x') rho^{e'}(y)]
#  x [f^{cbd} U^{ac}(w) rho^d(x) - f^{acd} U^{cb}(z) U^{de}(x) rho^e(x)]
for sk, kern in _k12:
    for s1, colr1, r1a, r1b in (
        (+1, [U("a", "bp", "xp"), U("cp", "dp", "y"), U("cp", "b2", "z")], ("bp", "xp"), ("dp", "y")),
        (-1, [U("dp", "b2", "z"), U("dp", "ep", "y"), U("a", "cp", "wp")], ("cp", "xp"), ("ep", "y")),
    ):
        for s2, colr2, r2 in (
            (+1, [F("c", "b2", "d"), U("a", "c", "w")], ("d", "x")),
            (-1, [F("a", "c", "d"), U("c", "b2", "z"), U("d", "e", "x")], ("e", "x")),
        ):
            add(-(1j / (8 * PI**5)) * sk * s1 * s2,
                kern, colr1 + colr2, [r1a, r1b, r2], True, "GIII_R1")

# ---- Row II  [3rho]:  -(i/8pi^5) K12 x
#  [U^{c'd'}(y) U^{ab'}(x') U^{c'b}(z) rho^{d'}(y) rho^{b'}(x') - U^{ab'}(x') rho^b(y) rho^{b'}(x')]
#  x [f^{cbd} U^{ac}(w) rho^d(x) - f^{acd} U^{cb}(z) U^{de}(x) rho^e(x)]
for sk, kern in _k12:
    for s1, colr1, rhos1 in (
        (+1, [U("cp", "dp", "y"), U("a", "bp", "xp"), U("cp", "b2", "z")], [("dp", "y"), ("bp", "xp")]),
        (-1, [U("a", "bp", "xp")], [("b2", "y"), ("bp", "xp")]),
    ):
        for s2, colr2, r2 in (
            (+1, [F("c", "b2", "d"), U("a", "c", "w")], ("d", "x")),
            (-1, [F("a", "c", "d"), U("c", "b2", "z"), U("d", "e", "x")], ("e", "x")),
        ):
            add(-(1j / (8 * PI**5)) * sk * s1 * s2,
                kern, colr1 + colr2, rhos1 + [r2], True, "GIII_R2")

# ---- Row III part 1  [3rho]  (G-III eq (75) blue extraction; rename y'->wp)
#  -(i f^{c'd'a}/16pi^5) K(wp-z).K(x-z) K(x'-wp).K(y-w) x
#  [U^{c'e'}(wp) U^{d'b}(z) rho^{e'}(x') - U^{d'b}(z) U^{c'e'}(x') rho^{e'}(x')]
#  x [U^{ac}(w) {rho^c(y),rho^b(x)} - U^{cb}(z) U^{ad}(y) U^{ce}(x) {rho^d(y),rho^e(x)}]
for s1, colr1 in ((+1, [U("cp", "ep", "wp"), U("dp", "b2", "z")]),
                  (-1, [U("dp", "b2", "z"), U("cp", "ep", "xp")])):
    for s2, colr2, r2a, r2b in (
        (+1, [U("a", "c", "w")], ("c", "y"), ("b2", "x")),
        (-1, [U("c", "b2", "z"), U("a", "d", "y"), U("c", "e", "x")], ("d", "y"), ("e", "x")),
    ):
        for rA, rB in ((r2a, r2b), (r2b, r2a)):
            add(-(1j / (16 * PI**5)) * s1 * s2,
                [K("wp", "z", "x", "z"), K("xp", "wp", "y", "w")],
                [F("cp", "dp", "a")] + colr1 + colr2,
                [("ep", "xp"), rA, rB], True, "GIII_R3p1")

# ---- Row III part 2  [2rho]  (G-III eq (89) blue extraction; rename y'->wp)
#  -(1/8pi^5) [ -K(x-w).K(x'-wp) K(z-w).K(z-wp)
#               + (1/2) K(x-z).K(wp-z) K(x-w).K(x'-wp) ] x
#  [f^{c'd'a} U^{c'e'}(wp) U^{d'b}(z) rho^{e'}(x') - f U^{d'b}(z) U^{c'e'}(x') rho^{e'}(x')]
#  x [f^{cbd} U^{ac}(w) rho^d(x) - f^{acd} U^{cb}(z) U^{de}(x) rho^e(x)]
_k_r3p2 = [(-1, [K("x", "w", "xp", "wp"), K("z", "w", "z", "wp")]),
           (+0.5, [K("x", "z", "wp", "z"), K("x", "w", "xp", "wp")])]
for sk, kern in _k_r3p2:
    for s1, colr1 in ((+1, [U("cp", "ep", "wp"), U("dp", "b2", "z")]),
                      (-1, [U("dp", "b2", "z"), U("cp", "ep", "xp")])):
        for s2, colr2, r2 in (
            (+1, [F("c", "b2", "d"), U("a", "c", "w")], ("d", "x")),
            (-1, [F("a", "c", "d"), U("c", "b2", "z"), U("d", "e", "x")], ("e", "x")),
        ):
            add(-(1 / (8 * PI**5)) * sk * s1 * s2,
                kern, [F("cp", "dp", "a")] + colr1 + colr2,
                [("ep", "xp"), r2], True, "GIII_R3p2")

# ---- Row V  [3rho]  (G-III eq (112) blue piece; rename y'->wp; f^{c'ca})
#  -(i f^{c'ca}/8pi^5) K(wp-z).K(x-z) K(x'-wp).K(y-w) x
#  [U^{c'e'}(wp) rho^{e'}(x') - U^{c'e'}(x') rho^{e'}(x')]
#  x [U^{ce}(x) U^{ad}(y) rho^e(x) rho^d(y) - U^{ad}(w) U^{ce}(x) rho^e(x) rho^d(y)]
for s1, p1 in ((+1, "wp"), (-1, "xp")):
    for s2, colr2 in ((+1, [U("c", "e", "x"), U("a", "d", "y")]),
                      (-1, [U("a", "d", "w"), U("c", "e", "x")])):
        add(-(1j / (8 * PI**5)) * s1 * s2,
            [K("wp", "z", "x", "z"), K("xp", "wp", "y", "w")],
            [F("cp", "c", "a"), U("cp", "ep", p1)] + colr2,
            [("ep", "xp"), ("e", "x"), ("d", "y")], True, "GIII_R5")

# ---- Row VI  [3rho]  (G-III eq (127) blue piece; rename y'->wp; f^{c'd'a})
#  -(i f^{c'd'a}/8pi^5) same kernel x
#  [U^{c'e'}(wp) U^{d'b}(z) rho^{e'}(x') - U^{d'b}(z) U^{c'e'}(x') rho^{e'}(x')]
#  x [U^{cb}(z) U^{ad}(y) U^{ce}(x) rho^d(y) rho^e(x) - U^{ad}(y) rho^d(y) rho^b(x)]
for s1, colr1 in ((+1, [U("cp", "ep", "wp"), U("dp", "b2", "z")]),
                  (-1, [U("dp", "b2", "z"), U("cp", "ep", "xp")])):
    for s2, colr2, rhos2 in (
        (+1, [U("c", "b2", "z"), U("a", "d", "y"), U("c", "e", "x")], [("d", "y"), ("e", "x")]),
        (-1, [U("a", "d", "y")], [("d", "y"), ("b2", "x")]),
    ):
        add(-(1j / (8 * PI**5)) * s1 * s2,
            [K("wp", "z", "x", "z"), K("xp", "wp", "y", "w")],
            [F("cp", "dp", "a")] + colr1 + colr2,
            [("ep", "xp")] + rhos2, True, "GIII_R6")
