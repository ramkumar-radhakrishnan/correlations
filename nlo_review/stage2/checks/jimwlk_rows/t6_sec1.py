# (1) section-1 lists of Evolution_JIMWLK__terms.tex versus the entries of Evolution_two_rho.tex
# entries are hand-transcribed in the variables of the lists (relabelling maps given in comments), units g^4/(16 pi^5)/(2pi)^3/k^+ log.
# UV bracket [2/eps + 2 gamma_E + log(r^2 mu^2/4)] -> (1/pi) int_z of the kernel given in the comment.
import sys, numpy as np
sys.path.insert(0, ".")
from nc import *
J = "9aa89920-Evolution_JIMWLK__terms.tex"; BJ = dict(align_blocks(J))
def d(p, a, b): return p(a) - p(b)
def Kd(u, v): return (u @ v) / ((u @ u) * (v @ v))
def inv(u): return 1.0 / (u @ u)
KK = lambda p: Kd(d(p, "x'", "w'"), d(p, "x", "w"))
FBU = r"[U^{ab'}(x')-U^{ab'}(w')]"
E = {}
# GI RII:  -(N_c/16pi^4) K.K [-2/eps-2g-log((x-w)^2..)] [U(x')rho-U(w')rho][U(w)rho(x)-U(x)rho(x)] ; bracket -> -(1/pi) int_z [1/(w-z)^2 - K(x-z).K(w-z)]
E[41] = (True, [(1.0, lambda p: KK(p) * (inv(d(p, "w", "z")) - Kd(d(p, "x", "z"), d(p, "w", "z"))),
                r"N_{c}[U^{ab'}(x')\rho^{b'}(x')-U^{ab'}(w')\rho^{b'}(x')][U^{ab}(w)\rho^{b}(x)-U^{ab}(x)\rho^{b}(x)]")])
# GI RV 2a: 1st term relabelled (z->w, y->x, x->z); 2nd term relabelled (x->w, y->x)
k2a = lambda p: KK(p) * (-inv(d(p, "w", "z")) + 0.5 * Kd(d(p, "w", "z"), d(p, "x", "z")))
E[66] = (True, [(-1.0, k2a, FBU + r"[-f^{adc}f^{hbe}U^{db}(z)U^{ce}(w)+N_{c}U^{ah}(z)]\rho^{b'}(x')\rho^{h}(x)"),
                (-1.0, k2a, FBU + r"[-f^{adc}f^{hbe}U^{db}(w)U^{ce}(z)+N_{c}U^{ah}(w)]\rho^{b'}(x')\rho^{h}(x)")])
# GI RV 2b UV: (N_c/8pi^4)[2/eps+..log((y-z)^2)] K(x'-w').K(y-z) [U(x')-U(w')][U(y)-U(z)] ; relabel (z->w, y->x) ; bracket -> (1/pi) int_z [1/(w-z)^2-K(x-z).K(w-z)]
kuv = lambda p: KK(p) * (inv(d(p, "w", "z")) - Kd(d(p, "x", "z"), d(p, "w", "z")))
E[93] = (True, [(2.0, kuv, r"N_{c}" + FBU + r"[U^{ah}(x)-U^{ah}(w)]\rho^{b'}(x')\rho^{h}(x)")])
# GI RV 2b finite, terms 1 (x->w, y->x) and 2 (z->w, y->x, x->z)
E[111] = (True, [(-1.0, k2a, r"N_{c}" + FBU + r"[U^{ah}(x)-U^{ah}(w)]\rho^{b'}(x')\rho^{h}(x)"),
                 (-1.0, k2a, r"N_{c}" + FBU + r"[U^{ah}(x)-U^{ah}(z)]\rho^{b'}(x')\rho^{h}(x)")])
# GI RV 2b finite, term 3 (z->w, y->x, x->z): kernel K(x'-w').K(x-w) (x-z).(x-w)/((x-z)^2 (w-z)^2)
k3 = lambda p: KK(p) * (d(p, "x", "z") @ d(p, "x", "w")) / ((d(p, "x", "z") @ d(p, "x", "z")) * (d(p, "w", "z") @ d(p, "w", "z")))
E[137] = (True, [(-2.0, k3, r"N_{c}" + FBU + r"[U^{ah}(x)-U^{ah}(w)]\rho^{b'}(x')\rho^{h}(x)")])
# GII RI (pointwise, no relabelling); full kernel as typed
def kg2(p):
    xw = d(p, "x", "w"); xpwp = d(p, "x'", "w'"); xz = d(p, "x", "z"); xpz = d(p, "x'", "z"); zw = d(p, "z", "w"); zwp = d(p, "z", "w'")
    den = (xpwp @ xpwp) * (xz @ xz) * (xw @ xw) - (xw @ xw) ** 2 * (xpz @ xpz)
    br = (xpwp @ xw) * (zw @ zwp) / ((zw @ zw) * (zwp @ zwp)) + (xpwp @ xw) * (xpz @ zw) / (2 * (zw @ zw) * (xpz @ xpz)) \
        + (xw @ xpwp) * (xz @ zwp) / (2 * (xz @ xz) * (zwp @ zwp)) + (xz @ xpz) * (xpwp @ xw) / (4 * (xz @ xz) * (xpz @ xpz))
    return br * ((xz @ xz) - (xw @ xw) * (xpz @ xpz) / (xpwp @ xpwp)) / den
E[157] = (False, [(2.0, kg2, r"[f^{c'bd'}f^{cbd}U^{ac'}(w')U^{ac}(w)\rho^{d'}(x')\rho^{d}(x)-f^{c'bd'}f^{acd}U^{ac'}(w')U^{cb}(z)U^{de}(x)\rho^{d'}(x')\rho^{e}(x)"
                             r"-f^{ac'd'}f^{cbd}U^{ac}(w)U^{c'b}(z)U^{d'e'}(x')\rho^{e'}(x')\rho^{d}(x)+N_{c}U^{de'}(x')U^{de}(x)\rho^{e'}(x')\rho^{e}(x)]")])
# GII RIV: (N_c/8pi^4) K.K [-d/eps - d gamma - d/2 log((y'-y)^2..) - 1] [U(y')-U(x')][U(y)-U(x)], d=2-eps -> -[2/eps+2g+log] ; (y->w, y'->w')
#          bracket -> (1/pi) int_z (w-w')^2/(2 (w-z)^2 (w'-z)^2)
ksym = lambda p: KK(p) * 0.5 * (d(p, "w", "w'") @ d(p, "w", "w'")) * inv(d(p, "w", "z")) * inv(d(p, "w'", "z"))
E[181] = (False, [(-2.0, ksym, r"N_{c}[U^{ce'}(w')-U^{ce'}(x')][U^{ce}(w)-U^{ce}(x)]\rho^{e'}(x')\rho^{e}(x)")])
# GIII RIII second part (y'->w')
k33 = lambda p: KK(p) * (Kd(d(p, "z", "w"), d(p, "w'", "z")) + 0.5 * Kd(d(p, "x", "z"), d(p, "w'", "z")))
E[205] = (True, [(-2.0, k33, r"[f^{c'd'a}U^{c'e'}(w')U^{d'b}(z)\rho^{e'}(x')-f^{c'd'a}U^{d'b}(z)U^{c'e'}(x')\rho^{e'}(x')]"
                             r"[f^{cbd}U^{ac}(w)\rho^{d}(x)-f^{acd}U^{cb}(z)U^{de}(x)\rho^{e}(x)]")])
NAMES = {41: "GI RII", 66: "GI RV 2a", 93: "GI RV 2b UV", 111: "GI RV 2b finite (terms 1+2)", 137: "GI RV 2b finite (term 3)",
         157: "GII RI", 181: "GII RIV", 205: "GIII RIII 2nd part"}
def entry(l):
    wcc, ents = E[l]; ms = []
    for c, kf, s in ents: ms += monos_from(s, coef=c, kfun=kf)
    return ms + (cc(ms) if wcc else [])
if __name__ == "__main__":
    Ws = [World(61), World(62)]
    for l in sorted(E):
        lst = list_monos_h(BJ[l], False)
        out = []
        for W in Ws:
            A = two_list(lst, W); Bm = two_list(entry(l), W); out.append((np.abs(A).max(), np.abs(A - Bm).max()))
        print("J line %3d %-28s (%2d list terms): |list| = %.3g ; max|list - entry| = %.2e / %.2e" % (l, NAMES[l], len(lst), out[0][0], out[0][1], out[1][1]))
    # GII RI kernel simplification
    W = Ws[0]; p = lambda nm: W.X[{"x": "P0", "x'": "P1"}.get(nm, nm)]
    simp = KK(p) * (Kd(d(p, "w", "z"), d(p, "w'", "z")) - 0.5 * Kd(d(p, "x'", "z"), d(p, "w", "z")) - 0.5 * Kd(d(p, "x", "z"), d(p, "w'", "z")) + 0.25 * Kd(d(p, "x", "z"), d(p, "x'", "z")))
    print("GII RI kernel as typed vs K.K[K(w-z).K(w'-z) - K(x'-z).K(w-z)/2 - K(x-z).K(w'-z)/2 + K(x-z).K(x'-z)/4]: diff %.1e" % abs(kg2(p) - simp))
