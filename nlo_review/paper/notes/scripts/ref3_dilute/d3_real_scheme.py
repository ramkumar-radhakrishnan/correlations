#!/usr/bin/env python3
"""
ref3_dilute / d3: rapidity scheme of the real emission.
Paper (plus-momentum cut-offs):   R_paper = int_0^1 dx/(x(1-x)) [A_gl/(qa^2 qb^2) - 1/(k^2 p^2)]   (+ ln(vee/Lambda)/(k^2 p^2))
BSV (hep-ph/0608154, Eq. eq:defbsquare, rapidity y2 of the second gluon, cut s_Lambda):
   R_BSV(s_L) = int dy2 [A_gl theta(s_L - s12)]/(qa^2 qb^2) - ln(s_L^2/(k^2 p^2)) /(k^2 p^2)
   (BSV's 1/2 in front of the MRK term belongs to their pair normalisation)
with y1 = 0 for the observed gluon, x = |k|/(|k|+|p| e^{y2}), s12 = Lambda^2/(x(1-x)) (massless pair, transverse
Lambda = (1-x)k - x p).  Since dx/(x(1-x)) = -dy2 exactly, R_paper = lim_{s_L->inf} R_BSV.  Checked numerically.
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import d1_qmrk_integrand as d1
from scipy.integrate import quad

def f(x, qa, k, p):
    qb = k + p - qa
    return d1.A_gluons(qa, qb, k, p, x) / ((qa @ qa) * (qb @ qb)) - 1 / ((k @ k) * (p @ p))

def R_paper(qa, k, p):
    g = lambda u: f(1 / (1 + np.exp(u)), qa, k, p)          # u = ln((1-x)/x), dx/(x(1-x)) = du
    a, e1 = quad(g, -30, 0, limit=400); b, e2 = quad(g, 0, 30, limit=400)
    return a + b

def R_bsv(qa, k, p, sL):
    nk, np_ = np.linalg.norm(k), np.linalg.norm(p)
    def h(y2):
        x = nk / (nk + np_ * np.exp(y2))
        La = (1 - x) * k - x * p
        s12 = La @ La / (x * (1 - x))
        qb = k + p - qa
        return d1.A_gluons(qa, qb, k, p, x) / ((qa @ qa) * (qb @ qb)) if s12 < sL else 0.0
    # region s12 < sL: bounded in y2
    ymax = np.log(sL / (nk * np_)) + 5
    val, err = quad(h, -ymax, ymax, limit=800, points=[0.0])
    # single-inclusive normalisation: A_gl/(qa^2 qb^2) -> 1/(k^2 p^2) at both ends, window width ln(sL^2/(k^2 p^2));
    # BSV's factor 1/2 refers to their pair normalisation (each gluon can be the jet: factor 2 in config (c)).
    return val - np.log(sL ** 2 / ((k @ k) * (p @ p))) / ((k @ k) * (p @ p))

if __name__ == '__main__':
    rng = np.random.default_rng(3)
    for t in range(3):
        qa, k, p = rng.normal(size=2), rng.normal(size=2), rng.normal(size=2)
        rp = R_paper(qa, k, p)
        line = f"point {t}: R_paper = {rp:.8f}   R_BSV(s_L) ="
        for sL in (1e2, 1e4, 1e6):
            line += f"  {R_bsv(qa, k, p, sL):.8f} (s_L={sL:.0e})"
        print(line)
