#!/usr/bin/env python3
"""dilute2 / e7: the soft region p -> 0, s = p+/k+ -> 0 at fixed u = p^2/(s k^2) = p^-/k^- of the virtual rows
(I.III, I.IV, resolved I.V; point quark, dilute), i.e. loop gluons with tiny transverse momentum whose light-cone energy
is comparable with that of the observed gluon.  At small s each row behaves as F_row(s,p) -> H_row(u,phi)/(s p^2).
A hard cut s > Lambda cuts this region at p^2 > Lambda u k^2; dimensional regularisation does not.  If
   hbar(u) = <H_tot(u,phi)>_phi - H_inf theta(u>1) - H_0 theta(u<1),   M_0 = PV int du/u hbar,   M_1 = int du/u ln(u) hbar,
then (for H_0 = H_inf, M_0 = 0)  [uncut, DR] - [hard cut] = M_1/(4 pi)  as a contribution to Group I/LO.
usage: e7_region_a.py <point name> [rho]"""
import sys, os, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v3c_fast as v3, v4_rowIV as v4, v5_rowV as v5
from v0_common import V0
from e0_points import PTS

def H_rows(u, phi, k, l, rho):
    k2 = k @ k; n = np.array([np.cos(phi), np.sin(phi)])
    r = min(rho, rho * u ** 0.75) if u < 1 else rho          # keep p << s|k| excluded: need rho << u (region b) and s << 1
    p = r * n; s = r * r / (u * k2)
    h3 = s * r * r * v3.Fvec(np.array([s]), p, k, l)[0]
    h4 = s * r * r * v4.Rdens(s, p, k, l)
    # Row V: soft daughter xi -> 0 with kappa = p (two soft ends: factor 2; the end xi -> 1 is obtained by kappa -> -kappa)
    v0 = V0(k, l); q = k + l
    def v5dens(xi, kap):
        Ws = v5.W_cases(xi, kap, k, l); val = v5.R_of(Ws, v0, Ws.keys())
        Mz2 = xi * (1 - xi) * (q @ q); Md2 = xi * (1 - xi) * k2; rr = kap @ kap
        wq = v5.W_UV_unavg(xi, kap, q) * rr / (rr + Mz2); wk = v5.W_UV_unavg(xi, kap, k) * rr / (rr + Md2)
        sub = 2 * (((v5.P['z1'] + v5.P['x1']) * (v0 @ wq)) + (v5.P['d1'] * (v0 @ wk))).real / (np.sqrt(2) * (v0 @ v0))
        return val - sub
    h5 = s * r * r * (v5dens(s, p) + v5dens(1 - s, -p))
    return h3, h4, h5

def profile(name, rho=1e-4, nphi=24, tmax=16.0, nt=160):
    P = PTS[name]; k = np.array([1.0, 0.0]); l = P['l']
    ts = (np.arange(nt) + 0.5) * (2 * tmax / nt) - tmax          # symmetric in t = ln u, t = 0 excluded (PV at u = 1)
    phis = 2 * np.pi * (np.arange(nphi) + 0.5) / nphi
    Hb = np.array([[np.mean([H_rows(np.exp(t), ph, k, l, rho)[i] for ph in phis]) for i in range(3)] for t in ts])
    return ts, Hb

if __name__ == '__main__':
    name = sys.argv[1]; rho = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4
    ts, Hb = profile(name, rho)
    dt = ts[1] - ts[0]; tot = Hb.sum(1)
    print(f"# {name}: rho={rho}; columns: ln u, <H_III>, <H_IV>, <H_V>, total")
    for t, h in list(zip(ts, Hb))[::8]:
        print(f"{t:8.3f} {h[0]: .6e} {h[1]: .6e} {h[2]: .6e} {h.sum(): .6e}")
    H0 = Hb[:6].mean(0); Hi = Hb[-6:].mean(0)
    print(f"limits u->0: III {H0[0]:.6f} IV {H0[1]:.6f} V {H0[2]:.6f} total {H0.sum():.6f};  u->inf: III {Hi[0]:.6f} IV {Hi[1]:.6f} V {Hi[2]:.6f} total {Hi.sum():.6f}")
    for lab, col in (('III', Hb[:, 0]), ('IV', Hb[:, 1]), ('V', Hb[:, 2]), ('total', tot)):
        h0, hi = col[:6].mean(), col[-6:].mean()
        hb = col - np.where(ts > 0, hi, h0)
        M0 = np.sum(hb) * dt; M1 = np.sum(hb * ts) * dt
        print(f"  {lab:5s}: H_0 = {h0: .6f}  H_inf = {hi: .6f}  M_0 = {M0: .6f}  M_1 = {M1: .6f}  M_1/(4pi) = {M1 / (4 * np.pi): .6f}")
