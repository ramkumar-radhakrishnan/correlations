#!/usr/bin/env python3
"""
ref3_dilute / v0_common: shared definitions for the dilute virtual (Group I) computation.
Point-like quark projectile rho^a(x) = t^a delta^2(x) (fundamental generators), dilute target
U = 1 + u, u^{ab}(l) = f^{aeb} alpha^e(l), g = 1, k^+ = 1 (all terms are homogeneous of degree -1/2 in k^+).

LO amplitude: A_LO^a_i = sqrt(2) f^{aeb} t^b V0_i alpha^e,  V0 = k/k^2 - q/q^2,  q = k + l.
A virtual term  sum_s T_s^{ae} W_{s,i} alpha^e  contributes to the ratio Group I / LO
    R = 2 Re sum_s proj_s (V0 . W_s) / (sqrt(2) |V0|^2),   proj_s = sum Tr(LOcol^dag T_s) / sum Tr(LOcol^dag LOcol)
(quark colour traced, target colour summed).  With g = 1, R is in units of g^2 = 4 pi alpha_s.
"""
import numpy as np, itertools

def gellmann():
    l = np.zeros((8, 3, 3), complex)
    l[0][0, 1] = l[0][1, 0] = 1
    l[1][0, 1] = -1j; l[1][1, 0] = 1j
    l[2][0, 0] = 1; l[2][1, 1] = -1
    l[3][0, 2] = l[3][2, 0] = 1
    l[4][0, 2] = -1j; l[4][2, 0] = 1j
    l[5][1, 2] = l[5][2, 1] = 1
    l[6][1, 2] = -1j; l[6][2, 1] = 1j
    l[7] = np.diag([1, 1, -2]) / np.sqrt(3)
    return l / 2

def su2():
    s = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)
    return s / 2

def group(N):
    t = gellmann() if N == 3 else su2()
    n = len(t)
    f = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            comm = t[a] @ t[b] - t[b] @ t[a]
            for c in range(n):
                f[a, b, c] = (-2j * np.trace(comm @ t[c])).real
    return f, t

class Colour:
    def __init__(self, N=3):
        self.N = N
        self.f, self.t = group(N)
        self.LO = np.einsum('aeb,bij->aeij', self.f, self.t)     # [a, e, quark i, quark j]
        self.norm = np.einsum('aeji,aeji->', self.LO.conj(), self.LO).real  # sum Tr(LO^dag LO)
    def proj(self, T):
        """T[a,e,i,j]: projection onto LO colour structure."""
        return np.einsum('aeji,aeji->', self.LO.conj(), T) / self.norm

def V0(k, l):
    q = k + l
    return k / (k @ k) - q / (q @ q)

def beta(kp, k, pp, p):
    """1-rho part of B2^{ba}_{ji}(k+,k;p+,p) = f^{abc} rho^c beta_{ji}; matrix [j,i], i <-> (kp,k), j <-> (pp,p). g=1."""
    k2, p2 = k @ k, p @ p
    s = p + k; s2 = s @ s
    t1 = 1j * np.outer(p, k) * np.sqrt(pp * kp) / (2 * (pp * k2 + kp * p2)) * (1 / (kp * p2) - 1 / (k2 * pp))
    w = p * kp - pp * k
    D2 = s2 * np.sqrt(pp * kp) * (pp + kp) * (p2 / pp + k2 / kp)
    br = (s @ w) * np.eye(2) / (pp + kp) - np.outer(s, w) / kp - np.outer(w, s) / pp
    t2 = -1j * br / D2
    t3 = -1j * (pp - kp) * np.eye(2) / (2 * np.sqrt(pp * kp) * (pp + kp) ** 2 * (k2 / kp + p2 / pp))
    return t1 + t2 + t3

def Amom(pp, p):
    """momentum-space emission coefficient without colour: A^b_j(p+,p) = -sqrt2 p_j t^b/(sqrt(p+) p^2) (point source)."""
    return -np.sqrt(2) * p / (np.sqrt(pp) * (p @ p))
