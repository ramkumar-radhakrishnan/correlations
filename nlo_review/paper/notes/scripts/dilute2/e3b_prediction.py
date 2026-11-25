#!/usr/bin/env python3
"""dilute2 / e3b: prediction for DD[C] if the CGC result agrees with the one-loop Lipatov vertex: DD[C] = DD[(abar/4) X] - DD[conv].
Written before the new values of C were available."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = [sys.argv[0], '.', '.']
import e3_assemble as e3
from e0_points import PTS
T = {}
for n, P in PTS.items():
    X = e3.abar / 4 * e3.X_bsv(P['qa'], P['qb']); conv = e3.abar / 2 * np.log(P['kn'] ** 2) * np.log((P['qa'] @ P['qa']) / (P['qb'] @ P['qb']))
    T[n] = (X, conv); print(f"{n}: (abar/4)X = {X:.6f}  conv = {conv:.6f}")
for a, a_ in (('a1', 'a2'), ('a1', 'a3'), ('a2', 'a3')):
    dd = lambda i: T[a + 'b1'][i] + T[a_ + 'b2'][i] - T[a + 'b2'][i] - T[a_ + 'b1'][i]
    print(f"({a},{a_})x(b1,b2): DD[X] = {dd(0):.5f}  DD[conv] = {dd(1):.5f}  predicted DD[C] = {dd(0) - dd(1):.5f}   (without conversion term: {dd(0):.5f})")
