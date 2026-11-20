# Singlet forms: printed brackets (singlet_data.py, the source of the TeX of the averaged section)
# against direct contraction of the colour tensors; trace forms and fundamental forms of G, Q_d.
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fe_colour import *
from s3_colour import STRUCT, randUs, E
from singlet_data import CUBIC, QUAD
out = []
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)
def ev(S, U, br):
    tot = 0
    for c, kind, args, nc in br:
        if kind == 'D': v = Dtr(S, U, *args)
        elif kind == 'F': v = Ftr(S, U, *args)
        elif kind == 'Q': v = Qtr(S, U, *args)
        elif kind == 'G': v = (-1j*Trd(S, U, args[0], args[1], args[2], args[0])).real
        elif kind == 'Qd': v = (1j*Trd(S, U, *args)).real
        tot += c*(S.N if nc else 1)*v
    return tot
SF = {s[0]: s for s in STRUCT}
gen = {'ord': {'12': 0.5j, '23': 0.5j, '13': -0.5j}, 'r{rr}': {'12': 1j, '13': -1j}, '{rr}r': {'23': 1j, '13': -1j}}
for N in [3, 4]:
    S = SUN(N); rng = np.random.default_rng(100+N)
    # trace forms
    w = 0
    for _ in range(4):
        U = randUs(S, rng)
        w = max(w, abs(Ftr(S, U, 'x', 'y', 'z')-Fdef(S, U, 'x', 'y', 'z')))
        w = max(w, abs((-1j*Trd(S, U, 'x', 'y', 'z', 'x')).real-Gdef(S, U, 'x', 'y', 'z')), abs((-1j*Trd(S, U, 'x', 'y', 'z', 'x')).imag))
        Qd = 1j*Trd(S, U, 'x', 'y', 'z', 'w')
        w = max(w, abs(Qd.imag))
        Qdd = np.einsum('hij,ki,jl,hmn,km,ln->', S.f, U['x'], U['y'], S.d, U['z'], U['w'])  # f^{hij}d^{hmn}[U^T(x)... ] check of Q_d by components
    pr('SU(%d): F trace = f f UUU, G = -i Tr{T U U D U U}(a,b;c,a) = f d UUU real: max dev %.1e' % (N, w))
    # fundamental forms
    w = 0
    for _ in range(4):
        V = {p: S.randV(rng) for p in ['a', 'b', 'c', 'e']}
        U = {p: S.adj(V[p]) for p in V}
        tr = np.trace
        Gf = np.imag(tr(V['a'].conj().T@V['b'])*tr(V['b'].conj().T@V['c'])*tr(V['c'].conj().T@V['a']))+np.imag(tr(V['a'].conj().T@V['b']@V['c'].conj().T@V['a']@V['b'].conj().T@V['c']))
        Ga = (-1j*Trd(S, U, 'a', 'b', 'c', 'a')).real
        A = V['a'].conj().T@V['b']; B = V['c'].conj().T@V['e']
        Qf = -np.imag(tr(A)*tr(B)*tr(A.conj().T@B.conj().T))-np.imag(tr(A@B@A.conj().T@B.conj().T))
        Qa = (1j*Trd(S, U, 'a', 'b', 'c', 'e')).real
        w = max(w, abs(Gf-Ga), abs(Qf-Qa))
    pr('SU(%d): fundamental forms of G and Q_d vs adjoint traces: max dev %.1e' % (N, w))
    # cubic rows
    wc = w3 = wcoef = 0
    for row, R in CUBIC.items():
        for lab, contacts, c3tex, c3, br3, oargs in R['structs']:
            name, fn, typ, pts = SF[lab]
            assert typ == R['typ'] and pts == R['pts'], (row, lab)
            for _ in range(3):
                U = randUs(S, rng)
                for pair, fr, to, ctex, cc, nargs, br in contacts:
                    U2 = dict(U); U2[fr] = U[to]
                    Phi = E('abc,abc->', fn(S, U2), S.f)
                    expected = R['A']*gen[typ][pair]
                    wcoef = max(wcoef, abs(expected-cc))
                    wc = max(wc, abs(Phi-ev(S, U2, br))/max(1, abs(Phi)))
                Cd = E('abc,abc->', fn(S, U), S.d)
                w3 = max(w3, abs(R['t']*R['A']*Cd-c3*ev(S, U, br3))/max(1, abs(Cd)))
    pr('SU(%d): cubic rows, 13 structures: contact coefficients A x (generic factor) vs printed: max dev %.1e' % (N, wcoef))
    pr('SU(%d): contact brackets (D,F,Q traces) vs C^{abc}f^{abc} at the coincident points: max rel dev %.1e' % (N, wc))
    pr('SU(%d): (3rho) lines t*A*C^{abc}d^{abc} vs printed coefficient x (G,Q_d trace) bracket: max rel dev %.1e' % (N, w3))
    # quadratic
    wq = {}
    for name, (pq, br) in QUAD.items():
        dev = 0
        for _ in range(3):
            U = randUs(S, rng); L = U['xp']-U['wp']
            if name == 'I.II N1': C = E('ab,ac->bc', L, U['w']-U['x'])
            elif name == 'I.V N2a': C = E('ab,af->bf', L, S.N*U['x'])-E('aB,adc,fbe,db,ce->Bf', L, S.f, S.f, U['x'], U['z'])
            elif name == 'I.V Nc[U(y)-U(x)]': C = E('ab,af->bf', L, U['y']-U['x'])
            elif name == 'II.I second':
                f = S.f
                C = (E('kbB,lbd,ak,al->Bd', f, f, U['wp'], U['w'])
                     - E('kbB,acd,ak,cb,dE->BE', f, f, U['wp'], U['z'], U['x'])
                     - E('akd,cbq,ac,kb,dB->Bq', f, f, U['w'], U['z'], U['xp'])
                     + S.N*E('dB,dE->BE', U['xp'], U['x']))
            elif name.startswith('III.III second'):
                v = 'yp' if "y'" in name else 'z'; f = S.f
                Lv = E('db,cB->cdbB', U['z'], U[v]-U['xp'])
                R2 = E('cbd,ac->abd', f, U['w'])-E('acd,cb,de->abe', f, U['z'], U['x'])
                C = E('cda,cdbB,abh->Bh', f, Lv, R2)
            elif name == 'pair C colour C_C':
                C = E('aB,ah->Bh', L, S.N*U['x'])-E('aB,adc,hbe,db,ce->Bh', L, S.f, S.f, U['z'], U['w'])
            dev = max(dev, abs(np.trace(C)-ev(S, U, br))/max(1, abs(np.trace(C))))
        wq[name] = dev
    pr('SU(%d): quadratic brackets vs C^{aa}: ' % N+', '.join('%s %.1e' % (k, v) for k, v in wq.items()))
open('s5_singlet.log', 'w').write('\n'.join(out)+'\n')
# check of the printed (3rho) form of Sec. E: (-i c3/A) [ sum_G c Tr_d(a,b;c,a) - sum_Qd c' Tr_d(a,b;c,e) ] = t A C.d / A = t C.d
for N in [3, 4]:
    S = SUN(N); rng = np.random.default_rng(200+N); w = 0
    for row, R in CUBIC.items():
        for lab, contacts, c3tex, c3, br3, oargs in R['structs']:
            fn = SF[lab][1]
            for _ in range(3):
                U = randUs(S, rng)
                tr = 0
                for cf, kind, args, nc in br3:
                    if kind == 'G': tr += cf*Trd(S, U, args[0], args[1], args[2], args[0])
                    else: tr += -cf*Trd(S, U, *args)
                printed = (-1j*c3/R['A'])*tr
                direct = R['t']*E('abc,abc->', fn(S, U), S.d)
                w = max(w, abs(printed-direct)/max(1, abs(direct)))
    pr('SU(%d): Sec. E (3rho) lines as printed [(-i c3/A) x Tr{T U U D U U} combination] vs t C^{abc}d^{abc}: max rel dev %.1e' % (N, w))
open('s5_singlet.log', 'w').write('\n'.join(out)+'\n')
