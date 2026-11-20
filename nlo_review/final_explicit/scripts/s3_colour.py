# Colour checks for Final_terms_explicit.tex:
# (a) exact cancellations / colour identities (tensor level, random SU(3), SU(4) Wilson lines)
# (b) operator-level symmetrisation of every cubic colour structure (SU(3) charges on 3 sites)
# (c) colour-singlet factors as printed (D,F,Q,G,Q_d traces) vs direct contraction; trace forms
import numpy as np, itertools, sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fe_colour import *
out = []
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s); sys.stdout.flush()
PTS = ['x', 'xp', 'y', 'yp', 'z', 'w', 'wp']
def randUs(S, rng): return {p: S.adj(S.randV(rng)) for p in PTS}
E = lambda *a: np.einsum(*a, optimize=True)
# ---------------- cubic colour tensors C^{abc} (charges a,b,c in the order of the type) --------------
def C_2b(S, U):   # I.II N_2b, rho(x'){rho(x),rho(y)}
    f = S.f; L = U['xp']-U['wp']
    return E('ga,hbc,gh->abc', L, f, U['w'])-E('ga,gde,db,ec->abc', L, f, U['x'], U['y'])
def C_I3(S, U, u):  # I.III, ordered x',x,y
    L = U['xp']-U['wp']
    return E('ga,ghd,db,hc->abc', L, S.f, U['x'], U[u]-U['y'])
def C_I5(S, U):   # I.V first part, rho(x'){rho(y),rho(y')}
    L = U['xp']-U['wp']; f = S.f
    return E('ga,gde,dc,eb->abc', L, f, U['x'], U['z'])-E('ga,gde,eb,dc->abc', L, f, U['y'], U['yp'])
def C_1a(S, U):   # II.I N_1a, {rho(x'),rho(y)} rho(x)
    f = S.f
    t1 = E('ga,kbc,gk->abc', U['wp'], f, U['w'])-E('ga,gkl,kb,lc->abc', U['wp'], f, U['z'], U['x'])
    t2 = E('ga,mh,mb,khc,gk->abc', U['xp'], U['z'], U['y'], f, U['w'])-E('ga,mh,mb,gkl,kh,lc->abc', U['xp'], U['z'], U['y'], f, U['z'], U['x'])
    return t1-t2
def C_1b(S, U):   # II.I N_1b, rho(x){rho(x'),rho(y)}
    f = S.f
    t1 = E('gb,kca,gk->abc', U['w'], f, U['wp'])-E('gb,gkl,kc,la->abc', U['w'], f, U['z'], U['x'])
    t2 = E('gb,mh,mc,kha,gk->abc', U['xp'], U['z'], U['y'], f, U['wp'])-E('gb,mh,mc,gkl,kh,la->abc', U['xp'], U['z'], U['y'], f, U['z'], U['x'])
    return t1-t2
def C_III2(S, U):  # III.II, ordered y,x',x
    f = S.f
    t1 = E('ma,gb,mh,khc,gk->abc', U['y'], U['xp'], U['z'], f, U['w'])-E('ma,gb,mh,gkl,kh,lc->abc', U['y'], U['xp'], U['z'], f, U['z'], U['x'])
    t2 = E('gb,kac,gk->abc', U['xp'], f, U['w'])-E('gb,gkl,ka,lc->abc', U['xp'], f, U['z'], U['x'])
    return t1-t2
def C_III3(S, U, v):  # III.III first, rho(x'){rho(y),rho(x)}
    return E('mnk,ma,nc,kb->abc', S.f, U[v]-U['xp'], U['z'], U['w'])-E('mnk,ma,kb,nc->abc', S.f, U[v]-U['xp'], U['y'], U['x'])
def C_III5(S, U, v):  # III.V, ordered x',x,y
    return E('mkl,ma,kb,lc->abc', S.f, U[v]-U['xp'], U['x'], U['y']-U['w'])
def C_III6(S, U, v):  # III.VI, ordered x',y,x
    return E('mnk,ma,kb,nc->abc', S.f, U[v]-U['xp'], U['y'], U['x']-U['z'])
# structures: name, tensor fn, type, points
STRUCT = [
    ('I.II 2b', C_2b, 'r{rr}', ('xp', 'x', 'y')),
    ('I.III C_D (u=z)', lambda S, U: C_I3(S, U, 'z'), 'ord', ('xp', 'x', 'y')),
    ('I.III C_A (u=w)', lambda S, U: C_I3(S, U, 'w'), 'ord', ('xp', 'x', 'y')),
    ('I.V first', C_I5, 'r{rr}', ('xp', 'y', 'yp')),
    ('II.I 1a', C_1a, '{rr}r', ('xp', 'y', 'x')),
    ('II.I 1b', C_1b, 'r{rr}', ('x', 'xp', 'y')),
    ('III.II', C_III2, 'ord', ('y', 'xp', 'x')),
    ("III.III L(y')", lambda S, U: C_III3(S, U, 'yp'), 'r{rr}', ('xp', 'y', 'x')),
    ('III.III L(z)', lambda S, U: C_III3(S, U, 'z'), 'r{rr}', ('xp', 'y', 'x')),
    ("III.V L(y')", lambda S, U: C_III5(S, U, 'yp'), 'ord', ('xp', 'x', 'y')),
    ('III.V L(z)', lambda S, U: C_III5(S, U, 'z'), 'ord', ('xp', 'x', 'y')),
    ("III.VI L(y')", lambda S, U: C_III6(S, U, 'yp'), 'ord', ('xp', 'y', 'x')),
    ('III.VI L(z)', lambda S, U: C_III6(S, U, 'z'), 'ord', ('xp', 'y', 'x')),
]
# ---------------- row colour operators, written as in the row equations (operator arrays) --------------
def rowops(S, U, R):  # R: dict point -> rho operator array (n, dim, dim)
    f = S.f; ac = lambda A, B: E('aij,bjk->abik', A, B)+E('bij,ajk->abik', B, A)   # {A^a,B^b}
    o3 = lambda T, A, B, C: E('abc,aij,bjk,ckl->il', T, A, B, C)
    L = E('ab,bij->aij', U['xp']-U['wp'], R['xp'])
    ops = {}
    # I.II 2b: L^a [f^{hbc}U^{ah}(w){rho^b(x),rho^c(y)} - f^{abc}U^{bd}(x)U^{ce}(y){rho^d(x),rho^e(y)}]
    AXY = ac(R['x'], R['y'])
    ops['I.II 2b'] = E('aij,hbc,ah,bcjk->ik', L, f, U['w'], AXY)-E('aij,abc,bd,ce,dejk->ik', L, f, U['x'], U['y'], AXY)
    for u, nm in [('z', 'I.III C_D (u=z)'), ('w', 'I.III C_A (u=w)')]:
        ops[nm] = E('abd,aij,de,bc,ejk,ckl->il', f, L, U['x'], U[u]-U['y'], R['x'], R['y'])
    AYYp = ac(R['y'], R['yp'])
    ops['I.V first'] = E('aij,adc,db,ce,ebjk->ik', L, f, U['x'], U['z'], AYYp)-E('aij,abc,cd,be,dejk->ik', L, f, U['y'], U['yp'], AYYp)
    AXpY = ac(R['xp'], R['y'])
    Lb = E('ac,cbij->abij', U['wp'], AXpY)-E('cb,ad,ce,deij->abij', U['z'], U['xp'], U['y'], AXpY)
    Rb = E('cbd,ac,dij->abij', f, U['w'], R['x'])-E('acd,cb,de,eij->abij', f, U['z'], U['x'], R['x'])
    ops['II.I 1a'] = E('abij,abjk->ik', Lb, Rb)
    Lb2 = E('cbd,ac,dij->abij', f, U['wp'], R['x'])-E('acd,cb,de,eij->abij', f, U['z'], U['x'], R['x'])
    Rb2 = E('ac,cbij->abij', U['w'], AXpY)-E('cb,ad,ce,deij->abij', U['z'], U['xp'], U['y'], AXpY)
    ops['II.I 1b'] = E('abij,abjk->ik', Lb2, Rb2)
    Lb3 = E('cd,ae,cb,dij,ejk->abik', U['y'], U['xp'], U['z'], R['y'], R['xp'])-E('ae,bij,ejk->abik', U['xp'], R['y'], R['xp'])
    ops['III.II'] = E('abij,abjk->ik', Lb3, Rb)
    AYX = ac(R['y'], R['x'])
    for v in ['yp', 'z']:
        tag = "(y')" if v == 'yp' else '(z)'
        Lv = E('db,ce,eij->cdbij', U['z'], U[v]-U['xp'], R['xp'])          # L^{c'd'b}(v)
        R1 = E('ac,cbij->abij', U['w'], AYX)-E('cb,ad,ce,deij->abij', U['z'], U['y'], U['x'], AYX)
        ops['III.III L'+tag] = E('cda,cdbij,abjk->ik', f, Lv, R1)
        LV = E('ce,eij->cij', U[v]-U['xp'], R['xp'])
        RV = E('ce,ad,eij,djk->acik', U['x'], U['y']-U['w'], R['x'], R['y'])
        ops['III.V L'+tag] = E('cqa,cij,aqjk->ik', f, LV, RV)
        # III.VI in the author's original form with U(z)U(z): [U(v)U^{d'b}(z) rho - U^{d'b}(z)U(x') rho][U^{cb}(z)U^{ad}(y)U^{ce}(x) rho^d(y) rho^e(x) - U^{ad}(y) rho^d(y) rho^b(x)]
        RVI = E('cb,ad,ce,dij,ejk->abik', U['z'], U['y'], U['x'], R['y'], R['x'])-E('ad,dij,bjk->abik', U['y'], R['y'], R['x'])
        ops['III.VI L'+tag] = E('cda,cdbij,abjk->ik', f, Lv, RVI)
    return ops
# ---------------- decomposition (symmetrised + contact + one-rho), generic identities ----------------
def decomp(S, T, typ, rp, site):
    # rp: rho arrays at p1,p2,p3 ; site: site labels ; returns operator of RHS
    f = S.f; A, B, Cr = rp; dim = A.shape[1]
    # {rho^a rho^b rho^c}_sym contracted with T
    def trip(T, X, Y, Z):   # sum_{abc} T^{abc} X^a Y^b Z^c  (X,Y,Z operator arrays with their own colour index)
        W = np.tensordot(T, Z, axes=([2], [0]))            # (a,b,k,l)
        W = np.einsum('bjk,abkl->ajl', Y, W, optimize=True)  # (a,j,l)
        return np.einsum('aij,ajl->il', X, W, optimize=True)
    ops = [A, B, Cr]; sym = 0
    for perm in itertools.permutations(range(3)):
        Tp = np.transpose(T, perm)   # T'^{p0 p1 p2} with indices following the operators
        sym = sym+trip(Tp, ops[perm[0]], ops[perm[1]], ops[perm[2]])
    sym = sym/6
    d12 = float(site[0] == site[1]); d23 = float(site[1] == site[2]); d13 = float(site[0] == site[2])
    s2 = lambda X, Y: 0.5*(E('rij,cjk->rcik', X, Y)+E('cij,rjk->rcik', Y, X))   # {X^r Y^c}_sym
    c12 = E('abc,abr,rcik->ik', T, f, s2(A, Cr))*d12      # f^{abr} delta(p1-p2) {rho^r(p1) rho^c(p3)}
    c23 = E('abc,bcr,arik->ik', T, f, s2(A, B))*d23       # f^{bcr} delta(p2-p3) {rho^a(p1) rho^r(p2)}
    c13 = E('abc,acr,rbik->ik', T, f, s2(A, B))*d13       # f^{acr} delta(p1-p3) {rho^r(p1) rho^b(p2)}
    dd = d12*d13
    if typ == 'ord':
        one = -1/6*(E('abc,abr,rcs,sij->ij', T, f, f, A)-E('abc,bcr,ras,sij->ij', T, f, f, A))*dd
        return sym+0.5j*(c12+c23+c13)+one
    if typ == 'r{rr}':
        one = -1/6*(E('abc,abr,rcs,sij->ij', T, f, f, A)+E('abc,acr,rbs,sij->ij', T, f, f, A))*dd
        return 2*sym+1j*(c12+c13)+one
    if typ == '{rr}r':
        one = 1/6*(E('abc,bcr,ras,sij->ij', T, f, f, A)+E('abc,acr,rbs,sij->ij', T, f, f, A))*dd
        return 2*sym+1j*(c23+c13)+one
def _trip(T, X, Y, Z):
    W = np.tensordot(T, Z, axes=([2], [0]))
    W = np.einsum('bjk,abkl->ajl', Y, W, optimize=True)
    return np.einsum('aij,ajl->il', X, W, optimize=True)
def ordered(T, typ, rp):
    A, B, C = rp
    o = _trip(T, A, B, C)
    if typ == 'ord': return o
    if typ == 'r{rr}': return o+_trip(np.transpose(T, (0, 2, 1)), A, C, B)
    if typ == '{rr}r': return o+_trip(np.transpose(T, (1, 0, 2)), B, A, C)
if __name__ == '__main__':
    # ===== (a) tensor identities =====
    for N in [3, 4]:
        S = SUN(N); rng = np.random.default_rng(10+N); f = S.f; worst = {}
        for trial in range(4):
            U = randUs(S, rng)
            # I.IV vs III.I : charges (x', y, x)
            L = U['xp']-U['wp']
            T_I4 = E('aB,cdq,ac,bd,be->Beq', L, f, U['w'], U['z'], U['y'])-E('aB,abc,be,cq->Beq', L, f, U['y'], U['x'])
            # III.I: [U^{ab'}(x')U^{c'd'}(y)U^{c'b}(z) rho^{b'}(x')rho^{d'}(y) - U^{d'b}(z)U^{d'e'}(y)U^{ac'}(w') rho^{c'}(x')rho^{e'}(y)] [f^{cbd}U^{ac}(w) rho^d(x) - f^{acd}U^{cb}(z)U^{de}(x) rho^e(x)]
            Lt = E('aB,cD,cb->aBbD', L, U['y'], U['z'])
            Rt = E('cbQ,ac->abQ', f, U['w'])-E('acd,cb,dQ->abQ', f, U['z'], U['x'])
            T_III1 = E('aBbD,abQ->BDQ', Lt, Rt)
            worst['I.IV + III.I (colour tensors equal)'] = max(worst.get('I.IV + III.I (colour tensors equal)', 0), np.abs(T_I4-T_III1).max()/np.abs(T_I4).max())
            # 11/3: Row V UV colour [U(x')-U(w')]^{ab'}[U^{af}(y)-U^{af}(x)] after x->w, y->x equals -[U(x')-U(w')][U(w)-U(x)]
            cV = E('ab,af->bf', L, U['y']-U['x']); cV_rel = E('ab,af->bf', L, U['x']-U['w']); cII = E('ab,af->bf', L, U['w']-U['x'])
            worst['11/3: Row V colour relabelled = -Row II colour'] = max(worst.get('11/3: Row V colour relabelled = -Row II colour', 0), np.abs(cV_rel+cII).max())
            # pair identities (y' -> w')
            Uq = dict(U); Uq['yp'] = U['wp']
            CA = C_I3(S, U, 'w'); CD = C_I3(S, U, 'z')
            worst['pair A: f L_V R_V = -C_A'] = max(worst.get('pair A: f L_V R_V = -C_A', 0), np.abs(C_III5(S, Uq, 'yp')+CA).max())
            # pair B: C_B = L^a[f^{adc}U^{db}(z)U^{ce}(w){rho^e(y),rho^b(x)} - f^{abc}U^{cd}(y)U^{be}(x){rho^d(y),rho^e(x)}], type rho(x'){rho(y),rho(x)}
            CB = E('ga,gdc,db,ce->aeb', L, f, U['z'], U['w'])-E('ga,gbc,cd,be->ade', L, f, U['y'], U['x'])
            worst['pair B: f L R_1 = -C_B'] = max(worst.get('pair B: f L R_1 = -C_B', 0), np.abs(C_III3(S, Uq, 'yp')+CB).max())
            # pair C (two-rho, rho(x') rho(x)): C_C^{b'h} = [U(x')-U(w')]^{ab'}[N U^{ah}(x) - f^{adc}f^{hbe}U^{db}(z)U^{ce}(w)]
            CC = E('aB,ah->Bh', L, N*U['x'])-E('aB,adc,hbe,db,ce->Bh', L, f, f, U['z'], U['w'])
            # f^{c'd'a} L^{c'd'b}(w') R_2^{ab}, R_2^{ab} = f^{cbd}U^{ac}(w) rho^d(x) - f^{acd}U^{cb}(z)U^{de}(x) rho^e(x)
            Lw = E('db,cB->cdbB', U['z'], U['wp']-U['xp'])
            R2 = E('cbd,ac->abd', f, U['w'])-E('acd,cb,de->abe', f, U['z'], U['x'])
            LR2 = E('cda,cdbB,abh->Bh', f, Lw, R2)
            worst['pair C: f L R_2 = -C_C'] = max(worst.get('pair C: f L R_2 = -C_C', 0), np.abs(LR2+CC).max())
            # pair D: f L R_VI|_{x<->y} = +C_D ; R_VI type rho(x') rho(y) rho(x): swap x<->y in Wilson lines and charge positions
            Us = dict(Uq); Us['x'], Us['y'] = U['y'], U['x']
            T6 = C_III6(S, Us, 'yp')            # charges (x', y->x, x->y) = (x', x, y) order after the swap
            worst['pair D: f L R_VI|_{x<->y} = +C_D'] = max(worst.get('pair D: f L R_VI|_{x<->y} = +C_D', 0), np.abs(T6-CD).max())
            Usw = dict(U); Usw['x'], Usw['y'] = U['y'], U['x']
            Tsw = C_2b(S, Usw).transpose(0, 2, 1)   # charges x<->y: rho(x'){rho(y),rho(x)} ; anticommutator symmetric in the last two
            T0 = C_2b(S, U); Tsym = 0.5*(T0+T0.transpose(0, 2, 1)); Tsws = 0.5*(Tsw+Tsw.transpose(0, 2, 1))
            worst['I.II 2b colour antisymmetric under x<->y'] = max(worst.get('I.II 2b colour antisymmetric under x<->y', 0), np.abs(Tsym+Tsws).max())
        for k, v in worst.items(): pr('(a) SU(%d) %s: max deviation %.1e' % (N, k, v))
    # ===== (b) operator level =====
    S = SUN(3); rng = np.random.default_rng(3); n = S.n; dim = 27
    I3 = np.eye(3)
    rho_site = []
    for s in range(3):
        mats = []
        for a in range(n):
            m = [I3, I3, I3]; m[s] = S.t[a]
            mats.append(np.kron(np.kron(m[0], m[1]), m[2]))
        rho_site.append(np.array(mats))
    worst = 0; worst_ctrl = 1e9; ncase = 0
    names = [s[0] for s in STRUCT]
    for variant in ['equal', 'independent']:
        U0 = randUs(S, rng)
        for name, Cfn, typ, pts in STRUCT:
            for sites in itertools.product(range(3), repeat=3):
                U = dict(U0)
                if variant == 'equal':
                    # Wilson lines equal at coinciding charge sites
                    for i in range(3):
                        for j in range(i):
                            if sites[i] == sites[j]: U[pts[i]] = U[pts[j]]
                site_of = {pts[i]: sites[i] for i in range(3)}
                R = {p: rho_site[site_of[p]] if p in site_of else None for p in PTS}
                ops = rowops(S, U, {p: (R[p] if R[p] is not None else rho_site[0]) for p in PTS})
                lhs = ops[name]
                T = Cfn(S, U)
                rp = [rho_site[s] for s in sites]
                rhs = decomp(S, T, typ, rp, sites)
                o2 = ordered(T, typ, rp)
                err = np.abs(lhs-rhs).max()/max(np.abs(lhs).max(), 1e-12)
                err2 = np.abs(lhs-o2).max()/max(np.abs(lhs).max(), 1e-12)
                worst = max(worst, err, err2); ncase += 1
    pr('(b) operator level: row colour operator == tensor C^{abc} ordered == sym + contact + one-rho, %d cases, max rel dev %.1e' % (ncase, worst))
    # control: flip the sign of the contact term in one case
    U = randUs(S, rng); name, Cfn, typ, pts = STRUCT[1]; sites = (0, 0, 1)
    T = Cfn(S, U); rp = [rho_site[s] for s in sites]
    good = ordered(T, typ, rp); A_, B_, C_ = rp
    f = S.f
    s2 = lambda X, Y: 0.5*(E('rij,cjk->rcik', X, Y)+E('cij,rjk->rcik', Y, X))
    sym = sum(E('abc,%sij,%sjk,%skl->il' % tuple('abc'[p] for p in perm), T, *[rp[p] for p in perm]) for perm in itertools.permutations(range(3)))/6
    wrong = sym-0.5j*E('abc,abr,rcik->ik', T, f, s2(A_, C_))
    pr('(b) control (sign of the contact term flipped, I.III C_D, x\'=x): rel dev %.2f' % (np.abs(good-wrong).max()/np.abs(good).max()))
    open('s3_colour.log', 'w').write('\n'.join(out)+'\n')
