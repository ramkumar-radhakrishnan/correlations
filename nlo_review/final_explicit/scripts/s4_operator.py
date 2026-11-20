# Operator-level check (faster version of s3 part (b)): for every cubic colour structure and every assignment
# of its three charge points to three sites (SU(3) fundamental charges, 27-dim space), the row's colour operator
# as written in the row equations == tensor C^{abc} with the charges in the row's order == (i)+(ii)+(iii).
import numpy as np, itertools, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fe_colour import SUN
import s3_colour as s3
from s3_colour import STRUCT, randUs, decomp, ordered, PTS
E = lambda *a: np.einsum(*a, optimize='greedy')
s3.E = E
def rowop(S, U, R, name):
    f = S.f
    ac = lambda A, B: E('aij,bjk->abik', A, B)+E('bij,ajk->abik', B, A)
    L = E('ab,bij->aij', U['xp']-U['wp'], R['xp'])
    if name == 'I.II 2b':
        AXY = ac(R['x'], R['y'])
        return E('aij,hbc,ah,bcjk->ik', L, f, U['w'], AXY)-E('aij,abc,bd,ce,dejk->ik', L, f, U['x'], U['y'], AXY)
    if name.startswith('I.III'):
        u = 'z' if 'u=z' in name else 'w'
        LU = E('abd,aij->bdij', f, L)
        return E('bdij,de,bc,ejk,ckl->il', LU, U['x'], U[u]-U['y'], R['x'], R['y'])
    if name == 'I.V first':
        AYYp = ac(R['y'], R['yp'])
        return E('aij,adc,db,ce,ebjk->ik', L, f, U['x'], U['z'], AYYp)-E('aij,abc,cd,be,dejk->ik', L, f, U['y'], U['yp'], AYYp)
    if name in ('II.I 1a', 'II.I 1b', 'III.II'):
        AXpY = ac(R['xp'], R['y'])
        Rb = E('cbd,ac,dij->abij', f, U['w'], R['x'])-E('acd,cb,de,eij->abij', f, U['z'], U['x'], R['x'])
        if name == 'II.I 1a':
            Lb = E('ac,cbij->abij', U['wp'], AXpY)-E('cb,ad,ce,deij->abij', U['z'], U['xp'], U['y'], AXpY)
            return E('abij,abjk->ik', Lb, Rb)
        if name == 'II.I 1b':
            Lb2 = E('cbd,ac,dij->abij', f, U['wp'], R['x'])-E('acd,cb,de,eij->abij', f, U['z'], U['x'], R['x'])
            Rb2 = E('ac,cbij->abij', U['w'], AXpY)-E('cb,ad,ce,deij->abij', U['z'], U['xp'], U['y'], AXpY)
            return E('abij,abjk->ik', Lb2, Rb2)
        Lb3 = E('cd,ae,cb,dij,ejk->abik', U['y'], U['xp'], U['z'], R['y'], R['xp'])-E('ae,bij,ejk->abik', U['xp'], R['y'], R['xp'])
        return E('abij,abjk->ik', Lb3, Rb)
    v = 'yp' if "(y')" in name else 'z'
    if name.startswith('III.III'):
        AYX = ac(R['y'], R['x'])
        Lv = E('db,ce,eij->cdbij', U['z'], U[v]-U['xp'], R['xp'])
        R1 = E('ac,cbij->abij', U['w'], AYX)-E('cb,ad,ce,deij->abij', U['z'], U['y'], U['x'], AYX)
        return E('cda,cdbij,abjk->ik', f, Lv, R1)
    if name.startswith('III.V '):
        LV = E('ce,eij->cij', U[v]-U['xp'], R['xp'])
        RV = E('ce,ad,eij,djk->acik', U['x'], U['y']-U['w'], R['x'], R['y'])
        return E('cqa,cij,aqjk->ik', f, LV, RV)
    if name.startswith('III.VI'):
        # the author's simplified colour operator: f^{c'd'a}[U^{c'e'}(v)-U^{c'e'}(x')]rho^{e'}(x') U^{ad}(y)[U^{d'e}(x)-U^{d'e}(z)] rho^d(y) rho^e(x)
        LV = E('ce,eij->cij', U[v]-U['xp'], R['xp'])
        Rt = E('ad,qe,dij,ejk->aqik', U['y'], U['x']-U['z'], R['y'], R['x'])
        simp = E('cqa,cij,aqjk->ik', f, LV, Rt)
        # the earlier (checked) form with U(z)U(z): [U(v)U^{d'b}(z) rho - U^{d'b}(z)U(x') rho][U^{cb}(z)U^{ad}(y)U^{ce}(x) rho^d(y)rho^e(x) - U^{ad}(y) rho^d(y) rho^b(x)]
        Lv = E('db,ce,eij->cdbij', U['z'], U[v]-U['xp'], R['xp'])
        RVI = E('cb,ad,ce,dij,ejk->abik', U['z'], U['y'], U['x'], R['y'], R['x'])-E('ad,dij,bjk->abik', U['y'], R['y'], R['x'])
        old = E('cda,cdbij,abjk->ik', f, Lv, RVI)
        return simp, old
    raise KeyError(name)
S = SUN(3); rng = np.random.default_rng(3); n = S.n
I3 = np.eye(3)
rho_site = []
for s in range(3):
    mats = []
    for a in range(n):
        m = [I3, I3, I3]; m[s] = S.t[a]
        mats.append(np.kron(np.kron(m[0], m[1]), m[2]))
    rho_site.append(np.array(mats))
worst = 0; worst6 = 0; ncase = 0; t0 = time.time()
for variant in ['equal', 'independent']:
    U0 = randUs(S, rng)
    for name, Cfn, typ, pts in STRUCT:
        for sites in itertools.product(range(3), repeat=3):
            U = dict(U0)
            if variant == 'equal':
                for i in range(3):
                    for j in range(i):
                        if sites[i] == sites[j]: U[pts[i]] = U[pts[j]]
            site_of = {pts[i]: sites[i] for i in range(3)}
            R = {p: rho_site[site_of.get(p, 0)] for p in PTS}
            lhs = rowop(S, U, R, name)
            if isinstance(lhs, tuple):
                worst6 = max(worst6, np.abs(lhs[0]-lhs[1]).max()/np.abs(lhs[1]).max()); lhs = lhs[0]
            T = Cfn(S, U); rp = [rho_site[s] for s in sites]
            rhs = decomp(S, T, typ, rp, sites); o2 = ordered(T, typ, rp)
            sc = max(np.abs(lhs).max(), 1.0)   # vanishing operators (e.g. I.II at x=y with equal U) compared absolutely
            worst = max(worst, np.abs(lhs-rhs).max()/sc, np.abs(lhs-o2).max()/sc); ncase += 1
    print(variant, 'done', time.time()-t0, flush=True)
s = ('operator level (SU(3) fundamental charges on 3 sites, 27 site assignments x 13 structures x 2 Wilson-line variants = %d cases): '
     'row colour operator as written == ordered tensor product == (i)+(ii)+(iii): max deviation (relative, or absolute for operators of norm < 1) %.1e; Group III Row VI author simplified form == checked form: %.1e' % (ncase, worst, worst6))
# negative control: drop the contact term (ii) for one coincident assignment
name, Cfn, typ, pts = STRUCT[1]; sites = (0, 0, 1); U = randUs(S, rng)
T = Cfn(S, U); rp = [rho_site[q] for q in sites]
full = decomp(S, T, typ, rp, sites); nocontact = decomp(S, T, typ, rp, (0, 1, 2))
s += '; control (contact terms omitted at x\'=x): rel dev %.2f' % (np.abs(full-nocontact).max()/np.abs(full).max())
print(s); open('s4_operator.log', 'w').write(s+'\n')
