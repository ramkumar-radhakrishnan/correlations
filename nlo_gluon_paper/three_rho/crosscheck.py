"""Independent checks against the region-T leading-log reference; results -> crosscheck.json.
All with symmetric adjoint Wilson lines (the convention U^{bc} = U^{cb} of the notes)."""
import itertools
import json
import numpy as np
import sympy as sp
import compare_ref as CR
import four_check as FC
import g1split as GS
from corrections import correction1_terms, correction2_terms, correction1_variant_terms


def add_ct(A, t, name):
    rl = []
    for _, p in t['word']:
        if p not in rl:
            rl.append(p)
    c = complex(sp.N(t['coef']))
    for perm in itertools.permutations(['p1', 'p2', 'p3']):
        m = dict(t['cmap']); m.update(dict(zip(rl, perm)))
        facs = [('U', u[1], u[2], m.get(u[3], u[3])) if u[0] == 'U' else u for u in t['colour']]
        facs += [('rho', i, m[p]) for i, p in t['word']]
        cv = A.col(facs)
        for kc, fs in t['monos']:
            parts, sign = [], 1
            for f in fs:
                if f[0] == 'KK':
                    k, s = CR.kk_key(*(m.get(l, l) for l in f[1:])); parts.append(k); sign *= s
                else:
                    parts.append(('Phi', m[f[2][0]], m['w']))
            A.data[tuple(sorted(parts, key=str))][name] += c * float(kc) * sign * cv / 6


def main(seeds=(17, 23, 31)):
    out = {'four_rho': [], 'AB_vs_N4': [], 'rows': {}, 'group1': [], 'group1_pieces': {}}
    for seed in seeds:
        out['four_rho'].append(FC.check_four(seed))
        A = CR.Accum(seed)
        for t in CR.ref_terms():
            A.add_ref(t)
        for u in CR.user_terms():
            A.add_user(*u[1:])
        for t in FC.pbw3_terms():
            A.add_ref(t, prefix='N4 ')
        keys = list(A.data.keys())
        vec = lambda p: np.array([sum(v for n, v in A.data[k].items() if p(n)) for k in keys])
        ab = vec(lambda n: n.startswith('USR ') and n.split()[1][0] in 'AB')
        n4 = vec(lambda n: n.startswith('N4 '))
        out['AB_vs_N4'].append((np.abs(ab - n4).max(), np.abs(n4).max()))
        n4r = lambda row: (lambda n: n.startswith('N4 4:%s:' % row))
        usr = lambda key: (lambda n: n.startswith('USR ') and n.split()[1].rstrip('*') == key)
        rows = [('G2-I', GS.any_of(usr('GII.I1'), usr('GII.I2'), n4r('G2-I'))),
                ('G2-II', n4r('G2-II')), ('G2-III', n4r('G2-III')),
                ('G3-I', GS.any_of(usr('GIII.I'), n4r('G3-I'))), ('G3-II', GS.any_of(usr('GIII.II'), n4r('G3-II'))),
                ('G3-III', usr('GIII.III')), ('G3-IV', n4r('G3-IV')),
                ('G3-V', usr('GIII.V')), ('G3-VI', usr('GIII.VI'))]
        for name, pu in rows:
            pr = lambda n, name=name: n in ('REF 3:' + name, 'REF 4:' + name)
            u, r = vec(pu), vec(pr)
            out['rows'].setdefault(name, []).append((np.abs(u - r).max(), np.abs(r).max()))
        # Group I as a whole, as written and corrected (both variants of correction 1)
        G = GS.build(seed)
        for t in correction1_terms():
            add_ct(G, t, 'FIX1 ')
        for t in correction1_variant_terms():
            add_ct(G, t, 'FIX1v ')
        for t in correction2_terms():
            add_ct(G, t, 'FIX2 ')
        keys = list(G.data.keys())
        vg = lambda p: np.array([sum(v for n, v in G.data[k].items() if p(n)) for k in keys])
        user = vg(lambda n: n.startswith('U:') or n.startswith('N4 '))
        f1, f1v, f2 = vg(lambda n: n.startswith('FIX1 ')), vg(lambda n: n.startswith('FIX1v ')), vg(lambda n: n.startswith('FIX2 '))
        ref = vg(lambda n: n.startswith('REF '))
        out['group1'].append(dict(written=np.abs(user - ref).max(), corrected=np.abs(user + f1 + f2 - ref).max(),
                                  corrected_variant=np.abs(user + f1v + f2 - ref).max(), scale=np.abs(ref).max()))
        # Group I pieces as written
        pieces = {'norm': (GS.any_of(GS.n4('L1', 1), GS.n4('L1', 2)), GS.ref('I_Nafter', 'I_Nbarbefore'))}
        for name, pu, pr in GS.GROUPS[1:]:
            pieces[name] = (pu, pr)
        for name, (pu, pr) in pieces.items():
            u, r = vg(pu), vg(pr)
            out['group1_pieces'].setdefault(name, []).append((np.abs(u - r).max(), np.abs(r).max()))
        # the same pieces after the corrections
        H = GS.build(seed)
        for t in correction1_terms():
            add_ct(H, t, 'C1:%d' % t['sub'])
        for t in correction2_terms():
            add_ct(H, t, '%s:%d' % ('C2B' if t['name'].startswith('fixB') else 'C2D', t['sub']))
        keys = list(H.data.keys())
        vh = lambda p: np.array([sum(v for n, v in H.data[k].items() if p(n)) for k in keys])
        sub = lambda tag, ks: (lambda n: n.startswith(tag + ':') and int(n.split(':')[1]) in ks)
        after = {
            'II_direct': (GS.any_of(GS.subs('GI.II', (1, 3)), GS.n4('L2', 1), sub('C1', (1, 3))), GS.ref('II_A3')),
            'III_direct': (GS.any_of(GS.subs('GI.III', (1, 3)), GS.n4('L4', 2), sub('C2B', (1, 3))), GS.ref('III')),
            'D_direct': (sub('C2D', (1, 3)), GS.ref('Ia')),
            'unitarity': (GS.any_of(GS.GROUPS[5][1], sub('C2B', (2, 4)), sub('C2D', (2, 4))), GS.ref('II_A3bar')),
        }
        for name, (pu, pr) in after.items():
            u, r = vh(pu), vh(pr)
            out.setdefault('group1_after', {}).setdefault(name, []).append((np.abs(u - r).max(), np.abs(r).max()))
    json.dump(out, open('crosscheck.json', 'w'), indent=1, default=float)
    return out


if __name__ == '__main__':
    o = main()
    print(json.dumps(o, indent=1, default=float)[:3000])
