"""Group I split into the 'direct' pieces and the unitarity combination, compared with the reference."""
import itertools
import numpy as np
import compare_ref as CR
import four_check as FC
from cancel3 import term_list

class Acc(CR.Accum):
    def add_user_fine(self, key_, cc, n, c, src, colour, word, cmap):
        rl = [p for _, p in word]
        name = 'U:%s.%d%s' % (key_, n + 1, '*' if cc else '')
        for perm in itertools.permutations(['p1', 'p2', 'p3']):
            m = dict(cmap); m.update(dict(zip(rl, perm)))
            facs = [('U', u[1], u[2], m[u[3]]) if u[0] == 'U' else u for u in colour]
            facs += [('rho', i, m[p]) for i, p in word]
            cv = self.col(facs)
            for kc, fs in CR.user_monomials(src):
                parts, sign = [], 1
                for f in fs:
                    if f[0] == 'KK':
                        k, s = CR.kk_key(*(m[l] for l in f[1:])); parts.append(k); sign *= s
                    else:
                        parts.append(('Phi', m[f[2][0]], m['w']))
                self.data[tuple(sorted(parts, key=str))][name] += c * kc * sign * cv / 6

def build(seed):
    A = Acc(seed)
    for t in CR.ref_terms():
        if t.tag.split(':')[1].startswith('G1'):
            A.add_ref(t)
    for u in term_list():
        if u[0].startswith('GI.'):
            A.add_user_fine(*u)
    for t in FC.pbw3_terms():
        if t.tag.startswith('4:G1'):
            A.add_ref(t, prefix='N4 ')
    return A

def subs(row, ks):
    return lambda n: n.startswith('U:%s.' % row) and int(n.split('.')[-1].rstrip('*')) in ks
def n4(line, k2):
    return lambda n: n.startswith('N4 4:G1:%s.' % line) and n.rstrip('*').endswith('.%d' % k2)
def ref(*rows):
    return lambda n: any(n == 'REF %s:G1-%s' % (w, r) for r in rows for w in '34')
def any_of(*ps):
    return lambda n: any(p(n) for p in ps)

GROUPS = [
 ('normalization (Row I)',            n4('L1', 1) , None),
 ('U A3          (Row II direct)',    any_of(subs('GI.II', (1, 3)), n4('L2', 1)), ref('II_A3')),
 ('B3bar^dag U A (Row III direct)',   any_of(subs('GI.III', (1, 3)), n4('L4', 2)), ref('III')),
 ('Abar^dag UU B2 (Row IV direct)',   any_of(subs('GI.IV', (1, 3)), n4('L3', 1)), ref('IV')),
 ('Cbar^dag UU B2 (Row V direct)',    any_of(subs('GI.V1', (1, 3)), subs('GI.V2', (1, 3))), ref('V')),
 ('unitarity combination for Abar(3)dag',
   any_of(subs('GI.II', (2, 4)), subs('GI.III', (2, 4)), subs('GI.IV', (2, 4)), subs('GI.V1', (2, 4)),
          subs('GI.V2', (2, 4)), n4('L2', 2), n4('L3', 2), n4('L4', 1)), ref('II_A3bar')),
]

if __name__ == '__main__':
    for seed in (17, 23):
        A = build(seed)
        keys = list(A.data.keys())
        vec = lambda p: np.array([sum(v for n, v in A.data[k].items() if p(n)) for k in keys])
        print('seed', seed)
        norm_user = vec(any_of(n4('L1', 1), n4('L1', 2)))
        norm_ref = vec(ref('I_Nafter', 'I_Nbarbefore'))
        print(f'   {"normalization (Row I)":42s} |notes - ref| = {np.abs(norm_user - norm_ref).max():.3e}  (|ref| {np.abs(norm_ref).max():.3f})')
        for name, pu, pr in GROUPS[1:]:
            u, r = vec(pu), vec(pr)
            print(f'   {name:42s} |notes - ref| = {np.abs(u - r).max():.3e}  (|notes| {np.abs(u).max():.3f}, |ref| {np.abs(r).max():.3f})')
        print(f'   {"row I prime (reference only)":42s} |ref| = {np.abs(vec(ref("Ia"))).max():.3f}')
