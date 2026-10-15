"""Colour networks of the three-rho terms and how they cancel.

Every term is f^{abc} L1^a L2^b L3^c, where each leg L is a chain of Wilson lines acting on one charge,
L = U(p_1) U(p_2) ... U(p_n) rho(q)  (p_1 next to f).  With your convention U^{bc} = U^{cb} the orientation of each U
is dropped.  The canonical form uses only exact identities:
  (i)  U(p) U(p) = 1 on a leg (U^T U = 1 with U = U^T);
  (ii) f^{abc} (U(p)A)^a (U(p)B)^b C^c = f^{abc} A^a B^b (U(p)C)^c  (adjoint U preserves f), applied when two legs
       start with the same U(p); it lowers the number of U's by one (or by three with (i));
  (iii) the legs are ordered by f^{abc} = -f^{bac}, ..., with the sign.
Each canonical network is checked against a numerical fingerprint (random symmetric U), so equal networks are
never split and different ones never merged.  Output: networks3.json."""
import json
import itertools
from collections import defaultdict, OrderedDict
import numpy as np
import env  # noqa
import numerics as nm
from core import NC
import parse3 as P

PTS = ['w', "w'", 'z', 'x', "x'", 'y']


def legs_of(e):
    """returns [(chain, q), ...] for the three f indices, chain = tuple of U points from f outwards"""
    (fa, fb, fc), = e['fs']
    rho = {i: p for i, p in e['rs']}
    Us = list(e['Us'])
    out = []
    for idx in (fa, fb, fc):
        chain = []
        cur = idx
        used = set()
        while cur not in rho:
            nxt = None
            for k, (i, j, p) in enumerate(Us):
                if k in used:
                    continue
                if i == cur:
                    nxt, used_k = (j, p), k
                    break
                if j == cur:
                    nxt, used_k = (i, p), k
                    break
            assert nxt is not None, (e['block'], e['n'], idx, cur)
            used.add(used_k)
            chain.append(nxt[1])
            cur = nxt[0]
        out.append((tuple(chain), rho[cur]))
    return out


def reduce_leg(chain):
    ch = list(chain)
    k = 0
    while k < len(ch) - 1:
        if ch[k] == ch[k + 1]:
            del ch[k:k + 2]
            k = max(k - 1, 0)
        else:
            k += 1
    return tuple(ch)


def canonical(legs):
    """legs: list of 3 (chain, q); returns (sign, legs) canonical"""
    legs = [(reduce_leg(c), q) for c, q in legs]
    changed = True
    while changed:
        changed = False
        for i, j in ((0, 1), (0, 2), (1, 2)):
            ci, cj = legs[i][0], legs[j][0]
            if ci and cj and ci[0] == cj[0]:
                p = ci[0]
                k = 3 - i - j
                legs[i] = (ci[1:], legs[i][1])
                legs[j] = (cj[1:], legs[j][1])
                legs[k] = (reduce_leg((p,) + legs[k][0]), legs[k][1])
                changed = True
                break
    key = lambda l: (len(l[0]), [PTS.index(p) for p in l[0]], PTS.index(l[1]))
    order = sorted(range(3), key=lambda m: key(legs[m]))
    perm_sign = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (1, 0, 2): -1, (0, 2, 1): -1, (2, 1, 0): -1}[tuple(order)]
    return perm_sign, tuple(legs[m] for m in order)


def leg_tex(l):
    chain, q = l
    us = ''.join(r'U_{%s}' % p for p in chain)
    return r'%s\rho_{%s}' % (us, q)


def net_tex(legs):
    return r'f\big[%s\big]' % ',\,'.join(leg_tex(l) for l in legs)


def fingerprint(legs, cfgs, f):
    out = []
    for U, r in cfgs:
        vs = []
        for chain, q in legs:
            v = r[q].copy()
            for p in reversed(chain):
                v = U[p] @ v
            vs.append(v)
        out.append(np.einsum('abc,a,b,c->', f, *vs))
    return np.array(out)


def analyse(E):
    ts, f = nm.structure_constants(NC)
    cfgs = []
    for s in (3, 11, 29):
        rng = np.random.default_rng(s)
        cfgs.append(({p: nm.random_adjoint(NC, ts, rng, True) for p in PTS}, {p: rng.normal(size=NC * NC - 1) for p in PTS}))
    res = OrderedDict()
    for cls in ('alpha', "alpha'", 'beta', "beta'"):
        nets = []           # list of dicts: legs, fp, terms
        for e in [e for e in E if e['cls'] == cls]:
            legs0 = legs_of(e)
            s, legs = canonical(legs0)
            fp_raw = fingerprint(legs0, cfgs, f)
            fp = fingerprint(legs, cfgs, f)
            assert np.allclose(s * fp, fp_raw), (e['block'], e['n'])
            coef = (e['pref'] / 1j) * e['sign'] * e['mult'] * s          # units of i x common prefactor
            for net in nets:
                if np.allclose(net['fp'], fp, atol=1e-10):
                    net['terms'].append((e, coef, legs0))
                    break
                if np.allclose(net['fp'], -fp, atol=1e-10):
                    net['terms'].append((e, -coef, legs0))
                    break
            else:
                nets.append(dict(legs=legs, fp=fp, terms=[(e, coef, legs0)]))
        res[cls] = nets
    return res


if __name__ == '__main__':
    E = P.parse()
    R = analyse(E)
    out = {}
    for cls, nets in R.items():
        print('==', cls, 'networks', len(nets))
        out[cls] = []
        for net in nets:
            tot = sum(c for _, c, _ in net['terms'])
            nUs = sorted(set(e['nU'] for e, _, _ in net['terms']))
            print('  %-60s sum %+.3f  nU %s  terms %s' % (net_tex(net['legs']), tot.real, nUs,
                  ['%s.%d(%+g)' % (e['block'], e['n'], c.real) for e, c, _ in net['terms']]))
            out[cls].append(dict(net=net_tex(net['legs']), sum=[tot.real, tot.imag], nU=nUs,
                                 terms=[dict(block=e['block'], n=e['n'], coef=[c.real, c.imag], nU=e['nU']) for e, c, _ in net['terms']]))
    json.dump(out, open('networks3.json', 'w'), indent=1)
