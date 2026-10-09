"""End-to-end check of the reordering with explicit charge operators.

Charges rho^a(s) live on three discrete sites s (the delta-function becomes a
Kronecker delta).  For every sub-term we build the operator
    O = sum_{sites} K(sites) C^{i1 i2 i3} (rho rho rho in the given order),
its fully symmetrized part S(O), and the two-rho pieces of process.py.
In a colour-singlet state the one-rho remainder has zero expectation value, so
    Tr[P (O - S(O))] = Tr[P (two-rho pieces)]
must hold, P the projector on the singlet subspace."""
import itertools
import string
import numpy as np
import sympy as sp
from engine import Nc
import numerics as nm
import process


def kron_list(ms):
    out = ms[0]
    for m in ms[1:]:
        out = np.kron(out, m)
    return out


def site_ops(N, nsite, per_site):
    """rho^a(s) operators; per_site = number of fundamental particles per site."""
    ts, f = nm.structure_constants(N)
    d1 = N
    dsite = d1 ** per_site
    # generators on one site
    gens = []
    for t in ts:
        g = np.zeros((dsite, dsite), complex)
        for k in range(per_site):
            ms = [np.eye(d1)] * per_site
            ms = list(ms)
            ms[k] = t
            g += kron_list(ms)
        gens.append(g)
    ops = {}
    for s in range(nsite):
        for a, g in enumerate(gens):
            ms = [np.eye(dsite)] * nsite
            ms = list(ms)
            ms[s] = g
            ops[(a, s)] = kron_list(ms)
    dim = dsite ** nsite
    # singlet projector: null space of total Casimir
    C = np.zeros((dim, dim), complex)
    for a in range(len(ts)):
        Qa = sum(ops[(a, s)] for s in range(nsite))
        C += Qa @ Qa
    w, v = np.linalg.eigh(C)
    null = v[:, np.abs(w) < 1e-9]
    return ts, f, ops, null, null.shape[1]


def colour_tensor(facs, open_idx, f, Um, ncol):
    """contract all colour factors except rho's, leaving open_idx open."""
    letters = {}
    pool = iter(string.ascii_letters)
    ops, arrs = [], []
    for fac in facs:
        if fac[0] == 'rho':
            continue
        idx = [fac[1], fac[2]] if fac[0] in ('U', 'd') else [fac[1], fac[2], fac[3]]
        sub = ''
        for x in idx:
            if x not in letters:
                letters[x] = next(pool)
            sub += letters[x]
        ops.append(sub)
        arrs.append(Um[fac[3]] if fac[0] == 'U' else (f if fac[0] == 'f' else np.eye(ncol)))
    out = ''
    for x in open_idx:
        if x not in letters:
            letters[x] = next(pool)
        out += letters[x]
    if not arrs:
        # only deltas between open indices (should not happen)
        raise RuntimeError
    return np.einsum(','.join(ops) + '->' + out, *arrs, optimize='greedy')


PATH3 = None
PATH2 = None

def run(N, per_site, seed, nsamp=1):
    global PATH3, PATH2
    ts, f, ops, Sv, nsing = site_ops(N, 3, per_site)
    Sc = Sv.conj()
    ncol = N * N - 1
    rng = np.random.default_rng(seed)
    res = process.process()
    D = Sv.shape[0]; ns = Sv.shape[1]
    C3 = np.zeros((ncol,)*3); Rd = np.zeros((ncol, D, D)); 
    PATH3 = np.einsum_path('abc,in,aij,bjk,ckl,ln->', C3, Sc, Rd, Rd, Rd, Sv, optimize='greedy')[0]
    PATH2 = np.einsum_path('ab,in,aij,bjk,kn->', np.zeros((ncol,ncol)), Sc, Rd, Rd, Sv, optimize='greedy')[0]
    worst = 0.0
    scale = 0.0
    for _ in range(nsamp):
        for sres in res:
            src = sres['src']
            rl = src.rho_labels
            fixed = {p: nm.random_adjoint(N, ts, rng, True) for p in src.fixed_labels}
            Usite = [nm.random_adjoint(N, ts, rng, True) for _ in range(3)]
            K = rng.normal(size=(3, 3, 3))
            for sub in sres['subs']:
                word = sub['word']
                facs = sub['colour']
                lhs = 0
                for sites in itertools.product(range(3), repeat=3):
                    lab2site = dict(zip(rl, sites))
                    Um = dict(fixed)
                    for lab, s in lab2site.items():
                        Um[lab] = Usite[s]
                    kval = K[tuple(lab2site[l] for l in src.kargs)]
                    C = colour_tensor(facs, [i for i, _ in word], f, Um, ncol)
                    # list of operator words
                    if src.form == 'plain':
                        words = [word]
                    elif src.form == 'A{BC}':
                        A, B, Cc = word
                        words = [[A, B, Cc], [A, Cc, B]]
                    else:
                        B, Cc, A = word
                        words = [[B, Cc, A], [Cc, B, A]]
                    Rs = {lab: np.array([ops[(a, lab2site[lab])] for a in range(ncol)]) for lab in rl}
                    for wd in words:
                        # colour tensor C is indexed in the order of `word`; map to wd order
                        perm_w = [ [i for i, _ in word].index(i) for i, _ in wd ]
                        Cw = np.transpose(C, perm_w)
                        R = [Rs[p] for _, p in wd]
                        full = np.einsum('abc,in,aij,bjk,ckl,ln->', Cw, Sc, R[0], R[1], R[2], Sv, optimize=PATH3)
                        sym = 0
                        for pm in itertools.permutations(range(3)):
                            Cp = np.transpose(Cw, pm)
                            sym = sym + np.einsum('abc,in,aij,bjk,ckl,ln->', Cp, Sc, R[pm[0]], R[pm[1]], R[pm[2]], Sv, optimize=PATH3)
                        lhs = lhs + sub['sign'] * kval * (full - sym / 6)
                rhs = 0
                for pc in sub['pieces']:
                    if not pc['facs_after']:
                        continue
                    old, new = pc['subst']
                    surv = [l for l in rl if l != old]
                    rhos = [x for x in pc['facs_after'] if x[0] == 'rho']
                    coeff = complex(sp.N(pc['coeff_after'].subs(Nc, N)))
                    for sites in itertools.product(range(3), repeat=2):
                        lab2site = dict(zip(surv, sites))
                        lab2site[old] = lab2site[new]
                        Um = {p: fixed[p] for p in fixed}
                        for lab in rl:
                            Um[lab] = Usite[lab2site[lab]]
                        kval = K[tuple(lab2site[l] for l in src.kargs)]
                        C = colour_tensor(pc['facs_after'], [r[1] for r in rhos], f, Um, ncol)
                        R1 = np.array([ops[(a, lab2site[rhos[0][2]])] for a in range(ncol)])
                        R2 = np.array([ops[(a, lab2site[rhos[1][2]])] for a in range(ncol)])
                        val = np.einsum('ab,in,aij,bjk,kn->', C, Sc, R1, R2, Sv, optimize=PATH2) + np.einsum('ab,in,bij,ajk,kn->', C, Sc, R2, R1, Sv, optimize=PATH2)
                        rhs = rhs + coeff * kval * val
                err = abs(lhs - rhs)
                worst = max(worst, err)
                scale = max(scale, abs(lhs))
    return nsing, worst, scale


if __name__ == '__main__':
    for N, per_site in ((2, 2), (3, 1)):
        nsing, worst, scale = run(N, per_site, seed=7)
        print(f'SU({N}), {per_site} particle(s)/site: singlets={nsing}, max |LHS-RHS| = {worst:.2e}, max |LHS| = {scale:.2e}')
