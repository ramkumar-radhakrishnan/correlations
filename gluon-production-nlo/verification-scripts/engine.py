"""Operator-level verification engine for the NLO gluon-production colour/kernel algebra.

A Term represents (coefficient) x e^{-ik.(P[ph1]-P[ph0])} x kernel(P) x colour x rho-product:
  coef  : complex number (all common factors 1/(2pi)^3 g^4 (1/k+) log stripped)
  phase : (v_plus, v_minus) -> e^{-i k.(P[v_plus]-P[v_minus])}; normally ("wp","w")
  kern  : list of kernel factors:
            ('K', a,b,c,d)   : K(a-b).K(c-d) = (a-b).(c-d)/((a-b)^2 (c-d)^2)
            ('inv2', a,b)    : 1/(a-b)^2
            ('fun', f)       : f(P) arbitrary scalar function of the position dict
  color : list of colour factors:
            ('f', i,j,k)     : f^{ijk}
            ('U', i,j, v)    : U^{ij}(P[v])
            ('delta', i,j)
          repeated index letters are contracted; indices on rhos are free.
  rhos  : ordered list [(idx, var), ...]  (operator product, left to right)
  Nc_pow: power of explicit N_c factor.

Evaluation: for a random configuration (sites with 2D coords + adjoint U matrix per site),
sum over all assignments of the term's variables to distinct sites ("injective sum"),
accumulating into a canonical (PBW) coefficient dictionary keyed by the multiset of
(site, colour-index-value) pairs of the rho product (rhos at distinct sites commute).

Reordering of rhos into canonical variable order is done SYMBOLICALLY first via
  [rho^a(x), rho^b(y)] = i f^{abc} delta^2(x-y) rho^c(x),
producing lower-rho contact terms with merged variables.
"""
import itertools, math, cmath
import numpy as np

# ---------------- group theory ----------------
def su2():
    d = 3
    f = np.zeros((d, d, d))
    for i, j, k in itertools.permutations(range(3)):
        f[i, j, k] = np.sign(np.linalg.det(np.eye(3)[[i, j, k]]))
    return d, f

def gellmann():
    l = np.zeros((8, 3, 3), dtype=complex)
    l[0][0, 1] = l[0][1, 0] = 1
    l[1][0, 1] = -1j; l[1][1, 0] = 1j
    l[2][0, 0] = 1; l[2][1, 1] = -1
    l[3][0, 2] = l[3][2, 0] = 1
    l[4][0, 2] = -1j; l[4][2, 0] = 1j
    l[5][1, 2] = l[5][2, 1] = 1
    l[6][1, 2] = -1j; l[6][2, 1] = 1j
    l[7][0, 0] = l[7][1, 1] = 1/math.sqrt(3); l[7][2, 2] = -2/math.sqrt(3)
    return l / 2

def su3():
    T = gellmann()
    d = 8
    f = np.zeros((d, d, d))
    for a in range(d):
        for b in range(d):
            comm = T[a] @ T[b] - T[b] @ T[a]
            for c in range(d):
                f[a, b, c] = np.real(-2j * np.trace(comm @ T[c]))
    return d, f, T

def adjU_su2(rng):
    # random SO(3)
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q

def adjU_su3(rng, T):
    A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    Q = Q @ np.diag(np.diag(R) / np.abs(np.diag(R)))
    Q = Q / np.linalg.det(Q) ** (1 / 3)
    d = 8
    U = np.zeros((d, d))
    for a in range(d):
        for b in range(d):
            U[a, b] = np.real(2 * np.trace(T[a] @ Q @ T[b] @ Q.conj().T))
    return U

# ---------------- term ----------------
class Term:
    __slots__ = ("coef", "phase", "kern", "color", "rhos", "Nc_pow", "tag")

    def __init__(self, coef, kern, color, rhos, phase=("wp", "w"), Nc_pow=0, tag=""):
        self.coef = complex(coef)
        self.phase = phase
        self.kern = kern
        self.color = color
        self.rhos = list(rhos)
        self.Nc_pow = Nc_pow
        self.tag = tag

    def vars(self):
        vs = set(self.phase)
        for kf in self.kern:
            if kf[0] == 'K':
                vs.update(kf[1:5])
            elif kf[0] == 'inv2':
                vs.update(kf[1:3])
            elif kf[0] == 'fun':
                vs.update(kf[2])  # ('fun', f, varlist)
        for cf in self.color:
            if cf[0] == 'U':
                vs.add(cf[3])
        for (_, v) in self.rhos:
            vs.add(v)
        return sorted(vs)

    def hc(self):
        """Hermitian conjugate: conjugate coefficient, swap phase, reverse rho order."""
        return Term(self.coef.conjugate(), self.kern, self.color,
                    list(reversed(self.rhos)), phase=(self.phase[1], self.phase[0]),
                    Nc_pow=self.Nc_pow, tag=self.tag + "_hc")

CANON = ["xp", "yp", "up", "x", "y", "z", "u", "v", "w", "wp"]  # rho ordering priority

def _rho_rank(var):
    return CANON.index(var) if var in CANON else 99

def substitute_var(term, old, new):
    def sub(v):
        return new if v == old else v
    kern = []
    for kf in term.kern:
        if kf[0] == 'K':
            kern.append(('K', sub(kf[1]), sub(kf[2]), sub(kf[3]), sub(kf[4])))
        elif kf[0] == 'inv2':
            kern.append(('inv2', sub(kf[1]), sub(kf[2])))
        else:
            kern.append((kf[0], kf[1], [sub(v) for v in kf[2]]))
    color = []
    for cf in term.color:
        if cf[0] == 'U':
            color.append(('U', cf[1], cf[2], sub(cf[3])))
        else:
            color.append(cf)
    rhos = [(i, sub(v)) for (i, v) in term.rhos]
    return Term(term.coef, kern, color, rhos,
                phase=(sub(term.phase[0]), sub(term.phase[1])),
                Nc_pow=term.Nc_pow, tag=term.tag)

_dummy_counter = [0]
def _newidx():
    _dummy_counter[0] += 1
    return f"D{_dummy_counter[0]}"

def canonical_order(term):
    """Reorder rhos into CANON order symbolically.
    Returns list of terms (original reordered + contact terms)."""
    out = []
    work = [term]
    while work:
        t = work.pop()
        rhos = t.rhos
        # bubble sort by variable rank
        swapped = False
        for i in range(len(rhos) - 1):
            (a, va), (b, vb) = rhos[i], rhos[i + 1]
            if _rho_rank(va) > _rho_rank(vb):
                # swap: rho^a(va) rho^b(vb) = rho^b(vb) rho^a(va) + [rho^a(va), rho^b(vb)]
                # commutator = i f^{abc} delta^2(va-vb) rho^c(va)
                new_rhos = rhos[:i] + [(b, vb), (a, va)] + rhos[i + 2:]
                t2 = Term(t.coef, t.kern, t.color, new_rhos, t.phase, t.Nc_pow, t.tag)
                work.append(t2)
                c = _newidx()
                contact_rhos = rhos[:i] + [(c, va)] + rhos[i + 2:]
                tc = Term(1j * t.coef, t.kern, t.color + [('f', a, b, c)],
                          contact_rhos, t.phase, t.Nc_pow, t.tag + "_ct")
                # merge vb -> va everywhere (delta^2(va-vb))
                tc = substitute_var(tc, vb, va)
                work.append(tc)
                swapped = True
                break
        if not swapped:
            out.append(t)
    return out

# ---------------- configuration ----------------
class Config:
    def __init__(self, nsites, group, seed, kvec=None):
        rng = np.random.default_rng(seed)
        self.nsites = nsites
        self.pos = rng.uniform(-1, 1, size=(nsites, 2))
        if group == "su2":
            self.d, self.f = su2()
            self.U = np.stack([adjU_su2(rng) for _ in range(nsites)])
            self.Nc = 2.0
        else:
            self.d, self.f, T = su3()
            self.U = np.stack([adjU_su3(rng, T) for _ in range(nsites)])
            self.Nc = 3.0
        self.kvec = np.array(kvec if kvec is not None else rng.uniform(-1, 1, 2))

# ---------------- evaluation ----------------
def kernel_value(term, P):
    val = 1.0
    for kf in term.kern:
        if kf[0] == 'K':
            a = P[kf[1]] - P[kf[2]]
            b = P[kf[3]] - P[kf[4]]
            val *= float(np.dot(a, b) / (np.dot(a, a) * np.dot(b, b)))
        elif kf[0] == 'inv2':
            a = P[kf[1]] - P[kf[2]]
            val *= 1.0 / float(np.dot(a, a))
        else:
            val *= kf[1](P)
    return val

def color_tensor(term, cfg, siteof):
    """Return tensor over the free rho indices (in order of term.rhos)."""
    # build einsum
    letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    idxmap = {}
    def L(ix):
        if ix not in idxmap:
            idxmap[ix] = letters[len(idxmap)]
        return idxmap[ix]
    ops, subs = [], []
    for cf in term.color:
        if cf[0] == 'f':
            ops.append(cfg.f); subs.append(L(cf[1]) + L(cf[2]) + L(cf[3]))
        elif cf[0] == 'U':
            ops.append(cfg.U[siteof[cf[3]]]); subs.append(L(cf[1]) + L(cf[2]))
        elif cf[0] == 'delta':
            ops.append(np.eye(cfg.d)); subs.append(L(cf[1]) + L(cf[2]))
    out = "".join(L(ix) for (ix, _) in term.rhos)
    if not ops:
        # pure rho structure with possible repeated indices among rhos
        expr = ",".join([ "" ]) # handled below
    einstr = ",".join(subs) + "->" + out if ops else "->" + out
    if ops:
        return np.einsum(einstr, *ops, optimize=True)
    else:
        # e.g. rho^a(x) rho^a(x') : indices repeat among rhos -> use identity tensors
        # build via delta network
        uniq = {}
        for (ix, _) in term.rhos:
            uniq.setdefault(ix, []).append(True)
        # represent: product over distinct indices of delta chains
        # simplest: tensor of ones with constraints
        dims = tuple(cfg.d for _ in term.rhos)
        t = np.zeros(dims)
        for multi in itertools.product(range(cfg.d), repeat=len(term.rhos)):
            vals = {}
            ok = True
            for (ix, _), v in zip(term.rhos, multi):
                if ix in vals and vals[ix] != v:
                    ok = False
                    break
                vals[ix] = v
            if ok:
                t[multi] = 1.0
        return t

def accumulate(terms, cfg, sector_filter=None):
    """Sum over injective assignments; returns dict: sector(nrho) -> dense complex
    array of shape (M*d,)*nrho, SYMMETRIZED over the rho slots (PBW at distinct
    sites; same-site ordering already canonical)."""
    M, d = cfg.nsites, cfg.d
    total = {}
    for t0 in terms:
        for t in canonical_order(t0):
            vs = t.vars()
            n = len(vs)
            if n > M:
                raise RuntimeError(f"not enough sites for {t.tag}: needs {n}")
            nr = len(t.rhos)
            if sector_filter is not None and nr not in sector_filter:
                continue
            if nr not in total:
                total[nr] = np.zeros((M * d,) * nr, dtype=complex)
            acc = total[nr]
            for perm in itertools.permutations(range(M), n):
                siteof = dict(zip(vs, perm))
                P = {v: cfg.pos[s] for v, s in siteof.items()}
                ph = cmath.exp(-1j * float(np.dot(cfg.kvec, P[t.phase[0]] - P[t.phase[1]])))
                sc = t.coef * (cfg.Nc ** t.Nc_pow) * ph * kernel_value(t, P)
                if sc == 0:
                    continue
                T = sc * color_tensor(t, cfg, siteof)
                rsites = [siteof[v] for (_, v) in t.rhos]
                sl = tuple(slice(s * d, (s + 1) * d) for s in rsites)
                acc[sl] += T
    # symmetrize each sector over rho slots
    out = {}
    for nr, acc in total.items():
        sym = np.zeros_like(acc)
        for p in itertools.permutations(range(nr)):
            sym += np.transpose(acc, p)
        out[nr] = sym / math.factorial(nr)
    return out

def compare(totA, totB, label, tol=1e-9):
    ok = True
    allsec = set(totA) | set(totB)
    for sec in sorted(allsec):
        A = totA.get(sec)
        B = totB.get(sec)
        if A is None:
            A = np.zeros_like(B)
        if B is None:
            B = np.zeros_like(A)
        scale = max(np.abs(A).max(), np.abs(B).max(), 1e-30)
        worst = np.abs(A - B).max()
        status = "PASS" if worst / scale < tol else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"  [{label}] sector {sec}rho: max|diff|/scale = {worst/scale:.3e}  ({status}; scale={scale:.3e})")
    return ok

def total_norm(tot):
    return {sec: float(np.abs(a).max()) for sec, a in tot.items()}
