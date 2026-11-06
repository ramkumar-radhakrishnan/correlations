"""Check: BFKL adjoint final (BFKL_corrected lines 606-633) vs its colour-singlet
trace form (lines 635-644), and vs the user's target LO-BFKL equation.

Method: substitute <rho^b(x') rho^c(x)> -> delta^{bc} mu2(x',x)/(Nc^2-1) (singlet
average, with mu2 symmetric); evaluate both forms summed over injective assignments
of the 5 coordinates to random points (kernels explicit, traces with random SU(N)
fundamental V's); compare.
Adjoint variables (x',x) map to singlet (y... we use mu2 symmetric so either naming).
"""
import itertools, math, cmath
import numpy as np
from texparse import parse_block
from blocks import block_body, FB, PI

def su_n_random(N, rng):
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    Q, R = np.linalg.qr(A)
    Q = Q @ np.diag(np.diag(R) / np.abs(np.diag(R)))
    return Q / np.linalg.det(Q) ** (1 / N)

def gens(N):
    """generalized Gell-Mann, normalized Tr(T^aT^b) = delta/2"""
    Ts = []
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), dtype=complex); m[i, j] = m[j, i] = 0.5
            Ts.append(m)
            m = np.zeros((N, N), dtype=complex); m[i, j] = -0.5j; m[j, i] = 0.5j
            Ts.append(m)
    for k in range(1, N):
        m = np.zeros((N, N), dtype=complex)
        for i in range(k):
            m[i, i] = 1
        m[k, k] = -k
        m = m / math.sqrt(2 * k * (k + 1))
        Ts.append(m)
    return Ts

class Cfg:
    def __init__(self, N, nsites, seed):
        rng = np.random.default_rng(seed)
        self.N = N
        self.T = gens(N)
        self.d = N * N - 1
        self.pos = rng.uniform(-1, 1, size=(nsites, 2))
        self.V = [su_n_random(N, rng) for _ in range(nsites)]
        # adjoint
        self.Uadj = np.zeros((nsites, self.d, self.d))
        for s in range(nsites):
            for a in range(self.d):
                for b in range(self.d):
                    self.Uadj[s, a, b] = np.real(2 * np.trace(self.T[a] @ self.V[s] @ self.T[b] @ self.V[s].conj().T))
        f = np.zeros((self.d,) * 3)
        for a in range(self.d):
            for b in range(self.d):
                comm = self.T[a] @ self.T[b] - self.T[b] @ self.T[a]
                for c in range(self.d):
                    f[a, b, c] = np.real(-2j * np.trace(comm @ self.T[c]))
        self.f = f
        self.kvec = rng.uniform(-1, 1, 2)
        rng2 = np.random.default_rng(seed + 1)
        self.mu2 = rng2.uniform(0.2, 1.0, size=(nsites, nsites))
        self.mu2 = 0.5 * (self.mu2 + self.mu2.T)
        self.nsites = nsites

def Kdot(a, b):
    return float(np.dot(a, b) / (np.dot(a, a) * np.dot(b, b)))

def Ksc(x, y, z):
    return float(np.dot(x - z, y - z) / (np.dot(x - z, x - z) * np.dot(y - z, y - z)))

def eval_adjoint(terms, cfg):
    """terms from texparse (2-rho); rho^b(p1) rho^c(p2) -> delta^{bc} mu2(p1,p2)/(d)"""
    tot = 0
    varsets = {}
    for t in terms:
        vs = tuple(t.vars())
        varsets.setdefault(vs, []).append(t)
    for vs, ts in varsets.items():
        for perm in itertools.permutations(range(cfg.nsites), len(vs)):
            siteof = dict(zip(vs, perm))
            P = {v: cfg.pos[s] for v, s in siteof.items()}
            for t in ts:
                ph = cmath.exp(-1j * float(np.dot(cfg.kvec, P[t.phase[0]] - P[t.phase[1]])))
                kv = 1.0
                for kf in t.kern:
                    if kf[0] == 'K':
                        kv *= Kdot(P[kf[1]] - P[kf[2]], P[kf[3]] - P[kf[4]])
                    elif kf[0] == 'inv2':
                        r = P[kf[1]] - P[kf[2]]
                        kv *= 1.0 / float(np.dot(r, r))
                # colour: contract with rho->delta mu2/d
                (i1, p1), (i2, p2) = t.rhos
                mu = cfg.mu2[siteof[p1], siteof[p2]] / cfg.d
                # build colour sum with i1 == i2 summed
                letters = "abcdefghijklmnop"
                idxmap = {}
                def L(ix):
                    if ix not in idxmap:
                        idxmap[ix] = letters[len(idxmap)]
                    return idxmap[ix]
                ops, subs = [], []
                for cfac in t.color:
                    if cfac[0] == 'f':
                        ops.append(cfg.f); subs.append(L(cfac[1]) + L(cfac[2]) + L(cfac[3]))
                    elif cfac[0] == 'U':
                        ops.append(cfg.Uadj[siteof[cfac[3]]]); subs.append(L(cfac[1]) + L(cfac[2]))
                    elif cfac[0] == 'delta':
                        ops.append(np.eye(cfg.d)); subs.append(L(cfac[1]) + L(cfac[2]))
                # contract rho indices together
                ops.append(np.eye(cfg.d)); subs.append(L(i1) + L(i2))
                cval = np.einsum(",".join(subs) + "->", *ops, optimize=True)
                Nc = cfg.N
                tot += t.coef * (Nc ** t.Nc_pow) * ph * kv * mu * cval
    return tot

def eval_singlet_bfkl(cfg):
    """user's singlet form, lines 636-643 (unit prefactor g^4/16pi^5/k+ L_B/(Nc^2-1) stripped
    consistently: we keep 1/(Nc^2-1) and compare absolute vs adjoint/16pi^5 etc."""
    N = cfg.N
    tot = 0
    vs = ["x", "y", "z", "w", "wp"]
    for perm in itertools.permutations(range(cfg.nsites), 5):
        siteof = dict(zip(vs, perm))
        P = {v: cfg.pos[s] for v, s in siteof.items()}
        V = {v: cfg.V[siteof[v]] for v in vs}
        ph = cmath.exp(-1j * float(np.dot(cfg.kvec, P["wp"] - P["w"])))
        mu = cfg.mu2[siteof["x"], siteof["y"]]
        x, y, z, w, wp = (P[v] for v in vs)
        Ud = lambda v: V[v].conj().T
        def tr(m):
            return np.trace(m)
        t1 = Kdot(x - w, y - wp) * Ksc(x, y, z) * 2 * N * tr((Ud("x") - Ud("w")) @ (V["y"] - V["wp"]))
        t2 = Kdot(x - w, y - wp) * Ksc(x, x, z) * (-N) * tr((Ud("x") - Ud("w")) @ (V["y"] - V["wp"]))
        t3 = Kdot(x - w, y - wp) * Ksc(y, y, z) * (-N) * tr((Ud("x") - Ud("w")) @ (V["y"] - V["wp"]))
        t4 = Ksc(w, wp, z) * Ksc(x, y, z) * 2 * N * tr((Ud("z") - Ud("w")) @ (V["z"] - V["wp"]))
        t5 = Kdot(x - wp, w - z) * Ksc(x, y, z) * 2 * N * tr((Ud("z") - Ud("w")) @ (V["x"] - V["wp"]))
        t6 = Kdot(x - wp, w - z) * Ksc(y, y, z) * (-N) * tr((Ud("z") - Ud("w")) @ (V["x"] - V["wp"]))
        t7 = Kdot(x - w, wp - z) * Ksc(x, y, z) * 2 * N * tr((Ud("x") - Ud("w")) @ (V["z"] - V["wp"]))
        t8 = Kdot(x - w, wp - z) * Ksc(y, y, z) * (-N) * tr((Ud("x") - Ud("w")) @ (V["z"] - V["wp"]))
        tot += ph * mu / (N * N - 1) * (t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8)
    return tot / (16 * PI**5)

if __name__ == "__main__":
    body = block_body(FB, 606, 633, marker="\\Bigg[")
    terms = parse_block(body, extra_coef=1 / (16 * PI**5), tag="BFKLfin")
    print("parsed adjoint-final terms:", len(terms))
    for N, seed in ((2, 21), (3, 22)):
        cfg = Cfg(N, 5, seed)
        A = eval_adjoint(terms, cfg)
        S = eval_singlet_bfkl(cfg)
        print(f"N={N}: adjoint={A:.6e}  singlet={S:.6e}  ratio={abs(A/S) if S else float('nan'):.6f}  "
              f"reldiff={abs(A-S)/max(abs(A),abs(S)):.3e}")
