"""ref3_beta / b1: light-cone vertices g->gg and g->qqbar from one Lagrangian.

Checks
 [1] the bracket of the paper's triple-gluon matrix element (app. hamiltonian, eq. m3)
     equals the covariant 3-gluon vertex contracted with light-cone polarisation vectors
     (eps^+=0, eps^- = eps.k/k^+), i.e. the paper's vertex has the standard normalisation
     g f^{abc} Gamma^{mu nu rho} / ((2pi)^{3/2} sqrt(2k+ 2p+ 2q+));
 [2] with the quark vertex g t^a ubar eps-slash v in the same normalisation, the ratio of the
     squared vertices (angular average, parent polarisation fixed) is
        n_f T_R sum|V_q|^2 / (1/2 N_c sum|V_g|^2) = (n_f T_R/N_c) [xi^2+(1-xi)^2] / K_gg(xi),
     K_gg = xi(1-xi) - 2 + 1/(xi(1-xi)) (the kernel of the paper's self-energy, app. B);
 [3] boost invariance: ubar eps-slash v depends only on the relative momentum (Galilean),
     and the eps^- piece is Q^l S/Q^+ with S = ubar gamma^+ v (quark) and S_g=(q^+-p^+) delta (gluon);
 [4] the cancellation of the instantaneous emission (m4) against the eps^- part of the
     emission rho -> g -> gg (m1, m3), with the paper's own matrix elements;
 [5] the spin/polarisation tensors T1^{lj,im}, T2^{im} entering the Row I.V kernel (used in b3).
Run: python3 -I b1_vertices.py
"""
import numpy as np
rng = np.random.default_rng(7)

# ---------------- Dirac algebra (Dirac representation), light-cone components -------------
I2 = np.eye(2); Z2 = np.zeros((2, 2))
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
g0 = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
def gk(s): return np.block([[Z2, s], [-s, Z2]])
g1, g2, g3 = gk(sx), gk(sy), gk(sz)
gp = (g0 + g3) / np.sqrt(2); gm = (g0 - g3) / np.sqrt(2)
gperp = [g1, g2]
def slash(a):  # a = dict with 'p' (+), 'm' (-), 't' (2-vector, upper comps)
    # a.gamma = a^+ gamma^- + a^- gamma^+ - a_perp.gamma_perp
    return a['p'] * gm + a['m'] * gp - a['t'][0] * g1 - a['t'][1] * g2
def dot(a, b):
    return a['p'] * b['m'] + a['m'] * b['p'] - a['t'] @ b['t']
def onshell(pp, pt):
    pt = np.asarray(pt, float)
    return {'p': pp, 'm': pt @ pt / (2 * pp), 't': pt}
# good-component basis: eigenvectors of Lambda_+ = gamma^- gamma^+ /2 with eigenvalue 1, helicity basis
Lp = gm @ gp / 2
w, V = np.linalg.eig(Lp)
good = V[:, np.isclose(w, 1)]
# diagonalise Sigma^3 inside the good subspace
S3 = 1j * g1 @ g2
h, U = np.linalg.eig(good.conj().T @ S3 @ good)
chi = good @ U
chi = chi / np.linalg.norm(chi, axis=0)
def spinor(p, k):
    # u_h(p) = pslash gamma^+ chi_h / norm ; massless: pslash u = 0
    s = slash(p) @ gp @ chi[:, k]
    return s
def bar(s): return s.conj() @ g0
# normalise: sum_h u ubar = pslash
p = onshell(0.7, [0.3, -0.2])
M = sum(np.outer(spinor(p, k), bar(spinor(p, k))) for k in range(2))
c = np.trace(M @ np.linalg.pinv(slash(p))) / 4  # crude
nrm = None
# determine normalisation constant from u-bar gamma^+ u = 2 p^+
val = bar(spinor(p, 0)) @ gp @ spinor(p, 0)
NRM = np.sqrt((2 * p['p']) / val.real)
def u(p, k): return NRM / np.sqrt(p['p']) * np.sqrt(p['p']) * spinor(p, k) / np.sqrt(p['p']) * np.sqrt(0.7) if False else spinor(p, k) * np.sqrt(2 * p['p'] / (bar(spinor(p, k)) @ gp @ spinor(p, k)).real)
v = u  # massless: same spinor space (helicity labels swapped), spin sums identical
M = sum(np.outer(u(p, k), bar(u(p, k))) for k in range(2))
print('[0] sum_h u ubar = pslash :', np.abs(M - slash(p)).max())

def eps(k, i):
    e = np.zeros(2); e[i] = 1.0
    return {'p': 0.0, 'm': e @ k['t'] / k['p'], 't': e}

# ---------------- [1] triple-gluon vertex vs m3 -------------------------------------------
def cov3(k1, e1, k2, e2, k3, e3):
    # Gamma^{mu nu rho} e1_mu e2_nu e3_rho, all momenta incoming
    def sub(a, b): return {'p': a['p'] - b['p'], 'm': a['m'] - b['m'], 't': a['t'] - b['t']}
    return dot(e1, e2) * dot(sub(k1, k2), e3) + dot(e2, e3) * dot(sub(k2, k3), e1) + dot(e3, e1) * dot(sub(k3, k1), e2)
def neg(a): return {'p': -a['p'], 'm': -a['m'], 't': -a['t']}
def m3_bracket(k, p, q, i, j, l):
    kt, pt, qt = k['t'], p['t'], q['t']
    d = np.eye(2)
    return ((pt[i] - qt[i] + (q['p'] - p['p']) / k['p'] * kt[i]) * d[j, l]
            + (kt[j] + qt[j] - (k['p'] + q['p']) / p['p'] * pt[j]) * d[i, l]
            + ((k['p'] + p['p']) / q['p'] * qt[l] - pt[l] - kt[l]) * d[i, j])
mx = 0; ratios = []
for trial in range(5):
    kp = 1.3; xi = rng.uniform(0.05, 0.95)
    Q = rng.normal(size=2); pr = rng.normal(size=2)
    k = onshell(kp, Q); pm = onshell(xi * kp, xi * Q + pr); qm = onshell((1 - xi) * kp, (1 - xi) * Q - pr)
    for i in range(2):
        for j in range(2):
            for l in range(2):
                a = cov3(k, eps(k, i), neg(pm), eps(pm, j), neg(qm), eps(qm, l))
                b = m3_bracket(k, pm, qm, i, j, l)
                if abs(b) > 1e-8: ratios.append(a / b)
                mx = max(mx, abs(a + b) if abs(a) > 0 else 0)
ratios = np.array(ratios)
print('[1] covariant 3g vertex / m3 bracket: min %.15f max %.15f' % (ratios.min(), ratios.max()))

# ---------------- [2] ratio of squared vertices ---------------------------------------------
def Vq(k, p1, p2, i, h1, h2):
    return bar(u(p1, h1)) @ slash(eps(k, i)) @ v(p2, h2)
def avg_sq(xi, quark, nphi=64, Q=np.array([0.0, 0.0])):
    kp = 1.0; s = 0.0
    for ph in np.linspace(0, 2 * np.pi, nphi, endpoint=False):
        pr = np.array([np.cos(ph), np.sin(ph)])
        k = onshell(kp, Q); p1 = onshell(xi * kp, xi * Q + pr); p2 = onshell((1 - xi) * kp, (1 - xi) * Q - pr)
        i = 0
        if quark:
            s += sum(abs(Vq(k, p1, p2, i, a, b)) ** 2 for a in range(2) for b in range(2))
        else:
            s += sum(abs(m3_bracket(k, p1, p2, i, j, l)) ** 2 for j in range(2) for l in range(2))
    return s / nphi
for xi in [0.1, 0.3, 0.5, 0.77]:
    Kgg = xi * (1 - xi) - 2 + 1 / (xi * (1 - xi))
    r = avg_sq(xi, True) / (0.5 * avg_sq(xi, False))       # per (T_R n_f)/(N_c)
    pred = (xi ** 2 + (1 - xi) ** 2) / Kgg
    rQ = avg_sq(xi, True, Q=np.array([0.4, -1.1])) / (0.5 * avg_sq(xi, False, Q=np.array([0.4, -1.1])))
    print('[2] xi=%.2f  ratio=%.15f  pred=%.15f  (Q!=0: %.15f)' % (xi, r, pred, rQ))

# ---------------- [3] boost invariance and the eps^- piece ----------------------------------
mx = 0; mxS = 0
for trial in range(4):
    xi = rng.uniform(0.1, 0.9); kp = 0.8; pr = rng.normal(size=2); Q = rng.normal(size=2)
    k0 = onshell(kp, [0, 0]); k = onshell(kp, Q)
    for i in range(2):
        for a in range(2):
            for b in range(2):
                A0 = Vq(k0, onshell(xi * kp, pr), onshell((1 - xi) * kp, -pr), i, a, b)
                A1 = Vq(k, onshell(xi * kp, xi * Q + pr), onshell((1 - xi) * kp, (1 - xi) * Q - pr), i, a, b)
                mx = max(mx, abs(A0 - A1))
                # transverse part at Q: -ubar gamma^i v ; eps^- part: (Q^i/Q^+) ubar gamma^+ v
                p1 = onshell(xi * kp, xi * Q + pr); p2 = onshell((1 - xi) * kp, (1 - xi) * Q - pr)
                T = -(bar(u(p1, a)) @ gperp[i] @ v(p2, b)) * (-1)  # gamma_perp with lower index sign handled below
                E = Q[i] / kp * (bar(u(p1, a)) @ gp @ v(p2, b))
                Tt = -(bar(u(p1, a)) @ gperp[i] @ v(p2, b))
                mxS = max(mxS, abs(A1 - (Tt + E)))
print('[3] quark vertex: |V(Q)-V(0)| max = %.2e ; |V - (transverse + Q^i S/Q^+)| max = %.2e' % (mx, mxS))
# ubar gamma^+ v = 2 sqrt(p1+ p2+) delta_{h,-hbar} (up to phases)
p1 = onshell(0.3, [0.2, 0.5]); p2 = onshell(0.5, [-0.1, 0.4])
Sm = np.array([[bar(u(p1, a)) @ gp @ v(p2, b) for b in range(2)] for a in range(2)])
print('    |ubar gamma^+ v| / (2 sqrt(p1+ p2+)) =\n', np.round(np.abs(Sm) / (2 * np.sqrt(0.15)), 12))

# ---------------- [4] instantaneous vs eps^- part with the paper's matrix elements ----------
# |Psi_gg^2> = 1/2 sum <gg|H_ggg|g><g|H_g|0>/(E_gg E_g), |Psi_gg^3> = -1/2 <gg|H_inst|0>/E_gg
# common factors g^2 f^{abc} rho^a stripped; colour structure identical (f^{abc} rho^a, b<->p, c<->q)
for trial in range(3):
    kp = 1.0; xi = rng.uniform(0.1, 0.9); pp, qp = xi * kp, (1 - xi) * kp
    Q = rng.normal(size=2); pr = rng.normal(size=2)
    pt, qt = xi * Q + pr, (1 - xi) * Q - pr
    Egg = pt @ pt / (2 * pp) + qt @ qt / (2 * qp); Eg = Q @ Q / (2 * kp)
    j = l = 0
    # eps^- part of m3: (q+ - p+) k^i/k^+ delta_jl ; m1 = g k^i rho/(4 pi^{3/2} k+^{3/2}); m3 prefactor i g f/(8 pi^{3/2} sqrt(k p q))
    em = sum((Q[i] / (4 * np.pi ** 1.5 * kp ** 1.5)) * (1j / (8 * np.pi ** 1.5 * np.sqrt(kp * pp * qp))) * (qp - pp) * Q[i] / kp for i in range(2))
    em = 0.5 * em / (Egg * Eg)
    inst = -0.5 * (1j * (qp - pp) / (2 * (2 * np.pi) ** 3 * np.sqrt(pp * qp) * (pp + qp) ** 2)) / Egg
    print('[4] eps^- part of emission = %.12e%+.12ej ; instantaneous = %.12e%+.12ej ; sum = %.1e' % (em.real, em.imag, inst.real, inst.imag, abs(em + inst)))

# ---------------- [5] tensors for the Row I.V kernel -----------------------------------------
# vertex at Q=0 linear in the relative momentum: V^l = ptilde^j v^{lj}; S/Q^+ = s
def quark_tensors(xi):
    kp = 1.0; k0 = onshell(kp, [0, 0])
    vlj = np.zeros((2, 2, 2, 2), complex)  # l, j, h1, h2
    for j in range(2):
        pr = np.zeros(2); pr[j] = 1.0
        for l in range(2):
            for a in range(2):
                for b in range(2):
                    vlj[l, j, a, b] = Vq(k0, onshell(xi * kp, pr), onshell((1 - xi) * kp, -pr), l, a, b)
    # linearity check
    pr = np.array([0.37, -1.3]); errl = 0
    for l in range(2):
        for a in range(2):
            for b in range(2):
                errl = max(errl, abs(Vq(k0, onshell(xi * kp, pr), onshell((1 - xi) * kp, -pr), l, a, b) - pr @ vlj[l, :, a, b]))
    # S (ubar gamma^+ v) depends on ptilde?  check
    s = np.zeros((2, 2), complex); errs = 0
    for a in range(2):
        for b in range(2):
            s[a, b] = bar(u(onshell(xi * kp, [1, 0]), a)) @ gp @ v(onshell((1 - xi) * kp, [-1, 0]), b) / kp
            errs = max(errs, abs(s[a, b] - bar(u(onshell(xi * kp, pr), a)) @ gp @ v(onshell((1 - xi) * kp, -pr), b) / kp))
    T1 = np.einsum('ljab,imab->ljim', vlj, vlj.conj())
    T2 = np.einsum('ab,imab->im', s, vlj.conj())
    return T1, T2, errl, errs
def gluon_tensors(xi, dperp=2):
    d = np.eye(2)
    v = np.zeros((2, 2, 2, 2))  # l, j(mom), a, b (daughter pols)
    for l in range(2):
        for j in range(2):
            for a in range(2):
                for b in range(2):
                    v[l, j, a, b] = 2 * (d[l, j] * d[a, b] - d[a, j] * d[l, b] / xi - d[b, j] * d[l, a] / (1 - xi))
    s = (1 - 2 * xi) * d
    T1 = np.einsum('ljab,imab->ljim', v, v)
    T2 = np.einsum('ab,imab->im', s, v)
    return T1, T2
if __name__ == '__main__':
    for xi in [0.2, 0.6]:
        T1, T2, el, es = quark_tensors(xi)
        A = xi * (1 - xi)
        # expected structure: T1 = c1 d^{lj}d^{im} + c2 d^{li}d^{jm} + c3 d^{lm}d^{ji} + c4 eps eps ...
        d = np.eye(2); e = np.array([[0, 1], [-1, 0]])
        basis = [np.einsum('lj,im->ljim', d, d), np.einsum('li,jm->ljim', d, d), np.einsum('lm,ji->ljim', d, d),
                 np.einsum('lj,im->ljim', e, e), np.einsum('li,jm->ljim', e, e), np.einsum('lm,ji->ljim', e, e)]
        B = np.array([b.ravel() for b in basis]).T
        coef, res, *_ = np.linalg.lstsq(B, T1.ravel(), rcond=None)
        print('[5] quark xi=%.2f linear err %.1e, S const err %.1e' % (xi, el, es))
        print('    T1 = sum c_n basis_n, c(dd,dd,dd,ee,ee,ee) * A =', np.round(coef * A, 12), ' resid', np.abs(B @ coef - T1.ravel()).max())
        print('    T2 =', np.round(T2, 12).tolist(), ' (expected prop. to (1-2xi) delta or eps)')
        # contractions used later
        print('    T1 contracted d^{lj}... :  T1^{lj,im} d^{li} d^{jm} * A =', np.round(np.einsum('ljim,li,jm->', T1, d, d).real * A, 12),
              '   (2D trace check: 2*(xi^2+(1-xi)^2)*2 =', 4 * (xi ** 2 + (1 - xi) ** 2), ')')
        T1g, T2g = gluon_tensors(xi)
        coefg, *_ = np.linalg.lstsq(B, T1g.ravel(), rcond=None)
        print('    gluon: T1 coef =', np.round(coefg, 12), '; T2 =', np.round(T2g, 12).tolist())
