"""Plus prescription for the UV-finite part of Group I Row V, part 2b (d^3N_{2b finite}), exact on Lambda < p+ < k+ - Lambda.

In xi = p+/k+ (so (1/k+) int dp+ = int dxi), with the common factor N0 = (1/(2pi)^3)(g^4/8pi^4)(1/2pi) and int_{x,y,z,x',w'}
stripped, the integrand is I(xi) = I_first(xi) + I_CT(xi).  The note writes
    int_d^{1-d} I = int_d^{1-d} [I - I0/xi - I1/(1-xi)] + (I0 + I1) ln((1-d)/d),     d = Lambda/k+,
with the closed-form endpoint coefficients I0, I1 of the note.  Checked here, at random transverse points, with a smooth
model of the colour operator O(x) = [U(x')-U(w')] N_c [U(y)-U(x)] rho rho (a smooth complex function of x):
  1. xi I(xi) -> I0 at xi -> 0 and (1-xi) I(xi) -> I1 at xi -> 1 (closed forms);
  2. the note's plus-prescribed integrand, written term by term, equals I - I0/xi - I1/(1-xi);
  3. plus form + endpoint terms = the original integral (quadrature);
  4. UV: the angular average of r^2 x (each piece) -> 0 as r = |x - z| -> 0, for the regular remainder, the two
     plus-subtracted pieces, the two endpoint groups, and the combined endpoint form;
  5. the combined endpoint form (p+ -> k+ endpoint relabelled x <-> z) has the same x,z integral as I0 + I1."""
import json
import numpy as np
import mpmath as mp

rng = np.random.default_rng(12)
DPERP = 2.0


def dot(a, b):
    return float(a @ b)


class Pt:
    def __init__(self, seed):
        r = np.random.default_rng(seed)
        self.x, self.y, self.z, self.xp, self.wp = (r.normal(size=2) * 0.7 for _ in range(5))
        self.k = r.normal(size=2) * 1.3
        self.kp = 1.7
        c = r.normal(size=6)
        # smooth colour model O(p) for the charge-y Wilson line at point p; O_x = O(x), O_z = O(z)
        self.O = lambda p, c=c: complex(np.cos(c[0] * p[0] + c[1] * p[1]) + c[2], np.sin(c[3] * p[0] - c[4] * p[1]) + c[5])


def Kf(a, b, c):
    """K(a,b,c) = (a-c).(b-c)/((a-c)^2 (b-c)^2)"""
    return dot(a - c, b - c) / (dot(a - c, a - c) * dot(b - c, b - c))


def first(P, xi, x=None):
    x = P.x if x is None else x
    y, z, A, kp = P.y, P.z, P.xp - P.wp, P.kp
    ph = np.exp(-1j * dot(P.k, P.wp - z)) * np.exp(-1j * dot(P.k, z - x) * xi)
    D = kp * (xi * dot(y - x, y - x) + (1 - xi) * dot(y - z, y - z))
    xz, yx, yz = x - z, y - x, y - z
    V = (z - x) / dot(z - x, z - x)
    a = A / dot(A, A)
    I2 = np.eye(2)
    B1 = DPERP * I2 / (2 * dot(xz, xz)) * (dot(yx, yx) - dot(yz, yz))
    B2 = np.outer(yx, xz) / dot(xz, xz) - np.outer(yx, yz) / (2 * dot(yz, yz))
    C3 = dot(yz, xz) / dot(xz, xz) + dot(yx, yz) / (2 * dot(yx, yx)) - (dot(yx, yx) - dot(yz, yz)) / (2 * dot(xz, xz))
    B4 = np.outer(yz, xz) / dot(xz, xz) + np.outer(yz, yx) / (2 * dot(yx, yx))      # [i, m]: (y-z)^i (x-z)^m, (y-x)^m (y-z)^i
    C5 = dot(yx, xz) / dot(xz, xz) - dot(yx, yz) / (2 * dot(yz, yz)) - (dot(yx, yx) - dot(yz, yz)) / (2 * dot(xz, xz))
    B6 = (np.outer(xz, yx) / dot(xz, xz) - np.outer(yz, yx) / (2 * dot(yz, yz)) + np.outer(xz, yz) / dot(xz, xz)
          + np.outer(yx, yz) / (2 * dot(yx, yx)))                                     # [i, m]
    con = lambda M: float(a @ M @ V)
    return dict(ph=ph, D=D, B1=con(B1), B2=con(B2), C3=con(I2) * C3, B4=con(B4), C5=con(I2) * C5, B6=con(B6), O=P.O(x))


def I_full(P, xi):
    f = first(P, xi)
    M = (xi * (1 - xi) * f['B1'] - xi / (1 - xi) * f['B2'] + (1 - xi) * f['C3'] - (1 - xi) / xi * f['B4']
         + xi * f['C5'] - f['B6'])
    return f['ph'] / f['D'] * M * f['O'] + I_ct(P, xi)


def ct_kernel(P):
    A, y, z, x = P.xp - P.wp, P.y, P.z, P.x
    return dot(A, y - z) / (dot(A, A) * dot(y - z, y - z)) * dot(y - x, y - z) / (dot(y - x, y - x) * dot(x - z, x - z))


def Phat(xi):
    return xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi


def I_ct(P, xi):
    return -np.exp(-1j * dot(P.k, P.wp - P.z)) / P.kp * Phat(xi) * ct_kernel(P) * P.O(P.z)


def I0_closed(P, x=None):
    x = P.x if x is None else x
    A, y, z, kp = P.xp - P.wp, P.y, P.z, P.kp
    pre = np.exp(-1j * dot(P.k, P.wp - z)) / kp * dot(A, y - z) / (dot(A, A) * dot(y - z, y - z))
    return pre * ((Kf(z, z, x) - 0.5 * Kf(y, z, x)) * P.O(x) - (Kf(z, z, x) - Kf(y, z, x)) * P.O(z))


def I1_closed(P, x=None):
    x = P.x if x is None else x
    A, y, z, kp = P.xp - P.wp, P.y, P.z, P.kp
    t1 = np.exp(-1j * dot(P.k, P.wp - x)) / kp * dot(A, y - x) / (dot(A, A) * dot(y - x, y - x)) * (Kf(x, x, z) - 0.5 * Kf(y, x, z)) * P.O(x)
    t2 = np.exp(-1j * dot(P.k, P.wp - z)) / kp * dot(A, y - z) / (dot(A, A) * dot(y - z, y - z)) * (Kf(z, z, x) - Kf(y, z, x)) * P.O(z)
    return t1 - t2


def note_plus_integrand(P, xi):
    """the plus-prescribed integrand exactly as written in the note"""
    f = first(P, xi)
    reg = (xi * (1 - xi) * f['B1'] + f['B2'] + (1 - xi) * f['C3'] + f['B4'] + xi * f['C5'] - f['B6']) * f['ph'] / f['D'] * f['O']
    F2 = lambda s: (lambda g: g['ph'] / g['D'] * g['B2'] * g['O'])(first(P, s))
    F4 = lambda s: (lambda g: g['ph'] / g['D'] * g['B4'] * g['O'])(first(P, s))
    plus = -(F4(xi) - F4(0.0)) / xi - (F2(xi) - F2(1.0)) / (1 - xi)
    ct = -np.exp(-1j * dot(P.k, P.wp - P.z)) / P.kp * (xi * (1 - xi) - 2) * ct_kernel(P) * P.O(P.z)
    return reg + plus + ct


def cquad(fun, a, b):
    re = mp.quad(lambda t: fun(float(t)).real, [a, (a + b) / 2, b])
    im = mp.quad(lambda t: fun(float(t)).imag, [a, (a + b) / 2, b])
    return complex(re, im)


# UV: angular average of r^2 * piece at x = z + r (cos t, sin t)
def ang(P, piece, r, n=400):
    ts = (np.arange(n) + 0.5) * 2 * np.pi / n
    vals = []
    for t in ts:
        x = P.z + r * np.array([np.cos(t), np.sin(t)])
        vals.append(piece(x))
    return r * r * np.mean(vals)


def with_x(P, x):
    Q = Pt.__new__(Pt)
    Q.__dict__.update(P.__dict__)
    Q.x = x
    return Q


# 5. combined endpoint form vs I0 + I1, integrated over x and z on a grid (y, x', w' fixed), with a Gaussian cutoff
def combined(P, x, z):
    A, y, kp = P.xp - P.wp, P.y, P.kp
    Ox, Oz = P.O(x), P.O(z)
    pre = np.exp(-1j * dot(P.k, P.wp - z)) / kp * dot(A, y - z) / (dot(A, A) * dot(y - z, y - z))
    return pre * (Kf(z, z, x) * (Ox - Oz) + Kf(y, z, x) * (1.5 * Oz - 0.5 * Ox))


def endpoints(P, x, z):
    Q = with_x(P, x)
    Q.z = z
    return I0_closed(Q) + I1_closed(Q)



if __name__ == '__main__':
    OUT = {}
    for s in (1, 2, 3):
        P = Pt(s)
        d = 0.11 / P.kp
        L = np.log((1 - d) / d)
        I0, I1 = I0_closed(P), I1_closed(P)
        lim0 = 1e-8 * I_full(P, 1e-8)
        lim1 = 1e-8 * I_full(P, 1 - 1e-8)
        xs = [0.07, 0.4, 0.83]
        d_note = max(abs(note_plus_integrand(P, t) - (I_full(P, t) - I0 / t - I1 / (1 - t))) for t in xs)
        orig = cquad(lambda t: I_full(P, t), d, 1 - d)
        plus = cquad(lambda t: note_plus_integrand(P, t), d, 1 - d) + (I0 + I1) * L
        OUT['config %d' % s] = dict(endpoint0_closed_vs_limit=abs(I0 - lim0) / abs(I0), endpoint1_closed_vs_limit=abs(I1 - lim1) / abs(I1),
                                    note_integrand_vs_definition=d_note, original=[orig.real, orig.imag], plus_form=[plus.real, plus.imag],
                                    rel_diff=abs(orig - plus) / abs(orig))
        print('config %d: endpoint closed forms rel. %.1e %.1e; note integrand %.1e; original %s, plus form %s, rel. diff %.1e' % (
            s, OUT['config %d' % s]['endpoint0_closed_vs_limit'], OUT['config %d' % s]['endpoint1_closed_vs_limit'], d_note, np.round(orig, 10), np.round(plus, 10), abs(orig - plus) / abs(orig)))


    P = Pt(1)
    xi0 = 0.37
    pieces = {
        '2b as written (fixed xi)': lambda x: (lambda Q: I_full(Q, xi0) - I_ct(Q, xi0))(with_x(P, x)),
        'rest (fixed xi)': lambda x: I_full(with_x(P, x), xi0),
        'plus-prescribed integrand (fixed xi)': lambda x: note_plus_integrand(with_x(P, x), xi0),
        'endpoint p+ -> 0 (I0)': lambda x: I0_closed(P, x),
        'endpoint p+ -> k+ (I1)': lambda x: I1_closed(P, x),
    }
    UV = {}
    for name, f in pieces.items():
        vals = [abs(ang(P, f, r)) for r in (1e-2, 1e-3, 1e-4)]
        UV[name] = vals
        print('%-40s r^2 <.> at r = 1e-2, 1e-3, 1e-4: %s' % (name, ['%.2e' % v for v in vals]))
    OUT['UV r^2 angular average'] = UV


    grid = np.linspace(-2.0, 2.0, 41)
    h = grid[1] - grid[0]
    pts = [np.array([a, b]) + 0.013 for a in grid for b in grid]
    w = lambda p: np.exp(-dot(p, p))                       # symmetric weight in x and z
    s1 = s2 = 0
    for xx in pts:
        for zz in pts:
            if dot(xx - zz, xx - zz) < 1e-6 or dot(P.y - xx, P.y - xx) < 1e-6 or dot(P.y - zz, P.y - zz) < 1e-6:
                continue
            ww = w(xx) * w(zz)
            s1 += ww * endpoints(P, xx, zz)
            s2 += ww * combined(P, xx, zz)
    OUT['combined endpoint form vs I0+I1 (grid, symmetric weight)'] = dict(I0_plus_I1=[s1.real, s1.imag], combined=[s2.real, s2.imag], rel_diff=abs(s1 - s2) / abs(s1))
    print('combined endpoint form: I0+I1 %s, combined %s, rel. diff %.1e' % (np.round(s1, 8), np.round(s2, 8), abs(s1 - s2) / abs(s1)))
    json.dump(OUT, open('plus2b_checks.json', 'w'), indent=1, default=float)
