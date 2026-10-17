"""Transcriptions of the displayed equations of Group_I.tex (integrands at fixed transverse points).

Each function returns the integrand of the row at fixed transverse points, WITH the prefactor 1/(2pi)^3 and the
phase e^{-ik(w'-w)} (or what replaces it after the w integration), and WITH the factor 1/(2pi) of dp+/2pi when the
p+ integral has been done (functions named *_int).  Before the p+ integral the functions return the integrand of
int dp+/2pi.  Charges are classical: {rho^a(x), rho^b(y)} = 2 rho^a(x) rho^b(y) at x != y.
Index placements are as written in the notes (they use U^{bc} = U^{cb}); colour sums are done literally."""
import numpy as np
from fp1 import F, NC, G, PI, I2, K, sq, PI3, phase

# ------------------------------------------------------------------ common pieces


def P(c):
    return c.pt


def conjcol(c):
    """[U^{ab'}(x') - U^{ab'}(w')] rho^{b'}(x')  ->  [a]"""
    return (c.U['xp'] - c.U['wp']) @ c.rho['xp']


def col3(c):
    """f^{abd} [U^{ab'}(x')rho^{b'}(x') - U^{ab'}(w')rho^{b'}(x')] [U^{de}(x)U^{bc}(z) rho^e(x)rho^c(y)
    - U^{de}(x)U^{bc}(y) rho^e(x)rho^c(y)]   (Rows III)"""
    U, r = c.U, c.rho
    L = conjcol(c)
    R = np.einsum('de,e,bc,c->db', U['x'], r['x'], U['z'] - U['y'], r['y'])
    return np.einsum('abd,a,db->', F, L, R)


def col4(c):
    """[U(x')-U(w')]^{ab'} rho^{b'}(x') [f^{cdc'} U^{ac}(w) U^{bd}(z) U^{be}(y) rho^e(y) rho^{c'}(x)
    - f^{abc} U^{be}(y) U^{cd}(x) rho^e(y) rho^d(x)]   (Row IV)"""
    U, r = c.U, c.rho
    L = conjcol(c)
    t1 = np.einsum('cdC,ac,bd,be,e,C->a', F, U['w'], U['z'], U['y'], r['y'], r['x'])
    t2 = np.einsum('abc,be,e,cd,d->a', F, U['y'], r['y'], U['x'], r['x'])
    return L @ (t1 - t2)


# ------------------------------------------------------------------ Row II, part 1 (one rho)
def r2_part1(c, kp, bracket):
    """line 46: 1/(2pi)^3 g^4 Nc/(32 pi^4) 1/k+ K(x'-w').K(x-w) [bracket] [U(x')-U(w')]rho(x') [U(w)-U(x)]rho(x)"""
    p = P(c)
    kern = K(p['xp'] - p['wp']) @ K(p['x'] - p['w'])
    R = (c.U['w'] - c.U['x']) @ c.rho['x']
    return phase(c) / PI3 * G**4 * NC / (32 * PI**4) / kp * kern * bracket * (conjcol(c) @ R)


# ------------------------------------------------------------------ Row III (lines 87-90): integrand of int dp+/2pi
def r3_kin(c):
    p = P(c)
    x, y, z, w, xp, wp = p['x'], p['y'], p['z'], p['w'], p['xp'], p['wp']
    return dict(x=x, y=y, z=z, w=w, xp=xp, wp=wp, A=sq(x - z), B=sq(x - w), KK1=K(xp - wp) @ K(y - z))


def r3_lines(c, kp, pp, Lam, strip=False):
    """the five lines of Eq. (90) as written: returns dict line -> value (0 outside its p+ range).
    strip=True multiplies by (p+ A - k+ B) (to remove the pole for the principal value)."""
    q = r3_kin(c)
    A_, B_ = q['A'], q['B']
    den = pp * A_ - kp * B_
    m = den if strip else 1.0
    col = col3(c)
    ph = phase(c) / PI3
    out = {}
    gt, lt = (kp + Lam < pp), (Lam < pp < kp - Lam)
    inst = 1j * G**4 / (8 * PI**4) * q['KK1'] * (pp + kp) / (kp - pp)**2 / den * m * col * ph
    out['L1'] = inst if gt else 0.0
    out['L2'] = inst if lt else 0.0
    kern3 = (K(q['xp'] - q['wp']) @ K(q['x'] - q['w'])) * (K(q['x'] - q['z']) @ K(q['y'] - q['z']))
    out['L3'] = (-1j * G**4 / (8 * PI**4) * kern3 * (kp * B_ + pp * A_) / den * m / (pp * kp) * col * ph) if (Lam < pp) else 0.0
    xw, xz, zw, yz = q['x'] - q['w'], q['x'] - q['z'], q['z'] - q['w'], q['y'] - q['z']
    Vv = kp * xw - pp * xz
    T = (np.einsum('ik,jm->ijkm', I2, I2) - kp / pp * np.einsum('jk,im->ijkm', I2, I2)
         - kp / (kp - pp) * np.einsum('ij,km->ijkm', I2, I2))
    if gt:
        t = np.einsum('i,j,k,m,ijkm->', K(q['xp'] - q['wp']), K(yz), K(zw), Vv, T)
        out['L4'] = 1j * G**4 / (4 * PI**4) / kp * t / (kp * B_ - pp * A_) * m / (kp - pp) * col * ph
    else:
        out['L4'] = 0.0
    if lt:
        t = np.einsum('i,j,k,m,ijkm->', K(q['xp'] - q['wp']), yz / sq(yz), zw, Vv, T) / sq(kp * xw - pp * xz)
        out['L5'] = -1j * G**4 / (4 * PI**4) * t / (kp * B_ - pp * A_) * m * pp / (kp - pp) * col * ph
    else:
        out['L5'] = 0.0
    return out


# ------------------------------------------------------------------ Row III after the p+ integral (Eqs. 94, 97, 100, 120)
def r3_int(c, kp, Lam, vee, absl=True):
    """every line of N1, N2, N3a, N3b as written; returns dict name -> value (includes 1/(2pi)^3 e^{..})"""
    q = r3_kin(c)
    A_, B_ = q['A'], q['B']
    x, y, z, w, xp, wp = q['x'], q['y'], q['z'], q['w'], q['xp'], q['wp']
    col = col3(c) * phase(c) / PI3
    c8, c16 = 1j * G**4 / (8 * PI**5), 1j * G**4 / (16 * PI**5)
    LvL = np.log(vee / Lam)
    KK1 = q['KK1']
    d = lambda a, b: a @ b
    xw, xz, zw, yz, xpw = x - w, x - z, z - w, y - z, xp - wp
    n2 = sq(xpw) * sq(yz)                       # (x'-w')^2 (y-z)^2
    o = {}
    # principal value: log of the modulus; absl=False takes the complex principal branch instead
    lg = (lambda v: np.log(abs(v))) if absl else (lambda v: np.log(complex(v)))
    # N1 (instantaneous p+ > k+, instantaneous p+ < k+, the 't4' entry)
    o['N1.1'] = c8 * KK1 / (A_ - B_) * (1 / Lam + 1 / (kp - vee)) * col
    o['N1.2'] = -c16 / kp * (LvL + np.log(1 - kp / vee)) * KK1 * (A_ + B_) / (A_ - B_)**2 * col
    o['N1.3'] = c16 / kp * KK1 * (B_ + A_) / (A_ - B_)**2 * lg((vee * A_ - kp * B_) / ((Lam + kp) * A_ - kp * B_)) * col
    o['N1.4'] = c8 * KK1 / (A_ - B_) * (1 / Lam - 1 / (kp - Lam)) * col
    o['N1.5'] = -c16 / kp * (np.log(vee / (kp - Lam)) - LvL) * KK1 * (A_ + B_) / (A_ - B_)**2 * col
    o['N1.6'] = c16 / kp * KK1 * (B_ + A_) / (A_ - B_)**2 * lg(((kp - Lam) * A_ - kp * B_) / (Lam * A_ - kp * B_)) * col
    k3 = d(xpw, xw) * d(xz, yz) / (sq(xpw) * B_ * A_ * sq(yz))
    o['N1.7'] = c16 / kp * LvL * k3 * col
    o['N1.8'] = -c8 / kp * k3 * lg((kp * B_ - vee * A_) / (kp * B_ - Lam * A_)) * col
    # N2 (non-instantaneous p+ > k+)
    Lv = LvL + np.log(1 - kp / vee)
    Lg = lg((B_ * kp - A_ * vee) / (B_ * kp - A_ * (kp + Lam)))
    o['N2.1'] = c8 * d(xpw, yz) / (n2 * (B_ - A_)) * (1 / (vee - kp) - 1 / Lam) * col
    o['N2.2'] = -c8 / kp * np.log(vee / (Lam + kp)) * d(xpw, xw) * d(yz, zw) / (n2 * sq(zw) * B_) * col
    o['N2.3red'] = -c8 / kp * Lv * d(xpw, zw) * d(yz, zw) / (n2 * sq(zw) * (B_ - A_)) * col
    o['N2.4red'] = +c8 / kp * Lv * d(xpw, zw) * d(yz, zw) / (n2 * sq(zw) * (B_ - A_)) * col
    o['N2.5'] = c8 / kp * Lv * d(xpw, yz) * d(zw, xz) * B_ / (n2 * sq(zw) * (B_ - A_)**2) * col
    o['N2.6'] = -c8 / kp * Lv * d(xpw, yz) * d(zw, xw) * A_ / (n2 * sq(zw) * (B_ - A_)**2) * col
    o['N2.7'] = c8 / kp * Lg * d(xpw, zw) * d(yz, xw) / (n2 * sq(zw) * (B_ - A_)) * col
    o['N2.8'] = -c8 / kp * Lg * d(xpw, xw) * d(yz, zw) * A_ / (n2 * sq(zw) * (B_ - A_) * B_) * col
    o['N2.9'] = c8 / kp * Lg * d(xpw, yz) * d(zw, xw) * A_ / (n2 * sq(zw) * (B_ - A_)**2) * col
    o['N2.10'] = -c8 / kp * Lg * d(xpw, zw) * d(yz, xz) * B_ / (n2 * sq(zw) * (B_ - A_) * A_) * col
    o['N2.11'] = c8 / kp * Lg * d(xpw, xz) * d(yz, zw) / (n2 * sq(zw) * (B_ - A_)) * col
    o['N2.12'] = -c8 / kp * Lg * d(xpw, yz) * d(zw, xz) * B_ / (n2 * sq(zw) * (B_ - A_)**2) * col
    # N3a (non-instantaneous p+ < k+)
    Lk = LvL + np.log((kp - Lam) / vee)
    Lq = lg((B_ * kp - A_ * (kp - Lam)) / (B_ * kp - A_ * Lam))
    Dl = B_ * A_ - d(xw, xz)**2
    Ac = (np.arctan(((kp - Lam) * A_ - kp * d(xw, xz)) / (kp * np.sqrt(Dl)))
          - np.arctan((Lam * A_ - kp * d(xw, xz)) / (kp * np.sqrt(Dl)))) / np.sqrt(Dl)
    o['N3a.1'] = c8 * (1 / Lam - 1 / (kp - Lam)) * d(xpw, yz) / (n2 * (B_ - A_)) * col
    o['N3a.2'] = -c8 / kp * Lk * d(xpw, yz) / (n2 * sq(zw)) * col
    o['N3a.3'] = -c8 / kp * Lk * d(xpw, yz) * d(zw, xw) * A_ / (n2 * sq(zw) * (B_ - A_)**2) * col
    o['N3a.4'] = c8 / kp * Lk * d(xpw, yz) * d(zw, xz) * B_ / (n2 * sq(zw) * (B_ - A_)**2) * col
    o['N3a.5'] = -c8 / kp * Lq * d(xpw, yz) * d(zw, xw) * A_ / (n2 * sq(zw) * (B_ - A_)**2) * col
    o['N3a.6'] = c8 / kp * Lq * d(xpw, yz) * d(zw, xz) * B_ / (n2 * sq(zw) * (B_ - A_)**2) * col
    o['N3a.7'] = -c8 / kp * Lq * d(xpw, zw) * d(yz, xw) / (n2 * sq(zw) * (B_ - A_)) * col
    o['N3a.8'] = c8 / kp * Lq * d(xpw, zw) * d(yz, xz) * B_ / (n2 * sq(zw) * (B_ - A_) * A_) * col
    o['N3a.9'] = c8 / kp * Lq * d(yz, zw) * d(xpw, xw) * A_ / (n2 * sq(zw) * (B_ - A_) * B_) * col
    o['N3a.10'] = -c8 / kp * Lq * d(yz, zw) * d(xpw, xz) / (n2 * sq(zw) * (B_ - A_)) * col
    o['N3a.11'] = c8 / kp * Ac * d(xpw, zw) * d(yz, xw) / (n2 * sq(zw)) * col
    o['N3a.12'] = c8 / kp * Ac * d(yz, zw) * d(xpw, xz) / (n2 * sq(zw)) * col
    o['N3a.13'] = -c8 / kp * Ac * d(xw, xz) * d(xpw, zw) * d(yz, xz) / (n2 * sq(zw) * A_) * col
    # N3b
    Lr = np.log(sq(kp * xw - (kp - Lam) * xz) / sq(kp * xw - Lam * xz))
    o['N3b.1'] = -c8 / kp * Ac * d(xw, xz) * d(yz, zw) * d(xpw, xw) / (n2 * sq(zw) * B_) * col
    o['N3b.2'] = -c16 / kp * Lr * d(xpw, zw) * d(yz, xz) / (n2 * sq(zw) * A_) * col
    o['N3b.3'] = c16 / kp * Lr * d(yz, zw) * d(xpw, xw) / (n2 * sq(zw) * B_) * col
    o['N3b.4'] = -c16 / kp * Lr * d(xpw, yz) / (n2 * sq(zw)) * col
    return o


# ------------------------------------------------------------------ Row IV (lines 134-142)
def r4_line(c, kp, pp):
    """Eq. (137): integrand of int_{Lambda}^{vee-k} dp+/2pi"""
    p = P(c)
    x, y, z, w, xp, wp = p['x'], p['y'], p['z'], p['w'], p['xp'], p['wp']
    xz, xw, zw = x - z, x - w, z - w
    D = pp * sq(xz) + kp * sq(xw)
    br = (I2 / (2 * sq(zw) * (pp + kp)) * (sq(xz) - sq(xw))
          + 1 / kp * (np.outer(zw, xz) / sq(zw) - np.outer(xw, xz) / (2 * sq(xw)))      # [i, j]: (x-z)^j (z-w)^i ...
          + 1 / pp * (np.outer(xw, zw) / sq(zw) + np.outer(xw, xz) / (2 * sq(xz))))
    t = np.einsum('i,j,ij->', K(xp - wp), K(y - z), br) / D
    return phase(c) / PI3 * 1j * G**4 / (4 * PI**4) * t * col4(c)


def r4_int(c, kp, Lam, vee):
    """Eq. (142): after the p+ integral (includes 1/(2pi)^3 e^{..})"""
    p = P(c)
    x, y, z, w, xp, wp = p['x'], p['y'], p['z'], p['w'], p['xp'], p['wp']
    xz, xw, zw = x - z, x - w, z - w
    A_, B_ = sq(xz), sq(xw)
    KK = np.outer(K(xp - wp), K(y - z))                       # [i, j]
    L1 = np.log(((vee - kp) * A_ + kp * B_) / (Lam * A_ + kp * B_))
    br1 = (np.outer(zw, xz) / (sq(zw) * A_) - np.outer(xw, zw) / (sq(zw) * B_) - np.outer(xw, xz) / (A_ * B_)
           + I2 / (2 * sq(zw)))
    t1 = 1j * G**4 / (8 * PI**5) / kp * np.sum(KK * br1) * L1
    t2 = -1j * G**4 / (16 * PI**5) / kp * (K(xp - wp) @ (y - z)) / (sq(y - z) * sq(zw)) * np.log(vee / (Lam + kp))
    br3 = np.outer(xw, zw) / (sq(zw) * B_) + np.outer(xw, xz) / (2 * A_ * B_)
    t3 = 1j * G**4 / (8 * PI**5) / kp * np.sum(KK * br3) * (np.log(vee / Lam) + np.log(1 - kp / vee))
    return phase(c) / PI3 * (t1 + t2 + t3) * col4(c), dict(t1=t1, t2=t2, t3=t3)


# ------------------------------------------------------------------ Row V, part 1 (lines 151-159)
def col5a(c):
    """[U(x')-U(w')]^{ab'}rho^{b'}(x') [f^{adc} U^{db}(x) U^{ce}(z) {rho^e(y), rho^b(y')}
    - f^{abc} U^{cd}(y) U^{be}(y') {rho^d(y), rho^e(y')}]"""
    U, r = c.U, c.rho
    t1 = 2 * np.einsum('adc,db,ce,e,b->a', F, U['x'], U['z'], r['y'], r['yp'])
    t2 = 2 * np.einsum('abc,cd,be,d,e->a', F, U['y'], U['yp'], r['y'], r['yp'])
    return conjcol(c) @ (t1 - t2)


def TCm(kp, pp):
    return (np.einsum('im,jk->ijkm', I2, I2) - kp / pp * np.einsum('jm,ik->ijkm', I2, I2)
            - kp / (kp - pp) * np.einsum('ij,km->ijkm', I2, I2))


def r5a_line(c, kp, pp):
    """Eq. (151): coefficient of delta^2(w - z + (p/k)(z-x)) in the integrand of int dp+/2pi (w placed by the caller)"""
    p = P(c)
    t = np.einsum('i,m,k,j,ijkm->', K(p['xp'] - p['wp']), K(p['z'] - p['x']), K(p['y'] - p['z']), K(p['yp'] - p['x']),
                  TCm(kp, pp))
    return phase(c) / PI3 * 1j * G**4 / (16 * PI**4) / kp**2 * t * col5a(c)


# ------------------------------------------------------------------ Row V, part 2 (lines 163-171), as in Group_V.tex
def brackets6(x, y, z):
    yx, yz, xz = y - x, y - z, x - z
    B1 = I2 * 2 / (2 * sq(xz)) * (sq(yx) - sq(yz))
    B2 = np.outer(yx, xz) / sq(xz) - np.outer(yx, yz) / (2 * sq(yz))
    B3 = I2 * (yz @ xz / sq(xz) + yx @ yz / (2 * sq(yx)) - (sq(yx) - sq(yz)) / (2 * sq(xz)))
    B4 = np.outer(yz, xz) / sq(xz) + np.outer(yz, yx) / (2 * sq(yx))
    B5 = I2 * (yx @ xz / sq(xz) - yx @ yz / (2 * sq(yz)) - (sq(yx) - sq(yz)) / (2 * sq(xz)))
    B6 = (np.outer(xz, yx) / sq(xz) - np.outer(yz, yx) / (2 * sq(yz)) + np.outer(xz, yz) / sq(xz)
          + np.outer(yx, yz) / (2 * sq(yx)))
    return [B1, B2, B3, B4, B5, B6]


def r5b_line(c, kp, pp):
    """Eq. (163) with d_perp = 2: coefficient of delta^2 in the integrand of int dp+/2pi"""
    p = P(c)
    x, y, z = p['x'], p['y'], p['z']
    Bs = brackets6(x, y, z)
    q = kp - pp
    br = (Bs[0] - kp**2 / q**2 * Bs[1] + kp / pp * Bs[2] - kp**2 / pp**2 * Bs[3] + kp / q * Bs[4] - kp**2 / (pp * q) * Bs[5])
    D = pp * sq(y - x) + q * sq(y - z)
    t = np.einsum('i,m,im->', K(p['xp'] - p['wp']), K(z - x), br) * pp * q / D
    U, r = c.U, c.rho
    colR = (np.einsum('adc,ebf,db,ce,f->a', F, F, U['x'], U['z'], r['y'])
            - np.einsum('abc,cbd,df,f->a', F, F, U['y'], r['y']))
    return phase(c) / PI3 * G**4 / (8 * PI**4) / kp**3 * t * (conjcol(c) @ colR)


# ------------------------------------------------------------------ Row V, part 1 after the w integral and the plus prescription
def r5a_pieces(c, kp):
    """phase phi(p) = e^{-ik(w'-z)} e^{-ik(z-x)p/k}, the three contractions T1 (d_im d_jk), T2 (d_jm d_ik),
    T3 (d_ij d_km), and the common factor C0 = 1/(2pi)^3 i g^4/(16 pi^4) x colour"""
    p = P(c)
    kv = c.kv
    phi = lambda pp: np.exp(-1j * kv @ (p['wp'] - p['z'])) * np.exp(-1j * kv @ (p['z'] - p['x']) * pp / kp)
    a, m, kk, j = K(p['xp'] - p['wp']), K(p['z'] - p['x']), K(p['y'] - p['z']), K(p['yp'] - p['x'])
    T1, T2, T3 = (a @ m) * (j @ kk), (j @ m) * (a @ kk), (a @ j) * (m @ kk)
    C0 = 1 / PI3 * 1j * G**4 / (16 * PI**4) * col5a(c)
    return phi, T1, T2, T3, C0


def r5a_plus(c, kp, Lam, lo, hi, quad):
    """Eq. (159) with the plus-prescribed integrals taken over (lo, hi): (0, k+) as written, (Lambda, k+-Lambda) exact.
    The endpoint terms are as written: ln((k+-Lambda)/Lambda) = log(vee/Lambda) + log((k+-Lambda)/vee)."""
    phi, T1, T2, T3, C0 = r5a_pieces(c, kp)
    reg1 = quad(lambda pp: phi(pp) * T1 - kp / pp * (phi(pp) - phi(0)) * T2, lo, hi) / (2 * PI) / kp**2
    reg2 = -quad(lambda q: (phi(kp - q) - phi(kp)) / q * T3, lo, hi) / (2 * PI) / kp
    L = np.log((kp - Lam) / Lam)
    end = -1 / (2 * PI) / kp * L * (phi(0) * T2 + phi(kp) * T3)
    return C0 * (reg1 + reg2 + end)


def r5a_exact(c, kp, Lam, quad):
    """the integral of Eq. (155) over Lambda < p+ < k+ - Lambda"""
    phi, T1, T2, T3, C0 = r5a_pieces(c, kp)
    f = lambda pp: phi(pp) * (T1 - kp / pp * T2 - kp / (kp - pp) * T3)
    return C0 * quad(f, Lam, kp - Lam) / (2 * PI) / kp**2


# ------------------------------------------------------------------ four-rho parts of Rows III and IV (not in Group_I.tex)
def r3_4rho(c, kp, pp):
    """(1/(2pi)^3) g^4/(8 pi^4 p+ k+) K(x'-w').K(x-w) K(y'-z).K(y-z) [U(x')-U(w')]^{ab'}rho^{b'}(x')
    {rho-bar^a(x), rho-bar^b(y')} [U(z)-U(y)]^{bc} rho^c(y)   (integrand of int_Lambda^vee dp+/2pi)"""
    p, U, r = P(c), c.U, c.rho
    kern = (K(p['xp'] - p['wp']) @ K(p['x'] - p['w'])) * (K(p['yp'] - p['z']) @ K(p['y'] - p['z']))
    col = conjcol(c) @ (2 * np.outer(U['x'] @ r['x'], U['yp'] @ r['yp'])) @ ((U['z'] - U['y']) @ r['y'])
    return phase(c) / PI3 * G**4 / (8 * PI**4 * pp * kp) * kern * col


def r4_4rho(c, kp, pp):
    """-(1/(2pi)^3) g^4/(8 pi^4 p+ k+) K(x'-w').K(y'-w) K(x-z).K(y-z) [U(x')-U(w')]^{ab'}rho^{b'}(x')
    [U^{ac}(w) U^{bd}(z) U^{be}(y)rho^e(y) {rho^c(y'), rho^d(x)} - U^{be}(y)rho^e(y) {rho-bar^a(y'), rho-bar^b(x)}]
    (integrand of int_Lambda^{vee-k+} dp+/2pi)"""
    p, U, r = P(c), c.U, c.rho
    kern = (K(p['xp'] - p['wp']) @ K(p['yp'] - p['w'])) * (K(p['x'] - p['z']) @ K(p['y'] - p['z']))
    Ry = U['y'] @ r['y']
    t1 = 2 * np.einsum('ac,bd,b,c,d->a', U['w'], U['z'], Ry, r['yp'], r['x'])
    t2 = 2 * np.einsum('a,b,b->a', U['yp'] @ r['yp'], U['x'] @ r['x'], Ry)
    return -phase(c) / PI3 * G**4 / (8 * PI**4 * pp * kp) * kern * (conjcol(c) @ (t1 - t2))


# ------------------------------------------------------------------ Row II, part 2 (line 53), integrand of int d^2p d^2k
def r2_part2(c, kp, kv, pv_, Bvec):
    """-(1/(2pi)^3) g^4/(64 pi^6) 1/k+ e^{ik(w-y)} e^{ip(y-x)} (x'-w')^i/((x'-w')^2 p^2 k^2) B^i
    [U(x')-U(w')]^{ab'}rho^{b'}(x') [U^{ah}(w) f^{hbc}{rho^b(x),rho^c(y)} - f^{abc} U^{bd}(x)U^{ce}(y){rho^d(x),rho^e(y)}]"""
    U, r, pt = c.U, c.rho, c.pt
    ph = np.exp(1j * kv @ (pt['w'] - pt['y'])) * np.exp(1j * pv_ @ (pt['y'] - pt['x']))
    kern = (K(pt['xp'] - pt['wp']) @ Bvec) / ((pv_ @ pv_) * (kv @ kv))
    R = (2 * np.einsum('ah,hbc,b,c->a', U['w'], F, r['x'], r['y'])
         - 2 * np.einsum('abc,bd,ce,d,e->a', F, U['x'], U['y'], r['x'], r['y']))
    return -phase(c) / PI3 * G**4 / (64 * PI**6) / kp * ph * kern * (conjcol(c) @ R)


# ------------------------------------------------------------------ Row V, part 2a: Eq. (1.24) and its plus prescription (1.26)
def r5_2a_parts(c, kp):
    p, U, r = P(c), c.U, c.rho
    x, y, z = p['x'], p['y'], p['z']
    Bs = brackets6(x, y, z)
    a = K(p['xp'] - p['wp'])
    m = K(z - x)
    S = [np.einsum('i,m,im->', a, m, B) for B in Bs]                   # K^i(x'-w') K^m(z-x) B_n^{im}
    phi = lambda pp: np.exp(-1j * c.kv @ (p['wp'] - z)) * np.exp(-1j * c.kv @ (z - x) * pp / kp)
    D = lambda pp: pp * sq(y - x) + (kp - pp) * sq(y - z)
    colR = (-np.einsum('adc,fbe,db,ce,f->a', F, F, U['x'], U['z'], r['y']) + NC * U['x'] @ r['y'])
    C0 = 1 / PI3 * G**4 / (8 * PI**4) / kp * (conjcol(c) @ colR)
    return S, phi, D, C0


def r5_2a_exact(c, kp, Lam, quad):
    """Eq. (1.24) integrated over Lambda < p+ < k+ - Lambda"""
    S, phi, D, C0 = r5_2a_parts(c, kp)
    def f(pp):
        xi = pp / kp
        co = [xi * (1 - xi), -xi / (1 - xi), 1 - xi, -(1 - xi) / xi, xi, -1.0]
        return phi(pp) / D(pp) * sum(cc * s for cc, s in zip(co, S))
    return C0 * quad(f, Lam, kp - Lam) / (2 * PI)


def r5_2a_plus(c, kp, Lam, lo, hi, quad):
    """Eq. (1.26) with the plus-prescribed integrals over (lo, hi)"""
    S, phi, D, C0 = r5_2a_parts(c, kp)
    def fA(pp):
        xi = pp / kp
        reg = phi(pp) / D(pp) * (xi * (1 - xi) * S[0] + (1 - xi) * S[2] + S[3] + xi * S[4] - S[5])
        sub = -kp / pp * (phi(pp) / D(pp) - phi(0) / D(0)) * S[3]
        return reg + sub
    def fB(q):
        pp = kp - q
        return -kp / q * (phi(pp) / D(pp) - phi(kp) / D(kp)) * S[1] + phi(pp) / D(pp) * S[1]
    L = np.log((kp - Lam) / Lam)
    end = -kp * L * (phi(0) / D(0) * S[3] + phi(kp) / D(kp) * S[1])
    return C0 * (quad(fA, lo, hi) + quad(fB, lo, hi) + end) / (2 * PI)
