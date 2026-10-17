"""The entries of B_3 in momentum space that were not re-derived before: the instantaneous entry (m5) and the two-rho
entries (m1 twice).  Each is rebuilt from the matrix elements of Calculation_main.tex and the energy denominators, and
compared with B_3 / (2 pi)^3 (the same factor as for B_2 and the non-instantaneous entries).

Two-rho: the two time orderings give c1 rho^b(-p) rho^a(k) + c2 rho^a(k) rho^b(-p)
        = (c1 + c2)/2 {rho^a(k), rho^b(-p)} + (c1 - c2)/2 [rho^b(-p), rho^a(k)],
with [rho^b(-p), rho^a(k)] = i f^{bac} rho^c(k-p) = -i f^{abc} rho^c(k-p)."""
import numpy as np

PI = np.pi
g = 1.3


def Eg(pp, p):
    return p @ p / (2 * pp)


def m1(pp, p):
    """<g^a_i(p)|H_g|0> without rho^a(-p): vector [i]"""
    return g * p / (4 * PI**1.5 * pp**1.5)


def inst_wf(kp, k, pp, p):
    """instantaneous one-gluon wave function: coefficient of f^{abc} rho^c(k-p) delta_ij |g^b_j(p)>, from m5"""
    m5 = 1j * g**2 * (kp + pp) / (2 * (2 * PI)**3 * np.sqrt(kp * pp) * (kp - pp)**2)   # <g(p)|H|g(k)>, f^{cab}rho^c
    return m5 / (Eg(kp, k) - Eg(pp, p))


def B3_inst(kp, k, pp, p):
    """B_3, first line (coefficient of f^{abc} rho^c(-p+k) delta_ij)"""
    return 1j * g**2 * np.sqrt(pp * kp) * (kp + pp) / (pp - kp)**2 / (pp * (k @ k) - kp * (p @ p))


def tworho(kp, k, pp, p):
    """c1, c2 [i, j]: c1 from |g(k)> -> |0> -> |g(p)>, c2 from |g(k)> -> |g(k) g(p)> -> |g(p)>"""
    em = m1(pp, p)                  # emission of p (rho^b(-p))
    ab = m1(kp, k)                  # absorption of k (rho^a(k)), the conjugate matrix element
    c1 = np.outer(ab, em) / (Eg(kp, k) * (Eg(kp, k) - Eg(pp, p)))
    c2 = np.outer(ab, em) / ((Eg(kp, k) - Eg(kp, k) - Eg(pp, p)) * (Eg(kp, k) - Eg(pp, p)))
    return c1, c2


def B3_two(kp, k, pp, p):
    """B_3, second line: coefficients of {rho^a(k), rho^b(-p)} and of i f^{abc} rho^c(k-p)  [i, j]"""
    base = -g**2 * np.outer(k, p) / ((k @ k) * (p @ p) * np.sqrt(kp * pp))
    anti = base
    comm = -base * ((p @ p) * kp + (k @ k) * pp) / ((p @ p) * kp - (k @ k) * pp)
    return anti, comm


def run(seed=7):
    rng = np.random.default_rng(seed)
    out = []
    for t in range(4):
        k, p = rng.normal(size=2), rng.normal(size=2)
        kp = 1.1
        for pp in (0.3 * kp, 2.4 * kp):
            r_inst = B3_inst(kp, k, pp, p) / inst_wf(kp, k, pp, p)
            c1, c2 = tworho(kp, k, pp, p)
            anti, comm = B3_two(kp, k, pp, p)
            wf_anti = (c1 + c2) / 2
            wf_comm = -(c1 - c2) / 2          # coefficient of i f^{abc} rho^c(k-p)
            out.append(dict(pp_over_kp=pp / kp, inst=complex(r_inst),
                            anti=float(np.mean(anti / wf_anti)), anti_spread=float(np.ptp(anti / wf_anti)),
                            comm=float(np.mean(comm / wf_comm)), comm_spread=float(np.ptp(comm / wf_comm))))
    return out


if __name__ == '__main__':
    for r in run():
        print(r, ' (8 pi^3 = %.8f)' % (8 * PI**3))
