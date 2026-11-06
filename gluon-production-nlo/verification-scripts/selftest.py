"""Engine self-tests."""
import numpy as np
from engine import *

def test_group(group):
    cfg = Config(4, group, seed=5)
    # adjoint invariance: f^{abc} U^{a a'} U^{b b'} U^{c c'} = f^{a'b'c'}
    U = cfg.U[0]
    lhs = np.einsum('abc,ax,by,cz->xyz', cfg.f, U, U, U)
    print(group, "f invariance:", np.max(np.abs(lhs - cfg.f)))
    print(group, "orthogonality:", np.max(np.abs(U @ U.T - np.eye(cfg.d))))
    # f^{acd} f^{bcd} = Nc delta^{ab}
    ff = np.einsum('acd,bcd->ab', cfg.f, cfg.f)
    print(group, "casimir:", np.max(np.abs(ff - cfg.Nc * np.eye(cfg.d))))

def test_reorder():
    # rho^a(x) rho^b(y) == rho^b(y) rho^a(x) + i f^{abc} delta(x-y) rho^c(x)
    # engine check: canonical_order of T_left should equal manual right side.
    cfg = Config(3, "su2", seed=7)
    kern = []
    t_left = Term(1.0, kern, [('U', 'a', 'p', 'x'), ('U', 'b', 'q', 'y')],
                  [('p', 'y'), ('q', 'x')], phase=('x', 'y'), tag="L")
    # canonical order puts x before y -> left ordering rho(y)rho(x) must reorder.
    tot_left = accumulate([t_left], cfg)
    # manual: rho(x) rho(y) ordering + contact
    t_r1 = Term(1.0, kern, [('U', 'a', 'p', 'x'), ('U', 'b', 'q', 'y')],
                [('q', 'x'), ('p', 'y')], phase=('x', 'y'), tag="R1")
    t_r2 = Term(1j, kern, [('U', 'a', 'p', 'x'), ('U', 'b', 'q', 'x'), ('f', 'p', 'q', 'r')],
                [('r', 'y')], phase=('x', 'x'), tag="R2")
    t_r2 = substitute_var(t_r2, 'y', 'x')
    tot_right = accumulate([t_r1, t_r2], cfg)
    ok = compare(tot_left, tot_right, "reorder")
    print("reorder test:", "PASS" if ok else "FAIL")

def test_hc():
    # (i f^{abc} U^{bd}(x) rho^d(x) rho^c(y)) + h.c. should equal
    # i f^{abc} U^{bd}(x) [rho^d(x), rho^c(y)] + ... just check hc engine runs & real total
    cfg = Config(3, "su3", seed=9)
    t = Term(0.3 + 0.1j, [('K', 'x', 'y', 'x', 'y')],
             [('f', 'a', 'b', 'c'), ('U', 'b', 'd', 'x'), ('U', 'c', 'e', 'y'), ('delta', 'a', 'a2'), ('U', 'a2', 'g', 'x')],
             [('d', 'x'), ('e', 'y'), ('g', 'x')], phase=('x', 'y'), tag="t")
    tot = accumulate([t, t.hc()], cfg)
    # total operator must be Hermitian: coefficient dict keys with sorted pairs;
    # hermiticity at the dict level: C[key] for reversed products... just check sum of
    # imaginary parts of diagonal-like keys is small? Simpler: full matrix check skipped.
    print("hc ran, sectors:", {k: len(v) for k, v in tot.items()})

if __name__ == "__main__":
    test_group("su2")
    test_group("su3")
    test_reorder()
    test_hc()
