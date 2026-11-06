"""Common utilities for numeric verification of p+ integration steps."""
import mpmath as mp
import random

mp.mp.dps = 30

def rvec(rng, scale=1.0):
    return (mp.mpf(rng.uniform(-scale, scale)), mp.mpf(rng.uniform(-scale, scale)))

def sub(u, v):
    return (u[0]-v[0], u[1]-v[1])

def dot(u, v):
    return u[0]*v[0] + u[1]*v[1]

def sq(u):
    return dot(u, u)

def quad(f, a, b):
    return mp.quad(f, [a, b])

def quad_pv(f, a, b, poles):
    """Principal value integral of f over [a,b] with simple poles inside."""
    poles = sorted([p for p in poles if a < p < b])
    if not poles:
        return mp.quad(f, [a, b])
    total = mp.mpf(0)
    pts = [a] + poles + [b]
    for i, p in enumerate(poles):
        lo, hi = pts[i], pts[i+2] if False else None
    # do it pole by pole: split interval into segments owned by each pole
    bounds = [a]
    for i in range(len(poles)-1):
        bounds.append((poles[i]+poles[i+1])/2)
    bounds.append(b)
    for p, lo, hi in zip(poles, bounds[:-1], bounds[1:]):
        d = min(p-lo, hi-p)
        # symmetric PV around pole over [p-d, p+d]
        def g(t, p=p):
            return f(p+t) + f(p-t)
        total += mp.quad(g, [mp.mpf('1e-12')*d, d])  # symmetric part; f(p+t)+f(p-t) is regular
        if p-d > lo:
            total += mp.quad(f, [lo, p-d])
        if hi > p+d:
            total += mp.quad(f, [p+d, hi])
    return total

def report(name, lhs, rhs, tol=mp.mpf('1e-6')):
    lhs, rhs = mp.mpf(lhs) if not isinstance(lhs, mp.mpc) else lhs, rhs
    diff = abs(lhs-rhs)
    scale = max(abs(lhs), abs(rhs), mp.mpf('1e-30'))
    ok = diff/scale < tol
    print(f"{'PASS' if ok else 'FAIL'}  {name}: lhs={mp.nstr(lhs,10)} rhs={mp.nstr(rhs,10)} reldiff={mp.nstr(diff/scale,3)}")
    return ok
