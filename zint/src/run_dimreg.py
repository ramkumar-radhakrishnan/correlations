"""Checks for the dim-reg / MSbar note.  -> ../RESULTS_dimreg.txt"""
import mpmath as mp
mp.mp.dps = 30
out=[]
def p(*a): out.append(" ".join(str(x) for x in a)); print(*a)
g = mp.euler

p("A. the d-dimensional transform, formula checked where the answer is known (d=3, a=1)")
d,a = 3,1
pref = mp.pi**(mp.mpf(d)/2)*2**(d-2*a)*mp.gamma((mp.mpf(d)-2*a)/2)/mp.gamma(a)
p(f"   pi^(d/2) 2^(d-2a) Gamma((d-2a)/2)/Gamma(a) = {mp.nstr(pref,14)}")
p(f"   known result  int d^3x e^(iqx)/x^2 = 2 pi^2/|q| ,  2 pi^2 = {mp.nstr(2*mp.pi**2,14)}")
p("")
p("B. the epsilon expansion:  pi Gamma(-eps) (q^2/(4 pi mu^2))^eps + pi/eps  ->  pi log(mubar^2/q^2)")
q2, mu2 = mp.mpf('2.7'), mp.mpf('1.3'); mubar2 = 4*mp.pi*mu2*mp.exp(-g)
p(f"   q^2 = {q2},  mu^2 = {mu2},  mubar^2 = 4 pi mu^2 e^-gamma = {mp.nstr(mubar2,12)}")
for e in ['1e-3','1e-4','1e-5','1e-6','1e-7']:
    eps = mp.mpf(e)
    v = mp.pi*mp.gamma(-eps)*(q2/(4*mp.pi*mu2))**eps + mp.pi/eps
    p(f"   eps={e:>5}:  Gamma-form + pi/eps = {mp.nstr(v,14):>16}   pi log(mubar^2/q^2) = {mp.nstr(mp.pi*mp.log(mubar2/q2),14)}")
p("")
p("C. the p+ moments   (lambda = 1e-9, l_k = log(1/lambda))")
lam = mp.mpf('1e-9'); lk = mp.log(1/lam)
C = lambda t: t*(1-t) + t/(1-t) + (1-t)/t
I0 = mp.quad(C,[lam,1-lam])
I1 = mp.quad(lambda t: C(t)*mp.log(1/(1-t)**2),[lam,1-lam])
Iodd = mp.quad(lambda t: C(t)*mp.log(t/(1-t)),[lam,1-lam])
Ic10 = mp.quad(lambda t: C(t)*mp.log(t*(1-t)),[lam,1-lam])
p(f"   int C_UV                = {mp.nstr(I0,14)}    2 l_k - 11/6            = {mp.nstr(2*lk-mp.mpf(11)/6,14)}")
p(f"   int C_UV log(1/xib^2)   = {mp.nstr(I1,14)}    l_k^2 + pi^2/3 - 67/18  = {mp.nstr(lk**2+mp.pi**2/3-mp.mpf(67)/18,14)}")
p(f"   int C_UV log(xi/xib)    = {mp.nstr(Iodd,8)}  (odd under xi <-> xibar, vanishes)")
p(f"   int C_UV log(xi xib)    = {mp.nstr(Ic10,14)}   [paper (C.10)];  sum with the row above: {mp.nstr(Ic10+I1,8)}")
open('/home/user/correlations/zint/RESULTS_dimreg.txt','w').write("\n".join(out)+"\n")
