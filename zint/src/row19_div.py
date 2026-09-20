"""Row 19: which of your six parts carries log(1/Lambda) and log(V), and how the
   would-be poles 1/(A-B), 1/(AB'-A'B) cancel WITHIN each part."""
import numpy as np
L=np.log; s=lambda v:v@v
from row19_lib import their_coef, exact_M
rng=np.random.default_rng(404); k=1.0
x,xp,z,w,wp=[rng.normal(size=2)*1.3 for _ in range(5)]
A,B,C,Ap,Bp,Cp=s(x-z),s(x-w),s(z-w),s(xp-z),s(xp-wp),s(z-wp)
print("="*80); print("1.  SENSITIVITY TO Lambda  (V fixed at 1e6)"); print("="*80)
print("   A=%.4f B=%.4f A'=%.4f B'=%.4f"%(A,B,Ap,Bp))
print("   %-8s %14s %14s %14s   %s"%("part","Lam=1e-4","Lam=1e-6","Lam=1e-8","d/dlog(1/Lam)  [predicted]"))
for nm,pred in (('23',0.0),('33',1/(k*B*Bp)),('22',0.0),('13',0.0),('12',0.0),('11_full',0.0)):
    v=[their_coef(A,B,C,Ap,Bp,Cp,k,l,1e6)[nm] for l in (1e-4,1e-6,1e-8)]
    sl=(v[2]-v[0])/L(1e4)
    if nm=='11_full': sl/= (2.0*(A-B)*(Ap-Bp)/(4*C*Cp))
    print("   %-8s %14.9f %14.9f %14.9f   %12.9f  [%.9f]"%(nm,*v,sl,pred))
print()
print("="*80); print("2.  SENSITIVITY TO V  (Lambda fixed at 1e-8)"); print("="*80)
print("   %-8s %14s %14s %14s   %s"%("part","V=1e5","V=1e7","V=1e9","d/dlogV  [predicted]"))
for nm,pred in (('23',0.0),('33',0.0),('22',1/(k*A*Ap)),('13',0.0),('12',0.0),('11_full',0.0)):
    v=[their_coef(A,B,C,Ap,Bp,Cp,k,1e-8,V)[nm] for V in (1e5,1e7,1e9)]
    sl=(v[2]-v[0])/L(1e4)
    if nm=='11_full': sl/= (2.0*(A-B)*(Ap-Bp)/(4*C*Cp))
    print("   %-8s %14.9f %14.9f %14.9f   %12.9f  [%.9f]"%(nm,*v,sl,pred))
print()
print("   => log(1/Lambda) lives ONLY in 2b, residue 1/[k+ (x-w)^2 (x'-w')^2].")
print("      log V        lives ONLY in 2c, residue 1/[k+ (x-z)^2 (x'-z)^2].")
print("      2d, 2e, 2f each contain explicit log V in their separate pieces, but the")
print("      coefficients cancel between the pieces: the PART is V-finite, the PIECES are not.")
print()
print("="*80); print("3.  THE 1/(A-B) AND 1/(AB'-A'B) POLES CANCEL WITHIN EACH PART"); print("="*80)
print("   (a) A -> B :")
print("   %-12s %14s %14s %14s %14s"%("A/B-1","2d","2e","2f/[t1t1']","roundoff?"))
for eps in (1e-2,1e-5,1e-8,1e-11,1e-13):
    ww=x+(w-x)*np.sqrt(A*(1-eps)/B); B2=s(x-ww); C2=s(z-ww)
    T=their_coef(A,B2,C2,Ap,Bp,Cp,k,1e-6,1e6); t11=2.0*(A-B2)*(Ap-Bp)/(4*C2*Cp)
    print("   %-12.0e %14.9f %14.9f %14.9f %14s"%(eps,T['13'],T['12'],T['11_full']/t11,
          "" if eps>1e-10 else "<- cancellation"))
print("   (b) A B' -> A' B :")
print("   %-12s %14s %14s %14s %14s"%("AB'/A'B-1","2a","2c","2b","2f/[t1t1']"))
for eps in (1e-2,1e-5,1e-8,1e-11,1e-13):
    wwp=xp+(wp-xp)*np.sqrt(Ap*B/A*(1+eps)/Bp); Bp2=s(xp-wwp); Cp2=s(z-wwp)
    T=their_coef(A,B,C,Ap,Bp2,Cp2,k,1e-6,1e6); t11=2.0*(A-B)*(Ap-Bp2)/(4*C*Cp2)
    print("   %-12.0e %14.9f %14.9f %14.9f %14.9f"%(eps,T['23'],T['22'],T['33'],T['11_full']/t11))
print()
print("   every part stays finite: the poles are removable and cancel INSIDE each part,")
print("   so you may integrate part by part.  But near those surfaces the cancellation is")
print("   between large numbers -- below ~1e-10 in the ratio, double precision loses it.")
print("   For numerics, expand analytically in a neighbourhood of A=B, A'=B', AB'=A'B.")
