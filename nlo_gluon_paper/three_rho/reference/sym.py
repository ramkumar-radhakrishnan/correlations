"""Symbolic terms:  c * prod U^{ij}(p) * prod f^{abc} * ordered prod rho^{a}(p) * kernel
   kernel = tuple of dot-pairs ((p1,q1),(p2,q2)) meaning  K(p1-q1).K(p2-q2),  K(X)=X/X^2.
   Positions: external 'w' (amp), 'wb' (= w', conj), 'z' (soft / p gluon); everything else is integrated."""
import numpy as np, itertools, sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from weyl import T as WT, weyl as wweyl
EXT=('w','wb','z')
class Term:
    __slots__=('c','Us','fs','rs','kern','tag')
    def __init__(s,c,Us,fs,rs,kern,tag=''):
        s.c=complex(c); s.Us=tuple(Us); s.fs=tuple(fs); s.rs=tuple(rs); s.kern=tuple(kern); s.tag=tag
    def __repr__(s): return f"{s.c:.3g} U{s.Us} f{s.fs} r{s.rs} K{s.kern} [{s.tag}]"
_ic=itertools.count()
def fresh(t,suffix):
    """rename all colour indices of t (append suffix) except the free ones 'a','h' """
    def r(i): return i if i in ('a','h') else i+suffix
    return Term(t.c,[(r(i),r(j),p) for i,j,p in t.Us],[tuple(r(x) for x in f) for f in t.fs],
                [(r(i),p) for i,p in t.rs],t.kern,t.tag)
def rename_pos(t,mp):
    g=lambda p: mp.get(p,p)
    return Term(t.c,[(i,j,g(p)) for i,j,p in t.Us],t.fs,[(i,g(p)) for i,p in t.rs],
                [((g(a),g(b)),(g(c),g(d))) for (a,b),(c,d) in t.kern],t.tag)
def conj(t):
    """hermitian conjugate: c*, reversed rho order (U, f real)"""
    return Term(np.conj(t.c),t.Us,t.fs,t.rs[::-1],t.kern,t.tag)
def mul(t1,t2,tag=''):
    return Term(t1.c*t2.c,t1.Us+t2.Us,t1.fs+t2.fs,t1.rs+t2.rs,t1.kern+t2.kern,tag)
# ---------------- numeric evaluation as operators on the valence space ----------------
class Num:
    def __init__(s,ctx,modelname={'w':'wp','wb':'w','z':'z'}):
        s.c=ctx; s.mn=modelname
        s.Rs={i:np.array([ctx.R[(i,a)] for a in range(8)]) for i in range(len(ctx.Xs))}
        from su3 import F
        s.F=F
    def pos(s,p,asg): return s.c.pos[s.mn[p]] if p in s.mn else s.c.pos[s.c.Xs[asg[p]]]
    def U(s,p,asg): return s.c.U[s.mn[p]] if p in s.mn else s.c.U[s.c.Xs[asg[p]]]
    def kern(s,t,asg):
        v=1.0
        for (a,b),(c_,d) in t.kern:
            A=np.subtract(s.pos(a,asg),s.pos(b,asg)); B=np.subtract(s.pos(c_,asg),s.pos(d,asg))
            v*=(A@B)/((A@A)*(B@B))
        return v
    def colour(s,t,asg,free=()):
        letters={}
        def L(i):
            if i not in letters: letters[i]=chr(ord('a')+len(letters)) if len(letters)<26 else chr(ord('A')+len(letters)-26)
            return letters[i]
        subs=[];ops=[]
        for i,j,p in t.Us: subs.append(L(i)+L(j)); ops.append(s.U(p,asg))
        for f in t.fs: subs.append(''.join(L(x) for x in f)); ops.append(s.F)
        slots=[]
        for k,(i,p) in enumerate(t.rs):
            sl=L(('slot',k)); subs.append(L(i)+sl); ops.append(np.eye(8)); slots.append(sl)
        out=''.join(slots)+''.join(L(x) for x in free)
        return np.einsum(','.join(subs)+'->'+out,*ops,optimize=True) if subs else np.array(1.0)
    def dummies(s,t):
        ps=[]
        for _,_,p in t.Us: ps.append(p)
        for _,p in t.rs: ps.append(p)
        for (a,b),(c_,d) in t.kern: ps+= [a,b,c_,d]
        out=[]
        for p in ps:
            if p not in EXT and p not in out: out.append(p)
        return out
    def op(s,t,symm=False):
        Vd=s.c.Vd; tot=np.zeros((Vd,Vd),complex); dm=s.dummies(t); n=len(t.rs)
        perms=list(itertools.permutations(range(n))) if symm else [tuple(range(n))]
        for img in itertools.product(range(len(s.c.Xs)),repeat=len(dm)):
            asg=dict(zip(dm,img)); kv=s.kern(t,asg)
            if kv==0: continue
            C=s.colour(t,asg)*t.c*kv
            sites=[asg[p] for _,p in t.rs]
            for pm in perms:
                Cp=np.transpose(C,pm) if n else C; st=[sites[k] for k in pm]
                if n==0: tot+=Cp*np.eye(Vd); continue
                acc=np.einsum('a...,aij->...ij',Cp,s.Rs[st[0]])
                for k in range(1,n):
                    acc=np.einsum('b...ij,bjk->...ik',acc,s.Rs[st[k]])
                tot+=acc/len(perms)
        return tot
