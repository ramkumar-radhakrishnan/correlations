# Minimal LaTeX-subset evaluator: turns an align-block right-hand side into an operator
#   sum_terms  T^{a1..an} rho^{a1}(p1)...rho^{an}(pn)   (ordered, non-commuting rho's)
import re, numpy as np, itertools
from scipy.linalg import expm
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
PTS=["x","x'","y","y'","z","w","w'","u"]
def setup(N=3,seed=0):
    rng=np.random.default_rng(seed)
    gens=[]
    for i in range(N):
        for j in range(i+1,N):
            m=np.zeros((N,N),complex); m[i,j]=m[j,i]=0.5; gens.append(m)
            m=np.zeros((N,N),complex); m[i,j]=-0.5j; m[j,i]=0.5j; gens.append(m)
    for k in range(1,N):
        m=np.zeros((N,N),complex)
        for i in range(k): m[i,i]=1
        m[k,k]=-k; gens.append(m/np.sqrt(2*k*(k+1)))
    t=np.array(gens); n=len(t)
    f=np.zeros((n,n,n))
    for a in range(n):
        for b in range(n):
            c=t[a]@t[b]-t[b]@t[a]
            for d in range(n): f[a,b,d]=np.real(-2j*np.trace(c@t[d]))
    def adj(pos):
        # Wilson line as a deterministic smooth function of the position (so U(y->x)=U(x))
        H=sum(np.sin(ca@pos+ph)*t[a] for a,(ca,ph) in enumerate(zip(CA,PH))); V=expm(1j*H)
        return np.real(np.array([[2*np.trace(t[a]@V@t[b]@V.conj().T) for b in range(n)] for a in range(n)]))
    CA=rng.normal(size=(n,2)); PH=rng.normal(size=n)
    env=dict(N=N,n=n,f=f,adj=adj,X={p:rng.normal(size=2) for p in PTS},
             const={"kp":rng.uniform(0.5,2),"Lam":rng.uniform(0.01,0.1),"vee":1e9,"pp":rng.uniform(3,5),
                    "eps":rng.uniform(0.1,1),"gam":0.5772,"mu":rng.uniform(0.5,2)},Ucache={})
    return env
def Uat(env,p):
    key=tuple(np.round(env["X"][p],14))
    if key not in env["Ucache"]: env["Ucache"][key]=env["adj"](env["X"][p])
    return env["Ucache"][key]
def clean(s):
    s=s.replace(r"^{\prime}","'").replace(r"^\prime","'")
    s=re.sub(r"\\mathbf\{([a-z])\}",r"\1",s)
    s=re.sub(r"\{\\color\{(red|orange)\}\s*\\alpha'?\}","",s)
    s=re.sub(r"\{\\color\{[a-z]+\}\s*([+-])",r"\1{",s)   # a coloured group that starts with a sign is a separate additive term
    s=re.sub(r"\\color\{[a-z]+\}","",s)
    s=re.sub(r"\\int_\{[^{}]*\}","",s)
    s=re.sub(r"e\^\{-ik[^{}]*\}","",s)
    for t in [r"\notag",r"\nonumber",r"\\",r"&",r"\,",r"\;",r"\!",r"\quad",r"\times",r"\left",r"\right",r"\Bigg",r"\Big",r"\bigg",r"\big",r"\int"]:
        s=s.replace(t," ")
    return s
class T:
    __slots__=("c","f","r")
    def __init__(s,c,f=(),r=()): s.c=c; s.f=list(f); s.r=list(r)
def collapse(v):
    if all(not t.r and not t.f for t in v): return [T(sum(t.c for t in v))]
    if all(not t.r and len(t.f)==1 and t.f[0][1]==("?",) for t in v): return [T(1,[(sum(t.c*t.f[0][0] for t in v),("?",))])]
    return v
def isscal(v): return len(v)==1 and not v[0].r and not v[0].f
def isvec(v): return len(v)==1 and not v[0].r and len(v[0].f)==1 and v[0].f[0][1]==("?",)
def mul(a,b): return [T(x.c*y.c,x.f+y.f,x.r+y.r) for x in a for y in b]
class P:
    def __init__(s,env,txt): s.env=env; s.s=clean(txt); s.p=0
    def ws(s):
        while s.p<len(s.s) and s.s[s.p].isspace(): s.p+=1
    def peek(s,t): s.ws(); return s.s.startswith(t,s.p)
    def eat(s,t):
        s.ws()
        if s.s.startswith(t,s.p): s.p+=len(t); return True
        return False
    def need(s,t):
        if not s.eat(t): raise SyntaxError("expected %r at ...%r"%(t,s.s[max(0,s.p-40):s.p+40]))
    def expr(s):
        out=[]; sign=1; first=True
        while True:
            s.ws()
            if s.eat("+"): sign=1
            elif s.eat("-"): sign=-1
            elif not first: break
            if s.atend(): 
                if first: return [T(0)]
                break
            t=s.term()
            for x in t: x.c*=sign
            out+=t; first=False; sign=1
            s.ws()
            if not (s.peek("+") or s.peek("-")): break
        return collapse(out)
    def atend(s):
        s.ws()
        return s.p>=len(s.s) or any(s.s.startswith(t,s.p) for t in [")","]","}","\\}",","])
    def term(s):
        v=[T(1)]
        while True:
            s.ws()
            if s.atend() or s.peek("+") or s.peek("-"): break
            f=s.factor()
            if s.eat(r"\cdot"):
                g=s.factor()
                if isvec(f) and isvec(g): f=[T(f[0].c*g[0].c*(f[0].f[0][0]@g[0].f[0][0]))]
                else: f=mul(f,g)
            v=mul(v,f)
        return v
    def idxs(s):
        s.ws(); m=re.compile(r"[\^_]\{([^{}]*)\}").match(s.s,s.p)
        if not m: raise SyntaxError("index group at %r"%s.s[s.p:s.p+30])
        s.p=m.end(); return tuple(re.findall(r"[a-z]'?",m.group(1)))
    def point(s):
        s.ws(); m=re.compile(r"\(\s*([a-z]'?)\s*\)?").match(s.s,s.p)   # tolerate a missing ')'
        if not m: raise SyntaxError("point at %r"%s.s[s.p:s.p+30])
        if not m.group(0).rstrip().endswith(")"): s.env.setdefault("typos",[]).append("missing ')' after %s"%s.s[max(0,s.p-12):m.end()])
        s.p=m.end(); return m.group(1)
    def group(s):
        s.ws()
        for o in "([{":
            if s.eat(o):
                v=s.expr(); s.ws()
                if not (s.eat(")") or s.eat("]") or s.eat("}")): raise SyntaxError("close at %r"%s.s[max(0,s.p-60):s.p+30])
                return v
        raise SyntaxError("group at %r"%s.s[s.p:s.p+30])
    def factor(s):
        v=s.atom()
        while True:
            s.ws()
            m=re.compile(r"\^\{([^{}]*)\}").match(s.s,s.p)
            if not m: break
            s.p=m.end(); e=m.group(1).strip()
            if e.isdigit():
                e=int(e)
                if isvec(v):
                    a=v[0].f[0][0]*v[0].c; assert e%2==0; v=[T((a@a)**(e//2))]
                else:
                    assert isscal(v),"power of non-scalar"; v=[T(v[0].c**e)]
            else:
                assert isvec(v) and re.fullmatch(r"[a-z]'?",e),("vector index",e)
                v=[T(v[0].c,[(v[0].f[0][0],(e,))])]
        return v
    def atom(s):
        env=s.env; C=env["const"]; s.ws()
        m=re.compile(r"\d+").match(s.s,s.p)
        if m: s.p=m.end(); return [T(float(m.group(0)))]
        if s.eat(r"\frac"):
            a=s.group(); b=s.group(); assert isscal(b),"non-scalar denominator"
            for t in a: t.c/=b[0].c
            return a
        if s.eat(r"\log"):
            pw=1
            if s.eat("^{2}"): pw=2
            if s.peek(r"\frac"): a=s.atom()
            else: a=s.group()
            assert isscal(a),"log of non-scalar"
            return [T(np.log(abs(a[0].c))**pw)]
        if s.eat(r"\mathcal{K}"):
            a=s.group(); assert isvec(a); vv=a[0].c*a[0].f[0][0]; return [T(1,[(vv/(vv@vv),("?",))])]
        if s.eat(r"\rho"):
            ix=s.idxs(); p=s.point(); assert len(ix)==1; return [T(1,[],[(ix[0],p)])]
        if s.eat(r"\delta"):
            ix=s.idxs(); assert len(ix)==2; return [T(1,[(np.eye(2),ix)])]
        for name,val in [(r"\pi",np.pi),(r"\vee",C["vee"]),(r"\Lambda",C["Lam"]),(r"\epsilon",C["eps"]),(r"\gamma",C["gam"]),
                         ("N_{c}",env["N"]),("k^{+2}",C["kp"]**2),("k^{+ 2}",C["kp"]**2),(r"d_{\perp}",2.0),("k^{+}",C["kp"]),("p^{+}",C["pp"]),(r"\mu^{2}_{\overline{MS}}",C["mu"])]:
            if s.eat(name): return [T(val)]
        if s.eat("f^"):
            s.p-=1; ix=s.idxs(); assert len(ix)==3,ix; return [T(1,[(env["f"],ix)])]
        if s.eat("U^"):
            s.p-=1; ix=s.idxs(); p=s.point(); assert len(ix)==2,ix; return [T(1,[(Uat(env,p),ix)])]
        if s.eat(r"\{"):
            a=s.expr()
            if s.eat(","):
                b=s.expr(); s.need(r"\}"); return mul(a,b)+mul(b,a)
            s.need(r"\}"); return a
        if s.peek("(") or s.peek("[") or s.peek("{"): return s.group()
        m=re.compile(r"([xyzwu])('?)").match(s.s,s.p)
        if m: s.p=m.end(); return [T(1,[(env["X"][m.group(0)].copy(),("?",))])]
        if s.eat("i"): return [T(1j)]
        if s.eat("g"): return [T(1.0)]
        raise SyntaxError("atom at ...%r"%s.s[max(0,s.p-50):s.p+50])
def evaluate(env,txt):
    """-> dict {(rho point names in order): tensor}"""
    p=P(env,txt); v=p.expr(); p.ws()
    if p.p<len(p.s): raise SyntaxError("trailing %r"%p.s[p.p:p.p+60])
    out={}
    L="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for t in v:
        if t.c==0: continue
        cnt={}
        # rho indices contracted directly with each other (e.g. rho^a(x) rho^a(x')): insert a colour delta
        ri=[i for i,_ in t.r]
        for i in set(ri):
            if ri.count(i)==2:
                k=[m for m,j in enumerate(ri) if j==i][1]; t.r[k]=(i+"#",t.r[k][1]); t.f.append((np.eye(env["n"]),(i,i+"#")))
        for _,ix in t.f:
            for i in ix: cnt[i]=cnt.get(i,0)+1
        for i,_ in t.r: cnt[i]=cnt.get(i,0)+1
        if env.get("ufilter") is not None:
            uids=set(id(a) for a in env["Ucache"].values())
            if sum(1 for a,_ in t.f if id(a) in uids)!=env["ufilter"]: continue
        bad={i:c for i,c in cnt.items() if c!=2}
        if bad: raise ValueError("index count != 2: %s  rhos=%s"%(bad,t.r))
        mp={i:L[k] for k,i in enumerate(cnt)}
        spec=",".join("".join(mp[i] for i in ix) for _,ix in t.f)+"->"+"".join(mp[i] for i,_ in t.r)
        val=t.c*np.einsum(spec,*[a for a,_ in t.f],optimize=True) if t.f else t.c
        key=tuple(p for _,p in t.r)
        out[key]=out.get(key,0)+val
    return out
def blocks(fn):
    s=open(UP+fn).read(); out={}
    for m in re.finditer(r"\\begin\{align\}(.*?)\\end\{align\}",s,re.S):
        out[s.count("\n",0,m.start())+1]=m.group(1)
    return out
def rhs(blk):
    return blk[blk.index("=")+1:]
def add(a,b,sb=1):
    out=dict(a)
    for k,v in b.items(): out[k]=out.get(k,0)+sb*v
    return out
def norm(a): return max([float(np.max(np.abs(v))) for v in a.values()]+[0.0])
