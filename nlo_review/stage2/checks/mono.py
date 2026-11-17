# Parser/evaluator for the monomial lists of Evolution_JIMWLK__terms.tex and BFKL.tex
#   each term:  coefficient * kernels * colour monomial with two rho's
import re, numpy as np
from fractions import Fraction
from scipy.linalg import expm
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
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
    def adj():
        H=sum(rng.normal()*t[a] for a in range(n)); V=expm(1j*H)
        return np.real(np.array([[2*np.trace(t[a]@V@t[b]@V.conj().T) for b in range(n)] for a in range(n)]))
    pts=["x","x'","y","y'","z","w","w'","u"]
    env=dict(N=N,n=n,f=f,U={p:adj() for p in pts},X={p:rng.normal(size=2) for p in pts},
             R=rng.normal(size=(n,n)))   # R[b,c] stands for rho^b(x') rho^c(x)  (generic, not assumed singlet)
    return env
def vec(env,expr):
    # expr like "x-w" or "x'-w'" or "-w'+x'"
    v=np.zeros(2)
    for sgn,name in re.findall(r"([+-]?)\s*([a-z]'?)",expr):
        v=v+(-1 if sgn=='-' else 1)*env["X"][name]
    return v
IDX=r"[a-z]'?"
def parse_term(env,s,rho_order=("x'","x")):
    """returns numeric value of one monomial; s is the LaTeX of the term (sign included)"""
    orig=s
    s=s.replace(r"\,","").replace(r"\notag","").replace(r"\\","").replace("&","").replace(r"\Bigg[","").replace(r"\Bigg]","").strip()
    s=s.replace(r"^{\prime}","'")
    val=1.0
    # kernels K(a).K(b)
    for m in list(re.finditer(r"\\mathcal\{K\}\(([^()]*)\)\\cdot\s*\\mathcal\{K\}\(([^()]*)\)",s)):
        a=vec(env,m.group(1)); b=vec(env,m.group(2)); val*= (a@b)/((a@a)*(b@b))
    s=re.sub(r"\\mathcal\{K\}\(([^()]*)\)\\cdot\s*\\mathcal\{K\}\(([^()]*)\)","",s)
    for m in list(re.finditer(r"\\frac\{1\}\{\(([^()]*)\)\^\{2\}\}",s)):
        a=vec(env,m.group(1)); val*=1.0/(a@a)
    s=re.sub(r"\\frac\{1\}\{\(([^()]*)\)\^\{2\}\}","",s)
    # colour
    ops=[]; sub=[]
    for m in re.finditer(r"f\^\{("+IDX+")("+IDX+")("+IDX+r")\}",s): ops.append(env["f"]); sub.append(m.groups())
    s=re.sub(r"f\^\{("+IDX+")("+IDX+")("+IDX+r")\}","",s)
    for m in re.finditer(r"U\^\{("+IDX+")("+IDX+r")\}\(([a-z]'?)\)",s): ops.append(env["U"][m.group(3)]); sub.append(m.groups()[:2])
    s=re.sub(r"U\^\{("+IDX+")("+IDX+r")\}\(([a-z]'?)\)","",s)
    rhos=re.findall(r"\\rho\^\{("+IDX+r")\}\(([a-z]'?)\)",s)
    s=re.sub(r"\\rho\^\{("+IDX+r")\}\(([a-z]'?)\)","",s)
    assert len(rhos)==2,("need two rho",orig)
    d=dict((p,i) for i,p in rhos)
    assert set(d)==set(rho_order),("rho points",rhos,orig)
    ops.append(env["R"]); sub.append((d[rho_order[0]],d[rho_order[1]]))
    if "N_{c}" in s: val*=env["N"]; s=s.replace("N_{c}","")
    # numeric coefficient
    c=s.replace(" ","")
    m=re.fullmatch(r"([+-]?)(?:\\frac\{(\d+)\}\{(\d+)\}|(\d+))?",c)
    assert m,("leftover",repr(c),orig)
    num=1.0
    if m.group(2): num=int(m.group(2))/int(m.group(3))
    elif m.group(4): num=int(m.group(4))
    if m.group(1)=='-': num=-num
    # einsum
    letters={}
    def L(ix):
        if ix not in letters: letters[ix]="abcdefghijklmnopqrstuvwxyz"[len(letters)]
        return letters[ix]
    spec=",".join("".join(L(i) for i in ss) for ss in sub)
    col=np.einsum(spec+"->",*ops)
    return num*val*col
def block_terms(blk):
    """split the content inside the outer \\Bigg[ ... \\Bigg] of an align block into monomial strings"""
    a=blk.index(r"\Bigg["); b=blk.rindex(r"\Bigg]")
    body=blk[a+len(r"\Bigg["):b]
    parts=re.split(r"\\notag\s*\\\\\s*&?",body)
    return [p.strip() for p in parts if p.strip()]
def blocks(fn):
    s=open(UP+fn).read(); out=[]
    for m in re.finditer(r"\\begin\{align\}(.*?)\\end\{align\}",s,re.S):
        out.append((s.count("\n",0,m.start())+1,m.group(1)))
    return out
def eval_block(env,blk):
    return sum(parse_term(env,t) for t in block_terms(blk))
