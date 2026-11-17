import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from symdec import *
rng=np.random.default_rng(0)
for n in (2,3,4):
    M=[rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)) for _ in range(n)]   # generic matrices: identity must hold in any associative algebra
    def val(t): return M[t] if isinstance(t,int) else val(t[0])@val(t[1])-val(t[1])@val(t[0])
    tot=0
    for c,trees in decomp(n):
        ms=[val(t) for t in trees]; s=0
        for perm in itertools.permutations(range(len(ms))):
            q=np.eye(5)
            for i in perm: q=q@ms[i]
            s=s+q
        tot=tot+float(c)*s/factorial(len(ms))
    q=np.eye(5)
    for m_ in M: q=q@m_
    print("n=%d  |product - sum of symmetrised terms| = %.2e"%(n,np.max(np.abs(tot-q))))
