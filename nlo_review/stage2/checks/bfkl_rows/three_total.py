import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from symdec import *
import re
np.seterr(all="ignore")
s=open("/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl/row2b.py").read()
eq=re.search(r'eq=r"""(.*?)"""',s,re.S).group(1)
B=blocks("9ea88190-BFKL.tex"); R=lambda *ls:[rhs(B[l]) for l in ls]
rows=[("4rho Na+Nb (l.34,40)",R(34,40),False),("GI RowII (z-form)",[eq],True),("GI RowIII",R(125,132),True),("GI RowIV",R(187),True),
      ("GII RowI p1",R(246,252),False),("GIII RowI",R(356),True),("GIII RowII",R(405),True),("GIII RowIII p1",R(459),True),
      ("GIII RowV",R(520),True),("GIII RowVI",R(546),True)]
m=lambda q: float(np.max(np.abs(q)))
if __name__=="__main__":
    for seed in (1,2):
        env=setup(3,seed); rng=np.random.default_rng(100+seed); p=list(rng.normal(size=(3,2)))
        v={nm:inv(env,tx,p,cc=cc) for nm,tx,cc in rows}
        tot=sum(v.values())
        print("seed",seed,"symmetrised three-rho coefficient at generic (p1,p2,p3;z,w,w'), units of common prefactor*log:")
        print("  ",{k:round(m(x),4) for k,x in v.items()})
        print("   TOTAL %.3e ; total without 4rho part %.3e ; total without GI RowII %.3e"%(m(tot),m(tot-v[rows[0][0]]),m(tot-v[rows[1][0]])))
        np.save("v3_%d.npy"%seed,np.array([v[k] for k,_,_ in rows]))
