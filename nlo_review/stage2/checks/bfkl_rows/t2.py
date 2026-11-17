import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from texeval import *
env=setup(3,1)
for s in [r"2", r"\frac{1}{(2\pi)^{3}}", r"(x-w)^{2}", r"\frac{(x-w)\cdot(y-z)}{(x-w)^{2}}", r"U^{ab}(x)\rho^{a}(x)\rho^{b}(y)", r"-\frac{1}{2}\,N_{c}\,\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\,\frac{1}{(x-z)^{2}}\,U^{ba}(w')\,U^{bc}(w)\,\rho^{a}(x')\,\rho^{c}(x)",r"\log{\left(\frac{\vee}{k^{+}}\right)}",r"\Bigg[ 2 U^{ab}(x)\rho^{a}(x)\rho^{b}(y)\Bigg]"]:
    r=evaluate(env,s); print(s[:40],{k:norm({k:v}) for k,v in r.items()})
