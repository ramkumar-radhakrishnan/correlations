"""Classical (symmetrized) three-rho values: engine terms (as LaTeX) vs the LaTeX first form minus J, per row."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_tex as ET
C, T3, R, J, U = ET.C, ET.T3, ET.R, ET.J, ET.U
out = {}
for name, (eng, first, jl, cc) in ET.SETS.items():
    ts = R.THREE[eng]()
    c = T3.configs(5)
    e = sum(t['c'] * 3 ** t.get('ncp', 0) * T3.sym_value(ET.tex(t), c, False) for t in ts)
    u = sum(T3.sym_value(C.prep(C.plus(U[l])), c, cc) for l in first) - sum(T3.sym_value(J[l], c, cc) for l in jl)
    out[name] = dict(engine=str(e), first_minus_J=str(u), rel=abs(e - u) / abs(u))
    print(name, out[name], flush=True)
json.dump(out, open(os.path.join(HERE, 'engine_tex3.json'), 'w'), indent=1)
