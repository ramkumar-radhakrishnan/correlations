"""The consolidated three-rho list ('All three rho term with c.c.') of the new Evolution_three_rho.tex:
identical to ../three_rho_v3/input/three_rho_all_v3.tex?  and does it cancel?  (re-runs the v3 parser and check)."""
import json, os, re, sys, importlib
HERE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.join(HERE, '..', 'three_rho_v3')
new = open(os.path.join(HERE, 'input', 'Evolution_three_rho.tex')).read()
sec = new[new.index(r'\section{All three'):new.index(r'\section{Cancellation of all three')]
old = open(os.path.join(V3, 'input', 'three_rho_all_v3.tex')).read()
norm = lambda s: re.sub(r'\s+', '', s)
same = norm(sec) == norm(old)
tmp = os.path.join(HERE, 'three_rho_all_new.tex')
open(tmp, 'w').write(sec)
sys.path.insert(0, V3)
os.chdir(V3)
import parse3
parse3.SRC = tmp
import check_cancel as CC
importlib.reload(CC)
out = dict(identical_to_v3=same)
print('consolidated list identical to the v3 list:', same)
json.dump(out, open(os.path.join(HERE, 'three_cancel.json'), 'w'), indent=1)
E = parse3.parse()
res = CC.summary(E, True)
out['n_terms'] = len(E)
out['cancellation_U_symmetric'] = res
print('terms', len(E), 'non-cancelling monomials per seed', [r['nonzero'] for r in res],
      'largest residual %.2e (largest single term %.2e)' % (max(r['max_residual'] for r in res), res[0]['max_single']))
json.dump(out, open(os.path.join(HERE, 'three_cancel.json'), 'w'), indent=1)
os.remove(tmp)
