"""Run each step once (seed 1) and report errors; usage: python3 quick.py Group_X [step index ...]"""
import sys
import time
import traceback
import check_all as CA
g = sys.argv[1]
idx = [int(a) for a in sys.argv[2:]] or range(len(CA.STEPS[g]))
for n in idx:
    name, kind, e1, e2, opt = CA.STEPS[g][n]
    t0 = time.time()
    try:
        r = CA.run_step(g, kind, e1, e2, opt, 1)
        print('%2d %-48s %s  (%.0fs)' % (n, name, r, time.time() - t0), flush=True)
    except Exception as e:
        print('%2d %-48s ERROR %s' % (n, name, str(e)[:300]), flush=True)
        traceback.print_exc(limit=3)
