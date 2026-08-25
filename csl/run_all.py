"""Regenerate RESULTS.txt."""
import subprocess, sys
from cases import report
from matrices import show
subprocess.run([sys.executable, "verify.py"], check=True)
for nq, N in [(1, 2), (1, 4), (2, 2), (2, 4)]:
    report(nq, N)
for nq, N in [(1, 2), (1, 4), (2, 2)]:
    show(nq, N)
subprocess.run([sys.executable, "spectra.py"], check=True)
