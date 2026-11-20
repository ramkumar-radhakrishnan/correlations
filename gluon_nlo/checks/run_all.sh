#!/bin/sh
# Re-runs every numerical check of the review; output -> RESULTS.txt
cd "$(dirname "$0")"
{
echo "== BFKL block (eqs 1.1-1.3) =="; python3 check_bfkl.py
echo; echo "== JIMWLK block (eqs 2.1-2.4) =="; python3 check_jimwlk.py
echo; echo "== Group II row I second part (7.3-7.8): singular surfaces =="; python3 probe_g2.py
echo; echo "== Group II (7.3-7.8): exchange symmetry =="; python3 sym_g2.py
echo; echo "== colour-factor symmetry, U->1, U->const =="; python3 misc_checks.py
} > RESULTS.txt 2>&1
cd "$(dirname "$0")"
{
echo; echo "== Group II Row IV derivation (DGLAP/RC origin) =="; python3 check_row4.py
echo; echo "== Group I row III (6.4)-(6.8): surfaces, collinear A0, imaginary parts =="; python3 probe_g1r3.py
echo; echo "== Group I row III: coincidence and large-distance limits =="; python3 probe_g1r3_limits.py
echo; echo "== Group I row III: per-term log slopes =="; python3 slopes_g1r3.py
echo; echo "== (7.1): limits =="; python3 probe_71.py
echo; echo "== Group II (7.3)-(7.8): limits =="; python3 probe_g2_limits.py
} >> RESULTS.txt 2>&1
