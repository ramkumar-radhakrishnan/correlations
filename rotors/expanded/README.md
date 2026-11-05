# Classical entanglement and entropy of two coupled quantum rotors — expanded draft

`rotor_entanglement_expanded.pdf` (31 pp., `elsarticle`, for *Annals of Physics*)
expands every calculation of the earlier draft step by step, at graduate level.

## Reproducing (Python 3.13: `numpy scipy matplotlib`)

```bash
python3.13 verify_all.py        # E2, Mathieu map, large-n asymptotics, deep libration
python3.13 verify2.py           # large-n check, fringe deficit at generic energies
python3.13 verify_integrals.py  # every closed-form integral
python3.13 bessel_check.py      # Wannier-Stark / Bessel solution, 1 - ln 2 up to gamma = 1e6
python3.13 sector_check.py      # even/odd-M Mathieu mapping, librating-state count
python3.13 figs.py && python3.13 figs_rot.py && python3.13 figs_rest.py
python3.13 table_checks.py      # Table B.1 rows -> table_checks.tex
pdflatex rotor_entanglement_expanded.tex   # x3
```

## Corrections to the earlier draft found while expanding

| where | issue |
|---|---|
| Eq. (38), large-n entropy | `ln(ħ²n/λI)` must be `ln(2ħ²n/λI)`; the draft was ~10 % low |
| §9.1, deep libration | intermediate `m_eff = I/2ħ²`, `ω_eff = √(2λħ²/I)` dimensionally inconsistent; redone with reduced moment `μ = I/2` (final result unchanged) |
| §3.4, second order | also requires `n ≠ ±2` (second-order exchange resonance) |
| abstract, §7.2 | the classical *entropy* is closed-form only in deep rotation; the classical *distribution* is closed-form everywhere |
| §8 | "agree to a few percent" refined: the absolute deficit tends to `1 − ln 2` |

## Results added by the expansion

- `E⁽²⁾ = λ²I / [2ħ²(n²−1)]`
- every sector `M` is exactly a Mathieu problem (periodic for even `M`, antiperiodic for odd), `q = 2λI/ħ² = 4α`; `N_lib ≈ (4/π)√q`
- classical normalisation via complete elliptic integrals; it equals the classical period
- the entropy maximum at the separatrix follows from a cusp in the support of the classical distribution
- deep rotation is a Wannier–Stark ladder: Schmidt weights `J_j²(γ)`, `γ = λI/ħ²n`, giving one universal crossover curve
- `S_q = S_cl − (1 − ln 2)` in the semiclassical limit (momentum-space WKB), approached as `γ^(−1/3)`; checked to `γ = 10⁶`
- deep libration: `S₀ = ¼ln(λI/ħ²) + ½ln(πe/√2) − 1/[8√(2λI/ħ²)]`

The bibliography was written from memory and must be checked before submission.
