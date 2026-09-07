"""
Parameterization convention test for pm.ZeroInflatedNegativeBinomial.
Verified 2026-09-07 on the estimation machine: zero share = 0.2473 at
psi=0.9, mu=5, alpha=1 (vs ~0.92 if psi were the zero-inflation prob).
Conclusion: psi = P(count process), NOT P(structural zero).
Inflate-equation gammas are therefore read on the ENGAGEMENT side.
"""
import pymc as pm
import numpy as np

draws = pm.draw(
    pm.ZeroInflatedNegativeBinomial.dist(psi=0.9, mu=5, alpha=1),
    draws=20000, random_seed=1,
)
zero_share = (draws == 0).mean()
print("zero share:", round(float(zero_share), 4))

# psi = P(count process): zeros come only from (1-psi) inflation plus the
# NegBin's own zeros -> ~0.25. If psi were P(structural zero) -> ~0.92.
assert zero_share < 0.5, "psi convention has changed - re-audit before trusting inflate signs"
print("PASS: psi is the count-process probability; inflate gammas read on the engagement side.")