# EII Lab

Measurement and estimation for enforcement data produced by the
institutions being measured.

---

## Why this is a lab and not a paper

On 2026-09-02 a parameter-ordering bug produced a frequentist
coefficient table in which every variable carried another variable's
estimate. `statsmodels` returns ZINB parameters as
`[inflate, count, alpha]`; the extraction code assumed
`[count, alpha, inflate]`. The error survived three weeks and a full
round of paper planning. It was caught on 2026-10-05 only because the
model was ported to R and the two implementations disagreed.

The lesson is not "check your indexing." It is that single-platform
work hides errors that cross-platform replication exposes in seconds,
and that a finding is only as good as the machinery that could have
falsified it. Everything below exists because of that day.

---

## Operating rules

1. **Two implementations or it isn't a result.** Every coefficient
   that enters a manuscript is produced by `pscl` and `brms` in R and
   agrees with the Python pipeline. Disagreement is a bug until proven
   otherwise.

2. **`validate.R` before anything leaves.** Run it before sending
   results, after any panel change, after any code change, after any
   package upgrade. A failure is not necessarily an error, but it must
   be explained before the number travels.

3. **Predictions and falsification conditions go in the script,
   written before the data is seen.** Every test script in this repo
   states what would confirm the hypothesis and what would kill it,
   above the code. This is the only defense against rebuilding an
   instrument until it gives the wanted sign.

4. **No number enters a manuscript from a chat log, a figure, or an
   old draft.** Only from a locked run of code in this repo.

5. **Nulls are results.** The PLR interaction is null. The layer
   substitution test is null. Both are in the repo with their
   verdicts.

---

## Verified baseline (2026-10-05)

Environment: R 4.5.3, pscl 1.5.9, brms 2.23.0, cmdstan 2.40.0, seed 42.
Library at `C:/R/library` (OneDrive blocks DLL loads).
brms must run from a terminal R session — Positron's console blocks
cmdstanr's `sink()`.

**Panel.** 942 country-years, 44 states, 1945–2025 frame.
**Estimation sample.** 535 complete cases, 32 states, 1970–2021.
49.9% zeros, variance/mean 1133.1.

**Count equation** (V-Dem primary; pscl and brms agree to MCMC noise):

| Variable | Freq | Bayes | P(dir.) |
|---|---|---|---|
| resource | 0.336 | 0.323 | 0.968 |
| friction | −0.412 | −0.442 | 0.977 |
| degradation | −1.083 | −1.051 | 0.995 |
| fiscal (cost) | −1.081 | −1.061 | 1.000 |
| trend | −0.340 | −0.310 | 0.854 |
| regime | 0.353 | 0.348 | 0.938 |
| gdp per capita | −0.669 | −0.659 | 1.000 |
| friction × fiscal | 0.057 | 0.038 | 0.582 |

**Zero-inflation.** `pscl` and `brms` both model P(structural zero);
positive = more non-enforcement. PyMC's `psi` is the opposite
convention — P(count process) — so every PyMC inflate sign reads
inverted. Verified across platforms 2026-10-05.

**Robust across everything** (two frameworks, three regime measures,
four priors, four DV weightings): fiscal cost negative, GDP per capita
negative, degradation negative, friction negative, interaction null.

---

## Open questions

**The regime break.** Frequentist period split finds both fiscal
mechanisms operating 1970–2002 (spend +1.58, cost −1.21) and neither
operating 2003–2021. Under Bayesian verification — see
`regime_break_bayes.R`. Not a finding until that run is read.

**Spend and cost are different things.** `cor = 0.042` over 1970–2021.
Variance profiles invert: spend SD 1.52 → 3.20, cost SD 0.70 → 0.16.
The original specification contained only cost, and cost stops moving
exactly when the theory says expansion begins.

**The panel ends before the test case.** Net interest was 1.49% of GDP
in 2021 and 3.15% in 2025. The constraint returned after the data
stopped. Extension to 2025 is the highest-value build available.

**Extensive vs intensive expansion.** The DV counts events, so it
measures states-entering-the-frame well and apparatus-thickening
poorly. A new reporting rule generates zero events and imposes
enormous cost. Scope condition, not flaw — but it should be stated.

**The fiscal divergence.** Residual from `net_interest ~ debt + yield`
clears all three diagnostic gates (trend r = 0.054, variance ratio
0.73, AR1 0.884) but the residual pattern is maturity structure, not
reporting distortion. A corrected benchmark needs Treasury average
interest rates, available only from 2001.

**IV-side triangulation does not exist.** The DV triangulates ten
sources. Every predictor is single-source: CBO for fiscal, Voeten for
friction, V-Dem/Polity for regime. The priors are diffuse because
nothing exists to make them informative. This is the measurement
asymmetry and it is the honest answer to "why not informative priors."

---

## Dead ends (do not revisit without new evidence)

**"MLE and Bayes disagree on corrupted data."** Was the indexing bug.
Both estimators agree on all eight count coefficients and, after the
psi correction, all five inflate coefficients.

**Sequential updating from the proxy panel.** The 25-country proxy EII
was built in an ephemeral container in May 2026 and no longer exists on
disk. Its posteriors cannot be verified and must not be used as priors.
The JCR-era proxy also found the interaction null (0.095, P(dir.) 0.688),
so "the proxy found it and the full panel lost it" was never true.

**Single-factor fiscal strain index.** Five candidate indicators do not
load on one factor. `net_interest` vs `debt_gdp` correlate at −0.010;
`yield10` vs `debt_gdp` at −0.773. Flow pressure and stock burden are
distinct and a one-factor model would discard the information in their
disagreement.

**Layer substitution under fiscal strain.** All four fiscal
coefficients null in 2003–2021 (p = 0.35 to 0.96). No asymmetry,
because the fiscal variable has no variance in that window.

---

## Files

| File | Purpose |
|---|---|
| `EII_pipeline.R` | Main estimation. 4 priors + elections + Polity + DV sensitivity |
| `validate.R` | Regression tests against the verified baseline |
| `fiscal_pull.R` | Five fiscal indicators from FRED, 1945–2025 |
| `divergence_diagnostic.R` | Three gates on the reporting-gap construct |
| `layer_decomposition.R` | Spectacle vs apparatus substitution test |
| `period_split.R` | Early/late split on the original specification |
| `spend_cost_split.R` | Spend and cost as separate mechanisms |
| `regime_break_bayes.R` | Bayesian verification of the period break |
