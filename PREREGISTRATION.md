# Pre-Registration: The 2022–2025 Extension

**Registered 2026-10-05, before any 2022–2025 enforcement data has been
coded.** No observation outside the existing 1945–2021 panel has been
examined. This document exists so the extension tests a standing
prediction rather than accommodating whatever the data turns out to show.

---

## Background

The panel's estimation sample runs 1970–2021. Within it, the fiscal
relationship is period-dependent:

| | 1970–2002 | 2003–2021 |
|---|---|---|
| Spend (net outlays, % GDP) | +1.004, P(dir.) 0.933 | +0.063, P(dir.) 0.637 |
| Cost (net interest, % GDP) | −0.728, P(dir.) 0.864 | −0.032, P(dir.) 0.599 |
| Raw SD of spend | 1.52 | 3.20 |
| Raw SD of cost | 0.70 | 0.16 |

Bayesian ZINB, tight priors (N(0, 1.5)), R-hat 1.0025, zero divergences.
Early sample 195 observations / 18 states / 77% zeros — thin, and the
early estimates are weakly supported rather than established.

The late-period nulls are not evidence of absence. Carrying cost barely
moved in that window: 0.16 percentage points of standard deviation
across nineteen years. A variable with no variance cannot produce a
coefficient.

**That changed after the panel ends.** Net interest as a share of GDP:

- 2021: 1.49%
- 2022: 1.83%
- 2023: 2.37%
- 2024: 3.00%
- 2025: 3.15%

Cost more than doubled in four years. For the first time since the
mid-1990s there is real variance in the fiscal constraint — and none of
it is in the data.

---

## The prediction: institutional lag

**Cost rises sharply 2022–2024 and enforcement intensity does not yet
respond.** The posterior for `cost_z` in the extended window will
remain centered near zero, P(dir.) below roughly 0.80, despite the
return of real variance in carrying cost.

**The mechanism is lag, not deafness.** An earlier formulation of this
prediction held that the architecture had permanently lost fiscal
sensitivity. That claim is not supported by the comparative record and
is withdrawn. The historical cases — Rome, Britain — show institutional
response to fiscal exhaustion running well behind the fiscal signal
itself. Commitments are contractual, personnel are tenured, obligations
are treaty-bound, and the apparatus continues operating on prior
authorization long after the conditions that authorized it have changed.
Retrenchment follows exhaustion by years, sometimes decades, and when it
comes it is abrupt rather than proportional.

The American case is predicted to exhibit the same structure. Carrying
cost returned in 2022. Response is not expected within three years.

This is a weaker claim than permanent insensitivity and it is the right
one: it is grounded in the comparative cases rather than asserted about
the present, and it generates a continuing test rather than a terminal
one.

---

## What falsifies it

**Cost coefficient returns negative with substantial posterior mass**
(P(dir.) > 0.90, HDI excluding zero) in the 2022–2024 window. The
architecture responds to fiscal signals within three years. The lag
mechanism is wrong, response is faster than the comparative cases
suggest, and the American case is not structurally analogous to the
historical ones on this dimension.

**Cost coefficient returns positive.** Enforcement rises with carrying
cost. Neither the lag mechanism nor the overextension thesis predicts
this. If it appears, both need revision rather than patching.

---

## What confirms it

Cost coefficient near zero in 2022–2024, with a tighter posterior than
2003–2021 produced — tighter because the IV now has variance to
estimate from. A null on a variable that moves is informative in a way
a null on a flat variable is not.

**Confirmation is provisional by construction.** A lag claim cannot be
confirmed by a single null; it predicts that the response arrives
later. The standing test is whether the cost coefficient turns negative
in subsequent extensions — 2027, 2029 — as the lag runs out. This
prediction therefore does not resolve within one paper. It is a
commitment the lab carries forward, and the timing of the turn, if it
comes, is itself the estimate of the lag.

---

## On what "collapse" means and why it is not the predicted outcome

The theory holds that coercive extraction architectures overextend and
reconfigure. It does not specify what the American architecture
reconfigures into, and that terminal state is treated here as unknown.

An outcome that cannot be specified cannot be predicted. "Collapse" is
therefore excluded as a registered prediction: reconfiguration into an
unspecified form is compatible with nearly any pattern in the data,
which makes it a frame rather than a hypothesis.

What the comparative cases do supply is the structural signature of
reconfiguration while it is underway — contraction at the periphery,
intensification at the core, substitution of cheap instruments for
expensive ones, and lag between fiscal exhaustion and institutional
response. Those are specific enough to test without knowing the
endpoint, and the lag prediction above is one of them.

**Scope condition, registered in advance.** The EII counts enforcement
events directed at states. If reconfiguration moves the architecture's
coercive output into modalities the index does not count — domestic
prosecution, regulatory conditioning, informal pressure, compliance
apparatus that generates no countable event — then a falling EII
registers the modality shift rather than a decline in coercive output.
The instrument would read reconfiguration as decline. This is a limit
of the measure, not of the theory, and any interpretation of a falling
EII must address it before claiming contraction.

---

## Secondary prediction: the friction interaction

`fric_x_cost` is positive and excludes zero in 1970–2002
(mean 0.979, HDI [0.025, 1.885], P(dir.) 0.985) and null in 2003–2021
(mean 0.045, P(dir.) 0.600). This is the only window in the entire
project where that interaction has been distinguishable from zero.

Two candidate mechanisms, registered without choosing between them:

**Adversary override.** Under bipolarity, architectural distance
signals adversary status. Enforcement against adversaries is
strategically required regardless of cost, so fiscal strain does not
deter it. Post-1991, high friction signals irrelevance rather than
threat, the override disappears, and friction becomes a pure cost term
(main effect falls from −1.27 to −0.33).

**Ally vetting.** Enforcement shifts from punishing adversaries to
conditioning partners — compliance requirements as the price of
access. This mechanism operates through apparatus rather than events
and would be largely invisible to an event-count DV.

**Prediction.** If adversary override is the mechanism, the interaction
should return positive in 2022–2025: great-power competition has
resumed and carrying cost has risen simultaneously, which is the
configuration that produced it in the first period. If ally vetting is
the mechanism, or if the early result was a coincidence of two
correlated time series, it should not return.

**Caveat registered in advance.** The early-period interaction sits on
195 observations with its HDI lower bound at 0.025 — one observation
from straddling zero. Collinearity between friction and cost in the
early window has not been checked and could produce this pattern
without any strategic logic. This is the weakest result in the project
and is registered as a question, not a finding.

---

## Principal threat to inference: political declination

The prediction above says enforcement will not respond to rising cost
because the architecture has lost fiscal sensitivity. A competing
explanation produces the identical null: enforcement is low in the
extension window because the executive is deliberately not enforcing —
declination policy, pardons, dropped and unfiled cases.

Fiscal insensitivity and political non-enforcement are observationally
equivalent in a model whose only fiscal covariate is carrying cost.
Both yield a flat coefficient. This is registered in advance as the
principal threat to inference, not raised after the fact.

Three features of the design bear on it.

**Timing is partly separable.** Carrying cost began rising in 2022,
under the prior administration; the current one took office in January
2025. If enforcement is already flat against rising cost over
2022–2024, that is evidence for fiscal insensitivity independent of
any current declination posture. 2025 is confounded.

*Decision registered in advance:* the primary test window is
**2022–2024**. 2025 is coded and reported but analyzed separately. The
prediction above is evaluated on the three clean years. If the clean
window and the full window disagree, both are reported and the
disagreement is the finding.

**Layers may diverge.** Political non-enforcement operates most
readily on Layer 3 — civil penalties against financial institutions
are discretionary, negotiated, and straightforward to simply not
pursue. UNSC designations (Layer 2) require multilateral process;
kinetic operations (Layer 4) are not a prosecutorial choice. A
declination signature should therefore appear as Layer 3 collapsing
while Layers 2 and 4 hold. A fiscal signature should appear across
layers or not at all.

*Registered test:* run the layer decomposition on the extension window.
Asymmetric collapse concentrated in Layer 3 favors declination.
Uniform flatness favors fiscal insensitivity.

**Declination is measurable, and measuring it reproduces the paper's
own problem.** DOJ and OFAC publish enforcement statistics; declination
rates and settlement counts are trackable and could enter as a
covariate. But those series are produced by the institutions doing the
declining — a single-source institutional measure of institutional
non-action, which is precisely the data-production problem this project
exists to address, recurring one level up. Any declination covariate
is therefore reported as a descriptive control, not as an identifying
strategy, and its own provenance is stated.

**What this threat does not excuse.** If the prediction fails — if cost
returns negative with substantial posterior mass — political
declination cannot be invoked to rescue it. A confound that could only
have produced the predicted result, and is unavailable to explain the
opposite one, is not a confound. It is an alibi. Registered here so
that it cannot be deployed that way later.

---

## Analysis plan

1. Code 2022–2025 using the existing four-layer protocol. No changes to
   source definitions, weighting, or inclusion rules. Any change is
   documented before coding begins.
2. Run `validate.R` on the extended panel. The 1945–2021 baseline
   coefficients must reproduce. If they move, the extension introduced
   an error.
3. **Primary test:** re-estimate the late period as 2003–2024 using the
   same specification, same priors, same seed. This is the window the
   prediction is evaluated on.
4. **Secondary:** re-estimate as 2003–2025 including the confounded
   year. Report both. Disagreement between them is itself reported.
5. Run the layer decomposition on 2022–2024 and on 2025 separately.
   Asymmetric Layer 3 collapse favors declination; uniform flatness
   favors lag.
6. Report the cost coefficient against the prediction above.
7. Report the interaction against the secondary prediction.
8. **Reconfiguration signature.** Two composition tests, reported
   whatever the cost coefficient shows, because they bear on the
   scope condition above:
   - *Instrument substitution.* Does the layer mix shift toward
     cheaper instruments — multilateral and financial up relative to
     kinetic — independent of the total?
   - *Target reallocation.* Do structural zeros reallocate across the
     panel rather than falling uniformly? The zero-inflation equation
     estimates this directly. Periphery contraction predicts sanctuary
     expanding among low-salience states while the core holds.

No specification search. No alternative windows. No dropping of years
that produce inconvenient results. If the specification needs to change,
that is a separate analysis reported as such.

---

## Power

Four years adds roughly 100–130 country-years to a 535-observation
sample. That is not enough to identify much on its own. The prediction
rests on those years supplying variance in a variable that had almost
none, not on sample size.

A null result on 19 flat years and a null result on 23 years with real
variance are different claims. Only the second is informative.

---

---

## Provenance

**2026-10-05, first draft.** Prediction stated as permanent fiscal
insensitivity ("the apparatus no longer reads the bill").

**2026-10-05, revised before commit.** Reformulated as institutional
lag grounded in the comparative cases. The insensitivity claim was
withdrawn as unsupported by the historical record and untestable within
one paper. Added the political-declination threat and the 2022–2024
primary window. Added the scope condition on modality shift and the
exclusion of "collapse" as a registered outcome. No 2022–2025
enforcement data had been examined at any point during drafting or
revision.

Both versions are preserved in this file's git history. The revision
preceded all data collection.

---

*Prediction: cost rises, enforcement does not respond within three years.*
*Mechanism: institutional lag, per the comparative cases.*
*Primary window 2022–2024. 2025 reported separately.*
*Standing test: does the coefficient turn negative in later extensions.*
*Registered 2026-10-05.*
