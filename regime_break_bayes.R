# ==============================================================================
# REGIME BREAK — BAYESIAN VERIFICATION
# ==============================================================================
# WHAT THIS TESTS
# The frequentist period split found that enforcement tracks fiscal
# conditions in 1970-2002 (spend +1.58, p = .053; cost -1.21, p = .070)
# and decouples in 2003-2021 (spend +0.02, p = .92; cost -0.05, p = .71).
#
# That result is NOT yet trustworthy:
#   - Early sample is 195 observations, 18 states, 76.9% zeros.
#   - The variance-covariance matrix returned NaNs: some standard
#     errors failed to compute, so the p-values are unreliable.
#   - Spend and cost correlate at 0.587 in the early window, so the
#     two coefficients are partly fighting each other.
#   - Effective N for spend and cost is 33 YEARS, not 195 country-years.
#
# Bayesian estimation does not fix the small sample. It does three
# things MLE could not here: it returns honest uncertainty instead of
# a failed Hessian, it regularizes coefficients that MLE lets run to
# the boundary on sparse data, and it gives a direct probability
# statement about the sign rather than a p-value from an SE that
# did not compute.
#
# WHAT WOULD FALSIFY THE REGIME BREAK
#   If the early-period posteriors for spend and cost straddle zero
#   with P(dir.) below ~0.90, the frequentist result was MLE noise on
#   a sparse sample and the regime-break story does not survive. In
#   that case the panel extension is not worth building on this basis.
#
#   If the early posteriors hold their signs with P(dir.) above ~0.95
#   while the late posteriors sit on zero, the break is real and the
#   extension to 2025 is the next step.
#
# HONEST EXPECTATION, STATED BEFORE RUNNING
#   Regularization will shrink the early coefficients toward zero.
#   Spend at +1.58 under MLE will come back smaller. The question is
#   whether it stays clearly positive, not whether it keeps its
#   magnitude. Expect P(dir.) in the 0.85-0.95 range if the effect is
#   real but modest. Anything above 0.99 on 33 time points should be
#   treated with suspicion, not celebration.
# ==============================================================================

.libPaths(c("C:/R/library", .libPaths()))
library(tidyverse)
library(brms)
library(posterior)
library(here)

options(brms.backend = "cmdstanr")
SEED <- 42
set.seed(SEED)

z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

hdi <- function(x, prob = 0.95) {
  x <- sort(x); n <- length(x); k <- floor(prob * n)
  w <- x[(k + 1):n] - x[1:(n - k)]
  i <- which.min(w)
  c(x[i], x[i + k])
}

# ---- Data ----
panel  <- read_csv(here("EIII_Panel_Definitive.csv"), show_col_types = FALSE)
fiscal <- read_csv(here("fiscal_indicators.csv"), show_col_types = FALSE)

panel <- panel |>
  left_join(fiscal |> select(year, net_outlays), by = "year") |>
  mutate(
    eiii_score  = as.integer(replace_na(eiii_score, 0)),
    degradation = as.integer(replace_na(degradation, 0)),
    across(c(friction, resource_density, fiscal_pressure,
             gdp_pc_log, trend, vdem_electoral, net_outlays), as.numeric)
  )

build <- function(data) {
  data |>
    drop_na(vdem_electoral, friction, resource_density,
            gdp_pc_log, fiscal_pressure, net_outlays) |>
    mutate(
      friction_z         = z_score(friction),
      resource_density_z = z_score(resource_density),
      cost_z             = z_score(fiscal_pressure),
      spend_z            = z_score(net_outlays),
      gdp_pc_log_z       = z_score(gdp_pc_log),
      trend_z            = z_score(trend),
      regime_z           = z_score(vdem_electoral),
      fric_x_cost        = friction_z * cost_z
    )
}

d_early <- build(panel |> filter(year <= 2002))
d_late  <- build(panel |> filter(year >= 2003))

cat("Early:", nrow(d_early), "obs,", n_distinct(d_early$country), "states,",
    n_distinct(d_early$year), "years,",
    round(100 * mean(d_early$eiii_score == 0), 1), "% zeros\n")
cat("Late: ", nrow(d_late), "obs,", n_distinct(d_late$country), "states,",
    n_distinct(d_late$year), "years,",
    round(100 * mean(d_late$eiii_score == 0), 1), "% zeros\n\n")

# ---- Formulas ----
f_count <- eiii_score ~ resource_density_z + friction_z + degradation +
                        cost_z + spend_z + trend_z + regime_z +
                        gdp_pc_log_z + fric_x_cost
f_zi <- ~ regime_z + friction_z + resource_density_z + degradation + trend_z

# ---- Priors ----
# Deliberately tighter than the main analysis. On 195 observations with
# 77% zeros, N(0, 2.5) lets coefficients wander to implausible
# magnitudes. N(0, 1.5) on a standardized scale still permits large
# effects (a 1-SD change moving the rate by a factor of e^1.5 = 4.5x)
# while preventing boundary solutions. This is regularization doing
# the job MLE could not do here.
pri <- c(
  prior(normal(0, 2.5), class = "Intercept"),
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 2.5), class = "Intercept", dpar = "zi"),
  prior(normal(0, 1.5), class = "b", dpar = "zi"),
  prior(exponential(1), class = "shape")
)

run <- function(d, label) {
  cat("\n--- Fitting", label, "---\n")
  fit <- brm(
    bf(f_count, zi = f_zi),
    data = d,
    family = zero_inflated_negbinomial(),
    prior = pri,
    chains = 4, iter = 3000, warmup = 1500, cores = 1,
    seed = SEED,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  dg <- summary(fit)
  rh <- max(rhat(fit), na.rm = TRUE)
  nd <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
  cat("  R-hat max:", round(rh, 4), "| Divergences:", nd, "\n")
  if (rh > 1.01) cat("  *** R-hat > 1.01 — do not interpret ***\n")
  if (nd > 10)  cat("  *** Many divergences — posterior unreliable ***\n")
  list(fit = fit, draws = as_draws_df(fit), rhat = rh, div = nd)
}

r_early <- run(d_early, "EARLY 1970-2002")
r_late  <- run(d_late,  "LATE 2003-2021")

saveRDS(list(early = r_early, late = r_late), here("regime_break_fits.rds"))

# ---- The comparison ----
vars <- c("resource_density_z", "friction_z", "degradation",
          "cost_z", "spend_z", "trend_z", "regime_z",
          "gdp_pc_log_z", "fric_x_cost")

summarize_period <- function(r, label) {
  cat("\n==================================================================\n")
  cat(label, "\n")
  cat("==================================================================\n")
  cat(sprintf("%-20s %8s %8s %9s %9s %9s\n",
              "Variable", "Mean", "SD", "HDI Lo", "HDI Hi", "P(dir.)"))
  cat(strrep("-", 70), "\n")
  out <- list()
  for (v in vars) {
    s <- r$draws[[paste0("b_", v)]]
    h <- hdi(s)
    pd <- max(mean(s > 0), mean(s < 0))
    st <- if (pd > 0.99) "**" else if (pd > 0.95) "*" else ""
    cat(sprintf("%-20s %8.3f %8.3f %9.3f %9.3f %9.3f %s\n",
                v, mean(s), sd(s), h[1], h[2], pd, st))
    out[[v]] <- tibble(variable = v, mean = mean(s), sd = sd(s),
                       hdi_lo = h[1], hdi_hi = h[2], p_dir = pd)
  }
  bind_rows(out)
}

t_early <- summarize_period(r_early, "EARLY 1970-2002 — POSTERIOR")
t_late  <- summarize_period(r_late,  "LATE 2003-2021 — POSTERIOR")

# ---- Side by side on the two fiscal mechanisms ----
cat("\n==================================================================\n")
cat("THE REGIME BREAK: SPEND AND COST ACROSS PERIODS\n")
cat("==================================================================\n")
cat(sprintf("\n%-10s %-8s %9s %9s %9s %9s\n",
            "Period", "Var", "Mean", "HDI Lo", "HDI Hi", "P(dir.)"))
cat(strrep("-", 60), "\n")
for (p in c("Early", "Late")) {
  r <- if (p == "Early") r_early else r_late
  for (v in c("spend_z", "cost_z")) {
    s <- r$draws[[paste0("b_", v)]]
    h <- hdi(s)
    cat(sprintf("%-10s %-8s %9.3f %9.3f %9.3f %9.3f\n",
                p, sub("_z$", "", v), mean(s), h[1], h[2],
                max(mean(s > 0), mean(s < 0))))
  }
}

# ---- Direct posterior comparison of the difference ----
# Does the spend effect actually DIFFER across periods, or are both
# just imprecise? Two independent fits, so this is a rough comparison
# of marginal posteriors, not a formal interaction test.
cat("\n--- Posterior difference, early minus late (rough) ---\n")
for (v in c("spend_z", "cost_z")) {
  e <- r_early$draws[[paste0("b_", v)]]
  l <- r_late$draws[[paste0("b_", v)]]
  n <- min(length(e), length(l))
  d <- sample(e, n) - sample(l, n)
  h <- hdi(d)
  cat(sprintf("%-8s diff = %7.3f  HDI [%.3f, %.3f]  P(diff>0) = %.3f\n",
              sub("_z$", "", v), mean(d), h[1], h[2], mean(d > 0)))
}
cat("\n(These are separate models, not a pooled interaction. Treat the\n")
cat(" difference as descriptive. A formal test requires one model with\n")
cat(" period interacted on every coefficient — more parameters than\n")
cat(" 195 early observations will support.)\n")

# ---- Verdict ----
cat("\n==================================================================\n")
cat("VERDICT\n")
cat("==================================================================\n")
sp_e <- r_early$draws$b_spend_z; co_e <- r_early$draws$b_cost_z
sp_l <- r_late$draws$b_spend_z;  co_l <- r_late$draws$b_cost_z
pd <- function(x) max(mean(x > 0), mean(x < 0))

early_live <- pd(sp_e) > 0.90 && mean(sp_e) > 0 && pd(co_e) > 0.90 && mean(co_e) < 0
late_dead  <- pd(sp_l) < 0.90 && pd(co_l) < 0.90

cat("Early spend P(dir.):", round(pd(sp_e), 3),
    "| sign:", ifelse(mean(sp_e) > 0, "+", "-"), "\n")
cat("Early cost  P(dir.):", round(pd(co_e), 3),
    "| sign:", ifelse(mean(co_e) > 0, "+", "-"), "\n")
cat("Late spend  P(dir.):", round(pd(sp_l), 3), "\n")
cat("Late cost   P(dir.):", round(pd(co_l), 3), "\n\n")

if (early_live && late_dead) {
  cat("REGIME BREAK SURVIVES. Both fiscal mechanisms operate in the\n")
  cat("early period with correct signs; neither operates in the late\n")
  cat("period. The panel extension to 2025 is worth building.\n")
} else if (!early_live) {
  cat("REGIME BREAK DOES NOT SURVIVE. The early-period effects do not\n")
  cat("hold under regularization. The frequentist result was MLE noise\n")
  cat("on a sparse sample. Do not build the extension on this basis.\n")
} else {
  cat("PARTIAL. Report both periods plainly; the break is not clean.\n")
}

write_csv(bind_rows(t_early |> mutate(period = "early"),
                    t_late  |> mutate(period = "late")),
          here("regime_break_posteriors.csv"))
cat("\nSaved: regime_break_posteriors.csv, regime_break_fits.rds\n")
