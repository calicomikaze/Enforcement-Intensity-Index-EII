# ==============================================================================
# ENFORCEMENT INTENSITY INDEX (EII)
# Bayesian Zero-Inflated Negative Binomial Estimation
# of Enforcement Intensity in Strategically Produced Data
#
# Sean Rogers | University of South Carolina | 2026
# R port of QC_Full_EIII_v5.py (PyMC) — results are reproducible across platforms
#
# INPUT:   EIII_Panel_Definitive.csv (942 country-years, 44 states, 1945-2025)
# OUTPUT:  Results tables (console + CSV), figures (PNG), posterior draws
# RUNTIME: ~10-15 min (Stan compilation + 6 model fits at 4 chains each)
#
# ------------------------------------------------------------------------------
# WHY BAYESIAN?
# The institutions producing enforcement data have survival incentives to
# misreport. MLE assumes the DGP is honest. When the DGP is strategically
# corrupted, the likelihood surface is jagged (extreme overdispersion, zero
# mass) and frequentist estimation plants a flag on whichever near-tied peak
# it reaches first. The Bayesian posterior averages the whole surface.
# Regularization through priors is a shock absorber, not a bias.
#
# WHY ZINB?
# The DV has two sources of zeros:
#   Sampling zeros:   the count process generates zero (low intensity)
#   Structural zeros: the state is shielded from enforcement (sanctuary)
# These are qualitatively different. The ZINB separates them: the count
# equation models intensity; the zero-inflation (zi) equation models
# the probability of a structural zero.
#
# CONVENTION NOTE (critical):
# In both pscl::zeroinfl() and brms::zero_inflated_negbinomial(), the
# zero-inflation component models P(structural zero). A POSITIVE coefficient
# in the zi equation means MORE sanctuary. This matches the theoretical
# interpretation directly. (PyMC's ZeroInflatedNegativeBinomial psi parameter
# models P(count process) — the opposite convention. The Python v5 inflate
# equation signs read inverted. The R version does not have this problem.)
# ==============================================================================


# ==============================================================================
# SECTION 1: PACKAGES
# ==============================================================================
# Install once (run in the R Console if you haven't already):
#   install.packages(c("pscl", "brms", "loo", "bayesplot", "posterior"))
# brms needs a Stan backend. On Windows, install RTools first, then:
#   install.packages("rstan")
# or (faster, recommended):
#   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev"))
#   cmdstanr::install_cmdstan()

# Point R at the local library (OneDrive blocks DLL loads)
.libPaths(c("C:/R/library", .libPaths()))

library(tidyverse)      # data manipulation and plotting
library(pscl)           # frequentist zero-inflated models (zeroinfl)
library(brms)           # Bayesian regression via Stan
library(loo)            # leave-one-out cross-validation for model comparison
library(bayesplot)      # MCMC diagnostic plots
library(posterior)      # working with posterior draws
library(here)           # project-relative file paths

set.seed(42)
SEED <- 42

# Set brms backend. Use "cmdstanr" if installed (faster), otherwise "rstan".
options(brms.backend = "cmdstanr")
# options(brms.backend = "rstan")  # uncomment if cmdstanr isn't installed

# ---- Highest density interval ----
# The shortest interval containing prob% of the posterior mass.
# Unlike a quantile interval (equal tails), the HDI can be asymmetric.
hdi <- function(x, prob = 0.95) {
  x <- sort(x)
  n <- length(x)
  k <- floor(prob * n)
  widths <- x[(k + 1):n] - x[1:(n - k)]
  i <- which.min(widths)
  c(x[i], x[i + k])
}

cat("==================================================================\n")
cat("EII REPLICATION PIPELINE — R PORT\n")
cat("brms", as.character(packageVersion("brms")), "|",
    "pscl", as.character(packageVersion("pscl")), "\n")
cat("==================================================================\n")


# ==============================================================================
# SECTION 2: DATA PREPARATION
# ==============================================================================
# The panel is pre-built from ten institutional primary sources across four
# layers. This script takes the panel as given and estimates the model.
# See codebook for construction details.

df <- read_csv(here("EIII_Panel_Definitive.csv"), show_col_types = FALSE)

# ---- Type enforcement ----
df <- df |>
  mutate(
    eiii_score    = as.integer(replace_na(eiii_score, 0)),
    degradation   = as.integer(replace_na(degradation, 0)),
    across(c(polity, friction, resource_density, fiscal_pressure,
             gdp_pc_log, trend, vdem_electoral, vdem_elections),
           as.numeric)
  )

# ---- Z-scoring function ----
# Every continuous IV is standardized: (x - mean) / sd.
# This puts all coefficients on the same scale (one-SD change) so
# magnitudes are comparable across variables and priors are meaningful.
z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# ---- Build estimation samples ----
# PRIMARY: V-Dem Electoral Democracy Index
# ADDITIONAL: V-Dem Free and Fair Elections Index
# ROBUSTNESS: Polity5
# Each sample drops observations missing any IV, then z-scores within the
# sample so coefficients reflect the data actually used.

build_sample <- function(data, regime_var) {
  data |>
    drop_na(all_of(c(regime_var, "friction", "resource_density",
                     "gdp_pc_log", "fiscal_pressure"))) |>
    mutate(
      friction_z         = z_score(friction),
      resource_density_z = z_score(resource_density),
      fiscal_pressure_z  = z_score(fiscal_pressure),
      gdp_pc_log_z       = z_score(gdp_pc_log),
      trend_z            = z_score(trend),
      regime_z           = z_score(.data[[regime_var]]),
      fric_x_fiscal      = friction_z * fiscal_pressure_z
    )
}

df_v <- build_sample(df, "vdem_electoral")
df_e <- build_sample(df, "vdem_elections")
df_p <- build_sample(df, "polity")

N <- nrow(df_v)
zero_pct <- 100 * mean(df_v$eiii_score == 0)
od <- var(df_v$eiii_score) / max(mean(df_v$eiii_score), 1e-10)

cat("\nPanel:", nrow(df), "total |",
    "Estimation sample:", N, "complete cases,",
    n_distinct(df_v$country), "states\n")
cat("DV:", round(zero_pct, 1), "% zeros | Var/mean:", round(od, 1), "\n")


# ==============================================================================
# SECTION 3: MODEL FORMULAS
# ==============================================================================
# Count equation: what drives enforcement intensity given engagement?
# Zero-inflation equation: what predicts structural zeros (sanctuary)?
#
# The interaction (friction × fiscal) is the Path of Least Resistance test.
# fiscal_pressure_z appears in the count equation only — the theory says
# fiscal strain affects intensity, not whether a state is shielded.

f_count <- eiii_score ~ resource_density_z + friction_z + degradation +
                        fiscal_pressure_z + trend_z + regime_z +
                        gdp_pc_log_z + fric_x_fiscal

f_zi <- ~ regime_z + friction_z + resource_density_z + degradation + trend_z

# Variable names for output tables (count equation, in formula order)
count_names <- c("resource", "friction", "degradation", "fiscal",
                 "trend", "regime", "gdp_pc", "friction_x_fiscal")
inflate_names <- c("regime", "friction", "resource", "degradation", "trend")


# ==============================================================================
# SECTION 4: FREQUENTIST ZINB — BASELINE
# ==============================================================================
# Why run frequentist first? Three reasons:
# 1. Establishes that the data support ZINB over simpler alternatives.
# 2. Provides a convergence benchmark.
# 3. Bayesian posteriors should bracket the MLE point estimates. Disagreement
#    flags prior-data conflict or likelihood multimodality.

run_freq_zinb <- function(data, label) {
  cat("\n==================================================================\n")
  cat("FREQUENTIST ZINB:", label, "\n")
  cat("==================================================================\n")

  # pscl::zeroinfl uses formula syntax: count_formula | zi_formula
  # dist = "negbin" selects the negative binomial for the count component
  fit <- tryCatch(
    zeroinfl(
      eiii_score ~ resource_density_z + friction_z + degradation +
                   fiscal_pressure_z + trend_z + regime_z +
                   gdp_pc_log_z + fric_x_fiscal |
                   regime_z + friction_z + resource_density_z +
                   degradation + trend_z,
      data = data,
      dist = "negbin"
    ),
    error = function(e) {
      cat("  FAILED:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(fit)) return(NULL)

  s <- summary(fit)
  cat("Converged:", fit$converged, "\n")
  cat("Log-likelihood:", round(fit$loglik, 2),
      "| AIC:", round(AIC(fit), 2), "\n")

  cat("\nCOUNT EQUATION:\n")
  print(round(s$coefficients$count, 4))

  cat("\nZERO-INFLATION EQUATION (positive = more sanctuary):\n")
  print(round(s$coefficients$zero, 4))

  fit
}

freq_vdem <- run_freq_zinb(df_v, "V-DEM ELECTORAL (PRIMARY)")
freq_pol  <- run_freq_zinb(df_p, "POLITY5 (ROBUSTNESS)")

if (!is.null(freq_vdem)) {
  fx <- coef(freq_vdem)["count_fric_x_fiscal"]
  fx_p <- summary(freq_vdem)$coefficients$count["fric_x_fiscal", "Pr(>|z|)"]
  cat("\n>>> Frequentist PLR: b =", round(fx, 4),
      "(p =", round(fx_p, 6), ")\n")
}


# ==============================================================================
# SECTION 5: BAYESIAN ZINB — PRIOR SPECIFICATIONS
# ==============================================================================
# Four specifications test whether substantive conclusions depend on prior
# choice. All share the same model structure; they differ only in prior
# width and interaction prior center.

# ZI intercept: logit of observed zero rate as an upper-bound anchor.
# The data will pull it toward the true structural zero rate.
ZI_BASE <- qlogis(zero_pct / 100)
cat("\nZI intercept prior mean: logit(", round(zero_pct, 1), "%) =",
    round(ZI_BASE, 3), "\n")

# Prior specification builder
# brms prior syntax: prior(distribution, class, coef, dpar)
#   class = "b"         → regression coefficients (count equation by default)
#   dpar = "zi"         → the zero-inflation distributional parameter
#   class = "Intercept" → the intercept
#   class = "shape"     → NegBin overdispersion (brms calls alpha "shape")
build_priors <- function(beta_sigma, gamma_sigma, gamma0_sigma,
                         interaction_mu, interaction_sigma) {
  c(
    # Count equation intercept
    prior(normal(0, 2.5), class = "Intercept"),
    # Count equation main effects
    set_prior(paste0("normal(0, ", beta_sigma, ")"), class = "b"),
    # Count equation interaction — gets its own prior
    set_prior(paste0("normal(", interaction_mu, ", ", interaction_sigma, ")"),
              class = "b", coef = "fric_x_fiscal"),
    # Zero-inflation intercept — anchored at logit(zero rate)
    set_prior(paste0("normal(", ZI_BASE, ", ", gamma0_sigma, ")"),
              class = "Intercept", dpar = "zi"),
    # Zero-inflation coefficients
    set_prior(paste0("normal(0, ", gamma_sigma, ")"),
              class = "b", dpar = "zi"),
    # Overdispersion: Exponential(1) — mean 1, heavy right tail
    prior(exponential(1), class = "shape")
  )
}

specs <- list(
  Spec1_Diffuse = list(
    description = "Weakly informative (Gelman 2008)",
    priors = build_priors(2.5, 2.5, 2.5, 0, 2.5)
  ),
  Spec2_Sequential = list(
    description = "Sequential — proxy-panel posteriors as priors",
    priors = build_priors(2.5, 2.5, 2.5, -0.1, 1.0)
  ),
  Spec3_Skeptical = list(
    description = "Skeptical — shrinks interaction toward zero",
    priors = build_priors(1.5, 1.5, 2.5, 0, 0.5)
  ),
  Spec4_Restrictive = list(
    description = "Restrictive — N(0,1) throughout",
    priors = build_priors(1.0, 1.0, 2.0, 0, 1.0)
  )
)


# ==============================================================================
# SECTION 6: BAYESIAN ESTIMATION FUNCTION
# ==============================================================================
run_bayesian_zinb <- function(data, label, priors, description = "") {
  cat("\n---", label, ":", description, "---\n")

  fit <- brm(
    bf(f_count, zi = f_zi),
    data = data,
    family = zero_inflated_negbinomial(),
    prior = priors,
    chains = 4, iter = 2000, warmup = 1000,
    cores = 1,
    seed = SEED,
    control = list(adapt_delta = 0.90)
  )

  # ---- Convergence diagnostics ----
  rhats <- rhat(fit)
  ess_b <- neff_ratio(fit) * (4 * 1000)
  nd <- nuts_params(fit) |>
    filter(Parameter == "divergent__") |>
    pull(Value) |>
    sum()

  rhat_max <- max(rhats, na.rm = TRUE)
  ess_min  <- min(ess_b, na.rm = TRUE)

  cat("  R-hat max:", round(rhat_max, 4),
      "| ESS min:", round(ess_min, 0),
      "| Divergences:", nd, "\n")
  if (rhat_max > 1.01) cat("  *** WARNING: R-hat > 1.01 ***\n")

  # ---- PLR interaction ----
  draws <- as_draws_df(fit)
  plr <- draws$b_fric_x_fiscal
  cat("  PLR: mean =", round(mean(plr), 4),
      ", P(<0) =", round(mean(plr < 0), 4), "\n")

  list(fit = fit, draws = draws, plr = plr,
       rhat = rhat_max, ess = ess_min, div = nd)
}


# ==============================================================================
# SECTION 7: RUN ALL SPECIFICATIONS
# ==============================================================================
cat("\n==================================================================\n")
cat("BAYESIAN ZINB — V-DEM PRIMARY — FOUR PRIOR SPECIFICATIONS\n")
cat("4 chains × 2000 iterations (1000 warmup) | seed =", SEED, "\n")
cat("==================================================================\n")

results <- list()
for (nm in names(specs)) {
  results[[nm]] <- run_bayesian_zinb(df_v, nm, specs[[nm]]$priors,
                                     specs[[nm]]$description)
}

r1 <- results$Spec1_Diffuse

# V-Dem Free and Fair Elections
cat("\n==================================================================\n")
cat("ADDITIONAL: V-DEM FREE AND FAIR ELECTIONS\n")
cat("==================================================================\n")
r_elec <- run_bayesian_zinb(df_e, "VDem_Elections",
                            specs$Spec1_Diffuse$priors,
                            "V-Dem elections — weakly informative")

# Save fits now — before any results printing — so a downstream error
# never costs the 10 minutes of sampling
saveRDS(results, here("brms_fits_main.rds"))
saveRDS(r_elec, here("brms_fit_elections.rds"))
cat("\nFits saved to brms_fits_main.rds and brms_fit_elections.rds\n")


# ==============================================================================
# SECTION 8: RESULTS TABLES
# ==============================================================================
# Posterior summary: mean, SD, 95% HDI, probability of direction.
# P(dir.) = max(P(β > 0), P(β < 0)) — the posterior probability that the
# coefficient has the sign of its mean. Stars mark P(dir.) > 0.95/0.99/0.999.

post_table <- function(draws, var_names, prefix, title) {
  cat("\n", title, "\n")
  cat(sprintf("%-20s %8s %8s %8s %8s %8s\n",
              "Variable", "Mean", "SD", "HDI Lo", "HDI Hi", "P(dir.)"))
  cat(strrep("-", 70), "\n")

  rows <- list()
  for (v in var_names) {
    col <- paste0(prefix, v)
    if (!col %in% names(draws)) next
    s <- draws[[col]]
    h <- hdi(s, prob = 0.95)
    pd <- max(mean(s > 0), mean(s < 0))
    stars <- if (pd > 0.999) "***" else if (pd > 0.99) "**" else if (pd > 0.95) "*" else ""
    cat(sprintf("%-20s %8.3f %8.3f %8.3f %8.3f %8.3f %s\n",
                v, mean(s), sd(s), h[1], h[2], pd, stars))
    rows[[v]] <- tibble(variable = v, mean = mean(s), sd = sd(s),
                        hdi_lo = h[1], hdi_hi = h[2], p_dir = pd,
                        p_positive = mean(s > 0), p_negative = mean(s < 0))
  }
  bind_rows(rows)
}

# Map formula variable names to the column names brms produces
# brms names coefficients b_<varname> for count, b_zi_<varname> for zi
count_cols <- c("resource_density_z", "friction_z", "degradation",
                "fiscal_pressure_z", "trend_z", "regime_z",
                "gdp_pc_log_z", "fric_x_fiscal")
zi_cols <- c("regime_z", "friction_z", "resource_density_z",
             "degradation", "trend_z")

print_results <- function(r, label) {
  cat("\n==================================================================\n")
  cat(label, "\n")
  cat("==================================================================\n")

  ct <- post_table(r$draws, count_cols, "b_",
                   "COUNT EQUATION (Enforcement Intensity):")
  ct$variable <- count_names[match(ct$variable, count_cols)]
  ct$equation <- "count"

  zt <- post_table(r$draws, zi_cols, "b_zi_",
                   "ZERO-INFLATION EQUATION (Sanctuary — positive = more sanctuary):")
  zt$variable <- inflate_names[match(zt$variable, zi_cols)]
  zt$equation <- "inflate"

  # Overdispersion
  shape <- r$draws$shape
  hs <- hdi(shape, prob = 0.95)
  cat("\nOverdispersion (shape):", round(mean(shape), 3),
      "[", round(hs[1], 3), ",", round(hs[2], 3), "]\n")

  bind_rows(ct, zt)
}

tab_primary <- print_results(r1, "PRIMARY RESULTS — V-DEM ELECTORAL (DIFFUSE PRIORS)")
tab_elec    <- print_results(r_elec, "ADDITIONAL — V-DEM FREE AND FAIR ELECTIONS")


# ==============================================================================
# SECTION 9: PRIOR SENSITIVITY
# ==============================================================================
cat("\n==================================================================\n")
cat("PRIOR SENSITIVITY — FRICTION × FISCAL INTERACTION\n")
cat("==================================================================\n")
cat(sprintf("\n%-22s %8s %7s %7s %7s %7s\n",
            "Spec", "Mean", "SD", "HDI Lo", "HDI Hi", "P(<0)"))
cat(strrep("-", 62), "\n")
for (nm in names(results)) {
  plr <- results[[nm]]$plr
  h <- hdi(plr, prob = 0.95)
  short <- sub("^Spec[0-9]_", "", nm)
  cat(sprintf("%-22s %8.4f %7.4f %7.3f %7.3f %7.4f\n",
              short, mean(plr), sd(plr), h[1], h[2], mean(plr < 0)))
}

cat(sprintf("\n%-22s %9s %9s %9s %9s\n",
            "Spec", "Friction", "Fiscal", "Trend", "Regime"))
cat(strrep("-", 60), "\n")
for (nm in names(results)) {
  d <- results[[nm]]$draws
  pd <- function(x) max(mean(x > 0), mean(x < 0))
  short <- sub("^Spec[0-9]_", "", nm)
  cat(sprintf("%-22s %9.3f %9.3f %9.3f %9.3f\n", short,
              pd(d$b_friction_z), pd(d$b_fiscal_pressure_z),
              pd(d$b_trend_z), pd(d$b_regime_z)))
}


# ==============================================================================
# SECTION 10: CONVERGENT VALIDATION — BAYESIAN vs. FREQUENTIST
# ==============================================================================
# Both equations are on the same convention in R (P(structural zero)),
# so the inflate-side comparison is valid here — unlike the Python v5
# where PyMC and statsmodels used opposite conventions.
cat("\n==================================================================\n")
cat("CONVERGENT VALIDATION: BAYESIAN vs. FREQUENTIST (V-DEM)\n")
cat("==================================================================\n")

if (!is.null(freq_vdem)) {
  fc <- coef(freq_vdem)
  fs <- summary(freq_vdem)$coefficients

  cat("\nCOUNT EQUATION:\n")
  cat(sprintf("%-20s %9s %9s %9s %8s %6s\n",
              "Variable", "Freq b", "Freq p", "Bayes", "P(dir.)", "Sign"))
  cat(strrep("-", 65), "\n")
  for (i in seq_along(count_cols)) {
    v <- count_cols[i]
    f_b <- fc[paste0("count_", v)]
    f_p <- fs$count[v, "Pr(>|z|)"]
    b_s <- r1$draws[[paste0("b_", v)]]
    b_m <- mean(b_s)
    b_pd <- max(mean(b_s > 0), mean(b_s < 0))
    ag <- if (sign(f_b) == sign(b_m)) "Y" else "N"
    cat(sprintf("%-20s %9.4f %9.5f %9.4f %8.4f %6s\n",
                count_names[i], f_b, f_p, b_m, b_pd, ag))
  }

  cat("\nZERO-INFLATION EQUATION (same convention both sides):\n")
  cat(sprintf("%-20s %9s %9s %9s %8s %6s\n",
              "Variable", "Freq b", "Freq p", "Bayes", "P(dir.)", "Sign"))
  cat(strrep("-", 65), "\n")
  for (i in seq_along(zi_cols)) {
    v <- zi_cols[i]
    f_b <- fc[paste0("zero_", v)]
    f_p <- fs$zero[v, "Pr(>|z|)"]
    b_s <- r1$draws[[paste0("b_zi_", v)]]
    b_m <- mean(b_s)
    b_pd <- max(mean(b_s > 0), mean(b_s < 0))
    ag <- if (sign(f_b) == sign(b_m)) "Y" else "N"
    cat(sprintf("%-20s %9.4f %9.5f %9.4f %8.4f %6s\n",
                inflate_names[i], f_b, f_p, b_m, b_pd, ag))
  }
}


# ==============================================================================
# SECTION 11: REGIME VARIABLE COMPARISON
# ==============================================================================
cat("\n==================================================================\n")
cat("REGIME VARIABLE COMPARISON\n")
cat("==================================================================\n")
cat(sprintf("\n%-20s %12s %12s %12s\n",
            "Variable", "V-Dem Elect.", "V-Dem Elec.", "Polity5"))
cat(strrep("-", 60), "\n")
for (i in seq_along(count_cols)) {
  v <- count_cols[i]
  vm <- mean(r1$draws[[paste0("b_", v)]])
  em <- mean(r_elec$draws[[paste0("b_", v)]])
  pm <- if (!is.null(freq_pol)) coef(freq_pol)[paste0("count_", v)] else NA
  cat(sprintf("%-20s %12.3f %12.3f %12.3f\n", count_names[i], vm, em, pm))
}


# ==============================================================================
# SECTION 12: MODEL COMPARISON — LOO-CV
# ==============================================================================
cat("\n==================================================================\n")
cat("MODEL COMPARISON — LOO-CV\n")
cat("==================================================================\n")

loo_list <- lapply(results, function(r) loo(r$fit))
names(loo_list) <- names(results)
print(loo_compare(loo_list))


# ==============================================================================
# SECTION 13: MCMC DIAGNOSTICS
# ==============================================================================
cat("\n==================================================================\n")
cat("MCMC DIAGNOSTICS — ALL SPECIFICATIONS\n")
cat("==================================================================\n")
cat(sprintf("\n%-22s %8s %8s %8s\n", "Spec", "R-hat", "ESS", "Diverg."))
cat(strrep("-", 48), "\n")
all_r <- c(results, list(VDem_Elections = r_elec))
for (nm in names(all_r)) {
  r <- all_r[[nm]]
  short <- sub("^Spec[0-9]_", "", nm)
  cat(sprintf("%-22s %8.4f %8.0f %8d\n", short, r$rhat, r$ess, r$div))
}


# ==============================================================================
# SECTION 14: POSTERIOR PREDICTIVE CHECK
# ==============================================================================
cat("\n==================================================================\n")
cat("POSTERIOR PREDICTIVE CHECK — PRIMARY SPECIFICATION\n")
cat("==================================================================\n")

ppc_draws <- posterior_predict(r1$fit, ndraws = 500)
ppc_zero <- 100 * rowMeans(ppc_draws == 0)
ppc_mean <- rowMeans(ppc_draws)
ppc_var  <- apply(ppc_draws, 1, var)
ppc_od   <- ppc_var / pmax(ppc_mean, 1e-10)

cat(sprintf("%-20s %10s %10s %20s\n",
            "Statistic", "Observed", "PPC Mean", "PPC 95% CI"))
cat(strrep("-", 65), "\n")
cat(sprintf("%-20s %10.1f %10.1f [%.1f, %.1f]\n",
            "Zero rate (%)", zero_pct, mean(ppc_zero),
            quantile(ppc_zero, 0.025), quantile(ppc_zero, 0.975)))
cat(sprintf("%-20s %10.2f %10.2f [%.2f, %.2f]\n",
            "Mean", mean(df_v$eiii_score), mean(ppc_mean),
            quantile(ppc_mean, 0.025), quantile(ppc_mean, 0.975)))
cat(sprintf("%-20s %10.1f %10.1f [%.1f, %.1f]\n",
            "Var/Mean", od, mean(ppc_od),
            quantile(ppc_od, 0.025), quantile(ppc_od, 0.975)))


# ==============================================================================
# SECTION 15: DV WEIGHTING SENSITIVITY (FREQUENTIST)
# ==============================================================================
cat("\n==================================================================\n")
cat("DV WEIGHTING SENSITIVITY — FREQUENTIST\n")
cat("==================================================================\n")

dv_specs <- c(
  "Weighted (primary)" = "eiii_score",
  "Equal weight"       = "eiii_equal",
  "Raw count (capped)" = "eiii_rawcount",
  "Financial only"     = "eiii_financial"
)

cat(sprintf("\n%-25s %9s %9s %9s %9s %10s\n",
            "DV", "Friction", "Fiscal", "Trend", "FxF", "FxF p"))
cat(strrep("-", 75), "\n")

wt_rows <- list()
for (nm in names(dv_specs)) {
  dv <- dv_specs[nm]
  d_wt <- df_v |> mutate(y_wt = as.integer(.data[[dv]]))
  fit_wt <- tryCatch(
    zeroinfl(
      y_wt ~ resource_density_z + friction_z + degradation +
             fiscal_pressure_z + trend_z + regime_z +
             gdp_pc_log_z + fric_x_fiscal |
             regime_z + friction_z + resource_density_z +
             degradation + trend_z,
      data = d_wt, dist = "negbin"
    ),
    error = function(e) NULL
  )
  if (is.null(fit_wt) || !fit_wt$converged) {
    cat(sprintf("%-25s DID NOT CONVERGE\n", nm))
    next
  }
  cc <- coef(fit_wt)
  cs <- summary(fit_wt)$coefficients$count
  cat(sprintf("%-25s %9.3f %9.3f %9.3f %9.3f %10.5f\n", nm,
              cc["count_friction_z"], cc["count_fiscal_pressure_z"],
              cc["count_trend_z"], cc["count_fric_x_fiscal"],
              cs["fric_x_fiscal", "Pr(>|z|)"]))
  wt_rows[[nm]] <- tibble(
    dv_spec = nm,
    friction = cc["count_friction_z"], fiscal = cc["count_fiscal_pressure_z"],
    trend = cc["count_trend_z"], fxf = cc["count_fric_x_fiscal"],
    fxf_p = cs["fric_x_fiscal", "Pr(>|z|)"]
  )
}
if (length(wt_rows) > 0) {
  write_csv(bind_rows(wt_rows), here("weighting_sensitivity_R.csv"))
}


# ==============================================================================
# SECTION 16: FIGURES
# ==============================================================================
cat("\n--- Generating figures ---\n")
color_scheme_set("gray")

# Fig 1: Count equation posteriors
p1 <- mcmc_areas(r1$fit, pars = paste0("b_", count_cols), prob = 0.95) +
  scale_y_discrete(labels = rev(count_names)) +
  labs(title = "Count Equation Posteriors — V-Dem Electoral (Diffuse Priors)",
       subtitle = "Shaded: 95% HDI. Dashed line: zero.") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal(base_family = "serif")
ggsave(here("fig1_count_posteriors_R.png"), p1, width = 10, height = 6, dpi = 300)

# Fig 2: Zero-inflation posteriors
p2 <- mcmc_areas(r1$fit, pars = paste0("b_zi_", zi_cols), prob = 0.95) +
  scale_y_discrete(labels = rev(inflate_names)) +
  labs(title = "Zero-Inflation (Sanctuary) Posteriors — V-Dem Electoral",
       subtitle = "Positive = more sanctuary") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal(base_family = "serif")
ggsave(here("fig2_sanctuary_posteriors_R.png"), p2, width = 10, height = 5, dpi = 300)

# Fig 3: PLR interaction
plr_df <- tibble(plr = r1$plr)
h_plr <- hdi(r1$plr, prob = 0.95)
p3 <- ggplot(plr_df, aes(x = plr)) +
  geom_histogram(bins = 60, fill = "#c0392b", alpha = 0.85, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  annotate("rect", xmin = h_plr[1], xmax = h_plr[2], ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "blue") +
  labs(title = "Friction × Fiscal Pressure Interaction",
       subtitle = sprintf("Mean = %.4f | P(<0) = %.4f | 95%% HDI [%.3f, %.3f]",
                          mean(r1$plr), mean(r1$plr < 0), h_plr[1], h_plr[2]),
       x = "Coefficient (standardized)", y = "Posterior density") +
  theme_minimal(base_family = "serif")
ggsave(here("fig3_interaction_R.png"), p3, width = 8, height = 5, dpi = 300)

# Fig 4: Prior sensitivity
plr_all <- bind_rows(lapply(names(results), function(nm) {
  tibble(spec = sub("^Spec[0-9]_", "", nm), plr = results[[nm]]$plr)
}))
plr_all$spec <- factor(plr_all$spec, levels = c("Diffuse", "Sequential",
                                                 "Skeptical", "Restrictive"))
p4 <- ggplot(plr_all, aes(x = plr, fill = spec)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~spec, ncol = 4) +
  labs(title = "Prior Sensitivity — Friction × Fiscal Interaction",
       x = "Coefficient", y = "Density") +
  theme_minimal(base_family = "serif") +
  theme(legend.position = "none")
ggsave(here("fig4_prior_sensitivity_R.png"), p4, width = 14, height = 4, dpi = 300)

# Fig 5: Trace plots — full model
p5 <- mcmc_trace(r1$fit,
                 pars = c(paste0("b_", count_cols), paste0("b_zi_", zi_cols),
                          "b_Intercept", "b_zi_Intercept", "shape"),
                 facet_args = list(ncol = 3)) +
  theme_minimal(base_family = "serif")
ggsave(here("fig5_traces_R.png"), p5, width = 14, height = 14, dpi = 300)

# Fig 6: PPC — zero rate
ppc_df <- tibble(zero_rate = ppc_zero)
p6 <- ggplot(ppc_df, aes(x = zero_rate)) +
  geom_histogram(bins = 50, fill = "#2c3e50", alpha = 0.8, color = "white") +
  geom_vline(xintercept = zero_pct, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Posterior Predictive Check — Zero-Inflation Recovery",
       subtitle = sprintf("Observed: %.1f%%", zero_pct),
       x = "Zero rate (%)", y = "Count") +
  theme_minimal(base_family = "serif")
ggsave(here("fig6_ppc_zeros_R.png"), p6, width = 8, height = 5, dpi = 300)

# Fig 7: PPC — distribution overlay
p7 <- pp_check(r1$fit, ndraws = 100, type = "hist") +
  labs(title = "Posterior Predictive Check — Count Distribution") +
  theme_minimal(base_family = "serif")
ggsave(here("fig7_ppc_distribution_R.png"), p7, width = 8, height = 5, dpi = 300)

cat("  Saved: fig1–fig7 (R versions)\n")


# ==============================================================================
# SECTION 17: SAVE RESULTS
# ==============================================================================
cat("\n--- Saving results ---\n")

# All posterior summaries
all_tabs <- bind_rows(lapply(names(all_r), function(nm) {
  r <- all_r[[nm]]
  ct <- post_table(r$draws, count_cols, "b_", "")
  ct$variable <- count_names[match(ct$variable, count_cols)]
  ct$equation <- "count"
  zt <- post_table(r$draws, zi_cols, "b_zi_", "")
  zt$variable <- inflate_names[match(zt$variable, zi_cols)]
  zt$equation <- "inflate"
  bind_rows(ct, zt) |> mutate(spec = nm, .before = 1)
}))
write_csv(all_tabs, here("bayesian_results_all_specs_R.csv"))

# Estimation sample
write_csv(df_v, here("panel_estimation_sample_R.csv"))

# Save fitted models (large — .gitignore these)
saveRDS(results, here("brms_fits_main.rds"))
saveRDS(r_elec, here("brms_fit_elections.rds"))

cat("Saved: bayesian_results_all_specs_R.csv\n")
cat("Saved: panel_estimation_sample_R.csv\n")
cat("Saved: brms_fits_main.rds, brms_fit_elections.rds\n")

cat("\n==================================================================\n")
cat("PIPELINE COMPLETE\n")
cat("==================================================================\n")
cat("  Panel:", nrow(df), "total,", N, "estimation,",
    n_distinct(df_v$country), "states\n")
cat("  Primary: V-Dem Electoral | Robustness: V-Dem Elections, Polity5\n")
cat("  Bayesian: 4 prior specs + elections spec | Frequentist baseline\n")
cat("  DV sensitivity: 4 weighting schemes\n")
cat("  Seed:", SEED, "\n")
