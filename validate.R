# ==============================================================================
# VALIDATION HARNESS
# ==============================================================================
# WHY THIS EXISTS
# On 2026-09-02 a parameter-ordering bug in the Python pipeline produced
# a frequentist coefficient table in which every variable was labeled
# with another variable's estimate. statsmodels' ZeroInflatedNegativeBinomialP
# returns parameters as [inflate, count, alpha]; the extraction code
# assumed [count, alpha, inflate]. The error survived three weeks and
# a full round of paper planning because nothing checked the output
# against an independent implementation. It was caught on 2026-10-05
# only because the model was ported to R.
#
# This script is the check that would have caught it in sixty seconds.
# It re-derives the core quantities from the panel and asserts they
# match values verified on 2026-10-05 against two independent
# implementations (pscl + brms in R, statsmodels + PyMC in Python).
#
# RUN THIS:
#   - before sending results to anyone
#   - after any change to the panel
#   - after any change to the estimation code
#   - after any package or R version upgrade
#
# A failure does not mean the new result is wrong. It means something
# changed and you must find out what before the number leaves the lab.
# ==============================================================================

.libPaths(c("C:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(tidyverse)
  library(pscl)
  library(here)
})

PASS <- 0; FAIL <- 0

check <- function(label, actual, expected, tol = 0.01) {
  ok <- !is.na(actual) && abs(actual - expected) <= tol
  if (ok) {
    PASS <<- PASS + 1
    cat(sprintf("  PASS  %-42s %10.4f\n", label, actual))
  } else {
    FAIL <<- FAIL + 1
    cat(sprintf("  FAIL  %-42s %10.4f  (expected %.4f +/- %.3f)\n",
                label, actual, expected, tol))
  }
}

check_exact <- function(label, actual, expected) {
  ok <- identical(as.numeric(actual), as.numeric(expected))
  if (ok) {
    PASS <<- PASS + 1
    cat(sprintf("  PASS  %-42s %10s\n", label, format(actual)))
  } else {
    FAIL <<- FAIL + 1
    cat(sprintf("  FAIL  %-42s %10s  (expected %s)\n",
                label, format(actual), format(expected)))
  }
}

z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

cat("==================================================================\n")
cat("EII VALIDATION HARNESS\n")
cat("Baseline verified 2026-10-05 (R: pscl 1.5.9 / brms 2.23.0)\n")
cat("==================================================================\n")

# ------------------------------------------------------------------
# 1. PANEL INTEGRITY
# ------------------------------------------------------------------
cat("\n[1] PANEL INTEGRITY\n")
df <- read_csv(here("EIII_Panel_Definitive.csv"), show_col_types = FALSE)

check_exact("panel rows", nrow(df), 942)
check_exact("panel columns", ncol(df), 24)
# 2026-10-05: the panel contains 43 states, not 44. Every planning
# document in this project asserted 44, including the Sep 7 note that
# declared the 43/44 discrepancy resolved in favor of 44. It was
# resolved in the wrong direction. The panel is the authority.
check_exact("states in frame", n_distinct(df$country), 43)
check_exact("first year", min(df$year), 1945)
check_exact("last year", max(df$year), 2025)

# DV must be integer-valued and non-negative: the ZINB requires it
dv_ok <- all(df$eiii_score >= 0, na.rm = TRUE) &&
         all(df$eiii_score == floor(df$eiii_score), na.rm = TRUE)
if (dv_ok) {
  PASS <- PASS + 1; cat("  PASS  DV is non-negative integer\n")
} else {
  FAIL <- FAIL + 1; cat("  FAIL  DV is NOT non-negative integer\n")
}

# Layers must sum to the DV, or the weighting scheme has drifted
layer_sum <- df$layer1_ofac_events + df$layer2_multilateral +
             df$layer3_agency_actions + df$layer4_kinetic
cat(sprintf("  INFO  rows where DV == sum(layers): %.1f%%\n",
            100 * mean(df$eiii_score == layer_sum, na.rm = TRUE)))

# ------------------------------------------------------------------
# 2. ESTIMATION SAMPLE
# ------------------------------------------------------------------
cat("\n[2] ESTIMATION SAMPLE (V-Dem complete cases)\n")
d <- df |>
  mutate(
    eiii_score  = as.integer(replace_na(eiii_score, 0)),
    degradation = as.integer(replace_na(degradation, 0))
  ) |>
  drop_na(vdem_electoral, friction, resource_density,
          gdp_pc_log, fiscal_pressure) |>
  mutate(
    friction_z         = z_score(friction),
    resource_density_z = z_score(resource_density),
    fiscal_pressure_z  = z_score(fiscal_pressure),
    gdp_pc_log_z       = z_score(gdp_pc_log),
    trend_z            = z_score(trend),
    regime_z           = z_score(vdem_electoral),
    fric_x_fiscal      = friction_z * fiscal_pressure_z
  )

check_exact("estimation N", nrow(d), 535)
check_exact("estimation states", n_distinct(d$country), 32)
check_exact("estimation first year", min(d$year), 1970)
check_exact("estimation last year", max(d$year), 2021)
check("zero rate (%)", 100 * mean(d$eiii_score == 0), 49.9, tol = 0.1)
check("variance / mean", var(d$eiii_score) / mean(d$eiii_score), 1133.1, tol = 1)

# Z-scoring must actually have worked
check("mean(friction_z)", mean(d$friction_z), 0, tol = 1e-8)
check("sd(friction_z)", sd(d$friction_z), 1, tol = 1e-8)

# ------------------------------------------------------------------
# 3. FREQUENTIST COEFFICIENTS
# ------------------------------------------------------------------
# These are the values the Python pipeline got wrong. Each is checked
# BY NAME against the pscl coefficient table, so a reordering in any
# future version of any package fails here instead of silently
# relabelling the results.
cat("\n[3] FREQUENTIST ZINB — V-DEM (count equation)\n")
f <- zeroinfl(
  eiii_score ~ resource_density_z + friction_z + degradation +
               fiscal_pressure_z + trend_z + regime_z +
               gdp_pc_log_z + fric_x_fiscal |
               regime_z + friction_z + resource_density_z +
               degradation + trend_z,
  data = d, dist = "negbin"
)

if (!f$converged) {
  FAIL <- FAIL + 1; cat("  FAIL  model did not converge\n")
} else {
  PASS <- PASS + 1; cat("  PASS  model converged\n")
}

cc <- summary(f)$coefficients$count
check("count: (Intercept)",        cc["(Intercept)", "Estimate"],         3.7784)
check("count: resource_density_z", cc["resource_density_z", "Estimate"],  0.3362)
check("count: friction_z",         cc["friction_z", "Estimate"],         -0.4123)
check("count: degradation",        cc["degradation", "Estimate"],        -1.0826)
check("count: fiscal_pressure_z",  cc["fiscal_pressure_z", "Estimate"],  -1.0808)
check("count: trend_z",            cc["trend_z", "Estimate"],            -0.3400)
check("count: regime_z",           cc["regime_z", "Estimate"],            0.3528)
check("count: gdp_pc_log_z",       cc["gdp_pc_log_z", "Estimate"],       -0.6689)
check("count: fric_x_fiscal",      cc["fric_x_fiscal", "Estimate"],       0.0570)
check("log-likelihood",            as.numeric(f$loglik),              -1527.59, tol = 0.5)

cat("\n[4] FREQUENTIST ZINB — V-DEM (zero-inflation equation)\n")
cat("    Convention: POSITIVE = more structural non-enforcement.\n")
cz <- summary(f)$coefficients$zero
check("zi: (Intercept)",        cz["(Intercept)", "Estimate"],        -3.0351)
check("zi: regime_z",           cz["regime_z", "Estimate"],           -0.1339)
check("zi: friction_z",         cz["friction_z", "Estimate"],          0.6993)
check("zi: resource_density_z", cz["resource_density_z", "Estimate"],  0.3183)
check("zi: degradation",        cz["degradation", "Estimate"],        -0.0461)
check("zi: trend_z",            cz["trend_z", "Estimate"],            -2.6823)

# ------------------------------------------------------------------
# 4. THE PLR NULL
# ------------------------------------------------------------------
# The interaction was reported for three weeks as b = -0.412, p = .055
# ("suggestive"). That was the friction main effect under the wrong
# label. The true value is null and must stay null.
cat("\n[5] PATH OF LEAST RESISTANCE — MUST REMAIN NULL\n")
plr_b <- cc["fric_x_fiscal", "Estimate"]
plr_p <- cc["fric_x_fiscal", "Pr(>|z|)"]
check("PLR coefficient", plr_b, 0.0570)
check("PLR p-value",     plr_p, 0.7579, tol = 0.02)
if (plr_p > 0.10) {
  PASS <- PASS + 1
  cat("  PASS  PLR remains statistically indistinguishable from zero\n")
} else {
  FAIL <- FAIL + 1
  cat("  FAIL  PLR is now significant — investigate before reporting\n")
}

# ------------------------------------------------------------------
# 5. REGIME ROBUSTNESS
# ------------------------------------------------------------------
cat("\n[6] POLITY5 ROBUSTNESS — SIGNS MUST MATCH V-DEM\n")
dp <- df |>
  mutate(
    eiii_score  = as.integer(replace_na(eiii_score, 0)),
    degradation = as.integer(replace_na(degradation, 0))
  ) |>
  drop_na(polity, friction, resource_density, gdp_pc_log, fiscal_pressure) |>
  mutate(
    friction_z         = z_score(friction),
    resource_density_z = z_score(resource_density),
    fiscal_pressure_z  = z_score(fiscal_pressure),
    gdp_pc_log_z       = z_score(gdp_pc_log),
    trend_z            = z_score(trend),
    regime_z           = z_score(polity),
    fric_x_fiscal      = friction_z * fiscal_pressure_z
  )

fp <- zeroinfl(
  eiii_score ~ resource_density_z + friction_z + degradation +
               fiscal_pressure_z + trend_z + regime_z +
               gdp_pc_log_z + fric_x_fiscal |
               regime_z + friction_z + resource_density_z +
               degradation + trend_z,
  data = dp, dist = "negbin"
)
cp <- summary(fp)$coefficients$count
sign_match <- all(sign(cp[rownames(cc), "Estimate"]) == sign(cc[, "Estimate"]))
if (sign_match) {
  PASS <- PASS + 1
  cat("  PASS  all count signs match between V-Dem and Polity\n")
} else {
  FAIL <- FAIL + 1
  cat("  FAIL  sign mismatch between regime measures\n")
  print(data.frame(vdem = cc[, "Estimate"], polity = cp[rownames(cc), "Estimate"]))
}

# ------------------------------------------------------------------
# 6. THE ROBUST SUBSTANTIVE FINDINGS
# ------------------------------------------------------------------
# These four survived every specification run on 2026-10-05: two
# estimation frameworks, three regime measures, four prior
# specifications, four DV weightings. If any flips sign, something
# fundamental has changed.
cat("\n[7] SUBSTANTIVE FINDINGS — SIGN STABILITY\n")
expect_neg <- c("fiscal_pressure_z", "gdp_pc_log_z", "degradation", "friction_z")
for (v in expect_neg) {
  b <- cc[v, "Estimate"]
  if (b < 0) {
    PASS <- PASS + 1; cat(sprintf("  PASS  %-22s negative (%.3f)\n", v, b))
  } else {
    FAIL <- FAIL + 1; cat(sprintf("  FAIL  %-22s POSITIVE (%.3f)\n", v, b))
  }
}

# ------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------
cat("\n==================================================================\n")
cat(sprintf("RESULT: %d passed, %d failed\n", PASS, FAIL))
cat("==================================================================\n")
if (FAIL == 0) {
  cat("All checks passed. The pipeline reproduces the 2026-10-05 baseline.\n")
} else {
  cat("ONE OR MORE CHECKS FAILED.\n")
  cat("Do not send results anywhere until you know why.\n")
  cat("A failure is not necessarily an error — the panel or the code may\n")
  cat("have changed deliberately. But the change must be identified and\n")
  cat("this file's expected values updated with a dated note.\n")
  quit(status = 1)
}
