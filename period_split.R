# ==============================================================================
# PERIOD SPLIT — IS THE FISCAL EFFECT CONDITIONAL ON THE CONSTRAINT BINDING?
# ==============================================================================
# MOTIVATION
# The pooled model finds fiscal pressure strongly negative on enforcement
# intensity (P(dir.) = 1.000). The 2003-2021 subsample finds nothing
# (p = 0.35 to 0.96 across all four layers). If the pooled effect is
# carried entirely by the earlier period, the finding is conditional,
# not general, and the paper must say so.
#
# THEORY
# Fiscal pressure can only suppress enforcement when the fiscal
# constraint actually binds. Carrying cost was high and rising
# 1979-1995 and again from 2022. Between 1996 and 2021 the stock grew
# while the reported cost stayed flat or fell: debt expanded without
# a binding budget constraint. An architecture facing no constraint
# has no reason to economize on enforcement.
#
# PREDICTION
#   Early period (1970-2002):  fiscal strongly negative.
#   Late period  (2003-2021):  fiscal null.
#
# FALSIFICATION
#   If fiscal is null in BOTH subperiods, the pooled result is an
#   artifact of pooling and the fiscal finding does not survive.
#   If fiscal is negative in BOTH, the effect is general and the
#   layer-test null was a power problem, not a period effect.
#
# POWER CAVEAT
# Splitting halves the sample. A null in a subperiod is weak evidence
# of absence. The script reports N and SEs for both periods so the
# reader can judge whether a null is informative or merely imprecise.
# ==============================================================================

.libPaths(c("C:/R/library", .libPaths()))
library(tidyverse)
library(pscl)
library(here)

set.seed(42)
z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

df <- read_csv(here("EIII_Panel_Definitive.csv"), show_col_types = FALSE)
df <- df |>
  mutate(
    eiii_score  = as.integer(replace_na(eiii_score, 0)),
    degradation = as.integer(replace_na(degradation, 0)),
    across(c(friction, resource_density, fiscal_pressure,
             gdp_pc_log, trend, vdem_electoral), as.numeric)
  )

# Z-score WITHIN each subsample: coefficients then describe variation
# available in that period, which is the quantity of interest.
build <- function(data) {
  data |>
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
}

d_all   <- build(df)
d_early <- build(df |> filter(year <= 2002))
d_late  <- build(df |> filter(year >= 2003))

fit <- function(d) {
  tryCatch(
    zeroinfl(
      eiii_score ~ resource_density_z + friction_z + degradation +
                   fiscal_pressure_z + trend_z + regime_z +
                   gdp_pc_log_z + fric_x_fiscal |
                   regime_z + friction_z + resource_density_z +
                   degradation + trend_z,
      data = d, dist = "negbin"
    ),
    error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL }
  )
}

cat("================= SAMPLES =================\n")
cat(sprintf("%-16s %6s %8s %10s %10s\n",
            "Period", "N", "States", "Years", "% Zero"))
cat(strrep("-", 54), "\n")
for (nm in c("Pooled", "Early", "Late")) {
  d <- switch(nm, Pooled = d_all, Early = d_early, Late = d_late)
  cat(sprintf("%-16s %6d %8d %10s %10.1f\n", nm, nrow(d),
              n_distinct(d$country),
              paste0(min(d$year), "-", max(d$year)),
              100 * mean(d$eiii_score == 0)))
}

# ---- Variance in the IV: the whole question ----
cat("\n========== FISCAL PRESSURE: RAW VARIATION ==========\n")
cat(sprintf("%-16s %10s %10s %10s %10s\n",
            "Period", "Mean", "SD", "Min", "Max"))
cat(strrep("-", 58), "\n")
for (nm in c("Pooled", "Early", "Late")) {
  d <- switch(nm, Pooled = d_all, Early = d_early, Late = d_late)
  x <- d$fiscal_pressure
  cat(sprintf("%-16s %10.3f %10.3f %10.3f %10.3f\n",
              nm, mean(x), sd(x), min(x), max(x)))
}
cat("\n(If the late period has little variation in fiscal pressure,\n")
cat(" a null coefficient there is uninformative, not evidence of\n")
cat(" absence. Read the SDs before reading the coefficients.)\n")

# ---- Fit all three ----
cat("\n================= FITTING =================\n")
f_all   <- fit(d_all);   cat("Pooled converged:", !is.null(f_all)  && f_all$converged, "\n")
f_early <- fit(d_early); cat("Early converged: ", !is.null(f_early) && f_early$converged, "\n")
f_late  <- fit(d_late);  cat("Late converged:  ", !is.null(f_late)  && f_late$converged, "\n")

# ---- The comparison ----
pull_coef <- function(f, v) {
  if (is.null(f)) return(rep(NA_real_, 3))
  s <- summary(f)$coefficients$count
  if (!v %in% rownames(s)) return(rep(NA_real_, 3))
  c(s[v, "Estimate"], s[v, "Std. Error"], s[v, "Pr(>|z|)"])
}

cat("\n==================================================================\n")
cat("FISCAL PRESSURE BY PERIOD\n")
cat("==================================================================\n")
cat(sprintf("\n%-16s %10s %10s %10s %8s\n",
            "Period", "Coef", "SE", "p", "N"))
cat(strrep("-", 58), "\n")
for (nm in c("Pooled", "Early", "Late")) {
  f <- switch(nm, Pooled = f_all, Early = f_early, Late = f_late)
  d <- switch(nm, Pooled = d_all, Early = d_early, Late = d_late)
  cf <- pull_coef(f, "fiscal_pressure_z")
  cat(sprintf("%-16s %10.4f %10.4f %10.5f %8d\n",
              nm, cf[1], cf[2], cf[3], nrow(d)))
}

# ---- Every coefficient, both periods ----
vars <- c("resource_density_z", "friction_z", "degradation",
          "fiscal_pressure_z", "trend_z", "regime_z",
          "gdp_pc_log_z", "fric_x_fiscal")
cat("\n==================================================================\n")
cat("ALL COUNT COEFFICIENTS BY PERIOD\n")
cat("==================================================================\n")
cat(sprintf("\n%-20s %9s %9s %9s %9s %9s %9s\n",
            "Variable", "Pool b", "Pool p", "Early b", "Early p",
            "Late b", "Late p"))
cat(strrep("-", 80), "\n")
for (v in vars) {
  a <- pull_coef(f_all, v); e <- pull_coef(f_early, v); l <- pull_coef(f_late, v)
  cat(sprintf("%-20s %9.3f %9.4f %9.3f %9.4f %9.3f %9.4f\n",
              v, a[1], a[3], e[1], e[3], l[1], l[3]))
}

# ---- Verdict ----
cat("\n==================================================================\n")
cat("VERDICT\n")
cat("==================================================================\n")
e <- pull_coef(f_early, "fiscal_pressure_z")
l <- pull_coef(f_late,  "fiscal_pressure_z")
sd_e <- sd(d_early$fiscal_pressure); sd_l <- sd(d_late$fiscal_pressure)

cat("Early fiscal:", round(e[1], 4), "(p =", round(e[3], 4),
    "), raw SD =", round(sd_e, 3), "\n")
cat("Late fiscal: ", round(l[1], 4), "(p =", round(l[3], 4),
    "), raw SD =", round(sd_l, 3), "\n\n")

if (!is.na(e[3]) && !is.na(l[3])) {
  if (e[3] < 0.05 && e[1] < 0 && l[3] > 0.10) {
    cat("CONDITIONAL. Fiscal strain suppresses enforcement when the\n")
    cat("constraint binds and does nothing when it does not. The\n")
    cat("pooled coefficient is an average across two regimes and\n")
    cat("should not be reported without the split.\n")
  } else if (e[3] > 0.10 && l[3] > 0.10) {
    cat("FALSIFIED. Fiscal is null in both subperiods. The pooled\n")
    cat("result is an artifact of pooling and does not survive.\n")
  } else if (e[3] < 0.05 && l[3] < 0.05 && e[1] < 0 && l[1] < 0) {
    cat("GENERAL. Fiscal negative in both periods; the layer-test\n")
    cat("null was a power problem, not a period effect.\n")
  } else {
    cat("AMBIGUOUS. Report both periods; do not force a reading.\n")
  }
  if (sd_l < 0.5 * sd_e) {
    cat("\nNOTE: late-period variation in fiscal pressure is less than\n")
    cat("half the early period's. A late null is weak evidence of\n")
    cat("absence — there is little variation to detect an effect from.\n")
  }
}

out <- tibble(
  period = c("Pooled", "Early", "Late"),
  years  = c(paste0(min(d_all$year), "-", max(d_all$year)),
             paste0(min(d_early$year), "-", max(d_early$year)),
             paste0(min(d_late$year), "-", max(d_late$year))),
  n      = c(nrow(d_all), nrow(d_early), nrow(d_late)),
  fiscal_b  = c(pull_coef(f_all, "fiscal_pressure_z")[1], e[1], l[1]),
  fiscal_se = c(pull_coef(f_all, "fiscal_pressure_z")[2], e[2], l[2]),
  fiscal_p  = c(pull_coef(f_all, "fiscal_pressure_z")[3], e[3], l[3]),
  fiscal_sd_raw = c(sd(d_all$fiscal_pressure), sd_e, sd_l)
)
write_csv(out, here("period_split_test.csv"))
cat("\nSaved: period_split_test.csv\n")
