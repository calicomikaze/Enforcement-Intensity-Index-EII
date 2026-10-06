# ==============================================================================
# SPEND AND COST AS SEPARATE MECHANISMS
# ==============================================================================
# THE RESPECIFICATION
# The pooled model contained one fiscal variable: net interest as % GDP.
# That is COST — what the architecture pays to carry accumulated
# obligations. It is not SPEND — what the architecture lays out to run
# itself.
#
# These are empirically distinct in this window: cor(spend, cost) = 0.042
# over 1970-2021. They are not two measures of one construct.
#
# Their variance profiles are inverted:
#                      Early (1970-2002)   Late (2003-2021)
#   Spend  (SD)              1.52                3.20
#   Cost   (SD)              0.70                0.16
#
# Cost stops moving exactly when the theory says expansion begins. The
# null fiscal result in the 2003-2021 layer test is therefore a dead
# variable, not an absent effect.
#
# PREDICTIONS
#   SPEND  positive on enforcement intensity. Unconstrained outlay
#          funds the apparatus; the architecture does more when it is
#          laying out more.
#   COST   negative. Carrying charges crowd out discretionary
#          enforcement.
#
# FALSIFICATION
#   If spend is null or negative, the expansion mechanism does not
#   operate through outlay and the "no constraint, maximum expansion"
#   argument loses its empirical leg.
#
# IDENTIFICATION CAVEAT — READ BEFORE INTERPRETING
#   Spend and cost are US-level annual series, constant within year
#   across all states. The effective N for these two coefficients is
#   the number of YEARS (52 pooled, 33 early, 19 late), not the number
#   of country-years. Conventional SEs are optimistic. The script
#   reports year counts alongside every coefficient so the reader can
#   judge. Treat these as associations across a few dozen time points,
#   not as precisely estimated panel effects.
# ==============================================================================

.libPaths(c("C:/R/library", .libPaths()))
library(tidyverse)
library(pscl)
library(here)

set.seed(42)
z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

panel  <- read_csv(here("EIII_Panel_Definitive.csv"), show_col_types = FALSE)
fiscal <- read_csv(here("fiscal_indicators.csv"), show_col_types = FALSE)

# ---- Merge spend onto the panel ----
# net_interest already in the panel as fiscal_pressure; bring in spend.
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
      cost_z             = z_score(fiscal_pressure),   # net interest, % GDP
      spend_z            = z_score(net_outlays),       # net outlays, % GDP
      gdp_pc_log_z       = z_score(gdp_pc_log),
      trend_z            = z_score(trend),
      regime_z           = z_score(vdem_electoral),
      fric_x_cost        = friction_z * cost_z
    )
}

d_all   <- build(panel)
d_early <- build(panel |> filter(year <= 2002))
d_late  <- build(panel |> filter(year >= 2003))

cat("================= SAMPLES =================\n")
cat(sprintf("%-10s %6s %7s %7s %12s %8s\n",
            "Period", "N", "States", "Years", "Span", "% Zero"))
cat(strrep("-", 56), "\n")
for (nm in c("Pooled", "Early", "Late")) {
  d <- switch(nm, Pooled = d_all, Early = d_early, Late = d_late)
  cat(sprintf("%-10s %6d %7d %7d %12s %8.1f\n", nm, nrow(d),
              n_distinct(d$country), n_distinct(d$year),
              paste0(min(d$year), "-", max(d$year)),
              100 * mean(d$eiii_score == 0)))
}

cat("\n========== SPEND AND COST: RAW VARIATION ==========\n")
cat(sprintf("%-10s %10s %8s %10s %8s %8s\n",
            "Period", "Spend SD", "range", "Cost SD", "range", "cor"))
cat(strrep("-", 60), "\n")
for (nm in c("Pooled", "Early", "Late")) {
  d <- switch(nm, Pooled = d_all, Early = d_early, Late = d_late)
  yr <- d |> distinct(year, net_outlays, fiscal_pressure)
  cat(sprintf("%-10s %10.3f %8.2f %10.3f %8.2f %8.3f\n", nm,
              sd(yr$net_outlays), diff(range(yr$net_outlays)),
              sd(yr$fiscal_pressure), diff(range(yr$fiscal_pressure)),
              cor(yr$net_outlays, yr$fiscal_pressure)))
}

# ---- Model: spend and cost entered separately ----
fit <- function(d) {
  tryCatch(
    zeroinfl(
      eiii_score ~ resource_density_z + friction_z + degradation +
                   cost_z + spend_z + trend_z + regime_z +
                   gdp_pc_log_z + fric_x_cost |
                   regime_z + friction_z + resource_density_z +
                   degradation + trend_z,
      data = d, dist = "negbin"
    ),
    error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL }
  )
}

cat("\n================= FITTING =================\n")
f_all   <- fit(d_all);   cat("Pooled:", if (!is.null(f_all))   f_all$converged   else FALSE, "\n")
f_early <- fit(d_early); cat("Early: ", if (!is.null(f_early)) f_early$converged else FALSE, "\n")
f_late  <- fit(d_late);  cat("Late:  ", if (!is.null(f_late))  f_late$converged  else FALSE, "\n")

pull_coef <- function(f, v) {
  if (is.null(f)) return(rep(NA_real_, 3))
  s <- summary(f)$coefficients$count
  if (!v %in% rownames(s)) return(rep(NA_real_, 3))
  c(s[v, "Estimate"], s[v, "Std. Error"], s[v, "Pr(>|z|)"])
}

# ---- THE TEST ----
cat("\n==================================================================\n")
cat("SPEND vs COST BY PERIOD\n")
cat("==================================================================\n")
cat("Prediction: spend POSITIVE, cost NEGATIVE.\n")
cat("N_years is the effective sample for these two coefficients.\n\n")
cat(sprintf("%-10s %-8s %10s %10s %10s %9s\n",
            "Period", "Var", "Coef", "SE", "p", "N_years"))
cat(strrep("-", 62), "\n")
for (nm in c("Pooled", "Early", "Late")) {
  f <- switch(nm, Pooled = f_all, Early = f_early, Late = f_late)
  d <- switch(nm, Pooled = d_all, Early = d_early, Late = d_late)
  ny <- n_distinct(d$year)
  for (v in c("spend_z", "cost_z")) {
    cf <- pull_coef(f, v)
    cat(sprintf("%-10s %-8s %10.4f %10.4f %10.5f %9d\n",
                nm, sub("_z$", "", v), cf[1], cf[2], cf[3], ny))
  }
}

# ---- Full table ----
vars <- c("resource_density_z", "friction_z", "degradation",
          "cost_z", "spend_z", "trend_z", "regime_z",
          "gdp_pc_log_z", "fric_x_cost")
cat("\n==================================================================\n")
cat("ALL COUNT COEFFICIENTS BY PERIOD\n")
cat("==================================================================\n")
cat(sprintf("\n%-20s %9s %8s %9s %8s %9s %8s\n",
            "Variable", "Pool b", "p", "Early b", "p", "Late b", "p"))
cat(strrep("-", 78), "\n")
for (v in vars) {
  a <- pull_coef(f_all, v); e <- pull_coef(f_early, v); l <- pull_coef(f_late, v)
  cat(sprintf("%-20s %9.3f %8.4f %9.3f %8.4f %9.3f %8.4f\n",
              v, a[1], a[3], e[1], e[3], l[1], l[3]))
}

# ---- Verdict ----
cat("\n==================================================================\n")
cat("VERDICT\n")
cat("==================================================================\n")
sp <- pull_coef(f_all, "spend_z"); co <- pull_coef(f_all, "cost_z")
cat("Pooled spend:", round(sp[1], 4), " p =", round(sp[3], 4), "\n")
cat("Pooled cost: ", round(co[1], 4), " p =", round(co[3], 4), "\n\n")

if (!is.na(sp[3]) && !is.na(co[3])) {
  if (sp[1] > 0 && sp[3] < 0.05 && co[1] < 0 && co[3] < 0.05) {
    cat("BOTH MECHANISMS SUPPORTED. Spend expands enforcement,\n")
    cat("cost contracts it. Pooling them into one fiscal variable\n")
    cat("was masking two opposing effects.\n")
  } else if (sp[1] > 0 && sp[3] < 0.05) {
    cat("SPEND SUPPORTED, cost not. The expansion mechanism operates\n")
    cat("through outlay; the carrying-cost mechanism is not detectable\n")
    cat("in this window.\n")
  } else if (co[1] < 0 && co[3] < 0.05) {
    cat("COST SUPPORTED, spend not. The original specification was\n")
    cat("adequate; separating spend adds nothing.\n")
  } else {
    cat("NEITHER. Report plainly. The expansion mechanism does not\n")
    cat("operate through federal outlay as measured here.\n")
  }
}

out <- bind_rows(lapply(c("Pooled", "Early", "Late"), function(nm) {
  f <- switch(nm, Pooled = f_all, Early = f_early, Late = f_late)
  d <- switch(nm, Pooled = d_all, Early = d_early, Late = d_late)
  s <- pull_coef(f, "spend_z"); c_ <- pull_coef(f, "cost_z")
  tibble(period = nm, n = nrow(d), n_years = n_distinct(d$year),
         spend_b = s[1], spend_se = s[2], spend_p = s[3],
         cost_b = c_[1], cost_se = c_[2], cost_p = c_[3])
}))
write_csv(out, here("spend_cost_split.csv"))
cat("\nSaved: spend_cost_split.csv\n")
