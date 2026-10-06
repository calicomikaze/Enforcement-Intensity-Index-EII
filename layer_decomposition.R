# ==============================================================================
# LAYER DECOMPOSITION — THE SUBSTITUTION TEST
# ==============================================================================
# THEORY
# The enforcement architecture produces two kinds of output.
#
#   SPECTACLE (Layers 1, 2, 4): executive designations, multilateral
#   listings, kinetic operations. Visible, announced, addressed to an
#   audience. Expensive: requires diplomatic capital, coalition
#   maintenance, military capacity. The state pays.
#
#   APPARATUS (Layer 3): civil penalties against financial institutions.
#   Negotiated, quietly announced, and the operative product is the
#   compliance infrastructure the settlement installs. Cheap to the
#   state: the penalized institutions absorb the cost of the apparatus
#   they are required to build.
#
# PREDICTION (asymmetric substitution)
#   Under fiscal strain the architecture cuts what it pays for and
#   preserves what others pay for. Fiscal pressure should be strongly
#   negative on Layers 1, 2, 4 and attenuated, null, or positive on
#   Layer 3.
#
# FALSIFICATION
#   If fiscal pressure is uniformly negative across all four layers,
#   there is no substitution. The architecture simply does less of
#   everything under strain — contraction without theater. The
#   spectacle/apparatus distinction would then be a description of
#   content, not a mechanism, and should be dropped from the theory.
#
# WINDOW
#   OFAC civil penalties (Layer 3) are published from 2003. Before
#   that, a Layer 3 zero is non-REPORTING, not non-enforcement — a
#   different generating process from the structural zeros the model
#   is built to estimate. All four layers are therefore restricted to
#   the common window so the comparison is like-for-like.
# ==============================================================================

.libPaths(c("C:/R/library", .libPaths()))
library(tidyverse)
library(pscl)
library(here)

set.seed(42)

z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

df <- read_csv(here("EIII_Panel_Definitive.csv"), show_col_types = FALSE)

layers <- c("layer1_ofac_events", "layer2_multilateral",
            "layer3_agency_actions", "layer4_kinetic")

df <- df |>
  mutate(
    across(all_of(layers), ~as.integer(replace_na(.x, 0))),
    degradation = as.integer(replace_na(degradation, 0)),
    across(c(polity, friction, resource_density, fiscal_pressure,
             gdp_pc_log, trend, vdem_electoral), as.numeric)
  )

# ---- Establish the common window empirically ----
# Do not assume 2003 — find the first year Layer 3 reports anything.
l3_first <- df |> filter(layer3_agency_actions > 0) |> pull(year) |> min()
cat("First year with any Layer 3 activity:", l3_first, "\n")

first_yr <- df |>
  group_by(year) |>
  summarize(l3 = sum(layer3_agency_actions), .groups = "drop") |>
  filter(l3 > 0) |>
  pull(year) |> min()

WINDOW_START <- first_yr
cat("Common window starts:", WINDOW_START, "\n\n")

# ---- Build the estimation sample on the common window ----
d <- df |>
  filter(year >= WINDOW_START) |>
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

cat("Window sample:", nrow(d), "country-years,",
    n_distinct(d$country), "states,",
    min(d$year), "-", max(d$year), "\n\n")

# ---- Layer descriptives on the common window ----
cat("=============== LAYER DESCRIPTIVES ===============\n")
cat(sprintf("%-24s %8s %8s %8s %10s\n",
            "Layer", "Mean", "Max", "% Zero", "Var/Mean"))
cat(strrep("-", 62), "\n")
for (L in layers) {
  x <- d[[L]]
  cat(sprintf("%-24s %8.2f %8d %8.1f %10.1f\n",
              L, mean(x), max(x), 100 * mean(x == 0),
              var(x) / max(mean(x), 1e-10)))
}

# ---- Fit the same specification to each layer ----
fit_layer <- function(dv) {
  d$y <- d[[dv]]
  tryCatch(
    zeroinfl(
      y ~ resource_density_z + friction_z + degradation +
          fiscal_pressure_z + trend_z + regime_z +
          gdp_pc_log_z + fric_x_fiscal |
          regime_z + friction_z + resource_density_z +
          degradation + trend_z,
      data = d, dist = "negbin"
    ),
    error = function(e) {
      message("  ", dv, " FAILED: ", conditionMessage(e))
      NULL
    }
  )
}

cat("\n=============== FITTING ===============\n")
fits <- list()
for (L in layers) {
  cat("Fitting", L, "... ")
  f <- fit_layer(L)
  if (!is.null(f)) {
    cat("converged:", f$converged, "\n")
    fits[[L]] <- f
  } else {
    cat("\n")
  }
}

# ---- THE TEST ----
cat("\n==================================================================\n")
cat("THE SUBSTITUTION TEST — FISCAL PRESSURE BY LAYER\n")
cat("==================================================================\n")
cat("\nPrediction: strongly negative on 1, 2, 4 (spectacle);\n")
cat("            attenuated/null/positive on 3 (apparatus).\n\n")
cat(sprintf("%-24s %10s %10s %10s %8s\n",
            "Layer", "Fiscal b", "SE", "p", "Type"))
cat(strrep("-", 66), "\n")

type_of <- c(layer1_ofac_events = "spectacle",
             layer2_multilateral = "spectacle",
             layer3_agency_actions = "APPARATUS",
             layer4_kinetic = "spectacle")

rows <- list()
for (L in layers) {
  if (is.null(fits[[L]])) next
  s <- summary(fits[[L]])$coefficients$count
  if (!"fiscal_pressure_z" %in% rownames(s)) next
  b  <- s["fiscal_pressure_z", "Estimate"]
  se <- s["fiscal_pressure_z", "Std. Error"]
  p  <- s["fiscal_pressure_z", "Pr(>|z|)"]
  cat(sprintf("%-24s %10.4f %10.4f %10.5f %8s\n", L, b, se, p, type_of[L]))
  rows[[L]] <- tibble(layer = L, type = type_of[L],
                      fiscal_b = b, fiscal_se = se, fiscal_p = p)
}

# ---- Full coefficient comparison across layers ----
cat("\n==================================================================\n")
cat("ALL COUNT COEFFICIENTS BY LAYER\n")
cat("==================================================================\n")
vars <- c("resource_density_z", "friction_z", "degradation",
          "fiscal_pressure_z", "trend_z", "regime_z",
          "gdp_pc_log_z", "fric_x_fiscal")
cat(sprintf("\n%-20s %11s %11s %11s %11s\n",
            "Variable", "L1 OFAC", "L2 Multi", "L3 Penalty", "L4 Kinetic"))
cat(strrep("-", 68), "\n")
for (v in vars) {
  vals <- sapply(layers, function(L) {
    if (is.null(fits[[L]])) return(NA_real_)
    s <- summary(fits[[L]])$coefficients$count
    if (!v %in% rownames(s)) return(NA_real_)
    s[v, "Estimate"]
  })
  cat(sprintf("%-20s %11.3f %11.3f %11.3f %11.3f\n", v,
              vals[1], vals[2], vals[3], vals[4]))
}

# ---- Verdict ----
cat("\n==================================================================\n")
cat("VERDICT\n")
cat("==================================================================\n")
res <- bind_rows(rows)
if (nrow(res) == 4) {
  spec <- res |> filter(type == "spectacle")
  app  <- res |> filter(type == "APPARATUS")
  cat("Spectacle layers, mean fiscal coefficient:",
      round(mean(spec$fiscal_b), 4), "\n")
  cat("Apparatus layer, fiscal coefficient:      ",
      round(app$fiscal_b, 4), "\n\n")
  if (all(spec$fiscal_b < 0) && app$fiscal_b >= min(spec$fiscal_b)) {
    cat("CONSISTENT with asymmetric substitution.\n")
    cat("The architecture cuts what it pays for and preserves\n")
    cat("what the penalized institutions pay for.\n")
  } else if (all(res$fiscal_b < 0) &&
             abs(app$fiscal_b - mean(spec$fiscal_b)) < 0.2) {
    cat("FALSIFIED. Fiscal strain suppresses all four layers alike.\n")
    cat("No substitution — contraction without theater. The\n")
    cat("spectacle/apparatus distinction describes content, not\n")
    cat("mechanism, and should be dropped from the theory.\n")
  } else {
    cat("MIXED. Report the pattern plainly; do not force a reading.\n")
  }
  write_csv(res, here("layer_substitution_test.csv"))
  cat("\nSaved: layer_substitution_test.csv\n")
} else {
  cat("Not all four layers estimated — verdict withheld.\n")
}
