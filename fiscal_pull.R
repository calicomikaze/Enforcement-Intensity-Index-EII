# ==============================================================================
# FISCAL STRAIN INDEX — DATA PULL
# ==============================================================================
# Pulls five indicators of U.S. fiscal crowding-out from FRED, 1945-2025.
#
# CONSTRUCT: the degree to which the enforcing state's capacity to fund
# discretionary activity is constrained by prior commitments. Not "is the
# deficit large" — a state can run a large deficit with cheap debt and face
# no constraint. The construct is crowding-out.
#
# INDICATORS (all annual, all public, no API key required):
#   FYOIGDA188S  Net interest outlays, % GDP      — prior commitments
#   FYFSGDA188S  Surplus/deficit, % GDP            — flow pressure
#   GFDEGDQ188S  Public debt, % GDP (quarterly)    — stock pressure
#   FYONGDA188S  Net outlays, % GDP                — total claim on output
#   DGS10        10-year Treasury yield (daily)    — price of carrying stock
#
# SIPRI military spending deliberately EXCLUDED: it is a choice variable,
# not a constraint, and proxies enforcement spending — including it risks
# measuring the outcome with the predictor.
#
# OUTPUT: fiscal_indicators.csv (year, five columns)
# ==============================================================================

.libPaths(c("C:/R/library", .libPaths()))
library(tidyverse)
library(here)

# FRED serves plain CSV at this endpoint with no key.
fred_get <- function(series_id) {
  url <- paste0("https://fred.stlouisfed.org/graph/fredgraph.csv?id=", series_id)
  message("Pulling ", series_id, " ...")
  d <- tryCatch(
    read_csv(url, show_col_types = FALSE),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(d)) return(NULL)
  names(d) <- c("date", "value")
  d |>
    mutate(
      date  = as.Date(date),
      value = suppressWarnings(as.numeric(value)),
      year  = as.integer(format(date, "%Y"))
    ) |>
    filter(!is.na(value))
}

# ---- Annual series: take as-is ----
net_interest <- fred_get("FYOIGDA188S") |>
  group_by(year) |> summarize(net_interest = mean(value), .groups = "drop")

deficit <- fred_get("FYFSGDA188S") |>
  group_by(year) |> summarize(deficit = mean(value), .groups = "drop")

net_outlays <- fred_get("FYONGDA188S") |>
  group_by(year) |> summarize(net_outlays = mean(value), .groups = "drop")

# ---- Quarterly / daily series: collapse to annual means ----
debt_gdp <- fred_get("GFDEGDQ188S") |>
  group_by(year) |> summarize(debt_gdp = mean(value), .groups = "drop")

yield10 <- fred_get("DGS10") |>
  group_by(year) |> summarize(yield10 = mean(value, na.rm = TRUE), .groups = "drop")

# ---- Merge on year ----
fiscal <- tibble(year = 1945:2025) |>
  left_join(net_interest, by = "year") |>
  left_join(deficit,      by = "year") |>
  left_join(net_outlays,  by = "year") |>
  left_join(debt_gdp,     by = "year") |>
  left_join(yield10,      by = "year")

# ---- Orient every indicator so HIGHER = MORE STRAIN ----
# net_interest: higher = more strain (already correct)
# deficit:      FRED reports deficits as negative, so flip the sign
# net_outlays:  higher = more strain (already correct)
# debt_gdp:     higher = more strain (already correct)
# yield10:      higher = more strain (already correct)
fiscal <- fiscal |> mutate(deficit = -deficit)

# ---- Coverage report ----
cat("\n================= COVERAGE =================\n")
cat(sprintf("%-16s %6s %8s %10s %10s\n", "Indicator", "N", "Missing", "Min", "Max"))
cat(strrep("-", 54), "\n")
for (v in c("net_interest", "deficit", "net_outlays", "debt_gdp", "yield10")) {
  x <- fiscal[[v]]
  cat(sprintf("%-16s %6d %8d %10.2f %10.2f\n",
              v, sum(!is.na(x)), sum(is.na(x)),
              min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
}

cat("\nFirst year with all five:",
    min(fiscal$year[complete.cases(fiscal)]), "\n")
cat("Complete-case years:", sum(complete.cases(fiscal)), "of 81\n")

# ---- Correlation structure ----
# A single-factor model assumes the indicators share one common cause.
# If correlations are weak or sign-inconsistent, the construct is not
# unidimensional and the model is wrong.
cat("\n============ CORRELATIONS (complete cases) ============\n")
print(round(cor(fiscal[, -1], use = "complete.obs"), 3))

write_csv(fiscal, here("fiscal_indicators.csv"))
cat("\nSaved: fiscal_indicators.csv\n")
