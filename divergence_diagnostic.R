# ==============================================================================
# FISCAL REPORTING DIVERGENCE — DIAGNOSTIC
# ==============================================================================
# CONSTRUCT: the gap between the REPORTED carrying cost of federal debt
# (net interest as % GDP — the number CBO and OMB publish, and the number
# in the enforcement model) and the cost IMPLIED by the obligation stock
# and its market price.
#
# When the stock rises and the reported cost does not, the published
# measure has stopped tracking the underlying position. That gap is
# institutional data production made visible — the same logic as
# Martinez's nightlights, with the bond market as the channel the
# reporting agency cannot touch.
#
# THIS SCRIPT DOES NOT BUILD THE MEASURE. It tests whether the measure
# is worth building, by answering one question: is the divergence
# distinguishable from a time trend? If the residual is just "later
# years," it is collinear with trend_z and contributes nothing.
#
# GATE: if |cor(divergence, year)| > 0.7, the construct is a trend in
# disguise and this line of work stops here.
# ==============================================================================

.libPaths(c("C:/R/library", .libPaths()))
library(tidyverse)
library(here)

fiscal <- read_csv(here("fiscal_indicators.csv"), show_col_types = FALSE)

# Restrict to the window the enforcement model actually uses
f <- fiscal |> filter(year >= 1966, year <= 2025, complete.cases(fiscal))

cat("Window:", min(f$year), "-", max(f$year), "|", nrow(f), "years\n\n")

# ---- The mechanical relationship ----
# Carrying cost should be a function of how much you owe and what it costs.
# net_interest ~ debt_gdp + yield10 is the accounting identity in
# reduced form. The residual is what the reported number does that the
# stock and the price do not explain.
m <- lm(net_interest ~ debt_gdp + yield10, data = f)

cat("=========== MECHANICAL MODEL ===========\n")
print(summary(m))

f$divergence <- resid(m)

# ---- GATE 1: is it a time trend? ----
r_trend <- cor(f$divergence, f$year)
cat("\n=========== GATE 1: TREND COLLINEARITY ===========\n")
cat("cor(divergence, year) =", round(r_trend, 3), "\n")
if (abs(r_trend) > 0.7) {
  cat("*** FAIL: divergence is a time trend in disguise. ***\n")
} else {
  cat("PASS: divergence carries information beyond the trend.\n")
}

# ---- GATE 2: does it have usable variance? ----
cat("\n=========== GATE 2: VARIANCE ===========\n")
cat("SD of divergence:", round(sd(f$divergence), 4), "\n")
cat("SD of net_interest:", round(sd(f$net_interest), 4), "\n")
cat("Ratio:", round(sd(f$divergence) / sd(f$net_interest), 3), "\n")
cat("(A ratio below ~0.15 means the stock and price explain nearly\n")
cat(" everything and there is no reporting gap to measure.)\n")

# ---- GATE 3: is it autocorrelated noise? ----
# A measure that is just last year's value plus noise is not a construct.
ar1 <- cor(f$divergence[-1], f$divergence[-nrow(f)])
cat("\n=========== GATE 3: AR(1) ===========\n")
cat("First-order autocorrelation:", round(ar1, 3), "\n")

# ---- Where does the reported number understate? ----
cat("\n=========== LARGEST NEGATIVE DIVERGENCES ===========\n")
cat("(reported carrying cost BELOW what stock and rates imply)\n\n")
f |>
  arrange(divergence) |>
  select(year, net_interest, debt_gdp, yield10, divergence) |>
  head(12) |>
  mutate(across(where(is.numeric), ~round(.x, 2))) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\n=========== LARGEST POSITIVE DIVERGENCES ===========\n")
cat("(reported carrying cost ABOVE what stock and rates imply)\n\n")
f |>
  arrange(desc(divergence)) |>
  select(year, net_interest, debt_gdp, yield10, divergence) |>
  head(12) |>
  mutate(across(where(is.numeric), ~round(.x, 2))) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ---- Plot ----
p <- ggplot(f, aes(x = year, y = divergence)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(
    title = "Fiscal Reporting Divergence, 1966-2025",
    subtitle = "Residual from net interest ~ debt stock + 10-year yield.\nNegative = reported carrying cost below what the position implies.",
    x = NULL, y = "Divergence (percentage points of GDP)"
  ) +
  theme_minimal(base_family = "serif")
ggsave(here("fig_divergence.png"), p, width = 9, height = 5, dpi = 300)

write_csv(f |> select(year, net_interest, debt_gdp, yield10, divergence),
          here("fiscal_divergence.csv"))

cat("\nSaved: fiscal_divergence.csv, fig_divergence.png\n")
