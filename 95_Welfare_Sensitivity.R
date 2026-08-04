################################################################################
# 95_Welfare_Sensitivity.R
#
# GAP 5 of 5  ->  MANUSCRIPT SECTION 6.7 "Welfare Sensitivity Analysis"
#
# Target values to reproduce:
#   (a) VSL RM 1.0m -> RM 128m ; RM 1.5m -> RM 193m ; RM 2.0m -> RM 257m
#   (b) Age-adjusted VSL (30-40% reduction) -> RM 116-135m
#   (c) Haze-season-only restriction         -> 89 excess deaths
#
# Coverage status before this script:
#   (a) DONE in revised_counterfactual.R (VSL_low / VSL_main / VSL_high)
#   (b) NOT CODED ANYWHERE
#   (c) NOT CODED ANYWHERE
#
# This script recomputes (a) from scratch so all three scenarios sit in one
# table with one set of coefficients, and adds (b) and (c).
#
# ALSO FIXES A BUG carried by revised_counterfactual.R line ~285, where
# `total_deaths = 47356` is hard-coded INSIDE group_by(year) %>% summarise().
# That gives every year the same denominator and makes the "% attributable"
# column meaningless. Here the denominator is computed per year from the data.
#
# Output: Output/Gap5_welfare_sensitivity.csv
################################################################################

source("90_Common_Setup.R")

cat("################################################################\n")
cat("# GAP 5: WELFARE SENSITIVITY (SECTION 6.7)\n")
cat("################################################################\n\n")

merged_daily <- prep_analysis_data()

POLL <- "PM2.5AVG"
d    <- complete_sample(merged_daily, POLL)

#===============================================================================
# 1. ESTIMATE BETA - FULL YEAR AND HAZE SEASON
#===============================================================================

cat("---------------------------------------------------------------\n")
cat("STEP 1: Coefficients\n")
cat("---------------------------------------------------------------\n")

# Full-sample beta (denominator for scenarios a and b)
r_full <- run_cf_iv(d, POLL, cluster = "twoway", label = "Full sample")
beta_full <- r_full$coefficient
se_full   <- r_full$se

# Haze-season beta. Section 4.4 reports 12.7% per 10 ug/m3 in June-September
# via an interaction model; here it is estimated on the haze-month subsample,
# which is the quantity Section 6.7 actually needs ("based on more precisely
# estimated coefficients").
d$haze_season <- as.integer(d$month %in% 6:9)
d_haze <- d %>% filter(haze_season == 1)

r_haze <- run_cf_iv(d_haze, POLL, cluster = "twoway",
                    label = "Haze season (Jun-Sep) only")
beta_haze <- r_haze$coefficient

cat(sprintf("  Full-sample beta : %.5f  (%.1f%% per 10 ug/m3)\n",
            beta_full, (exp(beta_full * 10) - 1) * 100))
cat(sprintf("  Haze-season beta : %.5f  (%.1f%% per 10 ug/m3)\n",
            beta_haze, (exp(beta_haze * 10) - 1) * 100))
cat("  Section 4.4 reports 12.7%% for haze season via interaction.\n\n")

#===============================================================================
# 2. FIRE-ATTRIBUTABLE POLLUTION AND EXCESS DEATHS
#===============================================================================

cat("---------------------------------------------------------------\n")
cat("STEP 2: Excess deaths\n")
cat("---------------------------------------------------------------\n")

fs <- r_full$first_stage
a1 <- coef(fs)[["FRP_u_wind_station"]]
a2 <- coef(fs)[["FRP_v_wind_station"]]

d$fire_pm25 <- pmax(a1 * d$FRP_u_wind_station + a2 * d$FRP_v_wind_station, 0)

#' Excess deaths = Deaths_obs * [1 - exp(-beta * Delta_fire)]   [manuscript S6.1]
#' Delta method CI on beta propagated through the same expression.
excess_deaths <- function(data, beta, beta_se, subset_idx = NULL) {
  if (!is.null(subset_idx)) data <- data[subset_idx, , drop = FALSE]
  
  ex     <- data$deaths * (1 - exp(-beta * data$fire_pm25))
  total  <- sum(ex, na.rm = TRUE)
  
  # d/dbeta of sum[D * (1 - exp(-b*F))] = sum[D * F * exp(-b*F)]
  grad   <- sum(data$deaths * data$fire_pm25 * exp(-beta * data$fire_pm25),
                na.rm = TRUE)
  se_tot <- abs(grad) * beta_se
  
  list(total = total,
       se    = se_tot,
       lo    = total - 1.96 * se_tot,
       hi    = total + 1.96 * se_tot,
       obs_deaths = sum(data$deaths, na.rm = TRUE),
       n = nrow(data))
}

ex_full <- excess_deaths(d, beta_full, se_full)
ex_haze <- excess_deaths(d, beta_haze, r_haze$se, which(d$haze_season == 1))

cat(sprintf("  Full year   : %.1f excess deaths  [95%% CI %.0f - %.0f]\n",
            ex_full$total, ex_full$lo, ex_full$hi))
cat(sprintf("  Haze season : %.1f excess deaths  [95%% CI %.0f - %.0f]\n",
            ex_haze$total, ex_haze$lo, ex_haze$hi))
cat("  Section 6.7 reports 89 for the haze-season restriction.\n\n")

#--- Per-year breakdown, with the CORRECT per-year denominator ----------------
by_year <- d %>%
  mutate(excess_i = deaths * (1 - exp(-beta_full * fire_pm25))) %>%
  group_by(year) %>%
  summarise(
    N_obs          = n(),
    Observed_deaths = sum(deaths, na.rm = TRUE),   # <- per year, NOT hard-coded
    Excess_deaths  = sum(excess_i, na.rm = TRUE),
    Mean_fire_pm25 = mean(fire_pm25, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Pct_attributable = 100 * Excess_deaths / Observed_deaths)

cat("  ---- Table 13 rebuilt with per-year denominators ----\n")
print(by_year %>% mutate(across(where(is.numeric), ~ round(.x, 2))))
cat(sprintf("\n  TOTAL observed deaths in this sample: %s\n",
            format(sum(by_year$Observed_deaths), big.mark = ",")))
cat(sprintf("  TOTAL excess deaths                : %.1f\n", sum(by_year$Excess_deaths)))
cat(sprintf("  Overall %% attributable             : %.2f%%\n\n",
            100 * sum(by_year$Excess_deaths) / sum(by_year$Observed_deaths)))

#===============================================================================
# 3. WELFARE SCENARIOS
#===============================================================================

cat("---------------------------------------------------------------\n")
cat("STEP 3: Welfare scenarios\n")
cat("---------------------------------------------------------------\n\n")

VSL_LOW  <- 1.0e6   # RM 1.0 million
VSL_MAIN <- 1.5e6   # RM 1.5 million (central, Viscusi & Masterman 2017)
VSL_HIGH <- 2.0e6   # RM 2.0 million

# Section 6.7: age adjustment "would reduce estimates by approximately 30-40%"
AGE_ADJ_LOW  <- 0.60   # 40% reduction
AGE_ADJ_HIGH <- 0.70   # 30% reduction

mk <- function(scenario, deaths, vsl, mult = 1, ci_lo = NA, ci_hi = NA, note = "") {
  data.frame(
    Scenario       = scenario,
    Excess_deaths  = deaths,
    VSL_RM         = vsl,
    Adjustment     = mult,
    Welfare_RM_mn  = deaths * vsl * mult / 1e6,
    CI_lo_RM_mn    = ifelse(is.na(ci_lo), NA, ci_lo * vsl * mult / 1e6),
    CI_hi_RM_mn    = ifelse(is.na(ci_hi), NA, ci_hi * vsl * mult / 1e6),
    Note           = note,
    stringsAsFactors = FALSE
  )
}

scen <- rbind(
  # (a) VSL sensitivity
  mk("Central: VSL RM 1.5m", ex_full$total, VSL_MAIN, 1,
     ex_full$lo, ex_full$hi, "Baseline reported in Table 14"),
  mk("Lower bound: VSL RM 1.0m", ex_full$total, VSL_LOW, 1,
     ex_full$lo, ex_full$hi, "Section 6.7 target: RM 128m"),
  mk("Upper bound: VSL RM 2.0m", ex_full$total, VSL_HIGH, 1,
     ex_full$lo, ex_full$hi, "Section 6.7 target: RM 257m"),
  
  # (b) Age-adjusted VSL - NOT PREVIOUSLY CODED
  mk("Age-adjusted, -40%", ex_full$total, VSL_MAIN, AGE_ADJ_LOW,
     NA, NA, "Section 6.7 target: RM 116m"),
  mk("Age-adjusted, -30%", ex_full$total, VSL_MAIN, AGE_ADJ_HIGH,
     NA, NA, "Section 6.7 target: RM 135m"),
  
  # (c) Haze-season restriction - NOT PREVIOUSLY CODED
  mk("Haze season only (Jun-Sep)", ex_haze$total, VSL_MAIN, 1,
     ex_haze$lo, ex_haze$hi, "Section 6.7 target: 89 deaths")
)

print(scen %>% mutate(across(c(Excess_deaths, Welfare_RM_mn,
                               CI_lo_RM_mn, CI_hi_RM_mn), ~ round(.x, 1))))

# USD conversion for the abstract
USD_RATE <- 4.3   # RM per USD, approximate 2020 average - VERIFY BEFORE USE
cat(sprintf("\n  Central estimate in USD: $%.1f million (at RM %.1f/USD)\n",
            scen$Welfare_RM_mn[1] / USD_RATE, USD_RATE))
cat("  Abstract claims USD 45 million against RM 193 million, implying a\n")
cat("  rate of RM 4.29/USD. Confirm the rate and year, and state it in the note.\n")

write.csv(scen, "Output/Gap5_welfare_sensitivity.csv", row.names = FALSE)
write.csv(by_year, "Output/Gap5_table13_rebuilt.csv", row.names = FALSE)
cat("\nSaved: Output/Gap5_welfare_sensitivity.csv\n")
cat("Saved: Output/Gap5_table13_rebuilt.csv\n")

#===============================================================================
# RECONCILIATION
#===============================================================================

cat("\n===============================================================\n")
cat("RECONCILIATION WITH SECTION 6.7 AND TABLE 13\n")
cat("===============================================================\n")

cat("\n*** THE DENOMINATOR PROBLEM - RESOLVE THIS BEFORE SUBMISSION ***\n")
cat("Table 1 reports 47,356 total deaths in the regression sample.\n")
cat("Table 13 reports observed deaths of 764 + 1,653 + 1,767 + 837 = 5,021,\n")
cat("and computes '2.6% of all deaths' as 128.4 / 5,021.\n")
cat("But the Table 13 note says the column is a subset of the 47,356.\n")
cat("Against 47,356 the share is 0.27%, not 2.6% - a factor of ten.\n")
cat("This number appears in Section 6.3 and drives the abstract's framing.\n\n")
cat(sprintf("This run computes observed deaths across the SAME station-days\n"))
cat(sprintf("used for the counterfactual: %s\n",
            format(sum(by_year$Observed_deaths), big.mark = ",")))
cat("Compare that against both 5,021 and 47,356 to find which subset Table 13\n")
cat("was actually built from, then fix the percentage to match.\n")

cat("\n*** ON THE HAZE-SEASON SCENARIO ***\n")
cat("Section 6.7 says the haze-only figure is 'lower than the full-year\n")
cat("estimate because we exclude non-haze-season deaths, but based on more\n")
cat("precisely estimated coefficients'. Note the haze-season beta here is\n")
cat("LARGER than the full-sample beta, so the reduction comes purely from\n")
cat("dropping eight months of exposure, not from a smaller coefficient.\n")
cat("Phrase it that way to avoid implying the coefficient fell.\n")

cat("\nGap 5 complete.\n")


