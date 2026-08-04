################################################################################
# 94_FireDays_Above5.R
#
# GAP 4 of 5  ->  MANUSCRIPT TABLE 12, FINAL COLUMN "Days > 5 ug/m3"
#
# Target values to reproduce:  412 / 687 / 2,341 / 198, total 3,638
#
# Why this script exists: mean fire PM2.5, max fire PM2.5 and fire share are
# all coded (0202_first.R, revised_counterfactual.R), but the count of
# station-days where fire-attributable pollution exceeds the WHO annual
# guideline of 5 ug/m3 is not computed anywhere.
#
# This script rebuilds ALL of Table 12 rather than just the missing column,
# so the whole table comes from one place and one first stage.
#
# CRITICAL DECISION - TRUNCATION AT ZERO:
#   fire_pm25 = a1*(FRP x u) + a2*(FRP x v) is negative whenever wind blows
#   pollution away from a station. 0202_first.R applies pmax(., 0);
#   revised_counterfactual.R does the same. Truncating changes the mean,
#   the share, AND the day count. This script reports BOTH so you can state
#   in the note which convention Table 12 uses. They are not interchangeable.
#
# Output: Output/Gap4_table12_fire_attributable.csv
################################################################################

source("90_Common_Setup.R")

cat("################################################################\n")
cat("# GAP 4: FIRE-ATTRIBUTABLE POLLUTION (TABLE 12)\n")
cat("################################################################\n\n")

merged_daily <- prep_analysis_data()

WHO_GUIDELINE <- 5   # ug/m3, WHO annual PM2.5 guideline

table12_all <- list()

for (poll in c("PM2.5AVG", "PM10AVG")) {
  
  cat("===============================================================\n")
  cat(poll, "\n")
  cat("===============================================================\n\n")
  
  d <- complete_sample(merged_daily, poll)
  
  #---------------------------------------------------------------------------
  # 1. FIRST STAGE -> fire-attributable pollution
  #---------------------------------------------------------------------------
  fs <- lm(as.formula(paste(
    poll, "~ FRP_u_wind_station + FRP_v_wind_station +",
    "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
    "factor(month) + factor(year) + factor(station_id)")), data = d)
  
  a1 <- coef(fs)[["FRP_u_wind_station"]]
  a2 <- coef(fs)[["FRP_v_wind_station"]]
  
  cat(sprintf("  First-stage coefficients: a1 (u-wind) = %.3e, a2 (v-wind) = %.3e\n",
              a1, a2))
  cat("  Manuscript Section 6.1 states a1 = 0.000226, a2 = -0.000443 for PM2.5\n\n")
  
  # Delta_Pollutant^fire = a1*(FRP x u) + a2*(FRP x v)      [manuscript S6.1]
  d$fire_raw  <- a1 * d$FRP_u_wind_station + a2 * d$FRP_v_wind_station
  d$fire_trunc <- pmax(d$fire_raw, 0)
  
  cat(sprintf("  Untruncated fire pollution: mean %.3f, min %.2f, max %.2f\n",
              mean(d$fire_raw), min(d$fire_raw), max(d$fire_raw)))
  cat(sprintf("  Truncated at zero         : mean %.3f, min %.2f, max %.2f\n",
              mean(d$fire_trunc), min(d$fire_trunc), max(d$fire_trunc)))
  cat(sprintf("  Share of station-days with NEGATIVE predicted fire pollution: %.1f%%\n\n",
              100 * mean(d$fire_raw < 0)))
  
  #---------------------------------------------------------------------------
  # 2. TABLE 12 BY YEAR - both truncation conventions
  #---------------------------------------------------------------------------
  build_t12 <- function(data, fire_col, tag) {
    data %>%
      group_by(year) %>%
      summarise(
        N_obs          = n(),
        Mean_pollutant = mean(.data[[poll]], na.rm = TRUE),
        Mean_fire      = mean(.data[[fire_col]], na.rm = TRUE),
        Max_fire       = max(.data[[fire_col]], na.rm = TRUE),
        # THE MISSING COLUMN:
        Days_above_5   = sum(.data[[fire_col]] > WHO_GUIDELINE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Fire_share_pct = 100 * Mean_fire / Mean_pollutant,
        Pct_days_above_5 = 100 * Days_above_5 / N_obs,
        Truncation = tag,
        Pollutant  = poll
      )
  }
  
  t12_trunc <- build_t12(d, "fire_trunc", "truncated at 0")
  t12_raw   <- build_t12(d, "fire_raw",   "untruncated")
  
  add_total <- function(tt, fire_col) {
    tot <- data.frame(
      year             = "Overall",
      N_obs            = nrow(d),
      Mean_pollutant   = mean(d[[poll]], na.rm = TRUE),
      Mean_fire        = mean(d[[fire_col]], na.rm = TRUE),
      Max_fire         = max(d[[fire_col]], na.rm = TRUE),
      Days_above_5     = sum(d[[fire_col]] > WHO_GUIDELINE, na.rm = TRUE),
      Fire_share_pct   = NA,
      Pct_days_above_5 = NA,
      Truncation       = unique(tt$Truncation),
      Pollutant        = poll,
      stringsAsFactors = FALSE
    )
    tot$Fire_share_pct   <- 100 * tot$Mean_fire / tot$Mean_pollutant
    tot$Pct_days_above_5 <- 100 * tot$Days_above_5 / tot$N_obs
    rbind(tt %>% mutate(year = as.character(year)), tot)
  }
  
  t12_trunc <- add_total(t12_trunc, "fire_trunc")
  t12_raw   <- add_total(t12_raw,   "fire_raw")
  
  cat("  ---- TABLE 12, truncated at zero (matches 0202_first.R convention) ----\n")
  print(t12_trunc %>%
          mutate(across(c(Mean_pollutant, Mean_fire, Max_fire,
                          Fire_share_pct, Pct_days_above_5), ~ round(.x, 2))) %>%
          dplyr::select(year, Mean_pollutant, Mean_fire, Fire_share_pct,
                        Max_fire, Days_above_5, Pct_days_above_5))
  
  cat("\n  ---- TABLE 12, untruncated ----\n")
  print(t12_raw %>%
          mutate(across(c(Mean_pollutant, Mean_fire, Max_fire,
                          Fire_share_pct, Pct_days_above_5), ~ round(.x, 2))) %>%
          dplyr::select(year, Mean_pollutant, Mean_fire, Fire_share_pct,
                        Max_fire, Days_above_5))
  cat("\n")
  
  table12_all[[paste0(poll, "_trunc")]] <- t12_trunc
  table12_all[[paste0(poll, "_raw")]]   <- t12_raw
}

#===============================================================================
# SAVE
#===============================================================================

out <- do.call(rbind, table12_all)
rownames(out) <- NULL
write.csv(out, "Output/Gap4_table12_fire_attributable.csv", row.names = FALSE)
cat("Saved: Output/Gap4_table12_fire_attributable.csv\n")

#===============================================================================
# RECONCILIATION
#===============================================================================

cat("\n===============================================================\n")
cat("RECONCILIATION WITH TABLE 12\n")
cat("===============================================================\n")
cat("Manuscript reports, for PM2.5:\n")
cat("   2017: 412   2018: 687   2019: 2,341   2020: 198   Total: 3,638\n")
cat("   and states this is 6.5%% of observations.\n\n")

chk <- table12_all[["PM2.5AVG_trunc"]]
tot <- chk$Days_above_5[chk$year == "Overall"]
cat(sprintf("This run (truncated)  : total = %d, i.e. %.1f%% of %s station-days\n",
            tot, 100 * tot / nrow(complete_sample(merged_daily, "PM2.5AVG")),
            format(nrow(complete_sample(merged_daily, "PM2.5AVG")), big.mark = ",")))

chk2 <- table12_all[["PM2.5AVG_raw"]]
cat(sprintf("This run (untruncated): total = %d\n",
            chk2$Days_above_5[chk2$year == "Overall"]))

cat("\nCHECK: 3,638 / 55,719 = 6.53%%, so the manuscript's 6.5%% is internally\n")
cat("consistent with its own count. If neither number above lands near 3,638,\n")
cat("the discrepancy is the truncation rule or the first-stage coefficients -\n")
cat("compare the a1/a2 printed above against the 0.000226 / -0.000443 quoted\n")
cat("in Section 6.1 before changing anything else.\n")

cat("\nGap 4 complete.\n")

