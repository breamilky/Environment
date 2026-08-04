#===============================================================================
# PM10 ROBUSTNESS ANALYSIS - Parallel to PM2.5 Analysis
# For JHE Resubmission
#===============================================================================

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("PM10 COMPREHENSIVE ROBUSTNESS ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

#-------------------------------------------------------------------------------
# SETUP
#-------------------------------------------------------------------------------

packages <- c("dplyr", "tidyr", "ggplot2", "zoo", "lubridate")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

set.seed(42)

# Create output directories
dir.create("Output/PM10_Robustness/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("Output/PM10_Robustness/Figures", recursive = TRUE, showWarnings = FALSE)

# Load data
if (file.exists("Output/RDS_files/script02_all_results.RData")) {
  load("Output/RDS_files/script02_all_results.RData")
  cat("✓ Loaded data from script02_all_results.RData\n")
} else {
  stop("No saved data found. Run Script 02 first.")
}

# Prepare data
analysis_data <- merged_daily_pm10 %>%
  mutate(
    date = as.Date(date),
    station_id = as.factor(LOCATION_clean),
    year = lubridate::year(date),
    month = lubridate::month(date)
  )

# Get PM10 baseline coefficient from results_pm10
if (exists("results_pm10")) {
  beta_pm10 <- results_pm10$coefficient
  se_pm10 <- results_pm10$se_twoway
  effect_pm10 <- results_pm10$effect_10ugm3
  n_obs_pm10 <- results_pm10$n_obs
} else {
  # Run baseline PM10 model
  cat("Running baseline PM10 IV-CF model...\n")
  
  pm10_data <- analysis_data %>% filter(!is.na(PM10AVG))
  
  fs_pm10 <- lm(PM10AVG ~ FRP_u_wind_station + FRP_v_wind_station +
                  RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
                  factor(month) + factor(year) + factor(station_id),
                data = pm10_data)
  pm10_data$cf_resid <- residuals(fs_pm10)
  
  ss_pm10 <- glm(deaths ~ PM10AVG + cf_resid +
                   RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
                   factor(month) + factor(year) + factor(station_id),
                 family = poisson(link = "log"), data = pm10_data)
  
  beta_pm10 <- coef(ss_pm10)["PM10AVG"]
  se_pm10 <- sqrt(vcov(ss_pm10)["PM10AVG", "PM10AVG"])
  effect_pm10 <- (exp(beta_pm10 * 10) - 1) * 100
  n_obs_pm10 <- nobs(ss_pm10)
}

cat("\n")
cat("PM10 Baseline Results:\n")
cat("  Coefficient: β =", sprintf("%.6f", beta_pm10), "\n")
cat("  SE:", sprintf("%.6f", se_pm10), "\n")
cat("  Effect per 10 µg/m³:", sprintf("%.1f%%", effect_pm10), "\n")
cat("  N:", format(n_obs_pm10, big.mark = ","), "\n\n")

#-------------------------------------------------------------------------------
# STANDARD IV-CF FUNCTION FOR PM10
#-------------------------------------------------------------------------------

run_iv_cf_pm10 <- function(data, instruments = c("FRP_u_wind_station", "FRP_v_wind_station")) {
  
  tryCatch({
    # Filter for PM10
    data <- data %>% filter(!is.na(PM10AVG))
    
    if (nrow(data) < 1000) stop("Insufficient data")
    
    # First stage
    fs_formula <- as.formula(paste(
      "PM10AVG ~", paste(instruments, collapse = " + "), "+",
      "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
      "factor(month) + factor(year) + factor(station_id)"
    ))
    fs <- lm(fs_formula, data = data)
    data$cf_resid <- residuals(fs)
    
    # Second stage
    ss <- glm(deaths ~ PM10AVG + cf_resid +
                RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
                factor(month) + factor(year) + factor(station_id),
              family = poisson(link = "log"), data = data)
    
    coef_val <- coef(ss)["PM10AVG"]
    se_val <- sqrt(vcov(ss)["PM10AVG", "PM10AVG"])
    
    list(
      coefficient = as.numeric(coef_val),
      se = as.numeric(se_val),
      p_value = as.numeric(2 * pnorm(-abs(coef_val / se_val))),
      effect_10 = as.numeric((exp(coef_val * 10) - 1) * 100),
      n_obs = nobs(ss),
      success = TRUE
    )
  }, error = function(e) {
    list(coefficient = NA, se = NA, p_value = NA, effect_10 = NA, n_obs = NA, success = FALSE)
  })
}

#===============================================================================
# SECTION 1: PLACEBO TESTS FOR PM10
#===============================================================================

cat("─── Section 1: PM10 Placebo Outcome Tests ───\n\n")

pm10_data <- analysis_data %>% filter(!is.na(PM10AVG))

placebo_pm10 <- list()

# 1.1 Actual deaths (reference)
placebo_pm10[["Actual"]] <- data.frame(
  Outcome = "Actual Deaths (Reference)",
  Coefficient = beta_pm10,
  SE = se_pm10,
  P_value = 2 * pnorm(-abs(beta_pm10 / se_pm10)),
  Result = "SIGNIFICANT"
)

cat("Actual deaths:   β =", sprintf("%.6f", beta_pm10), "✓\n")

# 1.2 Shuffled deaths (within station)
set.seed(12345)
pm10_data <- pm10_data %>%
  group_by(station_id) %>%
  mutate(deaths_shuffled = sample(deaths)) %>%
  ungroup()

temp_shuffled <- pm10_data
temp_shuffled$deaths_original <- temp_shuffled$deaths
temp_shuffled$deaths <- temp_shuffled$deaths_shuffled

shuffled_result <- run_iv_cf_pm10(temp_shuffled)

placebo_pm10[["Shuffled"]] <- data.frame(
  Outcome = "Shuffled Deaths",
  Coefficient = shuffled_result$coefficient,
  SE = shuffled_result$se,
  P_value = shuffled_result$p_value,
  Result = ifelse(shuffled_result$p_value > 0.05, "PASS", "CHECK")
)

cat("Shuffled deaths: β =", sprintf("%.6f", shuffled_result$coefficient), 
    ", p =", sprintf("%.4f", shuffled_result$p_value),
    ifelse(shuffled_result$p_value > 0.05, " ✓ PASS", " ⚠"), "\n")

# 1.3 Random deaths (global shuffle)
set.seed(54321)
pm10_data$deaths_random <- sample(pm10_data$deaths)

temp_random <- pm10_data
temp_random$deaths_original <- temp_random$deaths
temp_random$deaths <- temp_random$deaths_random

random_result <- run_iv_cf_pm10(temp_random)

placebo_pm10[["Random"]] <- data.frame(
  Outcome = "Random Deaths",
  Coefficient = random_result$coefficient,
  SE = random_result$se,
  P_value = random_result$p_value,
  Result = ifelse(random_result$p_value > 0.05, "PASS", "CHECK")
)

cat("Random deaths:   β =", sprintf("%.6f", random_result$coefficient), 
    ", p =", sprintf("%.4f", random_result$p_value),
    ifelse(random_result$p_value > 0.05, " ✓ PASS", " ⚠"), "\n")

# 1.4 Cross-year placebo (2019 pollution → 2018 deaths)
cross_year_pm10 <- tryCatch({
  data_2019 <- pm10_data %>%
    filter(year == 2019) %>%
    mutate(doy = yday(date)) %>%
    select(station_id, doy, PM10_2019 = PM10AVG,
           FRP_u_2019 = FRP_u_wind_station, FRP_v_2019 = FRP_v_wind_station,
           RH = RELATIVEHUMIDITYAVG, Temp = AmbientTemperatureAVG)
  
  data_2018 <- pm10_data %>%
    filter(year == 2018) %>%
    mutate(doy = yday(date)) %>%
    select(station_id, doy, deaths, month)
  
  cross_data <- inner_join(data_2019, data_2018, by = c("station_id", "doy"))
  
  if (nrow(cross_data) < 1000) stop("Insufficient data")
  
  fs <- lm(PM10_2019 ~ FRP_u_2019 + FRP_v_2019 + RH + Temp + factor(month) + factor(station_id), 
           data = cross_data)
  cross_data$cf_resid <- residuals(fs)
  
  ss <- glm(deaths ~ PM10_2019 + cf_resid + RH + Temp + factor(month) + factor(station_id),
            family = poisson, data = cross_data)
  
  coef_val <- coef(ss)["PM10_2019"]
  se_val <- sqrt(vcov(ss)["PM10_2019", "PM10_2019"])
  p_val <- 2 * pnorm(-abs(coef_val / se_val))
  
  list(coefficient = coef_val, se = se_val, p_value = p_val, success = TRUE)
}, error = function(e) list(coefficient = NA, se = NA, p_value = NA, success = FALSE))

if (cross_year_pm10$success) {
  placebo_pm10[["CrossYear"]] <- data.frame(
    Outcome = "Cross-Year (2019→2018)",
    Coefficient = cross_year_pm10$coefficient,
    SE = cross_year_pm10$se,
    P_value = cross_year_pm10$p_value,
    Result = ifelse(cross_year_pm10$p_value > 0.05, "PASS", "CHECK")
  )
  cat("Cross-year:      β =", sprintf("%.6f", cross_year_pm10$coefficient), 
      ", p =", sprintf("%.4f", cross_year_pm10$p_value),
      ifelse(cross_year_pm10$p_value > 0.05, " ✓ PASS", " ⚠"), "\n")
}

# Save placebo results
placebo_pm10_table <- do.call(rbind, placebo_pm10)
write.csv(placebo_pm10_table, "Output/PM10_Robustness/Tables/PM10_Placebo_Outcomes.csv", row.names = FALSE)
cat("\n✓ Saved: PM10_Placebo_Outcomes.csv\n\n")

#===============================================================================
# SECTION 2: PM10 DISTRIBUTED LAG MODEL (0-14 days)
#===============================================================================

cat("─── Section 2: PM10 Distributed Lag Model ───\n\n")

# Create lag variables for PM10
max_lag <- 14

pm10_lag_data <- pm10_data %>%
  arrange(station_id, date) %>%
  group_by(station_id)

for (k in 1:max_lag) {
  pm10_lag_data <- pm10_lag_data %>%
    mutate(
      !!paste0("PM10_lag", k) := lag(PM10AVG, k),
      !!paste0("FRP_u_lag", k) := lag(FRP_u_wind_station, k),
      !!paste0("FRP_v_lag", k) := lag(FRP_v_wind_station, k)
    )
}

pm10_lag_data <- ungroup(pm10_lag_data)

# Estimate each lag
pm10_lag_results <- data.frame(
  Lag = 0:max_lag,
  Coefficient = NA_real_,
  SE = NA_real_,
  P_value = NA_real_,
  Effect_10ugm3 = NA_real_,
  N_obs = NA_integer_
)

# Lag 0
pm10_lag_results[1, ] <- list(0, beta_pm10, se_pm10, 
                              2 * pnorm(-abs(beta_pm10 / se_pm10)), 
                              effect_pm10, n_obs_pm10)
cat("Lag  0: β =", sprintf("%.5f", beta_pm10), "*\n")

# Lags 1 to max_lag
for (k in 1:max_lag) {
  
  lag_var <- paste0("PM10_lag", k)
  iv_u <- paste0("FRP_u_lag", k)
  iv_v <- paste0("FRP_v_lag", k)
  
  result <- tryCatch({
    lag_data <- pm10_lag_data %>%
      filter(!is.na(!!sym(lag_var)), !is.na(!!sym(iv_u)), !is.na(!!sym(iv_v)))
    
    if (nrow(lag_data) < 1000) stop("Insufficient data")
    
    # First stage with lagged instruments
    fs_formula <- as.formula(paste(lag_var, "~", iv_u, "+", iv_v, "+",
                                   "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
                                   "factor(month) + factor(year) + factor(station_id)"))
    fs <- lm(fs_formula, data = lag_data)
    lag_data$cf_resid <- residuals(fs)
    
    # Second stage
    ss_formula <- as.formula(paste("deaths ~", lag_var, "+ cf_resid +",
                                   "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
                                   "factor(month) + factor(year) + factor(station_id)"))
    ss <- glm(ss_formula, family = poisson, data = lag_data)
    
    coef_val <- coef(ss)[lag_var]
    se_val <- sqrt(vcov(ss)[lag_var, lag_var])
    
    list(coef = coef_val, se = se_val, 
         p = 2*pnorm(-abs(coef_val/se_val)),
         effect = (exp(coef_val*10)-1)*100,
         n = nrow(lag_data))
  }, error = function(e) list(coef = NA, se = NA, p = NA, effect = NA, n = NA))
  
  pm10_lag_results[k+1, ] <- list(k, result$coef, result$se, result$p, result$effect, result$n)
  
  cat("Lag", sprintf("%2d", k), ": β =", sprintf("%.5f", result$coef))
  if (!is.na(result$p) && result$p < 0.05) cat(" *")
  cat("\n")
}

pm10_lag_results$CI_lower <- pm10_lag_results$Coefficient - 1.96 * pm10_lag_results$SE
pm10_lag_results$CI_upper <- pm10_lag_results$Coefficient + 1.96 * pm10_lag_results$SE

write.csv(pm10_lag_results, "Output/PM10_Robustness/Tables/PM10_Distributed_Lag_Model.csv", row.names = FALSE)

# Summary
peak_lag_pm10 <- pm10_lag_results$Lag[which.max(pm10_lag_results$Coefficient)]
n_sig_pm10 <- sum(pm10_lag_results$P_value < 0.05, na.rm = TRUE)

cat("\nPeak at lag", peak_lag_pm10, "days. Significant lags:", n_sig_pm10, "of", max_lag + 1, "\n")
cat("✓ Saved: PM10_Distributed_Lag_Model.csv\n\n")

#===============================================================================
# SECTION 3: PM10 CUMULATIVE EFFECTS (Moving Averages)
#===============================================================================

cat("─── Section 3: PM10 Cumulative Effects ───\n\n")

windows <- c(3, 7, 14)

pm10_cumulative <- data.frame(
  Window = integer(),
  Coefficient = numeric(),
  SE = numeric(),
  P_value = numeric(),
  Effect_10ugm3 = numeric(),
  N_obs = integer()
)

for (w in windows) {
  
  ma_var <- paste0("PM10_ma", w)
  
  # Create moving average
  temp_data <- pm10_data %>%
    arrange(station_id, date) %>%
    group_by(station_id) %>%
    mutate(!!ma_var := zoo::rollmean(PM10AVG, k = w, fill = NA, align = "right")) %>%
    ungroup() %>%
    filter(!is.na(!!sym(ma_var)))
  
  result <- tryCatch({
    fs <- lm(as.formula(paste(ma_var, "~ FRP_u_wind_station + FRP_v_wind_station +",
                              "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
                              "factor(month) + factor(year) + factor(station_id)")),
             data = temp_data)
    temp_data$cf_resid <- residuals(fs)
    
    ss <- glm(as.formula(paste("deaths ~", ma_var, "+ cf_resid +",
                               "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
                               "factor(month) + factor(year) + factor(station_id)")),
              family = poisson, data = temp_data)
    
    coef_val <- coef(ss)[ma_var]
    se_val <- sqrt(vcov(ss)[ma_var, ma_var])
    
    list(
      coefficient = coef_val,
      se = se_val,
      p_value = 2 * pnorm(-abs(coef_val / se_val)),
      effect_10 = (exp(coef_val * 10) - 1) * 100,
      n_obs = nobs(ss),
      success = TRUE
    )
  }, error = function(e) {
    list(coefficient = NA, se = NA, p_value = NA, effect_10 = NA, n_obs = NA, success = FALSE)
  })
  
  pm10_cumulative <- rbind(pm10_cumulative, data.frame(
    Window = w,
    Coefficient = result$coefficient,
    SE = result$se,
    P_value = result$p_value,
    Effect_10ugm3 = result$effect_10,
    N_obs = result$n_obs
  ))
  
  cat(sprintf("%2d-day MA: β = %.5f, effect = %.1f%%\n", w, result$coefficient, result$effect_10))
}

write.csv(pm10_cumulative, "Output/PM10_Robustness/Tables/PM10_Cumulative_Effects.csv", row.names = FALSE)
cat("✓ Saved: PM10_Cumulative_Effects.csv\n\n")

#===============================================================================
# SECTION 4: PM10 HETEROGENEITY (Season & Year)
#===============================================================================

cat("─── Section 4: PM10 Heterogeneity Analysis ───\n\n")

# 4.1 Season interaction (Haze vs Non-Haze)
cat("Season heterogeneity (interaction approach):\n")

pm10_data$haze_season <- ifelse(pm10_data$month %in% c(6,7,8,9), 1, 0)

season_het_pm10 <- tryCatch({
  # First stage with interaction
  fs <- lm(PM10AVG ~ FRP_u_wind_station + FRP_v_wind_station +
             I(FRP_u_wind_station * haze_season) + I(FRP_v_wind_station * haze_season) +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
             factor(month) + factor(year) + factor(station_id),
           data = pm10_data)
  pm10_data$cf_resid <- residuals(fs)
  
  # Second stage with interaction
  ss <- glm(deaths ~ PM10AVG * haze_season + cf_resid +
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
              factor(month) + factor(year) + factor(station_id),
            family = poisson, data = pm10_data)
  
  coef_base <- coef(ss)["PM10AVG"]
  coef_int <- coef(ss)["PM10AVG:haze_season"]
  se_int <- sqrt(vcov(ss)["PM10AVG:haze_season", "PM10AVG:haze_season"])
  
  list(
    effect_nonhaze = (exp(coef_base * 10) - 1) * 100,
    effect_haze = (exp((coef_base + coef_int) * 10) - 1) * 100,
    interaction_p = 2 * pnorm(-abs(coef_int / se_int)),
    success = TRUE
  )
}, error = function(e) list(success = FALSE))

if (season_het_pm10$success) {
  cat("  Non-haze (Oct-May):", sprintf("%.1f%%", season_het_pm10$effect_nonhaze), "\n")
  cat("  Haze (Jun-Sep):    ", sprintf("%.1f%%", season_het_pm10$effect_haze), "\n")
  cat("  Interaction p-value:", sprintf("%.4f", season_het_pm10$interaction_p), "\n")
  
  het_season_pm10 <- data.frame(
    Season = c("Non-Haze (Oct-May)", "Haze (Jun-Sep)"),
    Effect_Pct = c(season_het_pm10$effect_nonhaze, season_het_pm10$effect_haze),
    Interaction_P = c(NA, season_het_pm10$interaction_p)
  )
  write.csv(het_season_pm10, "Output/PM10_Robustness/Tables/PM10_Heterogeneity_Season.csv", row.names = FALSE)
}

# 4.2 Year heterogeneity
cat("\nYear heterogeneity (interaction approach, 2019 as reference):\n")

year_het_pm10 <- tryCatch({
  pm10_data$year_factor <- relevel(factor(pm10_data$year), ref = "2019")
  
  fs <- lm(PM10AVG ~ FRP_u_wind_station + FRP_v_wind_station +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
             factor(month) + year_factor + factor(station_id),
           data = pm10_data)
  pm10_data$cf_resid <- residuals(fs)
  
  ss <- glm(deaths ~ PM10AVG * year_factor + cf_resid +
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
              factor(month) + factor(station_id),
            family = poisson, data = pm10_data)
  
  V <- vcov(ss)
  coef_2019 <- as.numeric(coef(ss)["PM10AVG"])
  se_2019   <- sqrt(V["PM10AVG", "PM10AVG"])
  
  year_effects <- data.frame(
    Year             = "2019 (ref)",
    Coefficient      = coef_2019,
    SE               = se_2019,
    Effect_Pct       = (exp(coef_2019 * 10) - 1) * 100,
    CI_Lower         = coef_2019 - 1.96 * se_2019,
    CI_Upper         = coef_2019 + 1.96 * se_2019,
    Interaction_Coef = 0,
    Interaction_SE   = NA_real_,
    Interaction_P    = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Other years: effect_y = coef_2019 + coef_int, a linear combination.
  # Var(a + b) = Var(a) + Var(b) + 2*Cov(a, b). The covariance term is not
  # optional: the level and interaction coefficients are estimated from
  # overlapping data and are strongly negatively correlated.
  for (y in c("2017", "2018", "2020")) {
    int_name <- paste0("PM10AVG:year_factor", y)
    if (int_name %in% names(coef(ss))) {
      
      coef_int <- as.numeric(coef(ss)[int_name])
      coef_y   <- coef_2019 + coef_int
      
      se_y <- sqrt(V["PM10AVG", "PM10AVG"] +
                     V[int_name, int_name] +
                     2 * V["PM10AVG", int_name])
      
      se_int <- sqrt(V[int_name, int_name])
      p_int  <- 2 * pnorm(-abs(coef_int / se_int))
      
      year_effects <- rbind(year_effects, data.frame(
        Year             = y,
        Coefficient      = coef_y,
        SE               = se_y,
        Effect_Pct       = (exp(coef_y * 10) - 1) * 100,
        CI_Lower         = coef_y - 1.96 * se_y,
        CI_Upper         = coef_y + 1.96 * se_y,
        Interaction_Coef = coef_int,
        Interaction_SE   = se_int,
        Interaction_P    = p_int,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  year_effects <- year_effects[order(year_effects$Year), ]
  list(effects = year_effects, success = TRUE)
}, error = function(e) list(success = FALSE))

if (year_het_pm10$success) {
  for (i in 1:nrow(year_het_pm10$effects)) {
    cat("  ", year_het_pm10$effects$Year[i], ":", sprintf("%.1f%%", year_het_pm10$effects$Effect_Pct[i]), "\n")
  }
  write.csv(year_het_pm10$effects, "Output/PM10_Robustness/Tables/PM10_Heterogeneity_Year.csv", row.names = FALSE)
}

cat("  ", year_het_pm10$effects$Year[i], ":",
    sprintf("%.1f%%", year_het_pm10$effects$Effect_Pct[i]),
    "(SE", sprintf("%.5f", year_het_pm10$effects$SE[i]), 
    "| vs 2019 p =", sprintf("%.3f", year_het_pm10$effects$Interaction_P[i]), ")\n")

#===============================================================================
# SECTION 5: PM10 ALTERNATIVE IV SPECIFICATIONS
#===============================================================================

cat("─── Section 5: PM10 Alternative IV Specifications ───\n\n")

# Calculate F-stats for PM10
calc_f_stat_pm10 <- function(data, instruments) {
  data <- data %>% filter(!is.na(PM10AVG))
  
  fs_full <- lm(as.formula(paste("PM10AVG ~", paste(instruments, collapse = " + "), "+",
                                 "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
                                 "factor(month) + factor(year) + factor(station_id)")),
                data = data)
  fs_restr <- lm(PM10AVG ~ RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
                   factor(month) + factor(year) + factor(station_id), data = data)
  
  n <- nobs(fs_full)
  k <- length(coef(fs_full))
  p <- length(instruments)
  
  rss_full <- sum(resid(fs_full)^2)
  rss_restr <- sum(resid(fs_restr)^2)
  
  ((rss_restr - rss_full) / p) / (rss_full / (n - k))
}

pm10_alt_iv <- data.frame(
  Specification = character(),
  Coefficient = numeric(),
  SE = numeric(),
  First_stage_F = numeric(),
  Effect_10ugm3 = numeric()
)

# 5.1 Both instruments (main)
f_both_pm10 <- calc_f_stat_pm10(pm10_data, c("FRP_u_wind_station", "FRP_v_wind_station"))

pm10_alt_iv <- rbind(pm10_alt_iv, data.frame(
  Specification = "Both instruments (main)",
  Coefficient = beta_pm10,
  SE = se_pm10,
  First_stage_F = f_both_pm10,
  Effect_10ugm3 = effect_pm10
))

cat("Both IV:     β =", sprintf("%.6f", beta_pm10), ", F =", round(f_both_pm10, 0), "\n")

# 5.2 Only FRP × u_wind
result_u_pm10 <- run_iv_cf_pm10(pm10_data, instruments = "FRP_u_wind_station")
f_u_pm10 <- calc_f_stat_pm10(pm10_data, "FRP_u_wind_station")

pm10_alt_iv <- rbind(pm10_alt_iv, data.frame(
  Specification = "Only FRP × u_wind",
  Coefficient = result_u_pm10$coefficient,
  SE = result_u_pm10$se,
  First_stage_F = f_u_pm10,
  Effect_10ugm3 = result_u_pm10$effect_10
))

cat("u_wind only: β =", sprintf("%.6f", result_u_pm10$coefficient), ", F =", round(f_u_pm10, 0), "\n")

# 5.3 Only FRP × v_wind
result_v_pm10 <- run_iv_cf_pm10(pm10_data, instruments = "FRP_v_wind_station")
f_v_pm10 <- calc_f_stat_pm10(pm10_data, "FRP_v_wind_station")

pm10_alt_iv <- rbind(pm10_alt_iv, data.frame(
  Specification = "Only FRP × v_wind",
  Coefficient = result_v_pm10$coefficient,
  SE = result_v_pm10$se,
  First_stage_F = f_v_pm10,
  Effect_10ugm3 = result_v_pm10$effect_10
))

cat("v_wind only: β =", sprintf("%.6f", result_v_pm10$coefficient), ", F =", round(f_v_pm10, 0), "\n")

# 5.4 FRP only (no wind interaction)
result_frp_pm10 <- run_iv_cf_pm10(pm10_data, instruments = "total_frp")
f_frp_pm10 <- calc_f_stat_pm10(pm10_data, "total_frp")

pm10_alt_iv <- rbind(pm10_alt_iv, data.frame(
  Specification = "FRP only (no wind)",
  Coefficient = result_frp_pm10$coefficient,
  SE = result_frp_pm10$se,
  First_stage_F = f_frp_pm10,
  Effect_10ugm3 = result_frp_pm10$effect_10
))

cat("FRP only:    β =", sprintf("%.6f", result_frp_pm10$coefficient), ", F =", round(f_frp_pm10, 0), "\n")

write.csv(pm10_alt_iv, "Output/PM10_Robustness/Tables/PM10_Alternative_IV_Specifications.csv", row.names = FALSE)
cat("✓ Saved: PM10_Alternative_IV_Specifications.csv\n\n")

#===============================================================================
# SECTION 6: PM10 SENSITIVITY TO OUTLIERS
#===============================================================================

cat("─── Section 6: PM10 Sensitivity to Outliers ───\n\n")

pm10_sensitivity <- data.frame(
  Specification = character(),
  Coefficient = numeric(),
  SE = numeric(),
  Effect_10ugm3 = numeric(),
  N_obs = integer(),
  Pct_change = numeric()
)

# 6.1 Original
pm10_sensitivity <- rbind(pm10_sensitivity, data.frame(
  Specification = "Original",
  Coefficient = beta_pm10,
  SE = se_pm10,
  Effect_10ugm3 = effect_pm10,
  N_obs = n_obs_pm10,
  Pct_change = 0
))

cat("Original:        β =", sprintf("%.6f", beta_pm10), "\n")

# 6.2 Winsorized 1%/99%
p01_pm10 <- quantile(pm10_data$PM10AVG, 0.01, na.rm = TRUE)
p99_pm10 <- quantile(pm10_data$PM10AVG, 0.99, na.rm = TRUE)

wins_data_pm10 <- pm10_data %>%
  mutate(PM10AVG = pmin(pmax(PM10AVG, p01_pm10), p99_pm10))

wins_result_pm10 <- run_iv_cf_pm10(wins_data_pm10)

pm10_sensitivity <- rbind(pm10_sensitivity, data.frame(
  Specification = "Winsorized 1%/99%",
  Coefficient = wins_result_pm10$coefficient,
  SE = wins_result_pm10$se,
  Effect_10ugm3 = wins_result_pm10$effect_10,
  N_obs = wins_result_pm10$n_obs,
  Pct_change = (wins_result_pm10$coefficient / beta_pm10 - 1) * 100
))

cat("Winsorized:      β =", sprintf("%.6f", wins_result_pm10$coefficient),
    "(", sprintf("%+.1f%%", (wins_result_pm10$coefficient/beta_pm10-1)*100), ")\n")

# 6.3 Exclude PM10 > 150
excl150_result_pm10 <- run_iv_cf_pm10(pm10_data %>% filter(PM10AVG <= 150))

pm10_sensitivity <- rbind(pm10_sensitivity, data.frame(
  Specification = "Exclude PM10 > 150",
  Coefficient = excl150_result_pm10$coefficient,
  SE = excl150_result_pm10$se,
  Effect_10ugm3 = excl150_result_pm10$effect_10,
  N_obs = excl150_result_pm10$n_obs,
  Pct_change = (excl150_result_pm10$coefficient / beta_pm10 - 1) * 100
))

cat("Exclude >150:    β =", sprintf("%.6f", excl150_result_pm10$coefficient),
    "(", sprintf("%+.1f%%", (excl150_result_pm10$coefficient/beta_pm10-1)*100), ")\n")

# 6.4 Exclude PM10 > 100
excl100_result_pm10 <- run_iv_cf_pm10(pm10_data %>% filter(PM10AVG <= 100))

pm10_sensitivity <- rbind(pm10_sensitivity, data.frame(
  Specification = "Exclude PM10 > 100",
  Coefficient = excl100_result_pm10$coefficient,
  SE = excl100_result_pm10$se,
  Effect_10ugm3 = excl100_result_pm10$effect_10,
  N_obs = excl100_result_pm10$n_obs,
  Pct_change = (excl100_result_pm10$coefficient / beta_pm10 - 1) * 100
))

cat("Exclude >100:    β =", sprintf("%.6f", excl100_result_pm10$coefficient),
    "(", sprintf("%+.1f%%", (excl100_result_pm10$coefficient/beta_pm10-1)*100), ")\n")

write.csv(pm10_sensitivity, "Output/PM10_Robustness/Tables/PM10_Sensitivity_Outliers.csv", row.names = FALSE)
cat("✓ Saved: PM10_Sensitivity_Outliers.csv\n\n")

#===============================================================================
# SECTION 7: PM10 COUNTERFACTUAL ANALYSIS
#===============================================================================

cat("─── Section 7: PM10 Counterfactual Analysis ───\n\n")

# First stage for fire-attributable PM10
fs_pm10_cf <- lm(PM10AVG ~ FRP_u_wind_station + FRP_v_wind_station +
                   RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
                   factor(month) + factor(year) + factor(station_id),
                 data = pm10_data)

coef_u_pm10 <- coef(fs_pm10_cf)["FRP_u_wind_station"]
coef_v_pm10 <- coef(fs_pm10_cf)["FRP_v_wind_station"]

cat("  γ_u (FRP × u_wind):", sprintf("%.2e", coef_u_pm10), "\n")
cat("  γ_v (FRP × v_wind):", sprintf("%.2e", coef_v_pm10), "\n")

# Fire-attributable PM10
pm10_data$fire_pollution_pm10 <- pmax(
  coef_u_pm10 * pm10_data$FRP_u_wind_station + coef_v_pm10 * pm10_data$FRP_v_wind_station, 0)

mean_fire_pm10 <- mean(pm10_data$fire_pollution_pm10, na.rm = TRUE)
mean_total_pm10 <- mean(pm10_data$PM10AVG, na.rm = TRUE)
pct_from_fire_pm10 <- (mean_fire_pm10 / mean_total_pm10) * 100

cat("  Mean fire PM10:", sprintf("%.2f", mean_fire_pm10), "µg/m³\n")
cat("  Mean total PM10:", sprintf("%.2f", mean_total_pm10), "µg/m³\n")
cat("  Fire share:", sprintf("%.1f%%", pct_from_fire_pm10), "\n\n")

# Excess deaths
pm10_data$excess_pm10 <- pm10_data$deaths * (1 - exp(-beta_pm10 * pm10_data$fire_pollution_pm10))
pm10_data$excess_pm10[is.na(pm10_data$excess_pm10)] <- 0
pm10_data$excess_pm10[pm10_data$excess_pm10 < 0] <- 0

# Aggregate by year
pm10_cf_by_year <- pm10_data %>%
  group_by(year) %>%
  summarise(
    n_obs = n(),
    total_deaths = sum(deaths, na.rm = TRUE),
    excess_deaths = sum(excess_pm10, na.rm = TRUE),
    mean_fire_poll = mean(fire_pollution_pm10, na.rm = TRUE),
    mean_total_poll = mean(PM10AVG, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_attributable = (excess_deaths / total_deaths) * 100,
    fire_share = (mean_fire_poll / mean_total_poll) * 100
  )

# Totals
total_excess_pm10 <- sum(pm10_cf_by_year$excess_deaths)
total_deaths_pm10 <- sum(pm10_cf_by_year$total_deaths)
pct_attributable_pm10 <- (total_excess_pm10 / total_deaths_pm10) * 100

cat("Results by year:\n")
print(pm10_cf_by_year %>% 
        select(year, total_deaths, excess_deaths, pct_attributable) %>%
        mutate(excess_deaths = round(excess_deaths, 1),
               pct_attributable = round(pct_attributable, 2)))

cat("\n─────────────────────────────────────\n")
cat("TOTAL EXCESS DEATHS (PM10):", round(total_excess_pm10, 0), "\n")
cat("Percent attributable:", sprintf("%.2f%%", pct_attributable_pm10), "\n")

# Confidence intervals
deriv_pm10 <- sum(pm10_data$deaths * pm10_data$fire_pollution_pm10 * 
                    exp(-beta_pm10 * pm10_data$fire_pollution_pm10), na.rm = TRUE)
se_excess_pm10 <- abs(deriv_pm10) * se_pm10

ci_lower_pm10 <- max(0, total_excess_pm10 - 1.96 * se_excess_pm10)
ci_upper_pm10 <- total_excess_pm10 + 1.96 * se_excess_pm10

cat("95% CI: [", round(ci_lower_pm10, 0), ", ", round(ci_upper_pm10, 0), "]\n")

# Welfare
VSL <- 1500000
welfare_pm10 <- total_excess_pm10 * VSL / 1e6

cat("\nWelfare loss (central VSL): RM", round(welfare_pm10, 1), "million\n")

# Save
write.csv(pm10_cf_by_year, "Output/PM10_Robustness/Tables/PM10_Counterfactual_by_Year.csv", row.names = FALSE)

pm10_cf_summary <- data.frame(
  Metric = c("IV Coefficient", "SE", "Effect per 10 µg/m³", 
             "Mean fire PM10 (µg/m³)", "Fire share (%)", 
             "Total deaths", "Excess deaths", "95% CI lower", "95% CI upper",
             "Welfare loss (RM million)"),
  Value = c(beta_pm10, se_pm10, effect_pm10, 
            mean_fire_pm10, pct_from_fire_pm10,
            total_deaths_pm10, total_excess_pm10, ci_lower_pm10, ci_upper_pm10,
            welfare_pm10)
)
write.csv(pm10_cf_summary, "Output/PM10_Robustness/Tables/PM10_Counterfactual_Summary.csv", row.names = FALSE)

cat("✓ Saved: PM10_Counterfactual_by_Year.csv, PM10_Counterfactual_Summary.csv\n\n")

#===============================================================================
# SUMMARY
#===============================================================================

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("PM10 ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("Main Results:\n")
cat("  Coefficient: β =", sprintf("%.6f", beta_pm10), "\n")
cat("  Effect per 10 µg/m³:", sprintf("%.1f%%", effect_pm10), "\n")
cat("  First-stage F:", round(f_both_pm10, 0), "\n")
cat("  Peak lag:", peak_lag_pm10, "days\n")
cat("  Significant lags:", n_sig_pm10, "of 15\n")

if (season_het_pm10$success) {
  cat("\nHeterogeneity:\n")
  cat("  Haze season effect:", sprintf("%.1f%%", season_het_pm10$effect_haze), "\n")
  cat("  Non-haze effect:", sprintf("%.1f%%", season_het_pm10$effect_nonhaze), "\n")
  cat("  Interaction p:", sprintf("%.4f", season_het_pm10$interaction_p), "\n")
}

cat("\nCounterfactual:\n")
cat("  Excess deaths:", round(total_excess_pm10, 0), "\n")
cat("  95% CI: [", round(ci_lower_pm10, 0), ", ", round(ci_upper_pm10, 0), "]\n")
cat("  Welfare loss: RM", round(welfare_pm10, 1), "million\n")

cat("\nFiles saved to: Output/PM10_Robustness/Tables/\n")