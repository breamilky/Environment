################################################################################
#
#
#  Purpose: Complete robustness checks for IV-CF analysis of PM2.5 and mortality
#
#  Contents:
#    PART A: SETUP & DATA
#    PART B: BASELINE MODEL
#    PART C: PLACEBO OUTCOME TESTS
#    PART D: DISTRIBUTED LAG MODEL (0-14 days)
#    PART E: CUMULATIVE EFFECTS (Moving Averages)
#    PART F: HETEROGENEITY ANALYSIS (Interaction Approach)
#    PART G: SENSITIVITY TO OUTLIERS
#    PART H: ALTERNATIVE IV SPECIFICATIONS
#    PART I: FIGURES
#    PART J: SUMMARY
#  Date: February 2026
#
################################################################################

rm(list = ls())
gc()

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║       COMPREHENSIVE ROBUSTNESS ANALYSES - FINAL VERSION                  ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")


################################################################################
#                                                                              #
#  PART A: SETUP & DATA                                                        #
#                                                                              #
################################################################################

cat("================================================================\n")
cat("PART A: SETUP & DATA\n")
cat("================================================================\n\n")

# Load packages
required_packages <- c("dplyr", "tidyr", "ggplot2", "zoo", "lubridate", "broom")

cat("Loading packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

# Set seed for reproducibility
set.seed(42)

# Create output directories
dir.create("Output/Robustness/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("Output/Robustness/Figures", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
cat("\nLoading data...\n")

data_paths <- c(
  "Output/RDS_files/merged_daily_pm25.rds",
  "Data/Mid_process_data/merged_death_data_from_script_01.csv",
  "merged_death_data_from_script_01.csv"
)

data_loaded <- FALSE

# Try RDS first (faster)
if (file.exists("Output/RDS_files/merged_daily_pm25.rds")) {
  analysis_data <- readRDS("Output/RDS_files/merged_daily_pm25.rds")
  cat("  Loaded: Output/RDS_files/merged_daily_pm25.rds\n")
  data_loaded <- TRUE
}

# Fall back to CSV
if (!data_loaded) {
  for (path in data_paths) {
    if (file.exists(path)) {
      analysis_data <- read.csv(path, stringsAsFactors = FALSE)
      cat("  Loaded:", path, "\n")
      data_loaded <- TRUE
      break
    }
  }
}

if (!data_loaded) {
  stop("ERROR: Could not find data file!")
}

# -----------------------------------------------------------------------------
# Prepare Data
# -----------------------------------------------------------------------------
cat("\nPreparing data...\n")

analysis_data <- analysis_data %>%
  mutate(
    date = as.Date(date),
    station_id = as.factor(LOCATION_clean),
    year = lubridate::year(date),
    month = lubridate::month(date)
  )

# Create wind interaction instruments if not present
if (!"FRP_u_wind_station" %in% names(analysis_data)) {
  cat("  Creating wind interaction instruments...\n")
  analysis_data <- analysis_data %>%
    mutate(
      u_wind = WINDSPEEDAVG * sin(WINDDIRECTIONAVG * pi / 180),
      v_wind = WINDSPEEDAVG * cos(WINDDIRECTIONAVG * pi / 180),
      FRP_u_wind_station = total_frp * u_wind,
      FRP_v_wind_station = total_frp * v_wind
    )
}

# Filter to complete cases
analysis_data <- analysis_data %>%
  filter(
    !is.na(PM2.5AVG),
    !is.na(deaths),
    !is.na(FRP_u_wind_station),
    !is.na(FRP_v_wind_station),
    !is.na(RELATIVEHUMIDITYAVG),
    !is.na(AmbientTemperatureAVG)
  )

cat("\nData Summary:\n")
cat("  Observations:", format(nrow(analysis_data), big.mark = ","), "\n")
cat("  Stations:", length(unique(analysis_data$station_id)), "\n")
cat("  Date range:", as.character(min(analysis_data$date)), "to",
    as.character(max(analysis_data$date)), "\n")
cat("  Total deaths:", format(sum(analysis_data$deaths), big.mark = ","), "\n")


################################################################################
#                                                                              #
#  PART B: BASELINE MODEL                                                      #
#                                                                              #
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("PART B: BASELINE MODEL\n")
cat("================================================================\n\n")

# First stage
fs_baseline <- lm(
  PM2.5AVG ~ FRP_u_wind_station + FRP_v_wind_station +
    RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
    factor(month) + factor(year) + factor(station_id),
  data = analysis_data
)
analysis_data$cf_resid <- residuals(fs_baseline)

# Second stage
ss_baseline <- glm(
  deaths ~ PM2.5AVG + cf_resid +
    RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
    factor(month) + factor(year) + factor(station_id),
  family = poisson(link = "log"),
  data = analysis_data
)

beta_baseline <- coef(ss_baseline)["PM2.5AVG"]
se_baseline <- sqrt(vcov(ss_baseline)["PM2.5AVG", "PM2.5AVG"])
effect_baseline <- (exp(beta_baseline * 10) - 1) * 100

cat("Baseline IV-CF Result:\n")
cat("  β =", sprintf("%.6f", beta_baseline), "\n")
cat("  SE =", sprintf("%.6f", se_baseline), "\n")
cat("  Effect per 10 µg/m³:", sprintf("%.1f%%", effect_baseline), "\n")


################################################################################
#                                                                              #
#  PART C: PLACEBO OUTCOME TESTS                                               #
#                                                                              #
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("PART C: PLACEBO OUTCOME TESTS\n")
cat("================================================================\n\n")

cat("Logic: If our IV is valid, pollution should affect REAL deaths\n")
cat("       but NOT fake/shuffled deaths.\n\n")

placebo_results <- list()

# -----------------------------------------------------------------------------
# C.1 Actual Deaths (Reference)
# -----------------------------------------------------------------------------
cat("C.1 Actual Deaths (Reference):\n")

placebo_results[["Actual"]] <- data.frame(
  Outcome = "Actual Deaths (Reference)",
  Coefficient = beta_baseline,
  SE = se_baseline,
  P_value = 2 * pnorm(-abs(beta_baseline / se_baseline)),
  Significant = TRUE,
  Interpretation = "Should be significant"
)
cat("    β =", sprintf("%.6f", beta_baseline), "✓\n")

# -----------------------------------------------------------------------------
# C.2 Shuffled Deaths (within station)
# -----------------------------------------------------------------------------
cat("\nC.2 Shuffled Deaths (within station):\n")

set.seed(12345)
analysis_data <- analysis_data %>%
  group_by(station_id) %>%
  mutate(deaths_shuffled = sample(deaths)) %>%
  ungroup()

model_shuffled <- glm(
  deaths_shuffled ~ PM2.5AVG + cf_resid +
    RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
    factor(month) + factor(year) + factor(station_id),
  family = poisson(link = "log"),
  data = analysis_data
)

coef_shuffled <- coef(model_shuffled)["PM2.5AVG"]
se_shuffled <- sqrt(vcov(model_shuffled)["PM2.5AVG", "PM2.5AVG"])
p_shuffled <- 2 * pnorm(-abs(coef_shuffled / se_shuffled))

cat("    β =", sprintf("%.6f", coef_shuffled), 
    ", p =", sprintf("%.4f", p_shuffled))
if (p_shuffled > 0.05) cat(" ✓ PASS\n") else cat(" ⚠️\n")

placebo_results[["Shuffled"]] <- data.frame(
  Outcome = "Shuffled Deaths (within station)",
  Coefficient = coef_shuffled,
  SE = se_shuffled,
  P_value = p_shuffled,
  Significant = p_shuffled < 0.05,
  Interpretation = ifelse(p_shuffled > 0.05, "PASS", "CHECK")
)

# -----------------------------------------------------------------------------
# C.3 Random Deaths (global shuffle)
# -----------------------------------------------------------------------------
cat("\nC.3 Random Deaths (global shuffle):\n")

set.seed(54321)
analysis_data$deaths_random <- sample(analysis_data$deaths)

model_random <- glm(
  deaths_random ~ PM2.5AVG + cf_resid +
    RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
    factor(month) + factor(year) + factor(station_id),
  family = poisson(link = "log"),
  data = analysis_data
)

coef_random <- coef(model_random)["PM2.5AVG"]
se_random <- sqrt(vcov(model_random)["PM2.5AVG", "PM2.5AVG"])
p_random <- 2 * pnorm(-abs(coef_random / se_random))

cat("    β =", sprintf("%.6f", coef_random), 
    ", p =", sprintf("%.4f", p_random))
if (p_random > 0.05) cat(" ✓ PASS\n") else cat(" ⚠️\n")

placebo_results[["Random"]] <- data.frame(
  Outcome = "Random Deaths (global shuffle)",
  Coefficient = coef_random,
  SE = se_random,
  P_value = p_random,
  Significant = p_random < 0.05,
  Interpretation = ifelse(p_random > 0.05, "PASS", "CHECK")
)

# -----------------------------------------------------------------------------
# C.4 Cross-Year Placebo (2019 pollution → 2018 deaths)
# -----------------------------------------------------------------------------
cat("\nC.4 Cross-Year Placebo (2019 pollution → 2018 deaths):\n")

cross_year_result <- tryCatch({
  
  # Get 2019 pollution data
  data_2019 <- analysis_data %>%
    filter(year == 2019) %>%
    mutate(doy = lubridate::yday(date)) %>%
    dplyr::select(station_id, doy, PM2.5AVG_2019 = PM2.5AVG, 
                  cf_resid_2019 = cf_resid,
                  RELATIVEHUMIDITYAVG, AmbientTemperatureAVG)
  
  # Get 2018 deaths
  data_2018 <- analysis_data %>%
    filter(year == 2018) %>%
    mutate(doy = lubridate::yday(date)) %>%
    dplyr::select(station_id, doy, deaths_2018 = deaths, month)
  
  # Merge by station and day-of-year
  cross_data <- inner_join(data_2019, data_2018, by = c("station_id", "doy"))
  
  if (nrow(cross_data) < 1000) {
    cat("    Insufficient matched observations\n")
    return(NULL)
  }
  
  # Run placebo regression
  model_cross <- glm(
    deaths_2018 ~ PM2.5AVG_2019 + cf_resid_2019 +
      RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
      factor(month) + factor(station_id),
    family = poisson(link = "log"),
    data = cross_data
  )
  
  coef_cross <- coef(model_cross)["PM2.5AVG_2019"]
  se_cross <- sqrt(vcov(model_cross)["PM2.5AVG_2019", "PM2.5AVG_2019"])
  p_cross <- 2 * pnorm(-abs(coef_cross / se_cross))
  
  cat("    β =", sprintf("%.6f", coef_cross), 
      ", p =", sprintf("%.4f", p_cross))
  if (p_cross > 0.05) cat(" ✓ PASS\n") else cat(" ⚠️\n")
  
  data.frame(
    Outcome = "Cross-Year (2019 PM2.5 → 2018 deaths)",
    Coefficient = coef_cross,
    SE = se_cross,
    P_value = p_cross,
    Significant = p_cross < 0.05,
    Interpretation = ifelse(p_cross > 0.05, "PASS", "CHECK")
  )
  
}, error = function(e) {
  cat("    FAILED:", e$message, "\n")
  data.frame(
    Outcome = "Cross-Year (2019 PM2.5 → 2018 deaths)",
    Coefficient = NA, SE = NA, P_value = NA,
    Significant = NA, Interpretation = "Failed"
  )
})

if (!is.null(cross_year_result)) {
  placebo_results[["CrossYear"]] <- cross_year_result
}

# Save placebo results
placebo_table <- do.call(rbind, placebo_results)
write.csv(placebo_table, "Output/Robustness/Tables/Placebo_Outcomes.csv", row.names = FALSE)
cat("\nSaved: Output/Robustness/Tables/Placebo_Outcomes.csv\n")

# Summary
cat("\n┌─────────────────────────────────────────────────────────────────┐\n")
cat("│ PLACEBO OUTCOMES SUMMARY                                        │\n")
cat("├─────────────────────────────────────────────────────────────────┤\n")
n_pass <- sum(c(p_shuffled > 0.05, p_random > 0.05))
if (!is.null(cross_year_result) && !is.na(cross_year_result$P_value)) {
  n_pass <- n_pass + (cross_year_result$P_value > 0.05)
  n_total <- 3
} else {
  n_total <- 2
}
cat("│  Tests passed:", n_pass, "of", n_total, "                                       │\n")
if (n_pass >= 2) {
  cat("│  ✓ Supports causal interpretation                              │\n")
}
cat("└─────────────────────────────────────────────────────────────────┘\n")


################################################################################
#                                                                              #
#  PART D: DISTRIBUTED LAG MODEL (0-14 days)                                   #
#                                                                              #
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("PART D: DISTRIBUTED LAG MODEL (0-14 days)\n")
cat("================================================================\n\n")

cat("Testing individual lag effects with MATCHED instruments.\n")
cat("(Lag k pollution instrumented by lag k fire-wind)\n\n")

# Create all lag variables
max_lag <- 14

analysis_data <- analysis_data %>%
  arrange(station_id, date) %>%
  group_by(station_id)

for (k in 1:max_lag) {
  analysis_data <- analysis_data %>%
    mutate(
      !!paste0("PM25_lag", k) := dplyr::lag(PM2.5AVG, k),
      !!paste0("FRP_u_lag", k) := dplyr::lag(FRP_u_wind_station, k),
      !!paste0("FRP_v_lag", k) := dplyr::lag(FRP_v_wind_station, k)
    )
}

analysis_data <- ungroup(analysis_data)

# Estimate lag effects
lag_results <- list()

# Lag 0 (current day)
cat("Lag  0: β =", sprintf("%.5f", beta_baseline), 
    ", p =", sprintf("%.4f", 2 * pnorm(-abs(beta_baseline / se_baseline))), "\n")

lag_results[[1]] <- data.frame(
  Lag = 0,
  Coefficient = beta_baseline,
  SE = se_baseline,
  P_value = 2 * pnorm(-abs(beta_baseline / se_baseline)),
  CI_lower = beta_baseline - 1.96 * se_baseline,
  CI_upper = beta_baseline + 1.96 * se_baseline,
  Effect_10ugm3 = effect_baseline,
  N_obs = nrow(analysis_data)
)

# Lags 1 to max_lag
for (k in 1:max_lag) {
  
  lag_var <- paste0("PM25_lag", k)
  iv_u <- paste0("FRP_u_lag", k)
  iv_v <- paste0("FRP_v_lag", k)
  
  lag_results[[k + 1]] <- tryCatch({
    
    # Filter complete cases
    vars_needed <- c(lag_var, iv_u, iv_v, "deaths", 
                     "RELATIVEHUMIDITYAVG", "AmbientTemperatureAVG")
    lag_data <- analysis_data[complete.cases(analysis_data[, vars_needed]), ]
    
    if (nrow(lag_data) < 1000) {
      cat("Lag", sprintf("%2d", k), ": Insufficient data\n")
      return(data.frame(Lag = k, Coefficient = NA, SE = NA, P_value = NA,
                        CI_lower = NA, CI_upper = NA, Effect_10ugm3 = NA, 
                        N_obs = nrow(lag_data)))
    }
    
    # First stage with LAGGED instruments
    fs_formula <- as.formula(paste0(
      lag_var, " ~ ", iv_u, " + ", iv_v, " + ",
      "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
      "factor(month) + factor(year) + factor(station_id)"
    ))
    
    fs_lag <- lm(fs_formula, data = lag_data)
    lag_data$cf_resid_lag <- residuals(fs_lag)
    
    # Second stage
    ss_formula <- as.formula(paste0(
      "deaths ~ ", lag_var, " + cf_resid_lag + ",
      "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
      "factor(month) + factor(year) + factor(station_id)"
    ))
    
    ss_lag <- glm(ss_formula, family = poisson(link = "log"), data = lag_data)
    
    coef_val <- coef(ss_lag)[lag_var]
    se_val <- sqrt(vcov(ss_lag)[lag_var, lag_var])
    p_val <- 2 * pnorm(-abs(coef_val / se_val))
    
    cat("Lag", sprintf("%2d", k), ": β =", sprintf("%.5f", coef_val), 
        ", p =", sprintf("%.4f", p_val))
    if (p_val < 0.05) cat(" *")
    cat("\n")
    
    data.frame(
      Lag = k,
      Coefficient = coef_val,
      SE = se_val,
      P_value = p_val,
      CI_lower = coef_val - 1.96 * se_val,
      CI_upper = coef_val + 1.96 * se_val,
      Effect_10ugm3 = (exp(coef_val * 10) - 1) * 100,
      N_obs = nrow(lag_data)
    )
    
  }, error = function(e) {
    cat("Lag", sprintf("%2d", k), ": FAILED\n")
    data.frame(Lag = k, Coefficient = NA, SE = NA, P_value = NA,
               CI_lower = NA, CI_upper = NA, Effect_10ugm3 = NA, N_obs = NA)
  })
}

lag_table <- do.call(rbind, lag_results)
write.csv(lag_table, "Output/Robustness/Tables/Distributed_Lag_Model.csv", row.names = FALSE)
cat("\nSaved: Output/Robustness/Tables/Distributed_Lag_Model.csv\n")

# Summary
if (any(!is.na(lag_table$Coefficient))) {
  peak_idx <- which.max(lag_table$Coefficient)
  sig_lags <- sum(lag_table$P_value < 0.05, na.rm = TRUE)
  cat("\nPeak effect at Lag", lag_table$Lag[peak_idx], "days:",
      sprintf("%.1f%%", lag_table$Effect_10ugm3[peak_idx]), "per 10 µg/m³\n")
  cat("Significant lags:", sig_lags, "of", sum(!is.na(lag_table$P_value)), "\n")
}


################################################################################
#                                                                              #
#  PART E: CUMULATIVE EFFECTS (Moving Averages)                                #
#                                                                              #
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("PART E: CUMULATIVE EFFECTS (Moving Averages)\n")
cat("================================================================\n\n")

windows <- c(3, 7, 14)
cumulative_results <- list()

for (w in windows) {
  
  ma_var <- paste0("PM25_ma", w)
  
  # Create moving average
  analysis_data <- analysis_data %>%
    arrange(station_id, date) %>%
    group_by(station_id) %>%
    mutate(!!ma_var := zoo::rollmean(PM2.5AVG, k = w, fill = NA, align = "right")) %>%
    ungroup()
  
  cumulative_results[[as.character(w)]] <- tryCatch({
    
    ma_data <- analysis_data[!is.na(analysis_data[[ma_var]]), ]
    
    if (nrow(ma_data) < 1000) {
      cat("Window", sprintf("%2d", w), "days: Insufficient data\n")
      return(data.frame(Window = w, Coefficient = NA, SE = NA, P_value = NA,
                        CI_lower = NA, CI_upper = NA, Effect_10ugm3 = NA,
                        N_obs = nrow(ma_data)))
    }
    
    # First stage
    fs_formula <- as.formula(paste0(
      ma_var, " ~ FRP_u_wind_station + FRP_v_wind_station + ",
      "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
      "factor(month) + factor(year) + factor(station_id)"
    ))
    
    fs_ma <- lm(fs_formula, data = ma_data)
    ma_data$cf_resid_ma <- residuals(fs_ma)
    
    # Second stage
    ss_formula <- as.formula(paste0(
      "deaths ~ ", ma_var, " + cf_resid_ma + ",
      "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
      "factor(month) + factor(year) + factor(station_id)"
    ))
    
    ss_ma <- glm(ss_formula, family = poisson(link = "log"), data = ma_data)
    
    coef_val <- coef(ss_ma)[ma_var]
    se_val <- sqrt(vcov(ss_ma)[ma_var, ma_var])
    p_val <- 2 * pnorm(-abs(coef_val / se_val))
    effect_val <- (exp(coef_val * 10) - 1) * 100
    
    cat("Window", sprintf("%2d", w), "days: β =", sprintf("%.5f", coef_val),
        ", effect =", sprintf("%.1f%%", effect_val), "\n")
    
    data.frame(
      Window = w,
      Coefficient = coef_val,
      SE = se_val,
      P_value = p_val,
      CI_lower = coef_val - 1.96 * se_val,
      CI_upper = coef_val + 1.96 * se_val,
      Effect_10ugm3 = effect_val,
      N_obs = nrow(ma_data)
    )
    
  }, error = function(e) {
    cat("Window", sprintf("%2d", w), "days: FAILED\n")
    data.frame(Window = w, Coefficient = NA, SE = NA, P_value = NA,
               CI_lower = NA, CI_upper = NA, Effect_10ugm3 = NA, N_obs = NA)
  })
}

cumulative_table <- do.call(rbind, cumulative_results)
write.csv(cumulative_table, "Output/Robustness/Tables/Cumulative_Effects.csv", row.names = FALSE)
cat("\nSaved: Output/Robustness/Tables/Cumulative_Effects.csv\n")


################################################################################
#                                                                              #
#  PART F: HETEROGENEITY ANALYSIS (Interaction Approach)                       #
#                                                                              #
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("PART F: HETEROGENEITY ANALYSIS (Interaction Approach)\n")
cat("================================================================\n\n")

cat("Note: Using interaction terms for stable heterogeneity estimates.\n")
cat("      Separate subgroup regressions can be unstable with limited variation.\n\n")

het_results <- list()

# -----------------------------------------------------------------------------
# F.1 Season Interaction (Haze vs Non-Haze)
# -----------------------------------------------------------------------------
cat("--- F.1 Season Interaction (Haze vs Non-Haze) ---\n")

analysis_data$haze_season <- ifelse(
  analysis_data$month %in% c(6, 7, 8, 9), 1, 0
)

het_season <- tryCatch({
  
  # First stage with season-specific instruments
  fs_het <- lm(
    PM2.5AVG ~ FRP_u_wind_station + FRP_v_wind_station +
      I(FRP_u_wind_station * haze_season) + I(FRP_v_wind_station * haze_season) +
      RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
      factor(month) + factor(year) + factor(station_id),
    data = analysis_data
  )
  analysis_data$cf_resid_het <- residuals(fs_het)
  
  # Second stage with interaction
  ss_het <- glm(
    deaths ~ PM2.5AVG * haze_season + cf_resid_het +
      RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
      factor(month) + factor(year) + factor(station_id),
    family = poisson(link = "log"),
    data = analysis_data
  )
  
  coef_base <- coef(ss_het)["PM2.5AVG"]
  coef_interaction <- coef(ss_het)["PM2.5AVG:haze_season"]
  se_interaction <- sqrt(vcov(ss_het)["PM2.5AVG:haze_season", "PM2.5AVG:haze_season"])
  p_interaction <- 2 * pnorm(-abs(coef_interaction / se_interaction))
  
  effect_nonhaze <- coef_base
  effect_haze <- coef_base + coef_interaction
  
  cat("  Non-haze (Oct-May): β =", sprintf("%.5f", effect_nonhaze),
      ", effect =", sprintf("%.1f%%", (exp(effect_nonhaze * 10) - 1) * 100), "\n")
  cat("  Haze (Jun-Sep):     β =", sprintf("%.5f", effect_haze),
      ", effect =", sprintf("%.1f%%", (exp(effect_haze * 10) - 1) * 100), "\n")
  cat("  Interaction p-value:", sprintf("%.4f", p_interaction), "\n")
  
  data.frame(
    Analysis = "Season",
    Subgroup1 = "Non-Haze (Oct-May)",
    Subgroup2 = "Haze (Jun-Sep)",
    Effect1_Pct = (exp(effect_nonhaze * 10) - 1) * 100,
    Effect2_Pct = (exp(effect_haze * 10) - 1) * 100,
    Interaction_Coef = coef_interaction,
    Interaction_SE = se_interaction,
    Interaction_P = p_interaction
  )
  
}, error = function(e) {
  cat("  FAILED:", e$message, "\n")
  NULL
})

if (!is.null(het_season)) {
  het_results[["Season"]] <- het_season
}

# -----------------------------------------------------------------------------
# F.2 Year Interactions (Using 2019 as Reference)
# -----------------------------------------------------------------------------
cat("\n--- F.2 Year Effects (2019 as Reference) ---\n")

het_year <- tryCatch({
  
  analysis_data$year_factor <- factor(analysis_data$year)
  analysis_data$year_factor <- relevel(analysis_data$year_factor, ref = "2019")
  
  fs_het2 <- lm(
    PM2.5AVG ~ FRP_u_wind_station + FRP_v_wind_station +
      RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
      factor(month) + year_factor + factor(station_id),
    data = analysis_data
  )
  analysis_data$cf_resid_het2 <- residuals(fs_het2)
  
  ss_het2 <- glm(
    deaths ~ PM2.5AVG * year_factor + cf_resid_het2 +
      RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
      factor(month) + factor(station_id),
    family = poisson(link = "log"),
    data = analysis_data
  )
  
  coef_2019 <- coef(ss_het2)["PM2.5AVG"]
  se_2019 <- sqrt(vcov(ss_het2)["PM2.5AVG", "PM2.5AVG"])
  
  # Variance-covariance matrix, needed for the linear combinations below
  V <- vcov(ss_het2)
  
  year_effects <- NULL
  
  # 2019 (reference)
  year_effects <- rbind(year_effects, data.frame(
    Year             = "2019 (ref)",
    Coefficient      = as.numeric(coef_2019),
    SE               = as.numeric(se_2019),
    Effect_Pct       = (exp(coef_2019 * 10) - 1) * 100,
    CI_Lower         = as.numeric(coef_2019 - 1.96 * se_2019),
    CI_Upper         = as.numeric(coef_2019 + 1.96 * se_2019),
    Interaction_Coef = 0,
    Interaction_SE   = NA_real_,
    Interaction_P    = NA_real_,
    stringsAsFactors = FALSE
  ))
  
  cat("  2019 (reference): beta =", sprintf("%.5f", coef_2019),
      "(SE", sprintf("%.5f", se_2019), "), effect =",
      sprintf("%.1f%%", (exp(coef_2019 * 10) - 1) * 100), "\n")
  
  # Other years: effect_y = coef_2019 + coef_int, a linear combination.
  # Var(a + b) = Var(a) + Var(b) + 2*Cov(a, b). The covariance term is not
  # optional: the level and interaction coefficients are estimated from
  # overlapping data and are strongly negatively correlated, so
  # sqrt(Var(a) + Var(b)) would overstate the SE.
  for (y in c("2017", "2018", "2020")) {
    
    int_name <- paste0("PM2.5AVG:year_factor", y)
    
    if (int_name %in% names(coef(ss_het2))) {
      
      coef_int <- as.numeric(coef(ss_het2)[int_name])
      effect_y <- as.numeric(coef_2019) + coef_int
      
      se_y  <- sqrt(V["PM2.5AVG", "PM2.5AVG"] +
                      V[int_name, int_name] +
                      2 * V["PM2.5AVG", int_name])
      
      se_int <- sqrt(V[int_name, int_name])
      p_int  <- 2 * pnorm(-abs(coef_int / se_int))
      
      year_effects <- rbind(year_effects, data.frame(
        Year             = y,
        Coefficient      = effect_y,
        SE               = se_y,
        Effect_Pct       = (exp(effect_y * 10) - 1) * 100,
        CI_Lower         = effect_y - 1.96 * se_y,
        CI_Upper         = effect_y + 1.96 * se_y,
        Interaction_Coef = coef_int,
        Interaction_SE   = se_int,
        Interaction_P    = p_int,
        stringsAsFactors = FALSE
      ))
      
      cat("  ", y, ": beta =", sprintf("%.5f", effect_y),
          "(SE", sprintf("%.5f", se_y), "), effect =",
          sprintf("%.1f%%", (exp(effect_y * 10) - 1) * 100),
          "| vs 2019 p =", sprintf("%.3f", p_int), "\n")
    }
  }
  year_effects
  
}, error = function(e) {
  cat("  FAILED:", e$message, "\n")
  NULL
})

if (!is.null(het_year)) {
  het_results[["Year"]] <- het_year
}

# Save heterogeneity results
if (!is.null(het_season)) {
  write.csv(het_season, "Output/Robustness/Tables/Heterogeneity_Season.csv", row.names = FALSE)
  cat("\nSaved: Output/Robustness/Tables/Heterogeneity_Season.csv\n")
}

if (!is.null(het_year)) {
  write.csv(het_year, "Output/Robustness/Tables/Heterogeneity_Year.csv", row.names = FALSE)
  cat("Saved: Output/Robustness/Tables/Heterogeneity_Year.csv\n")
}


################################################################################
#                                                                              #
#  PART G: SENSITIVITY TO OUTLIERS                                             #
#                                                                              #
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("PART G: SENSITIVITY TO OUTLIERS\n")
cat("================================================================\n\n")

sensitivity_results <- list()

# G.1 Original (Reference)
sensitivity_results[["Original"]] <- data.frame(
  Specification = "Original",
  Coefficient = beta_baseline,
  SE = se_baseline,
  Effect_10ugm3 = effect_baseline,
  N_obs = nrow(analysis_data),
  Pct_change = 0
)
cat("G.1 Original: β =", sprintf("%.6f", beta_baseline), "\n")

# G.2 Winsorized 1%/99%
cat("\nG.2 Winsorized 1%/99%:\n")

p01 <- quantile(analysis_data$PM2.5AVG, 0.01, na.rm = TRUE)
p99 <- quantile(analysis_data$PM2.5AVG, 0.99, na.rm = TRUE)

wins_data <- analysis_data %>%
  mutate(PM2.5AVG_wins = pmin(pmax(PM2.5AVG, p01), p99))

sensitivity_results[["Winsorized"]] <- tryCatch({
  
  fs <- lm(PM2.5AVG_wins ~ FRP_u_wind_station + FRP_v_wind_station +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
             factor(month) + factor(year) + factor(station_id),
           data = wins_data)
  wins_data$cf_resid <- residuals(fs)
  
  ss <- glm(deaths ~ PM2.5AVG_wins + cf_resid +
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
              factor(month) + factor(year) + factor(station_id),
            family = poisson, data = wins_data)
  
  coef_val <- coef(ss)["PM2.5AVG_wins"]
  se_val <- sqrt(vcov(ss)["PM2.5AVG_wins", "PM2.5AVG_wins"])
  
  cat("    β =", sprintf("%.6f", coef_val), 
      "(", sprintf("%.1f%%", (coef_val/beta_baseline - 1) * 100), "change)\n")
  
  data.frame(
    Specification = "Winsorized 1%/99%",
    Coefficient = coef_val,
    SE = se_val,
    Effect_10ugm3 = (exp(coef_val * 10) - 1) * 100,
    N_obs = nrow(wins_data),
    Pct_change = (coef_val / beta_baseline - 1) * 100
  )
}, error = function(e) {
  cat("    FAILED\n")
  data.frame(Specification = "Winsorized 1%/99%", Coefficient = NA, SE = NA,
             Effect_10ugm3 = NA, N_obs = NA, Pct_change = NA)
})

# G.3 Exclude PM2.5 > 100
cat("\nG.3 Exclude PM2.5 > 100:\n")

excl_data <- analysis_data %>% filter(PM2.5AVG <= 100)

sensitivity_results[["Exclude_100"]] <- tryCatch({
  
  fs <- lm(PM2.5AVG ~ FRP_u_wind_station + FRP_v_wind_station +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
             factor(month) + factor(year) + factor(station_id),
           data = excl_data)
  excl_data$cf_resid <- residuals(fs)
  
  ss <- glm(deaths ~ PM2.5AVG + cf_resid +
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
              factor(month) + factor(year) + factor(station_id),
            family = poisson, data = excl_data)
  
  coef_val <- coef(ss)["PM2.5AVG"]
  se_val <- sqrt(vcov(ss)["PM2.5AVG", "PM2.5AVG"])
  
  cat("    β =", sprintf("%.6f", coef_val),
      ", n =", format(nrow(excl_data), big.mark = ","),
      "(", sprintf("%.1f%%", (coef_val/beta_baseline - 1) * 100), "change)\n")
  
  data.frame(
    Specification = "Exclude PM2.5 > 100",
    Coefficient = coef_val,
    SE = se_val,
    Effect_10ugm3 = (exp(coef_val * 10) - 1) * 100,
    N_obs = nrow(excl_data),
    Pct_change = (coef_val / beta_baseline - 1) * 100
  )
}, error = function(e) {
  cat("    FAILED\n")
  data.frame(Specification = "Exclude PM2.5 > 100", Coefficient = NA, SE = NA,
             Effect_10ugm3 = NA, N_obs = NA, Pct_change = NA)
})

# G.4 Exclude PM2.5 > 75
cat("\nG.4 Exclude PM2.5 > 75:\n")

excl75_data <- analysis_data %>% filter(PM2.5AVG <= 75)

sensitivity_results[["Exclude_75"]] <- tryCatch({
  
  fs <- lm(PM2.5AVG ~ FRP_u_wind_station + FRP_v_wind_station +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
             factor(month) + factor(year) + factor(station_id),
           data = excl75_data)
  excl75_data$cf_resid <- residuals(fs)
  
  ss <- glm(deaths ~ PM2.5AVG + cf_resid +
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
              factor(month) + factor(year) + factor(station_id),
            family = poisson, data = excl75_data)
  
  coef_val <- coef(ss)["PM2.5AVG"]
  se_val <- sqrt(vcov(ss)["PM2.5AVG", "PM2.5AVG"])
  
  cat("    β =", sprintf("%.6f", coef_val),
      ", n =", format(nrow(excl75_data), big.mark = ","),
      "(", sprintf("%.1f%%", (coef_val/beta_baseline - 1) * 100), "change)\n")
  
  data.frame(
    Specification = "Exclude PM2.5 > 75",
    Coefficient = coef_val,
    SE = se_val,
    Effect_10ugm3 = (exp(coef_val * 10) - 1) * 100,
    N_obs = nrow(excl75_data),
    Pct_change = (coef_val / beta_baseline - 1) * 100
  )
}, error = function(e) {
  cat("    FAILED\n")
  data.frame(Specification = "Exclude PM2.5 > 75", Coefficient = NA, SE = NA,
             Effect_10ugm3 = NA, N_obs = NA, Pct_change = NA)
})

sensitivity_table <- do.call(rbind, sensitivity_results)
write.csv(sensitivity_table, "Output/Robustness/Tables/Sensitivity_Outliers.csv", row.names = FALSE)
cat("\nSaved: Output/Robustness/Tables/Sensitivity_Outliers.csv\n")


################################################################################
#                                                                              #
#  PART H: ALTERNATIVE IV SPECIFICATIONS                                       #
#                                                                              #
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("PART H: ALTERNATIVE IV SPECIFICATIONS\n")
cat("================================================================\n\n")

alt_iv_results <- list()

# Calculate restricted model RSS for F-stat computation
fs_restricted <- lm(PM2.5AVG ~ RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
                      factor(month) + factor(year) + factor(station_id),
                    data = analysis_data)
rss_restricted <- sum(resid(fs_restricted)^2)
n <- nrow(analysis_data)

# H.1 Both instruments (baseline)
cat("H.1 Both instruments (main specification):\n")

fs_both <- lm(PM2.5AVG ~ FRP_u_wind_station + FRP_v_wind_station +
                RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
                factor(month) + factor(year) + factor(station_id),
              data = analysis_data)

rss_full <- sum(resid(fs_both)^2)
k <- length(coef(fs_both))
f_both <- ((rss_restricted - rss_full) / 2) / (rss_full / (n - k))

cat("    β =", sprintf("%.6f", beta_baseline), ", F =", round(f_both, 1), "\n")

alt_iv_results[["Both"]] <- data.frame(
  Specification = "Both instruments (main)",
  Coefficient = beta_baseline,
  SE = se_baseline,
  First_stage_F = f_both,
  Effect_10ugm3 = effect_baseline,
  N_obs = nrow(analysis_data)
)

# H.2 Only FRP × u_wind
cat("\nH.2 Only FRP × u_wind:\n")

alt_iv_results[["U_only"]] <- tryCatch({
  
  fs <- lm(PM2.5AVG ~ FRP_u_wind_station +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
             factor(month) + factor(year) + factor(station_id),
           data = analysis_data)
  
  rss_full_u <- sum(resid(fs)^2)
  k_u <- length(coef(fs))
  f_u <- ((rss_restricted - rss_full_u) / 1) / (rss_full_u / (n - k_u))
  
  analysis_data$cf_resid_u <- residuals(fs)
  
  ss <- glm(deaths ~ PM2.5AVG + cf_resid_u +
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
              factor(month) + factor(year) + factor(station_id),
            family = poisson, data = analysis_data)
  
  coef_val <- coef(ss)["PM2.5AVG"]
  se_val <- sqrt(vcov(ss)["PM2.5AVG", "PM2.5AVG"])
  
  cat("    β =", sprintf("%.6f", coef_val), ", F =", round(f_u, 1), "\n")
  
  data.frame(
    Specification = "Only FRP × u_wind",
    Coefficient = coef_val,
    SE = se_val,
    First_stage_F = f_u,
    Effect_10ugm3 = (exp(coef_val * 10) - 1) * 100,
    N_obs = nrow(analysis_data)
  )
}, error = function(e) {
  cat("    FAILED\n")
  data.frame(Specification = "Only FRP × u_wind", Coefficient = NA, SE = NA,
             First_stage_F = NA, Effect_10ugm3 = NA, N_obs = NA)
})

# H.3 Only FRP × v_wind
cat("\nH.3 Only FRP × v_wind:\n")

alt_iv_results[["V_only"]] <- tryCatch({
  
  fs <- lm(PM2.5AVG ~ FRP_v_wind_station +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
             factor(month) + factor(year) + factor(station_id),
           data = analysis_data)
  
  rss_full_v <- sum(resid(fs)^2)
  k_v <- length(coef(fs))
  f_v <- ((rss_restricted - rss_full_v) / 1) / (rss_full_v / (n - k_v))
  
  analysis_data$cf_resid_v <- residuals(fs)
  
  ss <- glm(deaths ~ PM2.5AVG + cf_resid_v +
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
              factor(month) + factor(year) + factor(station_id),
            family = poisson, data = analysis_data)
  
  coef_val <- coef(ss)["PM2.5AVG"]
  se_val <- sqrt(vcov(ss)["PM2.5AVG", "PM2.5AVG"])
  
  cat("    β =", sprintf("%.6f", coef_val), ", F =", round(f_v, 1), "\n")
  
  data.frame(
    Specification = "Only FRP × v_wind",
    Coefficient = coef_val,
    SE = se_val,
    First_stage_F = f_v,
    Effect_10ugm3 = (exp(coef_val * 10) - 1) * 100,
    N_obs = nrow(analysis_data)
  )
}, error = function(e) {
  cat("    FAILED\n")
  data.frame(Specification = "Only FRP × v_wind", Coefficient = NA, SE = NA,
             First_stage_F = NA, Effect_10ugm3 = NA, N_obs = NA)
})

# H.4 FRP only (no wind interaction)
cat("\nH.4 FRP only (no wind interaction):\n")

alt_iv_results[["FRP_only"]] <- tryCatch({
  
  fs <- lm(PM2.5AVG ~ total_frp +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
             factor(month) + factor(year) + factor(station_id),
           data = analysis_data)
  
  rss_full_frp <- sum(resid(fs)^2)
  k_frp <- length(coef(fs))
  f_frp <- ((rss_restricted - rss_full_frp) / 1) / (rss_full_frp / (n - k_frp))
  
  analysis_data$cf_resid_frp <- residuals(fs)
  
  ss <- glm(deaths ~ PM2.5AVG + cf_resid_frp +
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
              factor(month) + factor(year) + factor(station_id),
            family = poisson, data = analysis_data)
  
  coef_val <- coef(ss)["PM2.5AVG"]
  se_val <- sqrt(vcov(ss)["PM2.5AVG", "PM2.5AVG"])
  
  cat("    β =", sprintf("%.6f", coef_val), ", F =", round(f_frp, 1), "\n")
  
  data.frame(
    Specification = "FRP only (no wind)",
    Coefficient = coef_val,
    SE = se_val,
    First_stage_F = f_frp,
    Effect_10ugm3 = (exp(coef_val * 10) - 1) * 100,
    N_obs = nrow(analysis_data)
  )
}, error = function(e) {
  cat("    FAILED\n")
  data.frame(Specification = "FRP only (no wind)", Coefficient = NA, SE = NA,
             First_stage_F = NA, Effect_10ugm3 = NA, N_obs = NA)
})

alt_iv_table <- do.call(rbind, alt_iv_results)
write.csv(alt_iv_table, "Output/Robustness/Tables/Alternative_IV_Specifications.csv", row.names = FALSE)
cat("\nSaved: Output/Robustness/Tables/Alternative_IV_Specifications.csv\n")


################################################################################
#                                                                              #
#  PART I: FIGURES                                                             #
#                                                                              #
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("PART I: FIGURES\n")
cat("================================================================\n\n")

# -----------------------------------------------------------------------------
# Figure 1: Distributed Lag Structure
# -----------------------------------------------------------------------------
cat("Creating Figure 1: Distributed Lag Structure...\n")

lag_plot_data <- lag_table %>%
  filter(!is.na(Coefficient))

fig_lag <- ggplot(lag_plot_data, aes(x = Lag, y = Coefficient)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  scale_x_continuous(breaks = 0:14) +
  labs(
    title = "Distributed Lag Structure: PM2.5 Effect on Mortality",
    subtitle = "Each lag estimated separately with matched instruments",
    x = "Lag (days)",
    y = "IV-CF Coefficient",
    caption = "Shaded region = 95% CI. Each lag uses corresponding lagged instruments."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("Output/Robustness/Figures/Fig_Distributed_Lag.png", fig_lag,
       width = 10, height = 6, dpi = 300)
ggsave("Output/Robustness/Figures/Fig_Distributed_Lag.pdf", fig_lag,
       width = 10, height = 6)
cat("  Saved: Fig_Distributed_Lag.png\n")

# -----------------------------------------------------------------------------
# Figure 2: Coefficient Stability (Alternative IVs and Sensitivity)
# -----------------------------------------------------------------------------
cat("\nCreating Figure 2: Coefficient Stability...\n")

stability_data <- bind_rows(
  sensitivity_table %>%
    filter(!is.na(Coefficient)) %>%
    mutate(Type = "Sensitivity") %>%
    dplyr::select(Specification, Coefficient, SE, Type),
  
  alt_iv_table %>%
    filter(!is.na(Coefficient)) %>%
    mutate(Type = "Alternative IV") %>%
    dplyr::select(Specification, Coefficient, SE, Type)
) %>%
  mutate(
    CI_lower = Coefficient - 1.96 * SE,
    CI_upper = Coefficient + 1.96 * SE,
    Specification = factor(Specification, levels = rev(unique(Specification)))
  )

fig_stability <- ggplot(stability_data, aes(x = Coefficient, y = Specification, color = Type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = beta_baseline, linetype = "dotted", color = "red", alpha = 0.7) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Sensitivity" = "steelblue", "Alternative IV" = "darkorange")) +
  labs(
    title = "Coefficient Stability Across Specifications",
    subtitle = "Red dotted line = main estimate",
    x = "Coefficient on PM2.5",
    y = "",
    color = "Type",
    caption = "Error bars = 95% CI"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("Output/Robustness/Figures/Fig_Coefficient_Stability.png", fig_stability,
       width = 10, height = 7, dpi = 300)
ggsave("Output/Robustness/Figures/Fig_Coefficient_Stability.pdf", fig_stability,
       width = 10, height = 7)
cat("  Saved: Fig_Coefficient_Stability.png\n")

# -----------------------------------------------------------------------------
# Figure 3: Heterogeneity by Season
# -----------------------------------------------------------------------------
cat("\nCreating Figure 3: Heterogeneity by Season...\n")

if (!is.null(het_season)) {
  
  season_plot_data <- data.frame(
    Season = c("Haze\n(Jun-Sep)", "Non-Haze\n(Oct-May)"),
    Effect = c(het_season$Effect2_Pct, het_season$Effect1_Pct),
    Type = c("Haze", "Non-Haze")
  )
  
  fig_season <- ggplot(season_plot_data, aes(x = Season, y = Effect, fill = Type)) +
    geom_hline(yintercept = effect_baseline, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50") +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", Effect)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Haze" = "darkorange", "Non-Haze" = "steelblue")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      title = "Effect Heterogeneity by Season",
      subtitle = paste0("Red dashed line = pooled estimate (", sprintf("%.1f%%", effect_baseline), ")"),
      x = "",
      y = "Effect per 10 µg/m³ (%)",
      caption = paste0("Season interaction p-value: ", sprintf("%.3f", het_season$Interaction_P))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ggsave("Output/Robustness/Figures/Fig_Heterogeneity_Season.png", fig_season,
         width = 7, height = 6, dpi = 300)
  ggsave("Output/Robustness/Figures/Fig_Heterogeneity_Season.pdf", fig_season,
         width = 7, height = 6)
  cat("  Saved: Fig_Heterogeneity_Season.png\n")
}

# -----------------------------------------------------------------------------
# Figure 4: Heterogeneity by Year
# -----------------------------------------------------------------------------
cat("\nCreating Figure 4: Heterogeneity by Year...\n")

if (!is.null(het_year) && nrow(het_year) > 0) {
  
  year_plot_data <- het_year %>%
    filter(!is.na(Coefficient)) %>%
    mutate(Year = gsub(" \\(ref\\)", "", Year))
  
  fig_year <- ggplot(year_plot_data, aes(x = Year, y = Effect_Pct)) +
    geom_hline(yintercept = effect_baseline, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50") +
    geom_col(fill = "steelblue", alpha = 0.8, width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", Effect_Pct)), vjust = -0.5, size = 4) +
    scale_y_continuous(limits = c(-10, max(year_plot_data$Effect_Pct) * 1.2)) +
    labs(
      title = "Effect Heterogeneity by Year",
      subtitle = paste0("Red dashed line = pooled estimate (", sprintf("%.1f%%", effect_baseline), 
                        "); 2019 is reference year"),
      x = "Year",
      y = "Effect per 10 µg/m³ (%)",
      caption = "Note: Year effects estimated via interaction terms in pooled model."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9, color = "gray50")
    )
  
  ggsave("Output/Robustness/Figures/Fig_Heterogeneity_Year.png", fig_year,
         width = 8, height = 6, dpi = 300)
  ggsave("Output/Robustness/Figures/Fig_Heterogeneity_Year.pdf", fig_year,
         width = 8, height = 6)
  cat("  Saved: Fig_Heterogeneity_Year.png\n")
}


################################################################################
#                                                                              #
#  PART J: SUMMARY                                                             #
#                                                                              #
################################################################################

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                         ROBUSTNESS SUMMARY                               ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

cat("BASELINE RESULT:\n")
cat("  β =", sprintf("%.6f", beta_baseline), "\n")
cat("  Effect per 10 µg/m³:", sprintf("%.1f%%", effect_baseline), "\n\n")

cat("PLACEBO OUTCOMES:\n")
cat("  Shuffled deaths: p =", sprintf("%.4f", p_shuffled), 
    ifelse(p_shuffled > 0.05, " ✓ PASS", " ⚠️"), "\n")
cat("  Random deaths: p =", sprintf("%.4f", p_random),
    ifelse(p_random > 0.05, " ✓ PASS", " ⚠️"), "\n")
if (!is.null(cross_year_result) && !is.na(cross_year_result$P_value)) {
  cat("  Cross-year: p =", sprintf("%.4f", cross_year_result$P_value),
      ifelse(cross_year_result$P_value > 0.05, " ✓ PASS", " ⚠️"), "\n")
}
cat("\n")

cat("DISTRIBUTED LAG:\n")
if (any(!is.na(lag_table$Coefficient))) {
  peak_idx <- which.max(lag_table$Coefficient)
  sig_lags <- sum(lag_table$P_value < 0.05, na.rm = TRUE)
  cat("  Peak effect: Lag", lag_table$Lag[peak_idx], "days (",
      sprintf("%.1f%%", lag_table$Effect_10ugm3[peak_idx]), ")\n")
  cat("  Significant lags:", sig_lags, "of", sum(!is.na(lag_table$P_value)), "\n\n")
}

cat("CUMULATIVE EFFECTS:\n")
for (i in 1:nrow(cumulative_table)) {
  if (!is.na(cumulative_table$Coefficient[i])) {
    cat("  ", cumulative_table$Window[i], "-day MA:",
        sprintf("%.1f%%", cumulative_table$Effect_10ugm3[i]), "\n")
  }
}
cat("\n")

cat("HETEROGENEITY:\n")
if (!is.null(het_season)) {
  cat("  Non-haze season:", sprintf("%.1f%%", het_season$Effect1_Pct), "\n")
  cat("  Haze season:", sprintf("%.1f%%", het_season$Effect2_Pct), "\n")
  cat("  Interaction p-value:", sprintf("%.3f", het_season$Interaction_P), "\n\n")
}

cat("SENSITIVITY:\n")
if (nrow(sensitivity_table) > 0) {
  coef_range <- range(sensitivity_table$Coefficient, na.rm = TRUE)
  cat("  Coefficient range:", sprintf("%.4f", coef_range[1]), "to", 
      sprintf("%.4f", coef_range[2]), "\n\n")
}

cat("ALTERNATIVE IV:\n")
if (nrow(alt_iv_table) > 0) {
  coef_range_iv <- range(alt_iv_table$Coefficient, na.rm = TRUE)
  f_range <- range(alt_iv_table$First_stage_F, na.rm = TRUE)
  cat("  Coefficient range:", sprintf("%.4f", coef_range_iv[1]), "to",
      sprintf("%.4f", coef_range_iv[2]), "\n")
  cat("  F-stat range:", round(f_range[1], 0), "to", round(f_range[2], 0), "\n")
}

cat("\n")
cat("════════════════════════════════════════════════════════════════════════════\n")
cat("OUTPUT FILES CREATED:\n")
cat("════════════════════════════════════════════════════════════════════════════\n\n")

cat("Tables (Output/Robustness/Tables/):\n")
cat("  • Placebo_Outcomes.csv\n")
cat("  • Distributed_Lag_Model.csv\n")
cat("  • Cumulative_Effects.csv\n")
cat("  • Heterogeneity_Season.csv\n")
cat("  • Heterogeneity_Year.csv\n")
cat("  • Sensitivity_Outliers.csv\n")
cat("  • Alternative_IV_Specifications.csv\n")

cat("\nFigures (Output/Robustness/Figures/):\n")
cat("  • Fig_Distributed_Lag.png/pdf\n")
cat("  • Fig_Coefficient_Stability.png/pdf\n")
cat("  • Fig_Heterogeneity_Season.png/pdf\n")
cat("  • Fig_Heterogeneity_Year.png/pdf\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                      ROBUSTNESS ANALYSIS COMPLETE                        ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n")