print("Loading packages...")
packages <- c("ggplot2", "dplyr", "tidyr", "AER", "viridis", "gridExtra", 
              "ggridges", "cowplot", "boot", "sandwich", "lmtest", 
              "knitr", "kableExtra", "broom", "purrr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
print("packages loaded successfully")

print("Loading data...")
merged_daily <- read.csv('Data/Mid_process_data/merged_death_data_from_script_01.csv', header = T)
print("Data loaded successfully!")

merged_daily <- merged_daily %>%
  mutate(
    NOAVG_clean = ifelse(NOAVG < 0, NA, NOAVG),
    NOXAVG_clean = ifelse(NOXAVG < 0, NA, NOXAVG)
  )
merged_daily <- merged_daily %>%
  mutate(
    SO2_ugm3 = (SO2AVG * 64.066 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)),
    NO2_ugm3 = (NO2AVG * 46.0055 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)),
    NO_ugm3 = (NOAVG_clean * 30.01 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)),
    NOX_ugm3 = (NOXAVG_clean * 46.0055 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)), # using NO2 MW as reference
    O3_ugm3 = (O3AVG * 48 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)),
    CO_ugm3 = (COAVG * 28.01 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG))
  )

merged_daily_pm10 <- merged_daily %>%
  drop_na(PM10AVG, total_frp, FRP_u_wind, FRP_v_wind, RELATIVEHUMIDITYAVG, AmbientTemperatureAVG, year, month)

merged_daily_pm2.5 <- merged_daily %>%
  drop_na(PM2.5AVG, total_frp, FRP_u_wind, FRP_v_wind, RELATIVEHUMIDITYAVG, AmbientTemperatureAVG, year, month)

merged_daily_o3 <- merged_daily %>%
  drop_na(O3_ugm3, total_frp, FRP_u_wind, FRP_v_wind, RELATIVEHUMIDITYAVG, AmbientTemperatureAVG, year, month)

merged_daily_co <- merged_daily %>%
  drop_na(CO_ugm3, total_frp, FRP_u_wind, FRP_v_wind, RELATIVEHUMIDITYAVG, AmbientTemperatureAVG, year, month)

merged_daily_no <- merged_daily %>%
  drop_na(NO_ugm3, total_frp, FRP_u_wind, FRP_v_wind, RELATIVEHUMIDITYAVG, AmbientTemperatureAVG, year, month)

merged_daily_no2 <- merged_daily %>%
  drop_na(NO2_ugm3, total_frp, FRP_u_wind, FRP_v_wind, RELATIVEHUMIDITYAVG, AmbientTemperatureAVG, year, month)

merged_daily_nox <- merged_daily %>%
  drop_na(NOX_ugm3, total_frp, FRP_u_wind, FRP_v_wind, RELATIVEHUMIDITYAVG, AmbientTemperatureAVG, year, month)

merged_daily_so2 <- merged_daily %>%
  drop_na(SO2_ugm3, total_frp, FRP_u_wind, FRP_v_wind, RELATIVEHUMIDITYAVG, AmbientTemperatureAVG, year, month)

pollutants <- c("PM10AVG", 
                "PM2.5AVG", 
                "O3_ugm3", 
                "CO_ugm3", 
                "NO_ugm3", 
                "NO2_ugm3", 
                "NOX_ugm3", 
                "SO2_ugm3")

#------------------------------------------------------------------------
# First-Stage Diagnostics Test
#------------------------------------------------------------------------

# Function to calculate first-stage diagnostics
calculate_first_stage_diagnostics <- function(data, pollutant) {
  # Filter data
  filtered_data <- data %>%
    filter(!is.na(!!sym(pollutant)),
           !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
           !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))
  
  # First stage formula
  first_stage_formula <- as.formula(paste0(
    pollutant, " ~ total_frp + FRP_u_wind + FRP_v_wind +",
    " RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + factor(month) + factor(year)"
  ))
  
  # Run first stage regression
  first_stage_model <- tryCatch({
    lm(first_stage_formula, data = filtered_data)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(first_stage_model)) {
    return(data.frame(
      Pollutant = pollutant,
      F_statistic = NA,
      F_p_value = NA,
      R_squared = NA,
      Partial_R2 = NA,
      N = NA,
      total_frp_coef = NA,
      total_frp_p = NA,
      FRP_u_wind_coef = NA,
      FRP_u_wind_p = NA,
      FRP_v_wind_coef = NA,
      FRP_v_wind_p = NA
    ))
  }
  
  # Get coefficients and p-values for instruments
  model_summary <- summary(first_stage_model)
  coefs <- coef(model_summary)
  
  # Calculate F-statistic for instruments only - MODIFIED SECTION
  # Define the restricted model (without instruments)
  restricted_formula <- as.formula(paste0(
    pollutant, " ~ RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + factor(month) + factor(year)"
  ))
  
  # Run model without instruments
  restricted_model <- tryCatch({
    lm(restricted_formula, data = filtered_data)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(restricted_model)) {
    f_stat <- NA
    f_p <- NA
    partial_r2 <- NA
  } else {
    # Manually calculate F-statistic
    n <- nobs(first_stage_model)
    k <- length(coef(first_stage_model))
    p <- 3  # Number of instruments: total_frp, FRP_u_wind, FRP_v_wind
    rss_full <- sum(resid(first_stage_model)^2)
    rss_restricted <- sum(resid(restricted_model)^2)
    f_stat <- ((rss_restricted - rss_full)/p) / (rss_full/(n-k))
    f_p <- pf(f_stat, p, n-k, lower.tail = FALSE)
    
    # Calculate partial R-squared for instruments
    full_r2 <- summary(first_stage_model)$r.squared
    restricted_r2 <- summary(restricted_model)$r.squared
    partial_r2 <- (full_r2 - restricted_r2) / (1 - restricted_r2)
  }
  
  # Extract coefficients and p-values for instruments
  if ("total_frp" %in% rownames(coefs)) {
    total_frp_coef <- coefs["total_frp", "Estimate"]
    total_frp_p <- coefs["total_frp", "Pr(>|t|)"]
  } else {
    total_frp_coef <- NA
    total_frp_p <- NA
  }
  
  if ("FRP_u_wind" %in% rownames(coefs)) {
    FRP_u_wind_coef <- coefs["FRP_u_wind", "Estimate"]
    FRP_u_wind_p <- coefs["FRP_u_wind", "Pr(>|t|)"]
  } else {
    FRP_u_wind_coef <- NA
    FRP_u_wind_p <- NA
  }
  
  if ("FRP_v_wind" %in% rownames(coefs)) {
    FRP_v_wind_coef <- coefs["FRP_v_wind", "Estimate"]
    FRP_v_wind_p <- coefs["FRP_v_wind", "Pr(>|t|)"]
  } else {
    FRP_v_wind_coef <- NA
    FRP_v_wind_p <- NA
  }
  
  # Return diagnostics
  return(data.frame(
    Pollutant = pollutant,
    F_statistic = f_stat,
    F_p_value = f_p,
    R_squared = full_r2,
    Partial_R2 = partial_r2,
    N = nobs(first_stage_model),
    total_frp_coef = total_frp_coef,
    total_frp_p = total_frp_p,
    FRP_u_wind_coef = FRP_u_wind_coef,
    FRP_u_wind_p = FRP_u_wind_p,
    FRP_v_wind_coef = FRP_v_wind_coef,
    FRP_v_wind_p = FRP_v_wind_p
  ))
}

# Calculate first-stage diagnostics for all pollutants
first_stage_diagnostics_all <- lapply(pollutants, function(poll) {
  calculate_first_stage_diagnostics(merged_daily, poll)
})

# Combine results
first_stage_diagnostics_df <- bind_rows(first_stage_diagnostics_all)

# Save diagnostics
write.csv(first_stage_diagnostics_df, 
          "Output/FS_diagnostics.csv", 
          row.names = FALSE)

# Create a formatted table of first stage results
first_stage_table <- first_stage_diagnostics_df %>%
  select(Pollutant, F_statistic, F_p_value, Partial_R2, N) %>%
  mutate(
    F_statistic = round(F_statistic, 2),
    F_p_value = ifelse(F_p_value < 0.001, "<0.001", 
                       ifelse(F_p_value < 0.01, "<0.01",
                              ifelse(F_p_value < 0.05, "<0.05",
                                     round(F_p_value, 3)))),
    Partial_R2 = round(Partial_R2, 3),
    Weak_Instrument = ifelse(F_statistic < 10, "Yes", "No")
  )

# Print the table
print(kable(first_stage_table, 
            caption = "First-Stage Diagnostics for All Pollutants",
            format = "markdown"))

#------------------------------------------------------------------------
# AUTOCORRELATION TEST AND CUMBY-HUIZINGA IMPLEMENTATION
#------------------------------------------------------------------------

# Function that implements the 2-stage residual inclusion model
run_2sri_model <- function(data, pollutant, outcome = "deaths") {
  
  # Ensure we're using dplyr's filter
  model_data <- data %>% 
    dplyr::filter(!is.na(!!sym(outcome)), 
                  !is.na(!!sym(pollutant)), 
                  !is.na(total_frp), 
                  !is.na(FRP_u_wind), 
                  !is.na(FRP_v_wind),
                  !is.na(RELATIVEHUMIDITYAVG),
                  !is.na(AmbientTemperatureAVG))
  
  # First stage model with columns from your dataset
  first_stage_formula <- as.formula(paste0(
    pollutant, " ~ total_frp + FRP_u_wind + FRP_v_wind + 
    RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
    factor(month) + factor(year)"
  ))
  
  first_stage <- lm(first_stage_formula, data = model_data)
  
  # Calculate residuals for control function approach
  model_data$residuals <- stats::residuals(first_stage)
  
  # Second stage model with columns from your dataset
  second_stage_formula <- as.formula(paste0(
    outcome, " ~ ", pollutant, " + residuals + 
    RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
    factor(month) + factor(year)"
  ))
  
  # Assuming a Poisson model for count data
  second_stage <- glm(second_stage_formula, 
                      family = poisson(link = "log"), 
                      data = model_data)
  
  return(list(
    first_stage = first_stage,
    second_stage = second_stage,
    data = model_data
  ))
}

# Function to implement Cumby-Huizinga test
cumby_huizinga_test <- function(residuals, lags = 1:10) {
  # Remove NAs
  residuals <- residuals[!is.na(residuals)]
  n <- length(residuals)
  
  # Results dataframe
  results <- data.frame(
    lag = lags,
    statistic = NA_real_,
    p_value = NA_real_
  )
  
  for (i in seq_along(lags)) {
    lag <- lags[i]
    
    # Skip if lag is too large compared to sample size
    if (lag >= n - 1) {
      next
    }
    
    # Calculate sample autocorrelation at this lag
    r <- stats::acf(residuals, lag.max = lag, plot = FALSE)$acf[lag + 1]
    
    # Calculate Cumby-Huizinga statistic (simplified version)
    ch_stat <- n * r^2 / (1 + 2 * sum(stats::acf(residuals, lag.max = lag - 1, plot = FALSE)$acf[-1]^2))
    
    # P-value from chi-squared distribution with 1 df
    p_value <- 1 - stats::pchisq(ch_stat, df = 1)
    
    # Store results
    results$statistic[i] <- ch_stat
    results$p_value[i] <- p_value
  }
  
  return(results)
}

# Function to test for autocorrelation in both first-stage and second-stage models
test_autocorrelation <- function(data, pollutant) {
  # Run 2SRI model
  model_results <- run_2sri_model(data, pollutant)
  
  # Test for autocorrelation in first stage
  first_stage_autocorr <- cumby_huizinga_test(
    stats::residuals(model_results$first_stage), 
    lags = c(1, 4)
  )
  
  # Test for autocorrelation in second stage
  second_stage_autocorr <- cumby_huizinga_test(
    stats::residuals(model_results$second_stage), 
    lags = c(1, 4)
  )
  
  # Print results
  cat("\n--- Autocorrelation results for", pollutant, "---\n")
  cat("First stage:\n")
  print(first_stage_autocorr)
  cat("\nSecond stage:\n")
  print(second_stage_autocorr)
  
  return(list(
    first_stage = first_stage_autocorr,
    second_stage = second_stage_autocorr
  ))
}

all_results <- list()

summary_results <- data.frame(
  Pollutant = character(),
  FirstStage_Lag1_Stat = numeric(),
  FirstStage_Lag1_PValue = numeric(),
  FirstStage_Lag4_Stat = numeric(),
  FirstStage_Lag4_PValue = numeric(),
  SecondStage_Lag1_Stat = numeric(),
  SecondStage_Lag1_PValue = numeric(),
  SecondStage_Lag4_Stat = numeric(),
  SecondStage_Lag4_PValue = numeric(),
  stringsAsFactors = FALSE
)

for (pollutant in pollutants) {
  tryCatch({
    cat("\nTesting pollutant:", pollutant, "\n")
    
    # Run autocorrelation test
    all_results[[pollutant]] <- test_autocorrelation(merged_daily, pollutant)
    
    # Extract results for summary table
    summary_results <- rbind(summary_results, data.frame(
      Pollutant = pollutant,
      FirstStage_Lag1_Stat = all_results[[pollutant]]$first_stage$statistic[1],
      FirstStage_Lag1_PValue = all_results[[pollutant]]$first_stage$p_value[1],
      FirstStage_Lag4_Stat = all_results[[pollutant]]$first_stage$statistic[2],
      FirstStage_Lag4_PValue = all_results[[pollutant]]$first_stage$p_value[2],
      SecondStage_Lag1_Stat = all_results[[pollutant]]$second_stage$statistic[1],
      SecondStage_Lag1_PValue = all_results[[pollutant]]$second_stage$p_value[1],
      SecondStage_Lag4_Stat = all_results[[pollutant]]$second_stage$statistic[2],
      SecondStage_Lag4_PValue = all_results[[pollutant]]$second_stage$p_value[2]
    ))
    
    # Save detailed model results to RDS
    saveRDS(all_results[[pollutant]], file = paste0("Output/autocorrelation_", pollutant, ".rds"))
    
  }, error = function(e) {
    cat("Error processing", pollutant, ":", conditionMessage(e), "\n")
  })
}

# Save summary results to CSV
write.csv(summary_results, file = "Output/autocorrelation_summary.csv", row.names = FALSE)

# Create a formatted version of the summary table
formatted_summary <- summary_results %>%
  mutate(across(ends_with("Stat"), ~ format(., digits = 2))) %>%
  mutate(across(ends_with("PValue"), ~ format(., digits = 4)))


#------------------------------------------------------------------------
# FALSIFICATION TEST
#------------------------------------------------------------------------

# Function to create a simulated falsification outcome
create_falsification_outcome <- function(data) {
  # Set seed for reproducibility
  set.seed(123)
  
  # Create a copy of the data
  sim_data <- data
  
  # Create a random count variable similar to deaths but unrelated to pollution
  # Using a seasonal pattern plus noise
  sim_data$sim_chickenpox <- with(sim_data, {
    # Seasonal component
    seasonal <- 5 * sin(2 * pi * as.numeric(month) / 12)
    
    # Annual trend
    yearly <- 1 * as.numeric(year - min(year))
    
    # Random component (Poisson)
    random <- rpois(nrow(sim_data), lambda = 10)
    
    # Combine components and ensure non-negative
    pmax(0, round(10 + seasonal + yearly + random))
  })
  
  return(sim_data)
}

# Function to run falsification test
run_falsification_test <- function(data, pollutants) {
  # Create data with falsification outcome
  sim_data <- create_falsification_outcome(data)
  
  # Results dataframe
  results <- data.frame(
    Pollutant = character(),
    Coefficient = numeric(),
    SE = numeric(),
    P_value = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    N = integer(),
    stringsAsFactors = FALSE
  )
  
  # Run analysis for each pollutant
  for (pollutant in pollutants) {
    cat("Running falsification test for", pollutant, "\n")
    
    # Filter data
    filtered_data <- sim_data %>%
      dplyr::filter(!is.na(sim_chickenpox), !is.na(!!sym(pollutant)),
                    !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
                    !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))
    
    # First stage
    first_stage_formula <- as.formula(paste0(
      pollutant, " ~ total_frp + FRP_u_wind + FRP_v_wind + ",
      "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
      "factor(month) + factor(year)"
    ))
    
    first_stage_model <- tryCatch({
      lm(first_stage_formula, data = filtered_data)
    }, error = function(e) {
      cat("Error in first stage:", e$message, "\n")
      next
    })
    
    # Create predicted values and residuals
    predicted_name <- paste0("predicted_", sub("_ugm3|AVG", "", pollutant))
    residuals_name <- paste0("residuals_", sub("_ugm3|AVG", "", pollutant))
    
    filtered_data[[predicted_name]] <- predict(first_stage_model)
    filtered_data[[residuals_name]] <- residuals(first_stage_model)
    
    # Second stage
    second_stage_formula <- as.formula(paste0(
      "sim_chickenpox ~ ", predicted_name, " + ", residuals_name, " + ",
      "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
      "factor(month) + factor(year)"
    ))
    
    second_stage_model <- tryCatch({
      glm(second_stage_formula, 
          family = poisson(link = "log"), 
          data = filtered_data)
    }, error = function(e) {
      cat("Error in second stage:", e$message, "\n")
      next
    })
    
    # Calculate Newey-West standard errors
    nw_vcov <- NeweyWest(second_stage_model, lag = 4, prewhite = FALSE)
    robust_results <- coeftest(second_stage_model, vcov = nw_vcov)
    
    # Extract coefficient for predicted pollutant
    poll_index <- which(rownames(robust_results) == predicted_name)
    poll_coef <- robust_results[poll_index, "Estimate"]
    poll_se <- robust_results[poll_index, "Std. Error"]
    poll_p <- robust_results[poll_index, "Pr(>|z|)"]
    
    # Calculate confidence intervals
    ci_lower <- poll_coef - 1.96 * poll_se
    ci_upper <- poll_coef + 1.96 * poll_se
    
    # Add to results
    results <- rbind(results, data.frame(
      Pollutant = pollutant,
      Coefficient = poll_coef,
      SE = poll_se,
      P_value = poll_p,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      N = nobs(second_stage_model)
    ))
  }
  
  return(results)
}

# Run the falsification test
falsification_results <- run_falsification_test(merged_daily, pollutants)

# Print results
print(falsification_results)

# Save results to a CSV file
write.csv(falsification_results, 
          file = "Output/falsification_test_results.csv", 
          row.names = FALSE)

# Create a nicely formatted table
formatted_table <- falsification_results %>%
  mutate(
    Coefficient = sprintf("%.4f", Coefficient),
    SE = sprintf("%.4f", SE),
    P_value = sprintf("%.4f", P_value),
    CI = sprintf("[%.4f, %.4f]", CI_Lower, CI_Upper)
  ) %>%
  select(Pollutant, Coefficient, SE, P_value, CI, N)

print(formatted_table)
