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
# First-Stage Diagnostics for All Pollutants
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

#------------------------------------------------------------------------
# COUNTERFACTUAL ANALYSIS
#------------------------------------------------------------------------

# Create output directory if it doesn't exist
output_dir <- "Output"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to run two-stage residual inclusion model (2SRI)
run_2sri_model <- function(data, pollutant, outcome = "deaths") {
  # Filter data
  filtered_data <- data %>%
    filter(!is.na(!!sym(outcome)), !is.na(!!sym(pollutant)),
           !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
           !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))
  
  # First stage formula
  first_stage_formula <- as.formula(paste0(
    pollutant, " ~ total_frp + FRP_u_wind + FRP_v_wind + ",
    "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
    "factor(month) + factor(year)"
  ))
  
  # Run first stage regression
  first_stage_model <- lm(first_stage_formula, data = filtered_data)
  
  # Create predicted values and residuals
  predicted_name <- paste0("predicted_", sub("_ugm3|AVG", "", pollutant))
  residuals_name <- paste0("residuals_", sub("_ugm3|AVG", "", pollutant))
  
  filtered_data[[predicted_name]] <- predict(first_stage_model)
  filtered_data[[residuals_name]] <- residuals(first_stage_model)
  
  # Second stage formula with residual inclusion
  second_stage_formula <- as.formula(paste0(
    outcome, " ~ ", predicted_name, " + ", residuals_name, " + ",
    "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
    "factor(month) + factor(year)"
  ))
  
  # Run second stage - using GLM with Poisson for count data
  second_stage_model <- glm(
    second_stage_formula,
    family = poisson(link = "log"),
    data = filtered_data
  )
  
  # Extract coefficient for predicted pollutant
  coef_summary <- summary(second_stage_model)$coefficients
  poll_coef <- coef_summary[predicted_name, "Estimate"]
  poll_se <- coef_summary[predicted_name, "Std. Error"]
  poll_p <- coef_summary[predicted_name, "Pr(>|z|)"]
  
  # Return models and coefficient
  return(list(
    first_stage = first_stage_model,
    second_stage = second_stage_model,
    coefficient = poll_coef,
    se = poll_se,
    p_value = poll_p,
    data = filtered_data,
    predicted_name = predicted_name,
    residuals_name = residuals_name
  ))
}

# Function to calculate counterfactual deaths without fires
calculate_counterfactual_deaths <- function(data, pollutant, years = NULL) {
  cat("  Calculating counterfactual for", pollutant, "\n")
  
  # If years not provided, detect from data
  if (is.null(years)) {
    years <- sort(unique(data$year))
    cat("  Using detected years:", paste(years, collapse = ", "), "\n")
  }
  
  # Results dataframe
  results <- data.frame(
    Year = years,
    Percent_Increase = numeric(length(years)),
    Absolute_Increase = numeric(length(years))
  )
  
  # First, run the main model
  cat("  Running main 2SRI model...\n")
  model_result <- tryCatch({
    run_2sri_model(data, pollutant)
  }, error = function(e) {
    cat("  Error in main model:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(model_result)) {
    cat("  Failed to run main model. Skipping counterfactual analysis.\n")
    return(NULL)
  }
  
  # Extract coefficient for pollutant
  poll_coef <- model_result$coefficient
  cat("  Estimated coefficient:", poll_coef, "\n")
  
  # For each year, calculate the counterfactual
  for (i in seq_along(years)) {
    year <- years[i]
    cat("  Processing year", year, "...\n")
    
    # Filter data for this year
    year_data <- data %>%
      filter(year == !!year,
             !is.na(deaths), !is.na(!!sym(pollutant)),
             !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
             !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))
    
    if (nrow(year_data) == 0) {
      cat("  No data available for year", year, "\n")
      results$Percent_Increase[i] <- NA
      results$Absolute_Increase[i] <- NA
      next
    }
    
    cat("  Found", nrow(year_data), "observations for year", year, "\n")
    
    # First stage for this year
    first_stage_formula <- as.formula(paste0(
      pollutant, " ~ total_frp + FRP_u_wind + FRP_v_wind + ",
      "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
      "factor(month)"
    ))
    
    first_stage_model <- tryCatch({
      lm(first_stage_formula, data = year_data)
    }, error = function(e) {
      cat("  Error in year-specific model:", e$message, "\n")
      results$Percent_Increase[i] <- NA
      results$Absolute_Increase[i] <- NA
      return(NULL)
    })
    
    if (is.null(first_stage_model)) {
      next
    }
    
    # Create counterfactual data with no fires
    counterfactual_data <- year_data
    counterfactual_data$total_frp <- 0
    counterfactual_data$FRP_u_wind <- 0
    counterfactual_data$FRP_v_wind <- 0
    
    # Predict pollutant levels
    actual_pollution <- predict(first_stage_model, year_data)
    cf_pollution <- predict(first_stage_model, counterfactual_data)
    
    # Calculate difference in pollution
    pollution_diff <- mean(actual_pollution) - mean(cf_pollution)
    cat("  Estimated pollution reduction:", pollution_diff, "\n")
    
    # For Poisson models, effect is exp(beta*X) - 1
    effect_multiplier <- exp(poll_coef * pollution_diff) - 1
    
    # Calculate absolute and percentage increase
    mean_deaths <- mean(year_data$deaths)
    absolute_increase <- mean_deaths * effect_multiplier
    
    cat("  Estimated percent increase in mortality:", round(effect_multiplier * 100, 2), "%\n")
    cat("  Estimated absolute increase in mortality:", round(absolute_increase, 2), "\n")
    
    # Store results
    results$Percent_Increase[i] <- effect_multiplier * 100
    results$Absolute_Increase[i] <- absolute_increase
  }
  
  return(results)
}

# Function to run counterfactual analysis for all pollutants
run_all_counterfactuals <- function(data, pollutants, years = NULL) {
  cat("\n====== RUNNING COUNTERFACTUAL ANALYSIS ======\n")
  
  # If years not provided, detect from data
  if (is.null(years)) {
    years <- sort(unique(data$year))
    cat("Using years:", paste(years, collapse = ", "), "\n")
  }
  
  # Create results container
  counterfactual_results <- list()
  
  # For each pollutant
  for (pollutant in pollutants) {
    cat("\nRunning counterfactual for", pollutant, "\n")
    
    # Run counterfactual
    cf_result <- tryCatch({
      calculate_counterfactual_deaths(data, pollutant, years)
    }, error = function(e) {
      cat("Error:", e$message, "\n")
      return(NULL)
    })
    
    # Store results
    counterfactual_results[[pollutant]] <- cf_result
  }
  
  # Create summary dataframe
  summary_df <- data.frame(
    Pollutant = character(),
    Year = numeric(),
    Percent_Increase = numeric(),
    Absolute_Increase = numeric()
  )
  
  for (pollutant in names(counterfactual_results)) {
    result <- counterfactual_results[[pollutant]]
    
    if (!is.null(result)) {
      for (i in 1:nrow(result)) {
        if (!is.na(result$Percent_Increase[i])) {
          summary_df <- rbind(summary_df, data.frame(
            Pollutant = pollutant,
            Year = result$Year[i],
            Percent_Increase = result$Percent_Increase[i],
            Absolute_Increase = result$Absolute_Increase[i]
          ))
        }
      }
    }
  }
  
  # Save results to CSV
  write.csv(summary_df, file.path(output_dir, "counterfactual_summary.csv"), row.names = FALSE)
  
  # Create detailed results by pollutant and year
  detailed_results <- summary_df %>%
    mutate(
      Percent_Increase = round(Percent_Increase, 2),
      Absolute_Increase = round(Absolute_Increase, 3)
    )
  
  # Save detailed results
  write.csv(detailed_results, file.path(output_dir, "counterfactual_detailed.csv"), row.names = FALSE)
  
  # Create FIXED pivot tables that won't encounter the error
  # 1. First for percent increases
  percent_pivot <- summary_df %>%
    select(Pollutant, Year, Percent_Increase) %>%
    tidyr::pivot_wider(
      names_from = Year,
      values_from = Percent_Increase,
      names_prefix = "Percent_"
    )
  
  # 2. Then for absolute increases  
  absolute_pivot <- summary_df %>%
    select(Pollutant, Year, Absolute_Increase) %>%
    tidyr::pivot_wider(
      names_from = Year,
      values_from = Absolute_Increase,
      names_prefix = "Absolute_"
    )
  
  # 3. Join them together
  pivot_table <- dplyr::left_join(percent_pivot, absolute_pivot, by = "Pollutant")
  
  # Save pivot table
  write.csv(pivot_table, file.path(output_dir, "counterfactual_pivot.csv"), row.names = FALSE)
  
  # Create a clean summary focusing on key pollutants
  clean_summary <- summary_df %>%
    # Round for readability
    mutate(
      Percent_Increase = round(Percent_Increase, 2),
      Absolute_Increase = round(Absolute_Increase, 3)
    ) %>%
    # Filter to positive pollutants only for clearer findings
    filter(Pollutant %in% c("PM10AVG", "PM2.5AVG", "O3_ugm3")) %>%
    # Arrange by year and pollutant
    arrange(Year, Pollutant)
  
  # Save clean summary
  write.csv(clean_summary, file.path(output_dir, "key_pollutants_summary.csv"), row.names = FALSE)
  
  # Create visualizations
  
  # 1. Time series plot of percent increases
  if (nrow(summary_df) > 0) {
    # Time series of percent increases for all pollutants
    p1 <- ggplot(summary_df, aes(x = Year, y = Percent_Increase, color = Pollutant, group = Pollutant)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      theme_minimal() +
      labs(
        title = "Estimated Increase in Mortality Due to Fires Over Time",
        subtitle = "Percentage increase by pollutant and year",
        x = "Year",
        y = "Percent Increase in Mortality (%)"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black")
    
    # Save the plot
    ggsave(file.path(output_dir, "counterfactual_timeseries.png"), p1, width = 10, height = 6, dpi = 300)
    
    # 2. Focus on just positive pollutants
    if (nrow(clean_summary) > 0) {
      p2 <- ggplot(clean_summary, aes(x = Year, y = Percent_Increase, color = Pollutant, group = Pollutant)) +
        geom_line(size = 1.2) +
        geom_point(size = 3) +
        theme_minimal() +
        labs(
          title = "Estimated Mortality Increase from Key Pollutants",
          subtitle = "Particulate matter and ozone from fires",
          x = "Year",
          y = "Percent Increase in Mortality (%)"
        ) +
        theme(legend.position = "right")
      
      # Save the plot
      ggsave(file.path(output_dir, "key_pollutants_timeseries.png"), p2, width = 10, height = 6, dpi = 300)
    }
    
    # 3. Bar plot for pollutants by year
    p3 <- ggplot(summary_df, aes(x = Pollutant, y = Percent_Increase, fill = factor(Year))) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(
        title = "Mortality Impact by Pollutant and Year",
        x = "Pollutant",
        y = "Percent Increase in Mortality (%)",
        fill = "Year"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Save the plot
    ggsave(file.path(output_dir, "counterfactual_by_year.png"), p3, width = 12, height = 7, dpi = 300)
  }
  
  cat("\nCounterfactual analysis complete. Results saved to", output_dir, "directory.\n")
  
  return(list(
    by_pollutant = counterfactual_results,
    summary = summary_df,
    pivot_table = pivot_table,
    clean_summary = clean_summary
  ))
}

counterfactual_results <- run_all_counterfactuals(merged_daily, pollutants)
