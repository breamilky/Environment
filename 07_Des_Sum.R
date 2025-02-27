# Load required packages
library(dplyr)
library(ggplot2)
library(knitr)
library(stringr)
library(tidyr)
library(lubridate)

# Create output directory if it doesn't exist
if (!dir.exists("Output")) {
  dir.create("Output")
}
if (!dir.exists("Output/Plots")) {
  dir.create("Output/Plots")
}

# Generate comprehensive descriptive statistics
descriptive_analysis <- function(data) {
  # Print dataset dimensions
  cat("Dataset dimensions:", nrow(data), "rows,", ncol(data), "columns\n\n")
  
  # 1. Define variable categories
  pollutant_vars <- c("PM10AVG", "PM2.5AVG", "O3AVG", "COAVG", "NOAVG", "NO2AVG", "NOXAVG", "SO2AVG",
                      "O3_ugm3", "CO_ugm3", "NO_ugm3", "NO2_ugm3", "NOX_ugm3", "SO2_ugm3")
  
  instrument_vars <- c("total_frp", "FRP_u_wind", "FRP_v_wind", "u_wind", "v_wind", "wind_direction_rad")
  
  weather_vars <- c("RELATIVEHUMIDITYAVG", "AmbientTemperatureAVG", "WINDSPEEDAVG", 
                    "SOLARRADIATIONAVG", "WINDDIRECTIONAVG")
  
  time_vars <- c("date", "year", "month")
  
  location_vars <- c("STATIONID", "LOCATION", "LOCATION_clean")
  
  outcome_vars <- c("deaths")
  
  # Filter to variables actually in the dataset
  available_pollutants <- intersect(pollutant_vars, names(data))
  available_instruments <- intersect(instrument_vars, names(data))
  available_weather <- intersect(weather_vars, names(data))
  available_time <- intersect(time_vars, names(data))
  available_location <- intersect(location_vars, names(data))
  available_outcomes <- intersect(outcome_vars, names(data))
  
  # Create category list
  all_categories <- list(
    "Pollutants" = available_pollutants,
    "Instrumental Variables" = available_instruments,
    "Weather Variables" = available_weather,
    "Time Variables" = available_time,
    "Location Variables" = available_location,
    "Health Outcomes" = available_outcomes
  )
  
  # 2. Basic descriptive statistics for numeric variables
  stats_df <- data.frame(
    Category = character(),
    Variable = character(),
    N = integer(),
    Mean = numeric(),
    SD = numeric(),
    Variance = numeric(),
    Min = numeric(),
    Q1 = numeric(),
    Median = numeric(),
    Q3 = numeric(),
    Max = numeric(),
    Zeros = integer(),
    Pct_Zeros = numeric(),
    Missing = integer(),
    Pct_Missing = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cat_name in names(all_categories)) {
    cat_vars <- all_categories[[cat_name]]
    
    for (var in cat_vars) {
      # Skip if not numeric or integer
      if (!is.numeric(data[[var]]) && !is.integer(data[[var]])) {
        next
      }
      
      x <- data[[var]]
      n_obs <- sum(!is.na(x))
      n_missing <- sum(is.na(x))
      pct_missing <- 100 * n_missing / length(x)
      
      # Create a base row with all possible columns
      row <- data.frame(
        Category = cat_name,
        Variable = var,
        N = n_obs,
        Mean = mean(x, na.rm = TRUE),
        SD = sd(x, na.rm = TRUE),
        Variance = var(x, na.rm = TRUE),
        Min = min(x, na.rm = TRUE),
        Q1 = quantile(x, 0.25, na.rm = TRUE),
        Median = median(x, na.rm = TRUE),
        Q3 = quantile(x, 0.75, na.rm = TRUE),
        Max = max(x, na.rm = TRUE),
        Zeros = sum(x == 0, na.rm = TRUE),
        Pct_Zeros = 100 * sum(x == 0, na.rm = TRUE) / n_obs,
        Missing = n_missing,
        Pct_Missing = pct_missing,
        stringsAsFactors = FALSE
      )
      
      stats_df <- rbind(stats_df, row)
    }
  }
  
  # 3. Create formatted version for display/export
  formatted_stats <- stats_df %>%
    mutate(
      Variable = str_replace(Variable, "AVG$", ""),
      Variable = str_replace(Variable, "_ugm3$", ""),
      Mean_SD = sprintf("%.2f (%.2f)", Mean, SD),
      Range = sprintf("%.2f - %.2f", Min, Max),
      Median_IQR = sprintf("%.2f [%.2f, %.2f]", Median, Q1, Q3),
      Zeros_Pct = sprintf("%d (%.1f%%)", Zeros, Pct_Zeros),
      Missing_Pct = sprintf("%d (%.1f%%)", Missing, Pct_Missing)
    )
  
  # 4. Time series analysis if date variable exists
  ts_summary <- NULL
  if ("date" %in% names(data) && is.Date(data$date)) {
    # Calculate date range
    date_range <- range(data$date, na.rm = TRUE)
    study_duration <- as.numeric(difftime(date_range[2], date_range[1], units = "days"))
    
    # Break down by month and year
    data$year_month <- format(data$date, "%Y-%m")
    
    # Monthly summary of outcomes
    if (length(available_outcomes) > 0) {
      monthly_outcomes <- data %>%
        group_by(year_month) %>%
        summarize(across(all_of(available_outcomes), 
                         ~ sum(.x, na.rm = TRUE)), 
                  n = n(),
                  .groups = "drop")
    }
    
    # Monthly summary of pollutants
    if (length(available_pollutants) > 0) {
      monthly_pollutants <- data %>%
        group_by(year_month) %>%
        summarize(across(all_of(available_pollutants), 
                         ~ mean(.x, na.rm = TRUE)), 
                  n = n(),
                  .groups = "drop")
    }
    
    ts_summary <- list(
      date_range = date_range,
      study_duration = study_duration,
      monthly_outcomes = if(exists("monthly_outcomes")) monthly_outcomes else NULL,
      monthly_pollutants = if(exists("monthly_pollutants")) monthly_pollutants else NULL
    )
  }
  
  # 5. Regional analysis
  regional_summary <- NULL
  if (any(c("LOCATION", "LOCATION_clean") %in% names(data))) {
    location_var <- ifelse("LOCATION_clean" %in% names(data), "LOCATION_clean", "LOCATION")
    
    # Calculate pollutant means by region
    if (length(available_pollutants) > 0) {
      regional_pollutants <- data %>%
        group_by(across(all_of(location_var))) %>%
        summarize(across(all_of(available_pollutants), 
                         ~ mean(.x, na.rm = TRUE)), 
                  n = n(),
                  .groups = "drop")
    }
    
    # Calculate outcome sums by region
    if (length(available_outcomes) > 0) {
      regional_outcomes <- data %>%
        group_by(across(all_of(location_var))) %>%
        summarize(across(all_of(available_outcomes), 
                         ~ sum(.x, na.rm = TRUE)), 
                  n = n(),
                  .groups = "drop")
    }
    
    regional_summary <- list(
      regional_pollutants = if(exists("regional_pollutants")) regional_pollutants else NULL,
      regional_outcomes = if(exists("regional_outcomes")) regional_outcomes else NULL
    )
  }
  
  # 6. Correlation analysis for pollutants
  corr_matrix <- NULL
  if (length(available_pollutants) > 1) {
    corr_matrix <- data %>%
      select(all_of(available_pollutants)) %>%
      cor(use = "pairwise.complete.obs")
  }
  
  # 7. Instrument strength analysis
  instrument_analysis <- NULL
  iv_cols <- c("total_frp", "FRP_u_wind", "FRP_v_wind")
  if (all(iv_cols %in% names(data)) && length(available_pollutants) > 0) {
    # Create models to test instrumental variable strength
    iv_strength <- data.frame(
      Pollutant = character(),
      R_squared = numeric(),
      F_statistic = numeric(),
      P_value = numeric(),
      N = integer(),
      stringsAsFactors = FALSE
    )
    
    for (pollutant in available_pollutants) {
      tryCatch({
        model_data <- data %>%
          filter(!is.na(!!sym(pollutant)),
                 !is.na(total_frp),
                 !is.na(FRP_u_wind),
                 !is.na(FRP_v_wind))
        
        # First stage regression
        first_stage_formula <- as.formula(paste0(
          pollutant, " ~ total_frp + FRP_u_wind + FRP_v_wind"
        ))
        
        model <- lm(first_stage_formula, data = model_data)
        model_summary <- summary(model)
        
        # Extract key statistics
        iv_strength <- rbind(iv_strength, data.frame(
          Pollutant = pollutant,
          R_squared = model_summary$r.squared,
          F_statistic = model_summary$fstatistic[1],
          P_value = pf(model_summary$fstatistic[1], 
                       model_summary$fstatistic[2], 
                       model_summary$fstatistic[3], 
                       lower.tail = FALSE),
          N = nobs(model)
        ))
      }, error = function(e) {
        cat("Error analyzing instrument strength for", pollutant, ":", conditionMessage(e), "\n")
      })
    }
    
    instrument_analysis <- iv_strength
  }
  
  # Return all results
  return(list(
    basic_stats = stats_df,
    formatted_stats = formatted_stats,
    time_series = ts_summary,
    regional = regional_summary,
    correlations = corr_matrix,
    instrument_strength = instrument_analysis
  ))
}

# Run the analysis
tryCatch({
  results <- descriptive_analysis(merged_daily)
  
  # Save detailed results
  saveRDS(results, file = "Output/descriptive_analysis_full.rds")
  
  # Write formatted tables to CSV
  write.csv(results$formatted_stats, file = "Output/descriptive_statistics.csv", row.names = FALSE)
  
  if (!is.null(results$instrument_strength)) {
    write.csv(results$instrument_strength, file = "Output/instrument_strength.csv", row.names = FALSE)
  }
  
  if (!is.null(results$correlations)) {
    write.csv(results$correlations, file = "Output/pollutant_correlations.csv", row.names = FALSE)
  }
  
  # Create comprehensive report
  sink("Output/comprehensive_descriptive_report.txt")
  cat("COMPREHENSIVE DESCRIPTIVE ANALYSIS REPORT\n")
  cat("========================================\n\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Table 1: Basic Statistics
  cat("Table 1: Descriptive Statistics of Key Variables\n\n")
  
  # Create display table
  display_stats <- results$formatted_stats %>%
    select(Category, Variable, N, Mean_SD, Median_IQR, Range, Missing_Pct)
  
  # For outcome variables, include zero counts
  outcome_vars <- results$formatted_stats$Variable[results$formatted_stats$Category == "Health Outcomes"]
  if (length(outcome_vars) > 0) {
    display_stats_outcomes <- results$formatted_stats %>%
      filter(Category == "Health Outcomes") %>%
      select(Category, Variable, N, Mean_SD, Median_IQR, Range, Zeros_Pct, Missing_Pct)
    
    display_stats_others <- results$formatted_stats %>%
      filter(Category != "Health Outcomes") %>%
      select(Category, Variable, N, Mean_SD, Median_IQR, Range, Missing_Pct)
    
    # Print tables separately
    cat("Health Outcome Variables:\n")
    print(kable(display_stats_outcomes, 
                col.names = c("Category", "Variable", "N", "Mean (SD)", 
                              "Median [IQR]", "Range", "Zeros (%)", "Missing (%)")))
    
    cat("\n\nEnvironmental and Other Variables:\n")
    print(kable(display_stats_others, 
                col.names = c("Category", "Variable", "N", "Mean (SD)", 
                              "Median [IQR]", "Range", "Missing (%)")))
  } else {
    print(kable(display_stats, 
                col.names = c("Category", "Variable", "N", "Mean (SD)", 
                              "Median [IQR]", "Range", "Missing (%)")))
  }
  
  # Rest of report code remains the same
  # Table 2: Time Series Summary
  if (!is.null(results$time_series)) {
    cat("\n\nTable 2: Time Series Overview\n\n")
    cat("Study period:", format(results$time_series$date_range[1], "%Y-%m-%d"), 
        "to", format(results$time_series$date_range[2], "%Y-%m-%d"), "\n")
    cat("Duration:", results$time_series$study_duration, "days\n\n")
    
    if (!is.null(results$time_series$monthly_outcomes)) {
      cat("Monthly Health Outcomes:\n")
      print(head(results$time_series$monthly_outcomes, 10))
      if (nrow(results$time_series$monthly_outcomes) > 10) {
        cat("(showing first 10 of", nrow(results$time_series$monthly_outcomes), "months)\n")
      }
    }
  }
  
  # Table 3: Regional Analysis
  if (!is.null(results$regional) && !is.null(results$regional$regional_pollutants)) {
    cat("\n\nTable 3: Regional Pollution Levels\n\n")
    print(results$regional$regional_pollutants)
    
    if (!is.null(results$regional$regional_outcomes)) {
      cat("\nRegional Health Outcomes:\n")
      print(results$regional$regional_outcomes)
    }
  }
  
  # Table 4: Pollutant Correlations
  if (!is.null(results$correlations)) {
    cat("\n\nTable 4: Pollutant Correlation Matrix\n\n")
    print(round(results$correlations, 3))
  }
  
  # Table 5: Instrumental Variable Strength
  if (!is.null(results$instrument_strength)) {
    cat("\n\nTable 5: Instrumental Variable Strength\n\n")
    formatted_iv <- results$instrument_strength %>%
      mutate(
        R_squared = sprintf("%.3f", R_squared),
        F_statistic = sprintf("%.2f", F_statistic),
        P_value = sprintf("%.4f", P_value)
      )
    print(formatted_iv)
    
    cat("\nNote: F-statistics > 10 indicate strong instruments\n")
  }
  
  # Summary highlights
  cat("\n\nKey Findings:\n\n")
  
  # Death statistics
  death_vars <- results$basic_stats$Variable[results$basic_stats$Category == "Health Outcomes"]
  if (length(death_vars) > 0) {
    death_stats <- results$basic_stats[results$basic_stats$Variable %in% death_vars, ]
    cat("1. Mortality: Average of", round(death_stats$Mean, 3), 
        "deaths per observation, with", round(death_stats$Pct_Zeros, 1), 
        "% of observations having zero deaths\n")
  }
  
  # Pollutant levels
  pm25_var <- intersect(c("PM2.5AVG"), results$basic_stats$Variable)
  if (length(pm25_var) > 0) {
    pm25_stats <- results$basic_stats[results$basic_stats$Variable == pm25_var, ]
    cat("2. Fine particulate matter (PM2.5): Mean concentration of", 
        round(pm25_stats$Mean, 2), "μg/m³ (SD =", round(pm25_stats$SD, 2), ")\n")
  }
  
  pm10_var <- intersect(c("PM10AVG"), results$basic_stats$Variable)
  if (length(pm10_var) > 0) {
    pm10_stats <- results$basic_stats[results$basic_stats$Variable == pm10_var, ]
    cat("3. Coarse particulate matter (PM10): Mean concentration of", 
        round(pm10_stats$Mean, 2), "μg/m³ (SD =", round(pm10_stats$SD, 2), ")\n")
  }
  
  # Instrument strength
  if (!is.null(results$instrument_strength)) {
    strong_instruments <- results$instrument_strength$Pollutant[results$instrument_strength$F_statistic > 10]
    weak_instruments <- results$instrument_strength$Pollutant[results$instrument_strength$F_statistic <= 10]
    
    if (length(strong_instruments) > 0) {
      cat("4. Strong instrumental variables for:", paste(strong_instruments, collapse=", "), "\n")
    }
    
    if (length(weak_instruments) > 0) {
      cat("5. Weak instrumental variables for:", paste(weak_instruments, collapse=", "), "\n")
    }
  }
  
  # Missing data
  high_missing <- results$basic_stats$Variable[results$basic_stats$Pct_Missing > 20]
  if (length(high_missing) > 0) {
    cat("6. Variables with high missing data (>20%):", paste(high_missing, collapse=", "), "\n")
  }
  
  sink()
  
  # Generate key plots
  # Pollutant distributions
  for (pollutant in results$formatted_stats$Variable[results$formatted_stats$Category == "Pollutants"]) {
    # Get original variable name
    orig_var <- results$basic_stats$Variable[results$basic_stats$Variable == pollutant | 
                                               (results$basic_stats$Variable == paste0(pollutant, "AVG")) |
                                               (results$basic_stats$Variable == paste0(pollutant, "_ugm3"))]
    
    if (length(orig_var) > 0) {
      orig_var <- orig_var[1]
      
      tryCatch({
        p <- ggplot(merged_daily, aes_string(x = orig_var)) +
          geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
          labs(title = paste("Distribution of", pollutant),
               x = pollutant,
               y = "Frequency") +
          theme_minimal()
        
        ggsave(paste0("Output/Plots/", pollutant, "_distribution.png"), p, width = 8, height = 6)
      }, error = function(e) {
        cat("Error creating histogram for", pollutant, ":", conditionMessage(e), "\n")
      })
    }
  }
  
  # Health outcome distribution
  if (length(death_vars) > 0) {
    for (death_var in death_vars) {
      tryCatch({
        p <- ggplot(merged_daily, aes_string(x = death_var)) +
          geom_histogram(bins = max(15, min(30, max(merged_daily[[death_var]], na.rm = TRUE))), 
                         fill = "darkred", color = "black", alpha = 0.7) +
          labs(title = "Distribution of Death Counts",
               x = "Deaths",
               y = "Frequency") +
          theme_minimal()
        
        ggsave(paste0("Output/Plots/deaths_distribution.png"), p, width = 8, height = 6)
      }, error = function(e) {
        cat("Error creating histogram for deaths:", conditionMessage(e), "\n")
      })
    }
  }
  
  cat("\nComprehensive descriptive analysis complete. Results saved to the Output directory.\n")
  
}, error = function(e) {
  cat("Error in descriptive analysis:", conditionMessage(e), "\n")
  print(traceback())
})