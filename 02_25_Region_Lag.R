print("Loading packages...")
# Define the list of required packages
packages <- c("ggplot2", "dplyr", "tidyr", "viridis", "gridExtra", 
              "ggridges", "cowplot", "boot", "sandwich", "lmtest", 
              "knitr", "kableExtra", "broom", "purrr", "AER", 
              "stargazer", "lubridate")

# Loop through the packages: install if missing, then load them
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

# Removing implausible NO values 
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


# List of all pollutants to analyze
pollutants <- c("PM10AVG", 
                "PM2.5AVG", 
                "O3_ugm3", 
                "CO_ugm3", 
                "NO_ugm3", 
                "NO2_ugm3", 
                "NOX_ugm3", 
                "SO2_ugm3")

#------------------------------------------------------------------------
# 1. Regional Analysis for All Pollutants
#------------------------------------------------------------------------

# Function to run regional stratification IV analysis (same as your original)
run_regional_mortality_iv <- function(data, pollutant) {
  # Get unique regions (locations)
  regions <- unique(data$LOCATION_clean)
  
  # Run analysis for each region
  regional_results <- map_df(regions, function(region) {
    # Filter data for this region
    region_data <- data %>% 
      filter(LOCATION_clean == region,
             !is.na(deaths), !is.na(!!sym(pollutant)),
             !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
             !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))
    
    # Check if we have enough data
    if (nrow(region_data) < 30) {
      # Skip regions with too little data
      return(NULL)
    }
    
    # Formula for IV regression
    formula_str <- paste0(
      "deaths ~ ", pollutant, " + RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
      " factor(month) + factor(year) | total_frp + FRP_u_wind + FRP_v_wind +",
      " RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + factor(month) + factor(year)"
    )
    
    # Run IV regression with error handling
    model <- tryCatch({
      ivreg_model <- ivreg(as.formula(formula_str), data = region_data)
      
      # Calculate Newey-West standard errors
      nw_vcov <- NeweyWest(ivreg_model, lag = 4, prewhite = FALSE)
      nw_coef_test <- coeftest(ivreg_model, vcov = nw_vcov)
      
      # Extract coefficient information
      poll_index <- 2  # Assuming pollutant is the second coefficient
      poll_coef <- coef(ivreg_model)[poll_index]
      poll_se <- sqrt(diag(nw_vcov))[poll_index]
      poll_p <- nw_coef_test[poll_index, 4]
      
      # Return summary as tibble
      tibble(
        Region = region,
        Pollutant = pollutant,
        N = nobs(ivreg_model),
        Coefficient = poll_coef,
        SE = poll_se,
        P_value = poll_p,
        CI_Lower = poll_coef - 1.96 * poll_se,
        CI_Upper = poll_coef + 1.96 * poll_se
      )
    }, 
    error = function(e) {
      # Return NULL if error occurs
      return(NULL)
    })
    
    return(model)
  })
  
  return(regional_results)
}

# Run regional analysis for all pollutants
all_regional_results <- list()

for (pollutant in pollutants) {
  # Run regional analysis
  regional_results <- run_regional_mortality_iv(merged_daily, pollutant)
  
  # Store results if available
  if (!is.null(regional_results) && nrow(regional_results) > 0) {
    all_regional_results[[pollutant]] <- regional_results
    
    # Create forest plot for this pollutant
    regional_plot <- ggplot(regional_results, aes(x = reorder(Region, Coefficient), y = Coefficient)) +
      geom_point(size = 3, color = "darkblue") +
      geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "darkblue") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = paste("Regional Effects of", pollutant, "on Mortality"),
           subtitle = "Instrumental Variable Estimates by Location",
           x = "",
           y = "Effect Size (95% CI)") +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray80")
      )
    
    # Save the regional plot
    ggsave(paste0("Output/regional_mortality_effects_", gsub("[.]", "_", pollutant), ".png"), 
           regional_plot, width = 10, height = 8, dpi = 300)
    
    # Save regional results to CSV
    write.csv(regional_results, 
              paste0("Output/regional_mortality_results_", gsub("[.]", "_", pollutant), ".csv"), 
              row.names = FALSE)
  }
}

# Combine all regional results into one dataframe for comparison
combined_regional_results <- bind_rows(all_regional_results)

if (nrow(combined_regional_results) > 0) {
  # Save combined results
  write.csv(combined_regional_results, "Output/combined_regional_results.csv", row.names = FALSE)
  
  # Create a faceted plot comparing all pollutants
  combined_regional_plot <- ggplot(combined_regional_results, 
                                   aes(x = reorder(Region, Coefficient), y = Coefficient)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +
    facet_wrap(~ Pollutant, scales = "free_x") +
    labs(title = "Regional Effects of All Pollutants on Mortality",
         subtitle = "Instrumental Variable Estimates by Location",
         x = "",
         y = "Effect Size (95% CI)") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 7),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      panel.spacing = unit(1, "lines")
    )
  
  # Save combined plot
  ggsave("Output/combined_regional_effects.png", combined_regional_plot, 
         width = 14, height = 12, dpi = 300)
}

#------------------------------------------------------------------------
# 2. Lag Analysis for All Pollutants
#------------------------------------------------------------------------

# Function to analyze lag effects (same as your original)
analyze_lag_effects <- function(data, pollutant, max_lag = 7) {
  # Create lagged pollutant variables
  lag_data <- data
  
  for (lag in 1:max_lag) {
    # Create unique identifier to enable proper lagging
    lag_data <- lag_data %>%
      group_by(LOCATION_clean) %>%
      mutate(!!paste0(pollutant, "_lag", lag) := lag(!!sym(pollutant), n = lag)) %>%
      ungroup()
  }
  
  # Run IV for each lag
  lag_results <- tibble()
  
  # Current pollutant (lag 0)
  current_formula <- paste0(
    "deaths ~ ", pollutant, " + RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
    " factor(month) + factor(year) | total_frp + FRP_u_wind + FRP_v_wind +",
    " RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + factor(month) + factor(year)"
  )
  
  current_model <- tryCatch({
    # Filter data
    filtered_data <- lag_data %>%
      filter(!is.na(deaths), !is.na(!!sym(pollutant)),
             !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
             !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))
    
    ivreg(as.formula(current_formula), data = filtered_data)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(current_model)) {
    # Get Newey-West standard errors
    nw_vcov <- NeweyWest(current_model, lag = 4, prewhite = FALSE)
    nw_coef_test <- coeftest(current_model, vcov = nw_vcov)
    
    # Extract coefficient information
    coef_index <- 2  # Assuming pollutant is the second coefficient
    coef_val <- coef(current_model)[coef_index]
    se_val <- sqrt(diag(nw_vcov))[coef_index]
    p_val <- nw_coef_test[coef_index, 4]
    
    # Add to results
    lag_results <- bind_rows(lag_results, tibble(
      Lag = 0,
      Pollutant = pollutant,
      N = nobs(current_model),
      Coefficient = coef_val,
      SE = se_val,
      P_value = p_val,
      CI_Lower = coef_val - 1.96 * se_val,
      CI_Upper = coef_val + 1.96 * se_val
    ))
  }
  
  # Check all lags
  for (lag in 1:max_lag) {
    lag_var <- paste0(pollutant, "_lag", lag)
    
    lag_formula <- paste0(
      "deaths ~ ", lag_var, " + RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +",
      " factor(month) + factor(year) | total_frp + FRP_u_wind + FRP_v_wind +",
      " RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + factor(month) + factor(year)"
    )
    
    lag_model <- tryCatch({
      # Filter data
      filtered_data <- lag_data %>%
        filter(!is.na(deaths), !is.na(!!sym(lag_var)),
               !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
               !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))
      
      ivreg(as.formula(lag_formula), data = filtered_data)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(lag_model)) {
      # Get Newey-West standard errors
      nw_vcov <- NeweyWest(lag_model, lag = 4, prewhite = FALSE)
      nw_coef_test <- coeftest(lag_model, vcov = nw_vcov)
      
      # Extract coefficient information
      coef_index <- 2  # Assuming pollutant is the second coefficient
      coef_val <- coef(lag_model)[coef_index]
      se_val <- sqrt(diag(nw_vcov))[coef_index]
      p_val <- nw_coef_test[coef_index, 4]
      
      # Add to results
      lag_results <- bind_rows(lag_results, tibble(
        Lag = lag,
        Pollutant = pollutant,
        N = nobs(lag_model),
        Coefficient = coef_val,
        SE = se_val,
        P_value = p_val,
        CI_Lower = coef_val - 1.96 * se_val,
        CI_Upper = coef_val + 1.96 * se_val
      ))
    }
  }
  
  return(lag_results)
}

# Run lag analysis for all pollutants
all_lag_results <- list()

for (pollutant in pollutants) {
  # Run lag analysis
  lag_results <- analyze_lag_effects(merged_daily, pollutant, max_lag = 7)
  
  # Store results if available
  if (!is.null(lag_results) && nrow(lag_results) > 0) {
    all_lag_results[[pollutant]] <- lag_results
    
    # Create lag plot
    lag_plot <- ggplot(lag_results, aes(x = factor(Lag), y = Coefficient)) +
      geom_point(size = 3, color = "darkblue") +
      geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "darkblue") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Lag Effects of", pollutant, "on Mortality"),
           subtitle = "Instrumental Variable Estimates with 95% Confidence Intervals",
           x = "Lag (Days)",
           y = "Effect Size") +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray80")
      )
    
    # Save lag plot
    ggsave(paste0("Output/lag_effects_plot_", gsub("[.]", "_", pollutant), ".png"), 
           lag_plot, width = 10, height = 6, dpi = 300)
    
    # Save lag results
    write.csv(lag_results, 
              paste0("Output/lag_effects_results_", gsub("[.]", "_", pollutant), ".csv"), 
              row.names = FALSE)
  }
}

# Combine all lag results into one dataframe for comparison
combined_lag_results <- bind_rows(all_lag_results)

if (nrow(combined_lag_results) > 0) {
  # Save combined results
  write.csv(combined_lag_results, "Output/combined_lag_results.csv", row.names = FALSE)
  
  # Create a faceted plot comparing all pollutants
  combined_lag_plot <- ggplot(combined_lag_results, 
                              aes(x = factor(Lag), y = Coefficient, group = 1)) +
    geom_point(size = 2) +
    geom_line(linetype = "dashed") +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Pollutant, scales = "free_y") +
    labs(title = "Lag Effects of All Pollutants on Mortality",
         subtitle = "Instrumental Variable Estimates with 95% Confidence Intervals",
         x = "Lag (Days)",
         y = "Effect Size") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      panel.spacing = unit(1, "lines"),
      panel.grid.minor = element_blank()
    )
  
  # Save combined plot
  ggsave("Output/combined_lag_effects.png", combined_lag_plot, 
         width = 14, height = 10, dpi = 300)
}

#------------------------------------------------------------------------
# 3. Summary Table of Results
#------------------------------------------------------------------------

# Create a summary table for quick comparison
if (nrow(combined_regional_results) > 0) {
  # Summarize regional results by pollutant
  regional_summary <- combined_regional_results %>%
    group_by(Pollutant) %>%
    summarize(
      N_Regions = n(),
      Avg_Coefficient = mean(Coefficient, na.rm = TRUE),
      Min_Coefficient = min(Coefficient, na.rm = TRUE),
      Max_Coefficient = max(Coefficient, na.rm = TRUE),
      Significant_Positive = sum(CI_Lower > 0, na.rm = TRUE),
      Significant_Negative = sum(CI_Upper < 0, na.rm = TRUE),
      Non_Significant = sum(CI_Lower <= 0 & CI_Upper >= 0, na.rm = TRUE)
    )
  
  # Save summary
  write.csv(regional_summary, "Output/regional_analysis_summary.csv", row.names = FALSE)
}

if (nrow(combined_lag_results) > 0) {
  # Summarize lag results by pollutant
  lag_summary <- combined_lag_results %>%
    group_by(Pollutant) %>%
    summarize(
      Max_Effect_Lag = Lag[which.max(abs(Coefficient))],
      Max_Effect_Coef = Coefficient[which.max(abs(Coefficient))],
      Immediate_Effect = Coefficient[Lag == 0],
      Longest_Significant_Lag = ifelse(
        any(P_value < 0.05),
        max(Lag[P_value < 0.05]),
        NA
      )
    )
  
  # Save summary
  write.csv(lag_summary, "Output/lag_analysis_summary.csv", row.names = FALSE)
}

# Create a combined pollutant comparison table
if (exists("regional_summary") && exists("lag_summary")) {
  pollutant_comparison <- full_join(
    regional_summary, 
    lag_summary,
    by = "Pollutant"
  )
  
  write.csv(pollutant_comparison, "Output/pollutant_comparison_summary.csv", row.names = FALSE)
}

# Print completion message
cat("\nAnalysis completed for all pollutants.\n")
cat("Summary of outputs:\n")
cat("  - Individual regional plots and CSV files for each pollutant\n")
cat("  - Combined regional effects plot\n")
cat("  - Individual lag plots and CSV files for each pollutant\n")
cat("  - Combined lag effects plot\n")
cat("  - Summary tables of all analyses\n")
if (exists("plotly")) {
  cat("  - Interactive HTML visualizations\n")
}