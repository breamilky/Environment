#---------- ROBUSTNESS CHECKS ----------#

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

run_robustness_checks_all_pollutants <- function(data, pollutants) {
  # Create empty results dataframe
  all_results <- data.frame(
    Pollutant = character(),
    Specification = character(),
    Coefficient = numeric(),
    SE = numeric(),
    P_value = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    N = integer()
  )
  
  # Process each pollutant
  for (pollutant in pollutants) {
    cat("\nRunning robustness checks for", pollutant, "\n")
    
    # Different specifications
    specifications <- list(
      # Base model
      Base = list(
        first_stage = paste0(
          pollutant, " ~ total_frp + FRP_u_wind + FRP_v_wind + ",
          "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
          "factor(month) + factor(year)"
        ),
        second_stage = "deaths ~ predicted + residuals + RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + factor(month) + factor(year)"
      ),
      
      # Extended model with additional weather controls
      Extended = list(
        first_stage = paste0(
          pollutant, " ~ total_frp + FRP_u_wind + FRP_v_wind + ",
          "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
          "WINDSPEEDAVG + WINDDIRECTIONAVG + ",
          "factor(month) + factor(year)"
        ),
        second_stage = "deaths ~ predicted + residuals + RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + WINDSPEEDAVG + WINDDIRECTIONAVG + factor(month) + factor(year)"
      ),
      
      # Non-linear model with squared term
      Nonlinear = list(
        first_stage = paste0(
          pollutant, " ~ total_frp + I(total_frp^2) + FRP_u_wind + FRP_v_wind + ",
          "RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + ",
          "factor(month) + factor(year)"
        ),
        second_stage = "deaths ~ predicted + I(predicted^2) + residuals + RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + factor(month) + factor(year)"
      )
    )
    
    # Run each specification
    for (spec_name in names(specifications)) {
      spec <- specifications[[spec_name]]
      
      # Filter data appropriately for each specification
      if (spec_name == "Extended") {
        # For Extended model, include wind variables and filter out NAs
        filtered_data <- data %>%
          filter(!is.na(deaths), !is.na(!!sym(pollutant)),
                 !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
                 !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG),
                 !is.na(WINDSPEEDAVG), !is.na(WINDDIRECTIONAVG))
      } else {
        # For other models, just basic filtering
        filtered_data <- data %>%
          filter(!is.na(deaths), !is.na(!!sym(pollutant)),
                 !is.na(total_frp), !is.na(FRP_u_wind), !is.na(FRP_v_wind),
                 !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))
      }
      
      cat("  Specification:", spec_name, "- Observations:", nrow(filtered_data), "\n")
      
      # First stage
      first_stage_model <- tryCatch({
        lm(as.formula(spec$first_stage), data = filtered_data)
      }, error = function(e) {
        cat("  Error in first stage for", spec_name, ":", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(first_stage_model)) next
      
      # Create predicted values and residuals
      temp_data <- filtered_data
      temp_data$predicted <- predict(first_stage_model)
      temp_data$residuals <- residuals(first_stage_model)
      
      # Second stage
      second_stage_model <- tryCatch({
        glm(as.formula(spec$second_stage), 
            family = poisson(link = "log"), 
            data = temp_data)
      }, error = function(e) {
        cat("  Error in second stage for", spec_name, ":", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(second_stage_model)) next
      
      # Calculate Newey-West standard errors
      nw_vcov <- NeweyWest(second_stage_model, lag = 4, prewhite = FALSE)
      robust_results <- coeftest(second_stage_model, vcov = nw_vcov)
      
      # Extract coefficient for predicted pollutant
      poll_index <- which(rownames(robust_results) == "predicted")
      poll_coef <- robust_results[poll_index, "Estimate"]
      poll_se <- robust_results[poll_index, "Std. Error"]
      poll_p <- robust_results[poll_index, "Pr(>|z|)"]
      
      # Calculate confidence intervals
      ci_lower <- poll_coef - 1.96 * poll_se
      ci_upper <- poll_coef + 1.96 * poll_se
      
      # Add to results
      all_results <- rbind(all_results, data.frame(
        Pollutant = pollutant,
        Specification = spec_name,
        Coefficient = poll_coef,
        SE = poll_se,
        P_value = poll_p,
        CI_Lower = ci_lower,
        CI_Upper = ci_upper,
        N = nobs(second_stage_model)
      ))
      
      # For the non-linear model, also extract squared term
      if (spec_name == "Nonlinear") {
        squared_index <- which(rownames(robust_results) == "I(predicted^2)")
        if (length(squared_index) > 0) {
          sq_coef <- robust_results[squared_index, "Estimate"]
          sq_se <- robust_results[squared_index, "Std. Error"]
          sq_p <- robust_results[squared_index, "Pr(>|z|)"]
          
          # Add squared term to results
          all_results <- rbind(all_results, data.frame(
            Pollutant = paste0(pollutant, "^2"),
            Specification = "Nonlinear",
            Coefficient = sq_coef,
            SE = sq_se,
            P_value = sq_p,
            CI_Lower = sq_coef - 1.96 * sq_se,
            CI_Upper = sq_coef + 1.96 * sq_se,
            N = nobs(second_stage_model)
          ))
        }
      }
    }
  }
  
  return(all_results)
}

# Run the analysis for all pollutants
all_pollutants_results <- run_robustness_checks_all_pollutants(merged_daily, pollutants)

# Format and print the results
formatted_results <- all_pollutants_results %>%
  mutate(
    Coefficient = format(round(Coefficient, 5), nsmall = 5),
    SE = format(round(SE, 5), nsmall = 5),
    P_value = ifelse(P_value < 0.01, "<0.01",
                     ifelse(P_value < 0.05, "<0.05",
                            ifelse(P_value < 0.1, "<0.1",
                                   format(round(P_value, 3), nsmall = 3)))),
    CI = paste0("[", format(round(CI_Lower, 5), nsmall = 5), ", ", 
                format(round(CI_Upper, 5), nsmall = 5), "]")
  ) %>%
  select(Pollutant, Specification, Coefficient, SE, P_value, CI, N)

# Create a summary table that only shows significant results (p < 0.1)
significant_results <- all_pollutants_results %>%
  filter(P_value < 0.1) %>%
  mutate(
    Coefficient = format(round(Coefficient, 5), nsmall = 5),
    SE = format(round(SE, 5), nsmall = 5),
    P_value = ifelse(P_value < 0.01, "<0.01",
                     ifelse(P_value < 0.05, "<0.05",
                            ifelse(P_value < 0.1, "<0.1",
                                   format(round(P_value, 3), nsmall = 3)))),
    CI = paste0("[", format(round(CI_Lower, 5), nsmall = 5), ", ", 
                format(round(CI_Upper, 5), nsmall = 5), "]")
  ) %>%
  select(Pollutant, Specification, Coefficient, SE, P_value, CI)

# Print significant results
if (nrow(significant_results) > 0) {
  kable(significant_results, 
        caption = "Significant Results (p < 0.1) for All Pollutants")
} else {
  cat("No significant results found at p < 0.1 level.\n")
}


# Create output directory
output_dir <- "Output"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to create coefficient comparison plot
plot_coefficient_comparison <- function(results_df) {
  # Filter basic pollutants (not squared terms)
  basic_pollutants <- results_df %>%
    filter(!grepl("\\^2", Pollutant))
  
  # Create plot
  p <- ggplot(basic_pollutants, aes(x = Specification, y = Coefficient, fill = Specification)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
    facet_wrap(~ Pollutant, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "lightgrey"),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = "Comparison of Coefficients Across Model Specifications",
      subtitle = "Error bars show 95% confidence intervals",
      x = "Model Specification",
      y = "Coefficient (Effect on Mortality)"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  
  return(p)
}

# Function to create significant results plot
plot_significant_effects <- function(results_df, significance_threshold = 0.1) {
  # Filter significant results
  significant_results <- results_df %>%
    filter(P_value < significance_threshold & !grepl("\\^2", Pollutant))
  
  # Check if we have any significant results
  if(nrow(significant_results) == 0) {
    return(NULL)
  }
  
  # Create plot
  p <- ggplot(significant_results, 
              aes(x = reorder(paste(Pollutant, Specification, sep = "-"), Coefficient), 
                  y = Coefficient, fill = Pollutant)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste0("Significant Effects of Air Pollutants on Mortality (p < ", significance_threshold, ")"),
      subtitle = "Sorted by effect size",
      x = "",
      y = "Coefficient (Effect on Mortality)"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")
  
  return(p)
}

# Function to create consistency heatmap
plot_consistency_heatmap <- function(results_df) {
  # Create consistency data
  consistency_data <- results_df %>%
    filter(!grepl("\\^2", Pollutant)) %>%
    select(Pollutant, Specification, Coefficient) %>%
    mutate(
      Direction = ifelse(Coefficient > 0, "Positive", "Negative"),
      Magnitude = abs(Coefficient)
    )
  
  # Plot consistency heatmap
  p <- ggplot(consistency_data, aes(x = Specification, y = Pollutant, fill = Direction, alpha = Magnitude)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = c("Positive" = "#4CAF50", "Negative" = "#F44336")) +
    scale_alpha(range = c(0.2, 1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Consistency of Effect Direction Across Model Specifications",
      x = "Model Specification",
      y = "Pollutant",
      fill = "Effect Direction",
      alpha = "Magnitude"
    )
  
  return(p)
}

# Function to create significance heatmap
plot_significance_heatmap <- function(results_df) {
  # Create significance data
  significance_data <- results_df %>%
    filter(!grepl("\\^2", Pollutant)) %>%
    mutate(
      Significance = ifelse(P_value < 0.01, "p < 0.01",
                            ifelse(P_value < 0.05, "p < 0.05",
                                   ifelse(P_value < 0.1, "p < 0.1", "p ≥ 0.1"))),
      Significance = factor(Significance, levels = c("p < 0.01", "p < 0.05", "p < 0.1", "p ≥ 0.1"))
    )
  
  # Plot significance heatmap
  p <- ggplot(significance_data, aes(x = Specification, y = Pollutant, fill = Significance)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = c(
      "p < 0.01" = "#1a9641",
      "p < 0.05" = "#a6d96a",
      "p < 0.1" = "#ffffbf",
      "p ≥ 0.1" = "#f1f1f1"
    )) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Statistical Significance of Pollutant Effects Across Model Specifications",
      x = "Model Specification",
      y = "Pollutant",
      fill = "Significance Level"
    )
  
  return(p)
}

# Function to create and save summary tables
create_summary_tables <- function(results_df) {
  # Create basic summary of all results
  all_results <- results_df %>%
    mutate(
      P_value_formatted = ifelse(P_value < 0.01, "<0.01",
                                 ifelse(P_value < 0.05, "<0.05",
                                        ifelse(P_value < 0.1, "<0.1",
                                               round(P_value, 3)))),
      CI = paste0("[", round(CI_Lower, 5), ", ", round(CI_Upper, 5), "]"),
      Significant = P_value < 0.1
    ) %>%
    select(Pollutant, Specification, Coefficient, SE, P_value_formatted, CI, Significant, N)
  
  # Save to CSV
  write.csv(all_results, file.path(output_dir, "all_pollutant_effects.csv"), row.names = FALSE)
  
  # Create summary of significant results
  significant_results <- results_df %>%
    filter(P_value < 0.1) %>%
    mutate(
      P_value_formatted = ifelse(P_value < 0.01, "<0.01",
                                 ifelse(P_value < 0.05, "<0.05", "<0.1")),
      CI = paste0("[", round(CI_Lower, 5), ", ", round(CI_Upper, 5), "]"),
      Direction = ifelse(Coefficient > 0, "Positive", "Negative")
    ) %>%
    select(Pollutant, Specification, Coefficient, SE, P_value_formatted, CI, Direction, N)
  
  # Save to CSV
  write.csv(significant_results, file.path(output_dir, "significant_pollutant_effects.csv"), row.names = FALSE)
  
  # Return the summaries
  return(list(
    all_results = all_results,
    significant_results = significant_results
  ))
}

# Main function to generate all visualizations
visualize_pollution_effects <- function(results_df) {
  # Step 1: Create and save individual plots
  cat("Creating coefficient comparison plot...\n")
  p1 <- plot_coefficient_comparison(results_df)
  ggsave(file.path(output_dir, "coefficient_comparison.png"), p1, width = 12, height = 8, dpi = 300)
  
  cat("Creating significant effects plot...\n")
  p2 <- plot_significant_effects(results_df)
  if(!is.null(p2)) {
    ggsave(file.path(output_dir, "significant_effects.png"), p2, width = 10, height = 8, dpi = 300)
  } else {
    cat("No significant effects found to plot.\n")
  }
  
  cat("Creating consistency heatmap...\n")
  p4 <- plot_consistency_heatmap(results_df)
  ggsave(file.path(output_dir, "effect_consistency.png"), p4, width = 10, height = 8, dpi = 300)
  
  cat("Creating significance heatmap...\n")
  p5 <- plot_significance_heatmap(results_df)
  ggsave(file.path(output_dir, "significance_heatmap.png"), p5, width = 10, height = 8, dpi = 300)
  
  # Step 2: Create combined plot (if all plots are available)
  plots_list <- list()
  if(!is.null(p2)) plots_list <- c(plots_list, list(p2))
  plots_list <- c(plots_list, list(p5, p4))
  
  if(length(plots_list) > 0) {
    cat("Creating combined summary plot...\n")
    combined_plot <- gridExtra::grid.arrange(grobs = plots_list, ncol = 1, 
                                             top = "Air Pollution Effects on Mortality - Summary")
    ggsave(file.path(output_dir, "combined_summary.png"), combined_plot, width = 12, height = 16, dpi = 300)
  }
  
  # Step 3: Create and save summary tables
  cat("Creating summary tables...\n")
  tables <- create_summary_tables(results_df)
  
  cat("All visualizations and tables have been saved to the 'Output' directory.\n")
  return(list(
    individual_plots = list(
      coefficient_comparison = p1,
      significant_effects = p2,
      consistency_heatmap = p4,
      significance_heatmap = p5
    ),
    tables = tables
  ))
}

results <- run_robustness_checks_all_pollutants(merged_daily, pollutants)
outputs <- visualize_pollution_effects(results)