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


#---------- ROBUSTNESS CHECKS ----------#

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