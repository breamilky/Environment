print("Loading packages...")
# Define the list of required packages
packages <- c("ggplot2", "dplyr", "tidyr", "viridis", "gridExtra", 
              "ggridges", "cowplot", "boot", "sandwich", "lmtest", 
              "knitr", "kableExtra", "broom", "purrr")

# Remove duplicates (if any)
packages <- unique(packages)

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

### 01. Clean and convert pollutants

# Removing implausible NO values 
merged_daily <- merged_daily %>%
  mutate(
    NOAVG_clean = ifelse(NOAVG < 0, NA, NOAVG),
    NOXAVG_clean = ifelse(NOXAVG < 0, NA, NOXAVG)
  )

# Convert Unit: ppm to μg/m³
# (ppm * molecular weight * 1000) / (0.08205 * (273.15 + T))
# T = temperature in Celsius
# Molecular weights (g/mol):
# SO2 = 64.066
# NO2 = 46.0055
# NO = 30.01
# O3 = 48
# CO = 28.01
merged_daily <- merged_daily %>%
  mutate(
    SO2_ugm3 = (SO2AVG * 64.066 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)),
    NO2_ugm3 = (NO2AVG * 46.0055 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)),
    NO_ugm3 = (NOAVG_clean * 30.01 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)),
    NOX_ugm3 = (NOXAVG_clean * 46.0055 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)), # using NO2 MW as reference
    O3_ugm3 = (O3AVG * 48 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG)),
    CO_ugm3 = (COAVG * 28.01 * 1000) / (0.08205 * (273.15 + AmbientTemperatureAVG))
  )

### 02. Create subdataset for each individual pollutant

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


### 03. First stage analysis

# Create a list to store the models
first_stage_models <- list()

# First stage for PM10AVG
formula_pm10 <- "PM10AVG ~ total_frp + FRP_u_wind + FRP_v_wind + 
                RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
                factor(month) + factor(year)"
model_pm10 <- lm(as.formula(formula_pm10), data = merged_daily_pm10)
first_stage_models[[1]] <- list(
  pollutant = "PM10AVG",
  model = model_pm10,
  n_obs = nrow(merged_daily_pm10)
)

# First stage for PM2.5AVG
formula_pm25 <- "PM2.5AVG ~ total_frp + FRP_u_wind + FRP_v_wind + 
                RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
                factor(month) + factor(year)"
model_pm25 <- lm(as.formula(formula_pm25), data = merged_daily_pm2.5)
first_stage_models[[2]] <- list(
  pollutant = "PM2.5AVG",
  model = model_pm25,
  n_obs = nrow(merged_daily_pm2.5)
)

# First stage for SO2_ugm3
formula_so2 <- "SO2_ugm3 ~ total_frp + FRP_u_wind + FRP_v_wind + 
               RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
               factor(month) + factor(year)"
model_so2 <- lm(as.formula(formula_so2), data = merged_daily_so2)
first_stage_models[[3]] <- list(
  pollutant = "SO2_ugm3",
  model = model_so2,
  n_obs = nrow(merged_daily_so2)
)

# First stage for NO2_ugm3
formula_no2 <- "NO2_ugm3 ~ total_frp + FRP_u_wind + FRP_v_wind + 
               RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
               factor(month) + factor(year)"
model_no2 <- lm(as.formula(formula_no2), data = merged_daily_no2)
first_stage_models[[4]] <- list(
  pollutant = "NO2_ugm3",
  model = model_no2,
  n_obs = nrow(merged_daily_no2)
)

# First stage for NO_ugm3
formula_no <- "NO_ugm3 ~ total_frp + FRP_u_wind + FRP_v_wind + 
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
              factor(month) + factor(year)"
model_no <- lm(as.formula(formula_no), data = merged_daily_no)
first_stage_models[[5]] <- list(
  pollutant = "NO_ugm3",
  model = model_no,
  n_obs = nrow(merged_daily_no)
)

# First stage for NOX_ugm3
formula_nox <- "NOX_ugm3 ~ total_frp + FRP_u_wind + FRP_v_wind + 
               RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
               factor(month) + factor(year)"
model_nox <- lm(as.formula(formula_nox), data = merged_daily_nox)
first_stage_models[[6]] <- list(
  pollutant = "NOX_ugm3",
  model = model_nox,
  n_obs = nrow(merged_daily_nox)
)

# First stage for O3_ugm3
formula_o3 <- "O3_ugm3 ~ total_frp + FRP_u_wind + FRP_v_wind + 
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
             factor(month) + factor(year)"
model_o3 <- lm(as.formula(formula_o3), data = merged_daily_o3)
first_stage_models[[7]] <- list(
  pollutant = "O3_ugm3",
  model = model_o3,
  n_obs = nrow(merged_daily_o3)
)

# First stage for CO_ugm3
formula_co <- "CO_ugm3 ~ total_frp + FRP_u_wind + FRP_v_wind + 
              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
              factor(month) + factor(year)"
model_co <- lm(as.formula(formula_co), data = merged_daily_co)
first_stage_models[[8]] <- list(
  pollutant = "CO_ugm3",
  model = model_co,
  n_obs = nrow(merged_daily_co)
)

# Extract and combine results
fs_coef_sum <- map_df(first_stage_models, function(model_info) {
  # Get summary statistics
  model_summary <- summary(model_info$model)
  
  # Extract key statistics
  data.frame(
    Pollutant = model_info$pollutant,
    N = model_info$n_obs,
    R_squared = model_summary$r.squared,
    Adj_R_squared = model_summary$adj.r.squared,
    F_statistic = model_summary$fstatistic[1],
    p_value = pf(model_summary$fstatistic[1], 
                 model_summary$fstatistic[2], 
                 model_summary$fstatistic[3], 
                 lower.tail = FALSE)
  )
})

# Extract coefficients for fire variables (total_frp, FRP_u_wind, FRP_v_wind)
fs_fire_sum <- map_df(first_stage_models, function(model_info) {
  # Get coefficients and standard errors
  coef_data <- tidy(model_info$model)
  
  # Filter for fire-related variables
  fire_vars <- coef_data %>%
    filter(term %in% c("total_frp", "FRP_u_wind", "FRP_v_wind"))
  
  # Pivot to wide format
  result <- data.frame(
    Pollutant = model_info$pollutant,
    total_frp_coef = fire_vars$estimate[fire_vars$term == "total_frp"],
    total_frp_p = fire_vars$p.value[fire_vars$term == "total_frp"],
    u_wind_coef = fire_vars$estimate[fire_vars$term == "FRP_u_wind"],
    u_wind_p = fire_vars$p.value[fire_vars$term == "FRP_u_wind"],
    v_wind_coef = fire_vars$estimate[fire_vars$term == "FRP_v_wind"],
    v_wind_p = fire_vars$p.value[fire_vars$term == "FRP_v_wind"]
  )
  
  return(result)
})


# If you want to extract specific coefficients from all models
extract_coefficient <- function(models, coef_name) {
  map_df(models, function(model_info) {
    coefs <- coef(summary(model_info$model))
    if(coef_name %in% rownames(coefs)) {
      data.frame(
        Pollutant = model_info$pollutant,
        Estimate = coefs[coef_name, "Estimate"],
        Std_Error = coefs[coef_name, "Std. Error"],
        t_value = coefs[coef_name, "t value"],
        p_value = coefs[coef_name, "Pr(>|t|)"]
      )
    } else {
      data.frame(
        Pollutant = model_info$pollutant,
        Estimate = NA,
        Std_Error = NA,
        t_value = NA,
        p_value = NA
      )
    }
  })
}

# Example: Extract total_frp coefficient from all models
fs_ttfire_sum <- extract_coefficient(first_stage_models, "total_frp")


# First stage tests for multiple pollutants
first_stage_multi <- function(data, pollutants = c("PM10AVG", "SO2_ugm3", "NO2_ugm3", 
                                                   "NO_ugm3", "NOX_ugm3", "O3_ugm3", 
                                                   "CO_ugm3", "PM2.5AVG")) {
  results <- list()
  
  for(pollutant in pollutants) {
    formula <- as.formula(paste(pollutant, "~ total_frp + FRP_u_wind + FRP_v_wind + 
                              RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
                              year + month"))
    
    model <- lm(formula, data = data)
    
    # F-statistic
    f_stat <- summary(model)$fstatistic[1]
    
    # R-squared
    r_squared <- summary(model)$r.squared
    
    # Add predictions
    data[[paste0("predicted_", pollutant)]] <- predict(model, newdata = data)
    
    results[[pollutant]] <- list(
      model = model,
      f_stat = f_stat,
      r_squared = r_squared
    )
  }
  
  return(list(
    models = results,
    data = data
  ))
}

# Summary of first stage results
summarize_first_stage <- function(first_stage_results) {
  pollutants <- names(first_stage_results$models)
  summary_df <- data.frame(
    pollutant = pollutants,
    f_stat = sapply(pollutants, function(p) first_stage_results$models[[p]]$f_stat),
    r_squared = sapply(pollutants, function(p) first_stage_results$models[[p]]$r_squared)
  )
  return(summary_df)
}

fs_results <- first_stage_multi(merged_daily)
fs_test_sum <- summarize_first_stage(fs_results)

# Create a function to extract key statistics from each model
extract_model_summary <- function(model_name, model_summary) {
  # Get the R-squared, adjusted R-squared, and F-statistic
  r_squared <- model_summary$r.squared
  adj_r_squared <- model_summary$adj.r.squared
  f_stat <- model_summary$fstatistic[1]
  df1 <- model_summary$fstatistic[2]
  df2 <- model_summary$fstatistic[3]
  p_value <- pf(f_stat, df1, df2, lower.tail = FALSE)
  n_obs <- df2 + df1 + 1
  
  # Extract coefficients for key variables (fire variables)
  coefs <- coef(model_summary)
  
  # Extract coefficients and p-values for fire-related variables
  total_frp_coef <- coefs["total_frp", "Estimate"]
  total_frp_se <- coefs["total_frp", "Std. Error"]
  total_frp_p <- coefs["total_frp", "Pr(>|t|)"]
  
  u_wind_coef <- coefs["FRP_u_wind", "Estimate"]
  u_wind_se <- coefs["FRP_u_wind", "Std. Error"]
  u_wind_p <- coefs["FRP_u_wind", "Pr(>|t|)"]
  
  v_wind_coef <- coefs["FRP_v_wind", "Estimate"]
  v_wind_se <- coefs["FRP_v_wind", "Std. Error"]
  v_wind_p <- coefs["FRP_v_wind", "Pr(>|t|)"]
  
  # Create data frame with all model statistics
  data.frame(
    Pollutant = model_name,
    Observations = n_obs,
    R_squared = r_squared,
    Adj_R_squared = adj_r_squared,
    F_statistic = f_stat,
    total_frp_coef = total_frp_coef,
    total_frp_se = total_frp_se,
    total_frp_p = total_frp_p,
    u_wind_coef = u_wind_coef,
    u_wind_se = u_wind_se,
    u_wind_p = u_wind_p,
    v_wind_coef = v_wind_coef,
    v_wind_se = v_wind_se,
    v_wind_p = v_wind_p,
    stringsAsFactors = FALSE
  )
}

# Create summary table from all models
create_summary_table <- function() {
  # Extract information from all models
  pm10_summary <- extract_model_summary("PM10AVG", summary(model_pm10))
  pm25_summary <- extract_model_summary("PM2.5AVG", summary(model_pm25))
  so2_summary <- extract_model_summary("SO2_ugm3", summary(model_so2))
  no2_summary <- extract_model_summary("NO2_ugm3", summary(model_no2))
  no_summary <- extract_model_summary("NO_ugm3", summary(model_no))
  nox_summary <- extract_model_summary("NOX_ugm3", summary(model_nox))
  o3_summary <- extract_model_summary("O3_ugm3", summary(model_o3))
  co_summary <- extract_model_summary("CO_ugm3", summary(model_co))
  
  # Combine all summaries
  all_summaries <- rbind(
    pm10_summary,
    pm25_summary,
    so2_summary,
    no2_summary,
    no_summary,
    nox_summary,
    o3_summary,
    co_summary
  )
  
  return(all_summaries)
}

# Generate the full summary table
fs_result_sum <- create_summary_table()

# Function to format the table for publication
format_for_publication <- function(summary_df) {
  # Create a formatted table with coefficient and standard errors
  formatted_table <- data.frame(
    Pollutant = summary_df$Pollutant,
    Observations = format(summary_df$Observations, big.mark = ","),
    `R²` = sprintf("%.3f", summary_df$R_squared),
    `Adj. R²` = sprintf("%.3f", summary_df$Adj_R_squared),
    `F-stat` = sprintf("%.1f", summary_df$F_statistic),
    
    `Fire Intensity (total_frp)` = sprintf("%.6f%s\n(%.6f)", 
                                           summary_df$total_frp_coef,
                                           ifelse(summary_df$total_frp_p < 0.01, "***", 
                                                  ifelse(summary_df$total_frp_p < 0.05, "**", 
                                                         ifelse(summary_df$total_frp_p < 0.1, "*", ""))),
                                           summary_df$total_frp_se),
    
    `E-W Wind (FRP_u_wind)` = sprintf("%.6f%s\n(%.6f)", 
                                      summary_df$u_wind_coef,
                                      ifelse(summary_df$u_wind_p < 0.01, "***", 
                                             ifelse(summary_df$u_wind_p < 0.05, "**", 
                                                    ifelse(summary_df$u_wind_p < 0.1, "*", ""))),
                                      summary_df$u_wind_se),
    
    `N-S Wind (FRP_v_wind)` = sprintf("%.6f%s\n(%.6f)", 
                                      summary_df$v_wind_coef,
                                      ifelse(summary_df$v_wind_p < 0.01, "***", 
                                             ifelse(summary_df$v_wind_p < 0.05, "**", 
                                                    ifelse(summary_df$v_wind_p < 0.1, "*", ""))),
                                      summary_df$v_wind_se)
  )
  
  return(formatted_table)
}

# Create the publication-ready table
publication_table <- format_for_publication(fs_result_sum)

# Print as a formatted table (can be copied directly)
kable(publication_table, format = "latex", booktabs = TRUE, align = "lrrrccc",
      caption = "First-Stage Regression Results for Air Pollutants") %>%
  kable_styling(latex_options = c("striped", "scale_down")) %>%
  add_header_above(c(" " = 5, "Coefficient Estimates (Standard Errors)" = 3)) %>%
  footnote(general = "Note: Each column represents a separate regression. All regressions include controls for relative humidity, ambient temperature, month fixed effects, and year fixed effects. Standard errors are in parentheses.",
           symbol = c("* p<0.1; ** p<0.05; *** p<0.01"))

# For markdown output (useful for reports/presentations)
kable(publication_table, format = "markdown",
      caption = "First-Stage Regression Results for Air Pollutants")

### 04. Second stage analysis

# Define dataset and pollutant mapping
datasets <- list(
  list(data = "merged_daily_co", pollutant = "CO_ugm3"),
  list(data = "merged_daily_no", pollutant = "NO_ugm3"),
  list(data = "merged_daily_no2", pollutant = "NO2_ugm3"),
  list(data = "merged_daily_nox", pollutant = "NOX_ugm3"),
  list(data = "merged_daily_o3", pollutant = "O3_ugm3"),
  list(data = "merged_daily_pm10", pollutant = "PM10AVG"),
  list(data = "merged_daily_pm2.5", pollutant = "PM2.5AVG"),
  list(data = "merged_daily_so2", pollutant = "SO2_ugm3")
)

# Boostrapping for the predictions
block_bootstrap_pollutant <- function(data, pollutant_col, block_size = 30, R = 100) {
  n <- nrow(data)
  n_blocks <- ceiling(n / block_size)
  block_coefs <- numeric(R)
  
  predicted_name <- paste0("predicted_", sub("_ugm3|AVG", "", pollutant_col))
  residuals_name <- paste0("residuals_", sub("_ugm3|AVG", "", pollutant_col))
  
  for(i in 1:R) {
    # Generate block indices
    block_starts <- sample(1:(n - block_size + 1), n_blocks, replace = TRUE)
    indices <- unlist(lapply(block_starts, function(x) x:(x + block_size - 1)))
    indices <- indices[1:n]
    boot_data <- data[indices, ]
    
    # First stage
    first_stage <- lm(as.formula(paste0(pollutant_col, " ~ total_frp + FRP_u_wind + FRP_v_wind + 
                     RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
                     factor(month) + factor(year)")), data = boot_data)
    
    boot_data[[predicted_name]] <- predict(first_stage)
    boot_data[[residuals_name]] <- residuals(first_stage)
    
    # Second stage
    second_stage <- glm(as.formula(paste0("deaths ~ ", predicted_name, " + ", residuals_name, " + 
                     RELATIVEHUMIDITYAVG + AmbientTemperatureAVG + 
                     factor(month) + factor(year)")),
                        family = poisson(link = "log"),
                        data = boot_data)
    
    block_coefs[i] <- coef(second_stage)[predicted_name]
  }
  
  return(list(
    pollutant = sub("_ugm3|AVG", "", pollutant_col),
    mean = mean(block_coefs),
    se = sd(block_coefs),
    ci = quantile(block_coefs, c(0.025, 0.975)),
    coefs = block_coefs
  ))
}

# Run bootstraps
set.seed(123)
all_results <- list()

for(ds in datasets) {
  cat("\nProcessing", ds$data, "\n")
  data <- get(ds$data)
  all_results[[ds$pollutant]] <- block_bootstrap_pollutant(data, ds$pollutant)
}

# Create summary table
summary_df <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    Pollutant = x$pollutant,
    Mean = round(x$mean, 6),
    SE = round(x$se, 6),
    CI_Lower = round(x$ci[1], 6),
    CI_Upper = round(x$ci[2], 6)
  )
}))

### 05. Visualization
                             
# Create data
plot_data <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    Pollutant = x$pollutant,
    Coefficient = x$coefs
  )
}))

p <- ggplot(plot_data, aes(x = Pollutant, y = Coefficient)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Block Bootstrap Coefficients Distribution by Pollutant",
       x = "Pollutant",
       y = "Coefficient Value")
#------------------------------------------------------------------------
# 1. Forest Plot with Confidence Intervals
#------------------------------------------------------------------------
forest_plot <- ggplot(summary_df, aes(x = reorder(Pollutant, Mean), y = Mean)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Effect Estimates with 95% Confidence Intervals",
       subtitle = "Instrumental Variable Estimates Using Indonesian Fires",
       x = "",
       y = "Coefficient (95% CI)") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80")
  ) +
  # Add text labels with exact values
  geom_text(aes(label = sprintf("%.6f [%.6f, %.6f]", Mean, CI_Lower, CI_Upper)), 
            hjust = -0.1, size = 3)
#------------------------------------------------------------------------
# 2. Density Ridgeline Plot for Distribution Comparison
#------------------------------------------------------------------------

# Create ridgeline plot to show distribution of bootstrap coefficients
ridgeline_plot <- ggplot(plot_data, aes(x = Coefficient, y = reorder(Pollutant, Coefficient, median), fill = Pollutant)) +
  geom_density_ridges(scale = 0.9, alpha = 0.7, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis_d() +
  theme_ridges() +
  labs(title = "Distribution of Bootstrap Coefficients by Pollutant",
       subtitle = "Solid lines represent medians",
       x = "Coefficient Value",
       y = "") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold")
  )

#------------------------------------------------------------------------
# 3. Combined Violin and Box Plot
#------------------------------------------------------------------------

# Create violin plot combined with boxplot
violin_box_plot <- ggplot(plot_data, aes(x = reorder(Pollutant, Coefficient, median), y = Coefficient, fill = Pollutant)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  scale_fill_viridis_d() +
  coord_flip() +
  labs(title = "Distribution of Bootstrap Coefficients",
       subtitle = "Violin plot with embedded boxplot",
       x = "",
       y = "Coefficient Value") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80")
  ) +
  # Add reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

#------------------------------------------------------------------------
# 4. Histogram Facets by Pollutant
#------------------------------------------------------------------------

# Create histogram facets for each pollutant
histogram_facets <- ggplot(plot_data, aes(x = Coefficient, fill = Pollutant)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  facet_wrap(~ Pollutant, scales = "free_y") +
  scale_fill_viridis_d() +
  labs(title = "Histogram of Bootstrap Coefficients by Pollutant",
       x = "Coefficient Value",
       y = "Count") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90", color = NA)
  ) +
  # Add vertical line at zero and at the mean for each facet
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "red") +
  geom_vline(data = summary_df, aes(xintercept = Mean), linetype = "solid", color = "darkblue")

#------------------------------------------------------------------------
# 5. Bubble Plot with Size Representing Precision
#------------------------------------------------------------------------

# Create bubble plot where size represents precision (inverse of standard error)
bubble_plot <- ggplot(summary_df, aes(x = reorder(Pollutant, Mean), y = Mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(aes(size = 1/SE), alpha = 0.7, color = "darkblue") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "darkblue") +
  coord_flip() +
  labs(title = "Effect Estimates with Precision Indicators",
       subtitle = "Larger bubbles indicate higher precision (lower standard error)",
       x = "",
       y = "Coefficient Value",
       size = "Precision (1/SE)") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80")
  )

#------------------------------------------------------------------------
# 6. Coefficient Plot with Significance Color Coding
#------------------------------------------------------------------------

# Add significance indicator to summary_df
summary_df$Significant <- ifelse(summary_df$CI_Lower > 0 | summary_df$CI_Upper < 0, "Yes", "No")

# Create coefficient plot with significance color coding
significance_plot <- ggplot(summary_df, aes(x = reorder(Pollutant, Mean), y = Mean, color = Significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  scale_color_manual(values = c("Yes" = "darkred", "No" = "darkgray")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_flip() +
  labs(title = "Effect Estimates with Significance Indication",
       subtitle = "Red points indicate statistically significant effects",
       x = "",
       y = "Coefficient Value",
       color = "Statistically\nSignificant") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray80")
  )

#------------------------------------------------------------------------
# 7. Combined Visualizations for Publication
#------------------------------------------------------------------------

# Combine forest plot and ridgeline plot for a comprehensive view
combined_plot <- plot_grid(
  forest_plot + theme(plot.margin = margin(5, 5, 5, 5)),
  ridgeline_plot + theme(plot.margin = margin(5, 5, 5, 5)),
  labels = c("A", "B"),
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1.2)
)

print("Writing results...")
# Save results
write.csv(fs_result_sum, "Output/FS_result_sum.csv", row.names = FALSE)
write.csv(fs_test_sum, "Output/FS_test_sum.csv", row.names = FALSE)
write.csv(fs_coef_sum, "Output/FS_coef_sum.csv", row.names = FALSE)
write.csv(fs_fire_sum, "Output/FS_fire_sum.csv", row.names = FALSE)
write.csv(fs_ttfire_sum, "Output/FS_ttfire_sum.csv", row.names = FALSE)
write.csv(summary_df, "Output/SS_Boots_result.csv", row.names = FALSE)
                             
ggsave("Output/Plot_boots_pollutant.png", p, width = 12, height = 6, dpi = 300)
ggsave("Output/forest_plot.png", forest_plot, width = 10, height = 7, dpi = 300)
ggsave("Output/ridgeline_plot.png", ridgeline_plot, width = 10, height = 7, dpi = 300)
ggsave("Output/violin_box_plot.png", violin_box_plot, width = 10, height = 7, dpi = 300)
ggsave("Output/histogram_facets.png", histogram_facets, width = 12, height = 8, dpi = 300)
ggsave("Output/bubble_plot.png", bubble_plot, width = 10, height = 7, dpi = 300)
ggsave("Output/significance_plot.png", significance_plot, width = 10, height = 7, dpi = 300)
ggsave("Output/combined_plot.png", combined_plot, width = 12, height = 12, dpi = 300)

capture.output(summary(merged_daily[, c("PM10AVG", "PM2.5AVG", "SO2_ugm3", "NO2_ugm3", 
                                        "NO_ugm3", "NOX_ugm3", "O3_ugm3", "CO_ugm3")]), 
               file = "Output/Summary_pollutants.txt")

print("Results successfully saved!")
