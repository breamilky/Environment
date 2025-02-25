print("Loading packages...")
library(boot)
library(sandwich)
library(lmtest)
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(tidyr)
library(broom)
library(purrr)
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

# Capture the summary output to a file
capture.output(summary(merged_daily[, c("PM10AVG", "PM2.5AVG", "SO2_ugm3", "NO2_ugm3", 
                                        "NO_ugm3", "NOX_ugm3", "O3_ugm3", "CO_ugm3")]), 
               file = "Output/Summary_pollutants.txt")

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

# Create visualization
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

# Save results
write.csv(fs_result_sum, "Output/FS_result_sum.csv", row.names = FALSE)
write.csv(fs_test_sum, "Output/FS_test_sum.csv", row.names = FALSE)
write.csv(summary_df, "Output/SS_Boots_result.csv", row.names = FALSE)
ggsave("Output/Plot_boots_pollutant.png", p, width = 12, height = 6, dpi = 300)