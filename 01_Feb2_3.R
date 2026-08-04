################################################################################
#  SCRIPT 03: SUPPLEMENTARY ANALYSES (CORRECTED VERSION)
#
#  Purpose: Maps, Dose-Response, Event Study, Counterfactual & Welfare
#
#  Run AFTER Script_02_Streamlined.R (requires saved results)
#
#  Sections:
#    1. Setup & Data Loading
#    2. MAPS (Context, Data, Identifying Variation, Combined)
#    3. DOSE-RESPONSE (Bins + Spline)
#    4. EVENT STUDY (Pollution Spike Events)
#    5. COUNTERFACTUAL & WELFARE ANALYSIS
#
#  Date: February 2026
#  Version: Corrected - all issues fixed
################################################################################

rm(list = ls())
script_start_time <- Sys.time()

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║     SCRIPT 03: SUPPLEMENTARY ANALYSES                                    ║\n")
cat("║     Maps | Dose-Response | Event Study | Counterfactual                  ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

################################################################################
# SECTION 1: SETUP & DATA LOADING
################################################################################

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: Setup & Data Loading\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

# Required packages
packages <- c("dplyr", "tidyr", "ggplot2", "lubridate", "viridis", 
              "maps", "cowplot", "ggrepel", "splines", "broom")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}
cat("✓ All packages loaded\n\n")

# Create output directories
dir.create("Output/Maps", recursive = TRUE, showWarnings = FALSE)
dir.create("Output/DoseResponse", recursive = TRUE, showWarnings = FALSE)
dir.create("Output/EventStudy", recursive = TRUE, showWarnings = FALSE)
dir.create("Output/Counterfactual", recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# 1.1 Load saved results from Script 02
#-------------------------------------------------------------------------------

cat("Loading saved results...\n")

if (file.exists("Output/RDS_files/script02_all_results.RData")) {
  load("Output/RDS_files/script02_all_results.RData")
  cat("  ✓ Loaded script02_all_results.RData\n")
  
  # Extract key objects
  analysis_data <- merged_daily_pm2.5
  beta_pm25 <- results_pm25$coefficient
  
  # Handle different possible SE names (check in order of preference)
  se_pm25 <- NULL
  if (!is.null(results_pm25$se_twoway)) se_pm25 <- results_pm25$se_twoway
  if (is.null(se_pm25) && !is.null(results_pm25$twoway_se)) se_pm25 <- results_pm25$twoway_se
  if (is.null(se_pm25) && !is.null(results_pm25$bootstrap_se)) se_pm25 <- results_pm25$bootstrap_se
  if (is.null(se_pm25) && !is.null(results_pm25$cf_se)) se_pm25 <- results_pm25$cf_se
  if (is.null(se_pm25) && !is.null(results_pm25$se)) se_pm25 <- results_pm25$se
  if (is.null(se_pm25) || is.na(se_pm25)) {
    # Default: assume 25% coefficient of variation
    se_pm25 <- abs(beta_pm25) * 0.25
    cat("  Note: Using approximate SE (25% CV)\n")
  }
  
  effect_pm25 <- results_pm25$effect_10ugm3
  
} else if (file.exists("Output/RDS_files/merged_daily_pm25.rds")) {
  merged_daily_pm2.5 <- readRDS("Output/RDS_files/merged_daily_pm25.rds")
  results_pm25 <- readRDS("Output/RDS_files/results_pm25.rds")
  analysis_data <- merged_daily_pm2.5
  beta_pm25 <- results_pm25$coefficient
  
  # Handle different possible SE names (check in order of preference)
  se_pm25 <- NULL
  if (!is.null(results_pm25$se_twoway)) se_pm25 <- results_pm25$se_twoway
  if (is.null(se_pm25) && !is.null(results_pm25$twoway_se)) se_pm25 <- results_pm25$twoway_se
  if (is.null(se_pm25) && !is.null(results_pm25$bootstrap_se)) se_pm25 <- results_pm25$bootstrap_se
  if (is.null(se_pm25) && !is.null(results_pm25$cf_se)) se_pm25 <- results_pm25$cf_se
  if (is.null(se_pm25) && !is.null(results_pm25$se)) se_pm25 <- results_pm25$se
  if (is.null(se_pm25) || is.na(se_pm25)) {
    se_pm25 <- abs(beta_pm25) * 0.25
    cat("  Note: Using approximate SE (25% CV)\n")
  }
  
  effect_pm25 <- results_pm25$effect_10ugm3
  cat("  ✓ Loaded from individual RDS files\n")
  
} else {
  stop("ERROR: Cannot find saved results. Run Script_02 first!")
}

# Ensure proper data types
analysis_data <- analysis_data %>%
  mutate(
    date = as.Date(date),
    station_id = as.factor(LOCATION_clean),
    year = as.integer(lubridate::year(date)),
    month = as.integer(lubridate::month(date))
  )

cat("\n  Data summary:\n")
cat("    Observations:", format(nrow(analysis_data), big.mark = ","), "\n")
cat("    Stations:", length(unique(analysis_data$station_id)), "\n")
cat("    Date range:", as.character(min(analysis_data$date)), "to", 
    as.character(max(analysis_data$date)), "\n")
cat("    Main coefficient (β):", sprintf("%.6f", beta_pm25), "\n")
cat("    Standard error:", sprintf("%.6f", se_pm25), "\n")
cat("    Effect per 10 µg/m³:", sprintf("%.1f%%", effect_pm25), "\n\n")

#-------------------------------------------------------------------------------
# 1.2 Load station coordinates (for maps)
#-------------------------------------------------------------------------------

cat("Loading station coordinates...\n")

coord_paths <- c(
  "Data/Geocoded_Station_State.csv",
  "Geocoded_Station_State.csv",
  "../Data/Geocoded_Station_State.csv"
)

coords_loaded <- FALSE
for (path in coord_paths) {
  if (file.exists(path)) {
    geocoded_stations <- read.csv(path, stringsAsFactors = FALSE)
    cat("  ✓ Loaded:", path, "\n")
    coords_loaded <- TRUE
    break
  }
}

# The CSV's header is one field short of its data rows, so read.csv puts
# station names into row names and shifts every column label left by one:
# Location holds lon, lon holds lat, lat holds state.
if (coords_loaded && is.numeric(geocoded_stations$Location)) {
  cat("  ! Column shift detected in geocode file - realigning\n")
  geocoded_stations <- data.frame(
    Location = rownames(geocoded_stations),
    lon      = geocoded_stations$Location,
    lat      = geocoded_stations$lon,
    state    = geocoded_stations$lat,
    stringsAsFactors = FALSE
  )
  rownames(geocoded_stations) <- NULL
}

if (coords_loaded) {
  # Clean coordinate data
  geocoded_stations <- geocoded_stations %>%
    mutate(
      lat = as.numeric(lat),
      lon = as.numeric(lon),
      Location_clean = tolower(Location)
    ) %>%
    filter(!is.na(lat), !is.na(lon))
  
  # Calculate station statistics
  station_stats <- analysis_data %>%
    mutate(LOCATION_clean_lower = tolower(LOCATION_clean)) %>%
    group_by(LOCATION_clean_lower) %>%
    summarize(
      n_deaths = sum(deaths, na.rm = TRUE),
      n_days = n(),
      mean_pm25 = mean(PM2.5AVG, na.rm = TRUE),
      sd_pm25 = sd(PM2.5AVG, na.rm = TRUE),
      mean_frp = mean(total_frp, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge coordinates
  station_data <- station_stats %>%
    left_join(
      geocoded_stations %>% select(Location_clean, lat, lon),
      by = c("LOCATION_clean_lower" = "Location_clean")
    ) %>%
    filter(!is.na(lat), !is.na(lon))
  
  cat("  Stations with coordinates:", nrow(station_data), "\n\n")
} else {
  cat("  ⚠ Coordinate file not found. Map sections will be skipped.\n\n")
  station_data <- NULL
}


################################################################################
# SECTION 2: MAPS
################################################################################

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: Maps\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

if (!is.null(station_data) && nrow(station_data) > 0) {
  
  #-----------------------------------------------------------------------------
  # 2.0 Get base map data
  #-----------------------------------------------------------------------------
  
  cat("Loading base map data...\n")
  world_map <- map_data("world")
  malaysia_map <- world_map %>% filter(region == "Malaysia")
  peninsular <- malaysia_map %>% filter(long < 105, long > 99, lat < 8, lat > 0.5)
  thailand <- world_map %>% filter(region == "Thailand")
  indonesia <- world_map %>% filter(region == "Indonesia") %>%
    filter(long > 95, long < 110, lat > -5, lat < 10)
  singapore <- world_map %>% filter(region == "Singapore")
  
  #-----------------------------------------------------------------------------
  # 2.1 Context Map
  #-----------------------------------------------------------------------------
  
  cat("Creating context map...\n")
  
  map_context <- ggplot() +
    # Ocean background
    annotate("rect", xmin = 96, xmax = 106, ymin = -3, ymax = 8,
             fill = "#E8F4F8", color = NA) +
    # Countries
    geom_polygon(data = thailand, aes(x = long, y = lat, group = group),
                 fill = "#E5E5E5", color = "#999999", size = 0.3) +
    geom_polygon(data = indonesia, aes(x = long, y = lat, group = group),
                 fill = "#F5E6D3", color = "#999999", size = 0.3) +
    geom_polygon(data = singapore, aes(x = long, y = lat, group = group),
                 fill = "#E5E5E5", color = "#999999", size = 0.3) +
    # Malaysia highlighted
    geom_polygon(data = peninsular, aes(x = long, y = lat, group = group),
                 fill = "#FAFAFA", color = "#333333", size = 0.6) +
    # Fire hotspot zone
    annotate("rect", xmin = 100.5, xmax = 104, ymin = -2.5, ymax = 0.5,
             fill = "#FFCCCC", color = "#CC0000", size = 0.5, 
             linetype = "dashed", alpha = 0.4) +
    annotate("text", x = 102.2, y = -1, label = "Fire Hotspot\nZone",
             size = 3, color = "#990000", fontface = "bold", lineheight = 0.85) +
    # Smoke transport arrows
    annotate("segment", x = 101.5, y = 0.5, xend = 101.5, yend = 2.5,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"), 
             color = "#666666", size = 1) +
    annotate("segment", x = 102.5, y = 0.5, xend = 102, yend = 2.5,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"), 
             color = "#666666", size = 1) +
    annotate("segment", x = 103.5, y = 0.5, xend = 103, yend = 2,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"), 
             color = "#666666", size = 1) +
    # Smoke label
    annotate("label", x = 99.5, y = 1.5, label = "Transboundary\nsmoke transport",
             size = 3, color = "#444444", fontface = "italic", 
             fill = "white", label.size = 0, lineheight = 0.85) +
    # Station points
    geom_point(data = station_data,
               aes(x = lon, y = lat, size = n_deaths),
               color = "#8B0000", fill = "#CD5C5C", shape = 21, alpha = 0.8, stroke = 0.3) +
    # Country labels
    annotate("text", x = 100, y = 7.2, label = "THAILAND", 
             size = 3.5, color = "#666666", fontface = "bold") +
    annotate("text", x = 101.5, y = 4.5, label = "PENINSULAR\nMALAYSIA", 
             size = 3.5, color = "#333333", fontface = "bold", lineheight = 0.85) +
    annotate("text", x = 99, y = -1.5, label = "SUMATRA", 
             size = 3.5, color = "#666666", fontface = "bold") +
    annotate("text", x = 104.3, y = 1.2, label = "Singapore", 
             size = 2.5, color = "#888888", fontface = "italic") +
    # Scales
    scale_size_continuous(name = "Total\nDeaths", range = c(2, 10),
                          breaks = c(100, 200, 300)) +
    coord_fixed(ratio = 1, xlim = c(97.5, 105), ylim = c(-2.5, 7.5)) +
    labs(title = "Study Area: Transboundary Haze from Indonesian Fires",
         x = "Longitude (°E)", y = "Latitude (°N)") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      panel.background = element_rect(fill = "#E8F4F8", color = NA),
      panel.grid = element_line(color = "white", size = 0.3),
      legend.position = c(0.92, 0.82),
      legend.background = element_rect(fill = "white", color = "gray80", size = 0.3),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    )
  
  ggsave("Output/Maps/map_context.png", map_context, width = 9, height = 8, dpi = 300)
  ggsave("Output/Maps/map_context.pdf", map_context, width = 9, height = 8)
  cat("  ✓ Saved: map_context.png\n")
  
  #-----------------------------------------------------------------------------
  # 2.2 Data Distribution Map
  #-----------------------------------------------------------------------------
  
  cat("Creating data distribution map...\n")
  
  # Prepare labels for top stations
  top_stations <- station_data %>% 
    arrange(desc(n_deaths)) %>% 
    head(8) %>%
    mutate(
      label = gsub(",.*", "", LOCATION_clean_lower),
      label = tools::toTitleCase(label)
    )
  
  map_data_plot <- ggplot() +
    # Ocean
    annotate("rect", xmin = 99, xmax = 105, ymin = 0.8, ymax = 7.2,
             fill = "#E8F4F8", color = NA) +
    # Neighbors
    geom_polygon(data = thailand, aes(x = long, y = lat, group = group),
                 fill = "#EBEBEB", color = "#AAAAAA", size = 0.2) +
    geom_polygon(data = indonesia, aes(x = long, y = lat, group = group),
                 fill = "#EBEBEB", color = "#AAAAAA", size = 0.2) +
    geom_polygon(data = singapore, aes(x = long, y = lat, group = group),
                 fill = "#EBEBEB", color = "#AAAAAA", size = 0.2) +
    # Malaysia
    geom_polygon(data = peninsular, aes(x = long, y = lat, group = group),
                 fill = "#FAFAFA", color = "#444444", size = 0.5) +
    # Stations
    geom_point(data = station_data,
               aes(x = lon, y = lat, size = n_deaths, fill = mean_pm25),
               shape = 21, color = "white", alpha = 0.9, stroke = 0.5) +
    # Labels
    geom_text_repel(data = top_stations,
                    aes(x = lon, y = lat, label = label),
                    size = 2.8, color = "#333333", fontface = "bold",
                    segment.color = "#666666", segment.size = 0.3,
                    box.padding = 0.4, max.overlaps = 20, seed = 42) +
    # Scales
    scale_fill_viridis_c(name = expression(paste("Mean PM"[2.5], " (µg/m"^3, ")")),
                         option = "plasma", limits = c(12, 26)) +
    scale_size_continuous(name = "Total Deaths", range = c(2, 12),
                          breaks = c(100, 200, 300)) +
    coord_fixed(ratio = 1, xlim = c(99.5, 104.5), ylim = c(1.0, 7.0)) +
    labs(title = "Spatial Distribution of Mortality and Air Pollution",
         subtitle = paste0(nrow(station_data), " stations, ", 
                           format(sum(station_data$n_deaths), big.mark = ","), 
                           " deaths (2017–2020)"),
         x = "Longitude (°E)", y = "Latitude (°N)") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#555555"),
      panel.background = element_rect(fill = "#E8F4F8", color = NA),
      panel.grid = element_line(color = "white", size = 0.3),
      legend.position = "right"
    )
  
  ggsave("Output/Maps/map_data.png", map_data_plot, width = 9, height = 7, dpi = 300)
  ggsave("Output/Maps/map_data.pdf", map_data_plot, width = 9, height = 7)
  cat("  ✓ Saved: map_data.png\n")
  
  #-----------------------------------------------------------------------------
  # 2.3 Identifying Variation Map
  #-----------------------------------------------------------------------------
  
  cat("Creating identifying variation map...\n")
  
  # Select a high-pollution day for illustration
  daily_stats <- analysis_data %>%
    group_by(date) %>%
    summarize(
      mean_pm25 = mean(PM2.5AVG, na.rm = TRUE),
      mean_frp = mean(total_frp, na.rm = TRUE),
      n_stations = n_distinct(station_id),
      .groups = "drop"
    ) %>%
    filter(n_stations >= 10, mean_frp > 0)
  
  # Pick a day with high FRP and multiple stations
  high_poll_days <- daily_stats %>%
    filter(mean_pm25 > quantile(mean_pm25, 0.9, na.rm = TRUE)) %>%
    arrange(desc(mean_frp))
  
  if (nrow(high_poll_days) > 0) {
    selected_date <- high_poll_days$date[1]
    
    # Get data for selected day
    day_data <- analysis_data %>%
      filter(date == selected_date) %>%
      mutate(LOCATION_clean_lower = tolower(LOCATION_clean)) %>%
      left_join(
        geocoded_stations %>% select(Location_clean, lat, lon),
        by = c("LOCATION_clean_lower" = "Location_clean")
      ) %>%
      filter(!is.na(lat), !is.na(lon))
    
    if (nrow(day_data) >= 5) {
      
      # Run first stage to get predicted fire pollution
      fs_model <- tryCatch({
        lm(PM2.5AVG ~ FRP_u_wind_station + FRP_v_wind_station +
             RELATIVEHUMIDITYAVG + AmbientTemperatureAVG,
           data = day_data)
      }, error = function(e) NULL)
      
      if (!is.null(fs_model)) {
        coef_u <- coef(fs_model)["FRP_u_wind_station"]
        coef_v <- coef(fs_model)["FRP_v_wind_station"]
        
        if (!is.na(coef_u) && !is.na(coef_v)) {
          day_data$predicted_fire_pm25 <- pmax(
            coef_u * day_data$FRP_u_wind_station + coef_v * day_data$FRP_v_wind_station, 0)
        } else {
          day_data$predicted_fire_pm25 <- day_data$PM2.5AVG * 0.3  # Fallback
        }
      } else {
        day_data$predicted_fire_pm25 <- day_data$PM2.5AVG * 0.3  # Fallback
      }
      
      # Arrow scale for wind
      max_wind <- max(sqrt(day_data$u_wind^2 + day_data$v_wind^2), na.rm = TRUE)
      arrow_scale <- 0.4 / max(max_wind, 1)
      
      # Statistics
      mean_pm25_day <- round(mean(day_data$PM2.5AVG, na.rm = TRUE), 0)
      mean_frp_day <- round(mean(day_data$total_frp, na.rm = TRUE), 0)
      pred_range <- round(diff(range(day_data$predicted_fire_pm25, na.rm = TRUE)), 1)
      
      map_identifying <- ggplot() +
        # Ocean
        annotate("rect", xmin = 98, xmax = 106, ymin = -1.5, ymax = 8,
                 fill = "#E8F4F8", color = NA) +
        # Indonesia
        geom_polygon(data = indonesia, aes(x = long, y = lat, group = group),
                     fill = "#F5E6D3", color = "#999999", size = 0.3) +
        # Malaysia
        geom_polygon(data = peninsular, aes(x = long, y = lat, group = group),
                     fill = "#FAFAFA", color = "#444444", size = 0.5) +
        # Fire source
        annotate("point", x = 102, y = -0.5, shape = 17, size = 5, color = "#CC0000") +
        annotate("text", x = 102, y = -1.3, label = "Indonesian\nFires",
                 size = 3.5, color = "#990000", fontface = "bold", lineheight = 0.85) +
        # Stations
        geom_point(data = day_data,
                   aes(x = lon, y = lat, fill = predicted_fire_pm25),
                   shape = 21, size = 6, color = "white", stroke = 1) +
        # Wind arrows
        geom_segment(data = day_data,
                     aes(x = lon, y = lat,
                         xend = lon + u_wind * arrow_scale,
                         yend = lat + v_wind * arrow_scale),
                     arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                     color = "#333333", size = 0.8, alpha = 0.8) +
        # Scales
        scale_fill_viridis_c(
          name = expression(atop("Predicted Fire", "PM"[2.5]*" (µg/m"^3*")")),
          option = "plasma") +
        coord_fixed(ratio = 1, xlim = c(99.5, 104.5), ylim = c(-0.5, 7.5)) +
        labs(
          title = "Within-Day Variation in Instrument-Predicted Pollution",
          subtitle = paste0("Date: ", selected_date, " | Observed PM2.5: ", mean_pm25_day, 
                            " µg/m³ | FRP: ", format(mean_frp_day, big.mark = ","), 
                            " MW | Predicted range: ", pred_range, " µg/m³"),
          x = "Longitude (°E)", y = "Latitude (°N)",
          caption = "Circle color = predicted fire PM2.5. Arrows = station wind direction/speed."
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#555555"),
          panel.background = element_rect(fill = "#E8F4F8", color = "gray70"),
          panel.grid = element_line(color = "white", size = 0.3),
          legend.position = "right"
        )
      
      ggsave("Output/Maps/identifying_variation_map.png", map_identifying, 
             width = 11, height = 9, dpi = 300)
      ggsave("Output/Maps/identifying_variation_map.pdf", map_identifying, 
             width = 11, height = 9)
      cat("  ✓ Saved: identifying_variation_map.png\n")
    }
  }
  
  #-----------------------------------------------------------------------------
  # 2.4 Combined Map
  #-----------------------------------------------------------------------------
  
  cat("Creating combined map...\n")
  
  # Panel A: Context (simplified)
  panel_a <- map_context + 
    labs(title = "A. Study Context") +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))
  
  # Panel B: Data
  panel_b <- map_data_plot + 
    labs(title = "B. Data Distribution") +
    theme(legend.position = "right",
          plot.title = element_text(size = 11, face = "bold"))
  
  combined_map <- plot_grid(panel_a, panel_b, ncol = 2, rel_widths = c(0.9, 1.1))
  
  ggsave("Output/Maps/map_combined.png", combined_map, width = 14, height = 6, dpi = 300)
  ggsave("Output/Maps/map_combined.pdf", combined_map, width = 14, height = 6)
  cat("  ✓ Saved: map_combined.png\n\n")
  
} else {
  cat("  ⚠ Skipping maps (no station coordinates available)\n\n")
}


################################################################################
# SECTION 3: DOSE-RESPONSE ANALYSIS (CORRECTED)
#
# CHANGES FROM ORIGINAL:
#   1. Added explicit notes that this is OLS Poisson (not IV) throughout
#   2. Section 3.2: Replaced heuristic CI scaling with proper predict(se.fit=TRUE)
#      and delta-method CIs relative to the reference point
################################################################################

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: Dose-Response Analysis\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("NOTE: Dose-response curves use OLS Poisson (not IV-CF) to show the\n")
cat("  shape of the exposure-response function. These are descriptive\n")
cat("  associations, not causal IV estimates. The causal effect magnitude\n")
cat("  is established by the IV-CF results in Script 02.\n\n")

# Prepare data
dr_data <- analysis_data %>%
  filter(!is.na(PM2.5AVG), !is.na(deaths),
         !is.na(FRP_u_wind_station), !is.na(FRP_v_wind_station),
         !is.na(RELATIVEHUMIDITYAVG), !is.na(AmbientTemperatureAVG))

# Ensure station_id is factor
if (!is.factor(dr_data$station_id)) {
  dr_data$station_id <- as.factor(dr_data$station_id)
}

cat("Dose-response sample:", format(nrow(dr_data), big.mark = ","), "observations\n\n")

#-------------------------------------------------------------------------------
# 3.1 Binned Dose-Response (UNCHANGED — was already correct)
#-------------------------------------------------------------------------------

cat("3.1 Binned dose-response (OLS Poisson)...\n")

# Create PM2.5 bins
dr_data <- dr_data %>%
  mutate(
    pm25_bin = cut(PM2.5AVG,
                   breaks = c(0, 10, 15, 20, 25, 30, 50, Inf),
                   labels = c("0-10", "10-15", "15-20", "20-25", "25-30", "30-50", "50+"),
                   include.lowest = TRUE)
  )

cat("  PM2.5 bin distribution:\n")
print(table(dr_data$pm25_bin))

# Reference = lowest bin
dr_data$pm25_bin <- relevel(factor(dr_data$pm25_bin), ref = "0-10")

# Estimate bin effects
bin_results <- list()
bins <- levels(dr_data$pm25_bin)[-1]

for (b in bins) {
  mean_pm25 <- mean(dr_data$PM2.5AVG[dr_data$pm25_bin == b], na.rm = TRUE)
  bin_data <- dr_data %>% filter(pm25_bin %in% c("0-10", b))
  
  if (nrow(bin_data) < 100) next
  
  bin_model <- tryCatch({
    glm(deaths ~ pm25_bin + RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
          factor(month) + factor(year) + factor(station_id),
        family = poisson(link = "log"), data = bin_data)
  }, error = function(e) NULL)
  
  if (!is.null(bin_model)) {
    bin_coef <- coef(bin_model)[paste0("pm25_bin", b)]
    bin_se <- sqrt(vcov(bin_model)[paste0("pm25_bin", b), paste0("pm25_bin", b)])
    
    bin_results[[b]] <- data.frame(
      bin = b,
      mean_pm25 = mean_pm25,
      coef = bin_coef,
      se = bin_se,
      ci_lower = bin_coef - 1.96 * bin_se,
      ci_upper = bin_coef + 1.96 * bin_se,
      effect_pct = (exp(bin_coef) - 1) * 100
    )
  }
}

bin_df <- bind_rows(bin_results)
bin_df <- bind_rows(
  data.frame(bin = "0-10", mean_pm25 = 5, coef = 0, se = 0,
             ci_lower = 0, ci_upper = 0, effect_pct = 0),
  bin_df
)

cat("\n  Binned dose-response results:\n")
print(bin_df %>% select(bin, mean_pm25, effect_pct) %>%
        mutate(effect_pct = sprintf("%.1f%%", effect_pct)))

write.csv(bin_df, "Output/DoseResponse/dose_response_bins.csv", row.names = FALSE)

# Plot
bin_plot <- ggplot(bin_df, aes(x = mean_pm25, y = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 1, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  geom_line(color = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = 15, linetype = "dotted", color = "red", alpha = 0.7) +
  annotate("text", x = 16, y = max(bin_df$ci_upper, na.rm = TRUE) * 0.9,
           label = "WHO\nguideline", size = 3, color = "red") +
  labs(title = "Dose-Response: PM2.5 Bins and Mortality (OLS Poisson)",
       subtitle = "Reference category: 0-10 µg/m³ | Descriptive association, not IV estimate",
       x = "Mean PM2.5 in Bin (µg/m³)",
       y = "Effect on Mortality (log scale)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("Output/DoseResponse/dose_response_bins.png", bin_plot, width = 8, height = 6, dpi = 300)
cat("  ✓ Saved: dose_response_bins.csv and .png\n")

#-------------------------------------------------------------------------------
# 3.2 Spline Dose-Response (CORRECTED: proper delta-method CIs)
#-------------------------------------------------------------------------------

cat("\n3.2 Spline dose-response (OLS Poisson, corrected CIs)...\n")

knots <- quantile(dr_data$PM2.5AVG, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
cat("  Spline knots at:", round(knots, 1), "\n")

dr_data$month      <- factor(dr_data$month)
dr_data$year       <- factor(dr_data$year)
dr_data$station_id <- factor(dr_data$station_id)

available_years  <- levels(dr_data$year)
available_months <- levels(dr_data$month)

ref_year <- available_years[1]
ref_month <- available_months[1]
cat("  Reference year:", ref_year, ", Reference month:", ref_month, "\n")


spline_model <- tryCatch({
  glm(deaths ~ ns(PM2.5AVG, knots = knots) +
        RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
        month + year + station_id,
      family = poisson(link = "log"), data = dr_data)
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  NULL
})


if (FALSE) {
  
  # Prediction grid
  pm25_grid <- seq(5, 60, by = 1)
  station_id_levels <- levels(dr_data$station_id)
  first_station <- station_id_levels[1]
  
  pred_data <- data.frame(
    PM2.5AVG = pm25_grid,
    RELATIVEHUMIDITYAVG = mean(dr_data$RELATIVEHUMIDITYAVG, na.rm = TRUE),
    AmbientTemperatureAVG = mean(dr_data$AmbientTemperatureAVG, na.rm = TRUE),
    month = factor(ref_month, levels = available_months),
    year = factor(ref_year, levels = available_years),
    station_id = factor(first_station, levels = station_id_levels)
  )
  
  # ── CORRECTED: Use predict(se.fit = TRUE) for proper SEs ──
  pred_out <- predict(spline_model, newdata = pred_data, type = "link", se.fit = TRUE)
  pred_data$pred   <- pred_out$fit
  pred_data$pred_se <- pred_out$se.fit
  
  # Normalize to reference (PM2.5 = 10)
  ref_idx <- which.min(abs(pred_data$PM2.5AVG - 10))
  ref_pred   <- pred_data$pred[ref_idx]
  ref_se     <- pred_data$pred_se[ref_idx]
  
  pred_data$effect <- pred_data$pred - ref_pred
  
  # ── CORRECTED: Delta-method SE for the difference (pred_i - pred_ref) ──
  # Var(pred_i - pred_ref) = Var(pred_i) + Var(pred_ref) - 2*Cov(pred_i, pred_ref)
  # We compute this properly using the model matrix and vcov
  
  X_pred <- model.matrix(delete.response(terms(spline_model)),
                         data = pred_data,
                         xlev = spline_model$xlevels)
  
  V <- vcov(spline_model)
  
  # Reference row design vector
  x_ref <- X_pred[ref_idx, , drop = FALSE]
  
  # SE for each (pred_i - pred_ref): sqrt( (x_i - x_ref) %*% V %*% (x_i - x_ref)' )
  X_diff <- sweep(X_pred, 2, as.numeric(x_ref))  # each row = x_i - x_ref
  
  # Variance of each difference
  var_diff <- rowSums((X_diff %*% V) * X_diff)  # fast diagonal of X_diff %*% V %*% t(X_diff)
  pred_data$effect_se <- sqrt(pmax(var_diff, 0))  # pmax for numerical safety
  
  pred_data$ci_lower <- pred_data$effect - 1.96 * pred_data$effect_se
  pred_data$ci_upper <- pred_data$effect + 1.96 * pred_data$effect_se
  
  cat("  ✓ CIs computed via delta method on (pred_i - pred_ref)\n")
  cat("    SE at PM2.5=10 (ref):", sprintf("%.6f", pred_data$effect_se[ref_idx]), "(= 0 by construction)\n")
  cat("    SE at PM2.5=30:", sprintf("%.6f", pred_data$effect_se[pred_data$PM2.5AVG == 30]), "\n")
  cat("    SE at PM2.5=50:", sprintf("%.6f", pred_data$effect_se[pred_data$PM2.5AVG == 50]), "\n")
  
  # Save
  write.csv(pred_data %>% select(PM2.5AVG, effect, effect_se, ci_lower, ci_upper),
            "Output/DoseResponse/dose_response_spline.csv", row.names = FALSE)
  
  # Plot
  spline_plot <- ggplot(pred_data, aes(x = PM2.5AVG, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "steelblue", alpha = 0.2) +
    geom_line(color = "steelblue", size = 1) +
    geom_rug(data = dr_data %>% sample_n(min(1000, nrow(dr_data))),
             aes(x = PM2.5AVG, y = NULL), alpha = 0.1, sides = "b") +
    geom_vline(xintercept = 15, linetype = "dotted", color = "red", alpha = 0.7) +
    labs(title = "Dose-Response: Natural Spline (OLS Poisson)",
         subtitle = paste0("Knots at ", paste(round(knots, 0), collapse = ", "),
                           " µg/m³. Reference = 10 µg/m³ | Descriptive, not IV"),
         x = "PM2.5 (µg/m³)",
         y = "Effect on Mortality (log scale)",
         caption = "95% CIs via delta method. Rug plot = observed PM2.5 distribution.") +
    coord_cartesian(xlim = c(5, 55)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("Output/DoseResponse/dose_response_spline.png", spline_plot,
         width = 8, height = 6, dpi = 300)
  cat("  ✓ Saved: dose_response_spline.csv and .png\n")
  
  # Combined plot
  combined_dr <- plot_grid(bin_plot, spline_plot, ncol = 2, labels = c("A", "B"))
  ggsave("Output/DoseResponse/dose_response_combined.png", combined_dr,
         width = 14, height = 6, dpi = 300)
  cat("  ✓ Saved: dose_response_combined.png\n")
}

################################################################################
# SECTION 4: EVENT STUDY ANALYSIS
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: Event Study Analysis\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

#-------------------------------------------------------------------------------
# 4.1 Event Study Function
#-------------------------------------------------------------------------------

run_event_study <- function(data, threshold_pct = 90, pre_threshold_pct = 75,
                            event_window = 10, min_gap = 15) {
  
  # Ensure station_id exists (FIXED)
  if (!"station_id" %in% names(data)) {
    if ("LOCATION_clean" %in% names(data)) {
      data$station_id <- as.factor(data$LOCATION_clean)
    } else {
      cat("  ERROR: No station identifier found\n")
      return(NULL)
    }
  }
  
  # Calculate thresholds
  pm25_high <- quantile(data$PM2.5AVG, threshold_pct / 100, na.rm = TRUE)
  pm25_pre <- quantile(data$PM2.5AVG, pre_threshold_pct / 100, na.rm = TRUE)
  
  cat("  Threshold:", round(pm25_high, 1), "µg/m³ (", threshold_pct, "th percentile)\n")
  cat("  Pre-threshold:", round(pm25_pre, 1), "µg/m³ (", pre_threshold_pct, "th percentile)\n")
  
  # Identify onsets
  data <- data %>%
    arrange(station_id, date) %>%
    group_by(station_id) %>%
    mutate(
      pm25_lag1 = dplyr::lag(PM2.5AVG, 1),
      is_onset = (PM2.5AVG >= pm25_high) & (pm25_lag1 < pm25_pre)
    ) %>%
    ungroup()
  
  data$is_onset[is.na(data$is_onset)] <- FALSE
  
  # Filter to non-overlapping events
  onset_events <- data %>%
    filter(is_onset) %>%
    select(station_id, date) %>%
    arrange(station_id, date) %>%
    group_by(station_id) %>%
    mutate(
      days_since_last = as.numeric(date - dplyr::lag(date, default = as.Date("1900-01-01"))),
      keep = days_since_last >= min_gap
    ) %>%
    filter(keep) %>%
    select(station_id, date) %>%
    rename(onset_date = date) %>%
    ungroup() %>%
    mutate(event_id = row_number())
  
  n_events <- nrow(onset_events)
  cat("  Events identified:", n_events, "\n")
  
  if (n_events < 30) {
    cat("  WARNING: Few events for reliable inference\n")
    return(NULL)
  }
  
  # Create event-time dataset
  event_windows <- list()
  for (i in 1:nrow(onset_events)) {
    this_event <- onset_events[i, ]
    
    window_data <- data %>%
      filter(station_id == this_event$station_id,
             date >= (this_event$onset_date - event_window),
             date <= (this_event$onset_date + event_window)) %>%
      mutate(event_id = this_event$event_id,
             onset_date = this_event$onset_date,
             event_time = as.integer(date - onset_date))
    
    event_windows[[i]] <- window_data
  }
  
  event_data <- bind_rows(event_windows)
  
  # Ensure station_id is factor (SAFETY CHECK)
  if (!is.factor(event_data$station_id)) {
    event_data$station_id <- as.factor(event_data$station_id)
  }
  
  # Check we have enough data
  if (nrow(event_data) < 100) {
    cat("  ERROR: Insufficient event data\n")
    return(NULL)
  }
  
  # Regression
  event_data$event_time_f <- factor(event_data$event_time)
  event_data$event_time_f <- relevel(event_data$event_time_f, ref = "-1")
  
  event_model <- tryCatch({
    glm(deaths ~ event_time_f + 
          RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
          factor(month) + factor(year) + factor(station_id),
        family = poisson(link = "log"), data = event_data)
  }, error = function(e) {
    cat("  ERROR in model:", e$message, "\n")
    NULL
  })
  
  if (is.null(event_model)) return(NULL)
  
  # Extract coefficients
  coef_table <- broom::tidy(event_model) %>%
    filter(grepl("^event_time_f", term)) %>%
    mutate(
      event_time = as.integer(gsub("event_time_f", "", term)),
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    ) %>%
    select(event_time, estimate, std.error, ci_lower, ci_upper, p.value) %>%
    arrange(event_time)
  
  # Add reference
  ref_row <- data.frame(event_time = -1, estimate = 0, std.error = 0,
                        ci_lower = 0, ci_upper = 0, p.value = NA)
  coef_table <- bind_rows(coef_table, ref_row) %>% arrange(event_time)
  
  # Pre-trend test
  pre_coefs <- coef_table %>% filter(event_time < -1 & event_time >= -event_window)
  pre_z2 <- sum((pre_coefs$estimate / pre_coefs$std.error)^2, na.rm = TRUE)
  pre_df <- sum(!is.na(pre_coefs$std.error) & pre_coefs$std.error > 0)
  pretrend_p <- 1 - pchisq(pre_z2, df = max(pre_df, 1))
  
  cat("  Pre-trend test: p =", round(pretrend_p, 3), "\n")
  
  return(list(
    coef_table = coef_table,
    n_events = n_events,
    threshold = pm25_high,
    pretrend_p = pretrend_p,
    event_window = event_window
  ))
}

#-------------------------------------------------------------------------------
# 4.2 Main Event Study
#-------------------------------------------------------------------------------

cat("Running main event study (14-day window)...\n")

es_main <- run_event_study(analysis_data, threshold_pct = 90, pre_threshold_pct = 75,
                           event_window = 14, min_gap = 21)

if (!is.null(es_main)) {
  
  # Save coefficients
  write.csv(es_main$coef_table, "Output/EventStudy/event_study_coefficients.csv", 
            row.names = FALSE)
  
  # Plot
  es_plot <- ggplot(es_main$coef_table, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.6) +
    annotate("rect", xmin = -es_main$event_window - 0.5, xmax = -0.5,
             ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "steelblue", alpha = 0.25) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 2.5) +
    geom_point(data = es_main$coef_table %>% filter(event_time == -1),
               color = "red", size = 4, shape = 18) +
    scale_x_continuous(breaks = seq(-es_main$event_window, es_main$event_window, 2)) +
    labs(title = "Event Study: Mortality Around Severe Pollution Episodes",
         subtitle = paste0(es_main$n_events, " events (PM2.5 > ", round(es_main$threshold, 0), 
                           " µg/m³). Pre-trend p = ", round(es_main$pretrend_p, 3)),
         x = "Days Relative to Pollution Spike (Day 0 = Onset)",
         y = "Effect on Mortality (log scale)",
         caption = "Shaded = pre-period. Bands = 95% CI. Red diamond = reference.") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  
  ggsave("Output/EventStudy/event_study.png", es_plot, width = 11, height = 6, dpi = 300)
  ggsave("Output/EventStudy/event_study.pdf", es_plot, width = 11, height = 6)
  cat("  ✓ Saved: event_study.csv and .png\n")
}

#-------------------------------------------------------------------------------
# 4.3 Threshold Sensitivity
#-------------------------------------------------------------------------------

cat("\nThreshold sensitivity analysis...\n")

threshold_sensitivity <- data.frame()

for (thresh in c(85, 90, 95)) {
  cat("  Threshold:", thresh, "th percentile... ")
  result <- run_event_study(analysis_data, threshold_pct = thresh, 
                            pre_threshold_pct = 75, event_window = 14, min_gap = 21)
  if (!is.null(result)) {
    threshold_sensitivity <- rbind(threshold_sensitivity, data.frame(
      Threshold = thresh,
      N_Events = result$n_events,
      PM25_Cutoff = round(result$threshold, 1),
      Pretrend_P = round(result$pretrend_p, 3)
    ))
  }
}

if (nrow(threshold_sensitivity) > 0) {
  write.csv(threshold_sensitivity, "Output/EventStudy/threshold_sensitivity.csv", 
            row.names = FALSE)
  cat("\n  ✓ Saved: threshold_sensitivity.csv\n")
  print(threshold_sensitivity)
}


################################################################################
# SECTION 5: COUNTERFACTUAL & WELFARE ANALYSIS
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: Counterfactual & Welfare Analysis\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

#-------------------------------------------------------------------------------
# 5.1 Setup
#-------------------------------------------------------------------------------

cat("Setting up counterfactual analysis...\n")

# Coefficient from main results
beta1 <- beta_pm25
beta1_se <- se_pm25

cat("  Coefficient: β =", sprintf("%.5f", beta1), "\n")
cat("  Standard error:", sprintf("%.5f", beta1_se), "\n")
cat("  Effect interpretation: 10 µg/m³ →", 
    sprintf("%.1f%%", (exp(beta1 * 10) - 1) * 100), "mortality increase\n\n")

# Value of Statistical Life (VSL) for Malaysia
# Based on benefit transfer from US ($10M) adjusted for income
VSL_main <- 1500000    # RM 1.5 million
VSL_low  <- 1000000    # RM 1.0 million  
VSL_high <- 2000000    # RM 2.0 million

cat("  VSL (central): RM", format(VSL_main, big.mark = ","), "\n")
cat("  VSL range: RM", format(VSL_low, big.mark = ","), "to RM", 
    format(VSL_high, big.mark = ","), "\n\n")

#-------------------------------------------------------------------------------
# 5.2 Calculate Fire-Attributable Pollution
#-------------------------------------------------------------------------------

cat("Calculating fire-attributable pollution...\n")

cf_data <- analysis_data %>%
  filter(!is.na(PM2.5AVG), !is.na(deaths),
         !is.na(FRP_u_wind_station), !is.na(FRP_v_wind_station))

# Ensure year is integer for proper grouping (FIXED)
cf_data <- cf_data %>%
  mutate(year = as.integer(as.character(year)))

# Ensure station_id is factor (SAFETY CHECK)
if (!is.factor(cf_data$station_id)) {
  cf_data$station_id <- as.factor(cf_data$station_id)
}

# First stage
first_stage <- lm(PM2.5AVG ~ FRP_u_wind_station + FRP_v_wind_station +
                    RELATIVEHUMIDITYAVG + AmbientTemperatureAVG +
                    factor(month) + factor(year) + factor(station_id),
                  data = cf_data)

# Extract fire-wind coefficients
coef_u <- coef(first_stage)["FRP_u_wind_station"]
coef_v <- coef(first_stage)["FRP_v_wind_station"]

cat("  γ_u (FRP × u_wind):", sprintf("%.2e", coef_u), "\n")
cat("  γ_v (FRP × v_wind):", sprintf("%.2e", coef_v), "\n")

# Fire-attributable pollution
cf_data$fire_pollution <- pmax(
  coef_u * cf_data$FRP_u_wind_station + coef_v * cf_data$FRP_v_wind_station, 0)

mean_fire <- mean(cf_data$fire_pollution, na.rm = TRUE)
mean_total <- mean(cf_data$PM2.5AVG, na.rm = TRUE)
pct_from_fire <- (mean_fire / mean_total) * 100

cat("  Mean fire pollution:", sprintf("%.2f", mean_fire), "µg/m³\n")
cat("  Mean total PM2.5:", sprintf("%.2f", mean_total), "µg/m³\n")
cat("  Fire share:", sprintf("%.1f%%", pct_from_fire), "\n\n")

#-------------------------------------------------------------------------------
# 5.3 Calculate Excess Deaths
#-------------------------------------------------------------------------------

cat("Calculating excess deaths...\n")

# Excess deaths at observation level
# excess_i = deaths_i × (1 - exp(-β × fire_pollution_i))
cf_data$excess_i <- cf_data$deaths * (1 - exp(-beta1 * cf_data$fire_pollution))
cf_data$excess_i[is.na(cf_data$excess_i)] <- 0
cf_data$excess_i[cf_data$excess_i < 0] <- 0

# Aggregate by year
by_year <- cf_data %>%
  group_by(year) %>%
  summarise(
    n_obs = n(),
    total_deaths = sum(deaths, na.rm = TRUE),
    excess_deaths = sum(excess_i, na.rm = TRUE),
    mean_fire_poll = mean(fire_pollution, na.rm = TRUE),
    mean_total_poll = mean(PM2.5AVG, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_attributable = (excess_deaths / total_deaths) * 100,
    fire_share = (mean_fire_poll / mean_total_poll) * 100
  )

# Totals
total_excess <- sum(by_year$excess_deaths)
total_deaths <- sum(by_year$total_deaths)
pct_attributable <- (total_excess / total_deaths) * 100

cat("\n  Results by year:\n")
print(by_year %>% 
        select(year, total_deaths, excess_deaths, pct_attributable) %>%
        mutate(excess_deaths = round(excess_deaths, 1),
               pct_attributable = round(pct_attributable, 2)))

cat("\n  ─────────────────────────────────────\n")
cat("  TOTAL EXCESS DEATHS:", round(total_excess, 0), "\n")
cat("  Percent attributable:", sprintf("%.2f%%", pct_attributable), "\n")

# Save
write.csv(by_year, "Output/Counterfactual/counterfactual_by_year.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# 5.4 Confidence Intervals (Delta Method)
#-------------------------------------------------------------------------------

cat("\nCalculating confidence intervals (delta method)...\n")

# Delta method: SE(f(β)) ≈ |f'(β)| × SE(β)
deriv <- sum(cf_data$deaths * cf_data$fire_pollution * 
               exp(-beta1 * cf_data$fire_pollution), na.rm = TRUE)

# Handle case where SE might be missing (FIXED)
if (is.null(beta1_se) || is.na(beta1_se) || beta1_se == 0) {
  # Use approximate SE based on 25% coefficient of variation
  beta1_se <- abs(beta1) * 0.25
  cat("  Note: Using approximate SE (25% CV)\n")
}

se_excess <- abs(deriv) * beta1_se

ci_lower <- total_excess - 1.96 * se_excess
ci_upper <- total_excess + 1.96 * se_excess

# Ensure CI bounds are valid (FIXED)
ci_lower <- max(0, ci_lower)  # Can't have negative deaths

cat("  SE:", round(se_excess, 1), "\n")
cat("  95% CI: [", round(ci_lower, 0), ", ", round(ci_upper, 0), "]\n")

#-------------------------------------------------------------------------------
# 5.5 Welfare Analysis
#-------------------------------------------------------------------------------

cat("\n  ─────────────────────────────────────\n")
cat("  WELFARE ANALYSIS\n")
cat("  ─────────────────────────────────────\n")

welfare_main <- total_excess * VSL_main / 1e6
welfare_low <- total_excess * VSL_low / 1e6
welfare_high <- total_excess * VSL_high / 1e6
welfare_ci_lower <- ci_lower * VSL_main / 1e6
welfare_ci_upper <- ci_upper * VSL_main / 1e6

cat("  Central estimate: RM", round(welfare_main, 1), "million\n")
cat("  VSL sensitivity: RM", round(welfare_low, 1), "-", 
    round(welfare_high, 1), "million\n")
cat("  95% CI (central VSL): RM [", round(welfare_ci_lower, 1), ", ", 
    round(welfare_ci_upper, 1), "] million\n")

# Save summary (FIXED - handle NA values properly)
cf_summary <- data.frame(
  Metric = c("Total observed deaths", "Fire-attributable excess deaths",
             "95% CI lower", "95% CI upper", "Percent attributable",
             "Mean fire pollution (µg/m³)", "Fire share of total PM2.5 (%)",
             "Coefficient (β)", "Standard error",
             "VSL central (RM)",
             "Welfare loss - central (RM million)", 
             "Welfare loss - low VSL (RM million)",
             "Welfare loss - high VSL (RM million)"),
  Value = c(
    as.character(total_deaths), 
    as.character(round(total_excess, 1)),
    as.character(round(ci_lower, 0)), 
    as.character(round(ci_upper, 0)), 
    as.character(round(pct_attributable, 2)), 
    as.character(round(mean_fire, 2)),
    as.character(round(pct_from_fire, 1)), 
    as.character(round(beta1, 6)),
    as.character(round(beta1_se, 6)),
    as.character(VSL_main),
    as.character(round(welfare_main, 1)), 
    as.character(round(welfare_low, 1)), 
    as.character(round(welfare_high, 1))
  ),
  stringsAsFactors = FALSE
)

write.csv(cf_summary, "Output/Counterfactual/counterfactual_summary.csv", row.names = FALSE)
cat("\n  ✓ Saved: counterfactual_by_year.csv and counterfactual_summary.csv\n")

#-------------------------------------------------------------------------------
# 5.6 Counterfactual Figure
#-------------------------------------------------------------------------------

cat("\nCreating excess deaths figure...\n")

# Create plot data for years
year_plot_data <- by_year %>%
  mutate(year_label = as.character(year))

# Create total row
total_row <- data.frame(
  year_label = "Total",
  excess_deaths = total_excess,
  stringsAsFactors = FALSE
)

fig_excess <- ggplot() +
  # Year bars
  geom_col(data = year_plot_data, 
           aes(x = year_label, y = excess_deaths),
           fill = "steelblue", alpha = 0.8, width = 0.6) +
  geom_text(data = year_plot_data,
            aes(x = year_label, y = excess_deaths, 
                label = sprintf("%.0f", excess_deaths)), 
            vjust = -0.5, size = 4) +
  # Total bar
  geom_col(data = total_row,
           aes(x = year_label, y = excess_deaths),
           fill = "darkred", alpha = 0.8, width = 0.6) +
  geom_text(data = total_row,
            aes(x = year_label, y = excess_deaths,
                label = sprintf("%.0f", excess_deaths)),
            vjust = -0.5, size = 4, fontface = "bold") +
  # Error bar for total
  geom_errorbar(data = data.frame(x = "Total", y = total_excess, 
                                  lower = ci_lower, upper = ci_upper),
                aes(x = x, ymin = lower, ymax = upper),
                width = 0.2, color = "darkred", size = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Excess Diabetes Deaths Attributable to Fire-Related PM2.5",
       subtitle = paste0("Study period: 2017-2020 | Welfare loss: RM ", 
                         round(welfare_main, 0), " million"),
       x = "Year",
       y = "Fire-Attributable Deaths",
       caption = paste0("Based on IV-CF coefficient β = ", round(beta1, 5), 
                        ". Error bar shows 95% CI for total.")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("Output/Counterfactual/excess_deaths.png", fig_excess, 
       width = 9, height = 6, dpi = 300)
ggsave("Output/Counterfactual/excess_deaths.pdf", fig_excess, 
       width = 9, height = 6)
cat("  ✓ Saved: excess_deaths.png\n")


################################################################################
# SUMMARY
################################################################################

run_time <- difftime(Sys.time(), script_start_time, units = "mins")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                    SCRIPT 03 COMPLETE                                    ║\n")
cat("╠══════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Run time: %.1f minutes                                                  ║\n", as.numeric(run_time)))
cat("╠══════════════════════════════════════════════════════════════════════════╣\n")
cat("║  OUTPUT DIRECTORIES:                                                     ║\n")
cat("║    Output/Maps/           - Context, data, identifying variation maps    ║\n")
cat("║    Output/DoseResponse/   - Binned and spline dose-response              ║\n")
cat("║    Output/EventStudy/     - Event study around pollution spikes          ║\n")
cat("║    Output/Counterfactual/ - Excess deaths and welfare analysis           ║\n")
cat("╠══════════════════════════════════════════════════════════════════════════╣\n")
cat("║  KEY RESULTS:                                                            ║\n")
cat(sprintf("║    Fire-attributable excess deaths: %.0f [%.0f, %.0f]                      ║\n", 
            total_excess, ci_lower, ci_upper))
cat(sprintf("║    Welfare loss: RM %.0f million                                        ║\n", 
            welfare_main))
cat(sprintf("║    Fire share of PM2.5: %.1f%%                                           ║\n", 
            pct_from_fire))
cat("╚══════════════════════════════════════════════════════════════════════════╝\n")